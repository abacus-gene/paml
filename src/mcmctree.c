/* mcmctree.c

   Markov chain Monte Carlo on trees (Bayesian phylogenetic analysis)

                   Ziheng YANG, since December 2002

   cc -o mcmctree -O3 mcmctree.c tools.c -lm
   cc -o infinitesites -DINFINITESITES -O3 mcmctree.c tools.c -lm
   cl -O2 -W3 -D_CRT_SECURE_NO_DEPRECATE mcmctree.c tools.c
   cl -FeInfiniteSites.exe -DINFINITESITES -W3 -D_CRT_SECURE_NO_DEPRECATE -O2 mcmctree.c tools.c

   mcmctree <ControlFileName>

   InfiniteSites <ControlFileName>
   FixedDsClock1.txt or FixedDsClock23.txt are hard-coded file names
*/

/*
#define INFINITESITES
*/

#include "paml.h"


#define NS            500
#define NBRANCH      (NS*2-2)
#define NNODE        (NS*2-1)
#define MAXNSONS      3
#define NGENE         8000          /* used for gnodes[NGENE] */
#define LSPNAME       100
#define NCODE         64
#define NCATG         50
#define MaxNFossils   200


double(*rndSymmetrical)(void);



extern int noisy, NFunCall;
extern char BASEs[];
extern double PjumpOptimum;

int GetOptions(char *ctlf);
int ReadTreeSeqs(FILE*fout);
int ProcessNodeAnnotation();
int ReadBlengthGH(char infile[]);
int GenerateBlengthGH(char infile[]);
int GetMem(void);
void FreeMem(void);
int UseLocus(int locus, int copycondP, int setmodel, int setSeqName);
int AcceptRejectLocus(int locus, int accept);
void switchconPin(void);
int SetPGene(int igene, int _pi, int _UVRoot, int _alpha, double xcom[]);
int DownSptreeSetTime(int inode);
void getSinvDetS(double space[]);
int GetInitials(void);
int GenerateGtree(int locus);
int printGtree(int printBlength);
int SetParameters(double x[]);
int ConditionalPNode(int inode, int igene, double x[]);
double lnpData(double lnpDi[]);
double lnpD_locus(int locus);
double lnpD_locus_Approx(int locus);
double lnptNCgiventC(void);
double lnptC(void);
double lnptCalibrationDensity(double t, int fossil, double p[7]);
int SetupPriorTimesFossilErrors(void);
double lnpriorTimesBDS_Approach1(void);
double lnpriorTimesTipDate(void);
double lnpriorTimes(void);
double lnpriorRates(void);
void copySptree(void);
void printStree(void);
double InfinitesitesClock(double *FixedDs);
double Infinitesites(FILE *fout);
int collectx(FILE* fout, double x[]);
int MCMC(FILE* fout);
int LabelOldCondP(int spnode);
int UpdateTimes(double *lnL, double steplength[], char accept[]);
int UpdateTimesClock23(double *lnL, double steplength[], char accept[]);
int UpdateParaRates(double *lnL, double steplength[], char accept[], double space[]);
int UpdateRates(double *lnL, double steplength[], char accept[]);
int UpdateParameters(double *lnL, double steplength[], char accept[]);
int mixing(double *lnL, double steplength, char *accept);
int mixingTipDate(double *lnL, double steplength, char *accept);
int mixingCladeStretch(double *lnL, double steplength, char *accept);
int UpdatePFossilErrors(double steplength, char *accept);
int getPfossilerr(double postEFossil[], double nround);
int DescriptiveStatisticsSimpleMCMCTREE(FILE *fout, char infile[], int nbin);

struct CommonInfo {
   unsigned char *z[NS];
   char *spname[NS];
   char seqf[2048], outf[2048], treef[2048], daafile[2048], mcmcf[2048], inBVf[2048], checkpointf[2048];
   char oldconP[NNODE];       /* update conP for node? (0 yes; 1 no) */
   int seqtype, ns, ls, ngene, posG[2], lgene[3], *pose, npatt, readpattern;
   int np, ncode, ntime, nrate, nrgene, nalpha, npi, ncatG, print, verbose, checkpoint;
   int cleandata, ndata;
   int model, clock, fix_kappa, fix_alpha, fix_rgene, Mgene;
   int method, icode, codonf, aaDist, NSsites;
   double *fpatt, kappa, alpha, TipDate, TipDate_TimeUnit;
   double rgene[NGENE], piG[NGENE][NCODE];  /* not used */
   double(*plfun)(double x[], int np), freqK[NCATG], rK[NCATG], *conP, *fhK;
   double pi[NCODE];
   int curconP;                       /* curconP = 0 or 1 */
   size_t sconP;
   double *conPin[2], space[100000];  /* change space[] to dynamic memory? */
   int conPSiteClass, NnodeScale;
   char *nodeScale;          /* nScale[ns-1] for interior nodes */
   double *nodeScaleF;       /* nScaleF[npatt] for scale factors */
}  com;

#if(0)

struct NODE { /* ipop is node number in species tree */
   int father, nson, sons[2], ibranch, ipop;
   double branch, age, label, label2, *conP;
   char fossil, *name, *annotation;
};

struct TREE {
   int  ns, rooted, root;            /* number of tips and root index */
   struct NODE *nodes;
} tree;

#else

struct TREEB {
   int  nbranch, nnode, root, branches[NBRANCH][2];
}  tree;
struct TREEN { /* ipop is node number in species tree */
   int father, nson, sons[2], ibranch, ipop;
   double branch, age, label, label2, *conP;
   char fossil, *name, *annotation;
}  *nodes, **gnodes, nodes_t[2 * NS - 1];
/* nodes_t[] is working space.  nodes is a pointer and moves around.
gnodes[] holds the gene trees, subtrees constructed from the master species
tree.  Each locus has a full set of rates (rates) for all branches on the
master tree, held in stree.nodes[].rates.  Branch lengths in the gene
tree are calculated by using those rates and the divergence times.

gnodes[][].label in the gene tree is used to store branch lengths estimated
under no clock when mcmc.usedata=2 (normal approximation to likelihood).
*/

#endif


struct SPECIESTREE {
   int nspecies, nbranch, nnode, root, nfossil, duplication;
   double RootAge[4];
   struct TREESPN {
      char name[LSPNAME + 1], fossil, usefossil;    /* fossil: 0, 1(L), 2(U), 3(B), 4(G) */
      int father, nson, sons[2], label;
      double age, pfossil[7];                       /* parameters in fossil distribution */
      double *rates;                                /* log rates for loci */
   } nodes[2 * NS - 1];
}    stree;
/* all trees are binary & rooted, with ancestors unknown. */


struct DATA { /* locus-specific data and tree information */
   int ns[NGENE], ls[NGENE], npatt[NGENE], ngene, lgene[NGENE];
   int root[NGENE + 1], conP_offset[NGENE];
   int priortime, priorrate;
   char datatype[NGENE], *z[NGENE][NS], cleandata[NGENE];
   double *zmorph[NGENE][NS * 2 - 1], *Rmorph[NGENE];
   double zpopvar[NGENE], ldetRm[NGENE]; /* MdR */
   double *fpatt[NGENE], lnpT, lnpR, lnpDi[NGENE], pi[NGENE][NCODE];
   double kappa[NGENE], alpha[NGENE];
   double BDS[4];  /* parameters in the birth-death-sampling model */
   double kappagamma[2], alphagamma[2];
   double pfossilerror[3], /* (p_beta, q_beta, NminCorrect) */ Pfossilerr, *CcomFossilErr;
   int    rgeneprior;         /* 0: gamma-Dirichlet; 1: conditional iid */
   double rgene[NGENE + 1], sigma2[NGENE + 1], rgenepara[3], sigma2para[3];
   double *blMLE[NGENE], *Gradient[NGENE], *Hessian[NGENE];
   int    transform;
}  data;

struct MCMCPARAMETERS {
   int burnin, nsample, sampfreq, usedata, saveconP, print;
   int nsteplength, steplengthOffset[8];
   char  *accept;
   double *steplength, *Pjump;
}  mcmc; /* control parameters */


char *models[] = { "JC69", "K80", "F81", "F84", "HKY85", "T92", "TN93", "REV" };
enum { JC69, K80, F81, F84, HKY85, T92, TN93, REV } MODELS;
enum { BASE, AA, CODON, MORPHC } DATATYPE;

int nR = 4;
double PMat[16], Cijk[64], Root[4];
double _rateSite = 1, OldAge = 999;
int LASTROUND = 0, BayesEB, testlnL = 0, NPMat = 0; /* no use for this */
int debug = 0;

double BFbeta = 0;

/* for stree.nodes[].fossil: lower, upper, bounds, gamma, inverse-gamma */
enum { LOWER_F = 1, UPPER_F, BOUND_F, GAMMA_F, SKEWN_F, SKEWT_F, S2N_F } FOSSIL_FLAGS;
char *fossils[] = { " ", "L", "U", "B", "G", "SN", "ST", "S2N" };
int npfossils[] = { 0,   4,   2,   4,   2,   3,    4,     7 };
char *clockstr[] = { "", "Global clock", "Independent rates", "Autocorrelated rates" };
enum { SQRT_B = 1, LOG_B, ARCSIN_B } B_Transforms;
char *Btransform[] = { "", "square root", "logarithm", "arcsin" };

#define MCMCTREE  1
#include "treesub.c"


int main(int argc, char *argv[])
{
   char ctlf[2048] = "mcmctree.ctl";
   int i, j, k = 4;
   FILE  *fout;

#if(0)
   printf("Select the proposal kernel\n");
   printf("  0: uniform\n  1: Triangle\n  2: Laplace\n  3: Normal\n  4: Bactrian\n  5: BactrianTriangle\n  6: BactrianLaplace\n");
   scanf("%d", &k);
#endif
   if (k == 0) { rndSymmetrical = rnduM0V1;             PjumpOptimum = 0.4; }
   if (k == 1) { rndSymmetrical = rndTriangle;          PjumpOptimum = 0.4; }
   if (k == 2) { rndSymmetrical = rndLaplace;           PjumpOptimum = 0.4; }
   if (k == 3) { rndSymmetrical = rndNormal;            PjumpOptimum = 0.4; }
   if (k == 4) { rndSymmetrical = rndBactrian;          PjumpOptimum = 0.3; }
   if (k == 5) { rndSymmetrical = rndBactrianTriangle;  PjumpOptimum = 0.3; }
   if (k == 6) { rndSymmetrical = rndBactrianLaplace;   PjumpOptimum = 0.3; }

   data.priortime = 0;  /* 0: BDS;        1: beta */
   data.priorrate = 0;  /* 0: LogNormal;  1: gamma */

   noisy = 3;
   com.alpha = 0.;     com.ncatG = 1;
   com.ncode = 4;      com.clock = 1;

   printf("MCMCTREE in %s\n", pamlVerStr);
   if (argc > 2 && !strcmp(argv[argc - 1], "--stdout-no-buf"))
      setvbuf(stdout, NULL, _IONBF, 0);
   if (argc > 1)
      strncpy(ctlf, argv[1], 127);

   data.BDS[0] = 1;    data.BDS[1] = 1;  data.BDS[2] = 0;
   strcpy(com.outf, "out");
   strcpy(com.mcmcf, "mcmc.txt");

   starttimer();
   GetOptions(ctlf);

   fout = gfopen(com.outf, "w");
   fprintf(fout, "MCMCTREE (%s) %s\n", pamlVerStr, com.seqf);

   ReadTreeSeqs(fout);

   if (data.pfossilerror && (data.pfossilerror[2]<0 || data.pfossilerror[2]>stree.nfossil))
      error2("nMinCorrect for fossil errors is out of range.");

   if (mcmc.usedata == 1) {
      if (com.seqtype != 0) error2("usedata = 1 for nucleotides only");
      if (com.alpha == 0)
         com.plfun = lfun;
      else {
         if (com.ncatG<2 || com.ncatG>NCATG) error2("ncatG");
         com.plfun = lfundG;
      }
      if (com.model > HKY85)  error2("Only HKY or simpler models are allowed.");
      if (com.model == JC69 || com.model == F81) { com.fix_kappa = 1; com.kappa = 1; }
   }
   else if (mcmc.usedata == 3) {
      GenerateBlengthGH("out.BV");  /* this is used so that in.BV is not overwritten */
      exit(0);
   }

   /* Do we want RootAge constraint at the root if (com.clock==1)? */
   if (com.clock != 1 && stree.RootAge[1] <= 0 && stree.nodes[stree.root].fossil <= LOWER_F)
      error2("set RootAge in control file when there is no upper bound on root");
   if (data.pfossilerror[0] == 0 && !stree.nodes[stree.root].fossil) {
      if (stree.RootAge[1] <= 0)
         error2("set RootAge in control file when there is no upper bound on root");

      stree.nodes[stree.root].fossil = (stree.RootAge[0] > 0 ? BOUND_F : UPPER_F);
      for (i = 0; i < 4; i++)
         stree.nodes[stree.root].pfossil[i] = stree.RootAge[i];
   }
   if (com.TipDate_TimeUnit == 0)
      printf("\nFossil calibration information used.\n");
   for (i = stree.nspecies; i < stree.nspecies * 2 - 1; i++) {
      if ((k = stree.nodes[i].fossil) == 0) continue;
      printf("Node %3d: %3s ( ", i + 1, fossils[k]);
      for (j = 0; j < npfossils[k]; j++) {
         printf("%7.4f", stree.nodes[i].pfossil[j + (k == UPPER_F)]);
         printf("%s", (j == npfossils[k] - 1 ? " )\n" : ", "));
      }
   }

#if(defined INFINITESITES)
   Infinitesites(fout);
#else
   if (mcmc.print < 0) {
      mcmc.print *= -1;
      fputs("\nSpecies tree for FigTree.", fout);
      copySptree();
      FPN(fout); OutTreeN(fout, 1, PrNodeNum); FPN(fout);
      DescriptiveStatisticsSimpleMCMCTREE(fout, com.mcmcf, 1);
   }
   else
      MCMC(fout);
#endif
   fclose(fout);
   exit(0);
}


int GetMem(void)
{
   /* This allocates memory for conditional probabilities (conP).
      gnodes[locus] is not allocated here but in GetGtree().

      Conditional probabilities for internal nodes are com.conPin[2], allocated
      according to data.ns[locus] and data.npatt[locus] at all loci.  Two copies
      of the space are allocated, hence the [2].  The copy used for the current
      gene trees is com.conPin[com.curconP] while the copy for proposed gene trees
      is com.conPin[!com.curconP].  data.conP_offset[locus] marks the starting
      position in conPin[] for each locus.

      Memory arrangement if(com.conPSiteClass=1):
      ncode*npatt for each node, by node, by iclass, by locus
   */
   int g = data.ngene, g1 = g + (data.rgeneprior == 1);
   int locus, j, k, s = stree.nspecies, s1, sG = 1, sfhK = 0;
   double *conP, *rates;

   /* get mem for conP (internal nodes) for sequence loci.  Morph loci have npatt = 0.
   */
   if (mcmc.usedata == 1) {
      if (!com.fix_alpha && mcmc.saveconP) {
         com.conPSiteClass = 1;  sG = com.ncatG;
      }
      data.conP_offset[0] = 0;
      for (locus = 0, com.sconP = 0; locus < g; locus++) {
         s1 = com.ncode * data.npatt[locus];
         com.sconP += sG*s1*(data.ns[locus] - 1) * sizeof(double);
         if (locus < g - 1)
            data.conP_offset[locus + 1] = data.conP_offset[locus] + sG*s1*(data.ns[locus] - 1);
      }

      if ((com.conPin[0] = (double*)malloc(2 * com.sconP)) == NULL)
         error2("oom conP");

      com.conPin[1] = com.conPin[0] + com.sconP / sizeof(double);
      printf("\n%zu bytes for conP\n", 2 * com.sconP);

      /* set gnodes[locus][].conP for tips and internal nodes */
      com.curconP = 0;
      for (locus = 0; locus < g; locus++) {
         conP = com.conPin[0] + data.conP_offset[locus];
         for (j = data.ns[locus]; j < data.ns[locus] * 2 - 1; j++, conP += com.ncode*data.npatt[locus])
            gnodes[locus][j].conP = conP;
         if (!data.cleandata[locus]) {
            /* Is this call to UseLocus still needed?  Ziheng 28/12/2009 */
            UseLocus(locus, 0, mcmc.usedata, 0);
         }
      }

      if (!com.fix_alpha) {
         for (locus = 0; locus < g; locus++)
            sfhK = max2(sfhK, data.npatt[locus]);
         sfhK *= com.ncatG * sizeof(double);
         if ((com.fhK = (double*)realloc(com.fhK, sfhK)) == NULL) error2("oom");
      }
   }
   else if (mcmc.usedata == 2) { /* allocate data.Gradient & data.Hessian */
      for (locus = 0, k = 0; locus < data.ngene; locus++)
         k += (2 * data.ns[locus] - 1)*(2 * data.ns[locus] - 1 + 2);
      if ((com.conPin[0] = (double*)malloc(k * sizeof(double))) == NULL)
         error2("oom g & H");
      for (j = 0; j < k; j++)  com.conPin[0][j] = -1;
      for (locus = 0, j = 0; locus < data.ngene; locus++) {
         data.blMLE[locus] = com.conPin[0] + j;
         data.Gradient[locus] = com.conPin[0] + j + (2 * data.ns[locus] - 1);
         data.Hessian[locus] = com.conPin[0] + j + (2 * data.ns[locus] - 1) * 2;
         j += (2 * data.ns[locus] - 1)*(2 * data.ns[locus] - 1 + 2);
      }
   }

   if (com.clock > 1) {  /* space for rates */
      s1 = (s * 2 - 1)*g * sizeof(double);
      if (noisy) printf("%d bytes for rates.\n", s1);
      if ((rates = (double*)malloc(s1)) == NULL) error2("oom for rates");
      for (j = 0; j < s * 2 - 1; j++)
         stree.nodes[j].rates = rates + g*j;
   }

#if(!defined INFINITESITES)
   /* setting up space and offset for MCMC steplengths (steplength) */
   mcmc.nsteplength = (s - 1);                   /* t (node ages)  */
   mcmc.nsteplength += g1 + g1*(com.clock > 1);    /* mu[] & sigma2[] */
   if (com.clock > 1)      mcmc.nsteplength += g*(s * 2 - 2);   /* branch rates */
   if (mcmc.usedata == 1)  mcmc.nsteplength += g*(!com.fix_kappa + !com.fix_alpha);   /* subst paras */
   mcmc.nsteplength += 1 + (data.pfossilerror[0] > 0);    /* mixing & fossilerror */

   mcmc.steplength = (double*)malloc(mcmc.nsteplength * 3 * sizeof(double));
   if (mcmc.steplength == NULL)  error2("oom steplength");
   memset(mcmc.steplength, 0, mcmc.nsteplength * 3 * sizeof(double));
   mcmc.Pjump = mcmc.steplength + mcmc.nsteplength;
   mcmc.accept = (char*)(mcmc.Pjump + mcmc.nsteplength);
   mcmc.steplengthOffset[0] = 0;                 /* t */
   mcmc.steplengthOffset[1] = s - 1;             /* mu[] & sigma2[] */
   mcmc.steplengthOffset[2] = mcmc.steplengthOffset[1] + g1 + g1*(com.clock > 1);              /* branch rates */
   mcmc.steplengthOffset[3] = mcmc.steplengthOffset[2] + (com.clock > 1)*g*(s * 2 - 2);        /* subst paras */
   mcmc.steplengthOffset[4] = mcmc.steplengthOffset[3] + (mcmc.usedata == 1)*g*(!com.fix_kappa + !com.fix_alpha);  /* mixing */
   mcmc.steplengthOffset[5] = mcmc.steplengthOffset[4] + 1;                                    /* fossilerror */
#endif

   return(0);
}

void FreeMem(void)
{
   int locus, j;

   for (locus = 0; locus < data.ngene; locus++)
      free(gnodes[locus]);
   free(gnodes);
   if (mcmc.usedata)
      free(com.conPin[0]);
   if (mcmc.usedata == 1) {
      for (locus = 0; locus < data.ngene; locus++) {
         if (data.datatype[locus] == MORPHC) {
            free(data.zmorph[locus][0]);
            free(data.Rmorph[locus]);
         }
         else {
            for (j = 0; j < data.ns[locus]; j++)
               free(data.z[locus][j]);
            free(data.fpatt[locus]);
         }
      }
   }
   if (com.clock > 1)
      free(stree.nodes[0].rates);

   if (mcmc.usedata == 1 && com.alpha)
      free(com.fhK);

   free(mcmc.steplength);
}


int SaveMCMCstate(char *filename, int ir, double lnpR, double lnL)
{
/* This saves stree, the gene trees, and the MCMC step lengths.
*/
   FILE *f = (FILE*)fopen(filename, "wb");
   int i,j, s = stree.nspecies, g = data.ngene, g1 = g + (data.rgeneprior == 1);
   double *tmp = (double*)malloc((s*2-1)*sizeof(double));

   printf("\nSaving MCMC state to file %s...\n", filename);
   if (f == NULL) error2("file open error");
   if (tmp == NULL) error2("oom in SaveMCMCState()");
   fwrite(&ir, sizeof(int), 1, f);
   for (i = 0; i < s - 1; i++) tmp[i] = stree.nodes[s + i].age;
   fwrite(tmp, (s - 1) * sizeof(double), 1, f);
   if (com.clock > 1)  /* rates for clock2 & clock3 */
      fwrite(stree.nodes[0].rates, (s * 2 - 1)*g*sizeof(double), 1, f);
   for (i = 0; i < g; i++) {
      for (j = 0; j < data.ns[i] - 1; j++)
         tmp[j] = gnodes[i][data.ns[i]+j].age;
      fwrite(tmp, (data.ns[i] - 1)*sizeof(double), 1, f);
   }
   
   fwrite(data.rgene, g1*sizeof(double), 1, f);
   if (com.clock > 1)
      fwrite(data.sigma2, g1*sizeof(double), 1, f);
   if (mcmc.usedata == 1 && !com.fix_kappa)
      fwrite(data.kappa, g*sizeof(double), 1, f);
   if (mcmc.usedata == 1 && !com.fix_alpha)
      fwrite(data.alpha, g*sizeof(double), 1, f);

   fwrite(mcmc.steplength, mcmc.nsteplength * sizeof(double), 1, f);
   fwrite(&lnpR, sizeof(double), 1, f);
   fwrite(&lnL, sizeof(double), 1, f);
   fclose(f);
   free(tmp);
   printf("lnpR = %.6f lnL = %.6f\n", lnpR, lnL);
   return(0);
}

int ReadMCMCstate(char *filename)
{
/* This reads the MCMC state saved earlier.
*/
   FILE *f = (FILE*)fopen(filename, "rb");
   int ir, i, j, s = stree.nspecies, g = data.ngene, g1 = g + (data.rgeneprior == 1);
   double lnpR = 0, lnL = 0, *tmp = (double*)malloc((s * 2 - 1) * sizeof(double));

   printf("\nReading MCMC state from %s...\n", filename);
   if (f == NULL) error2("file open error");
   if (tmp == NULL) error2("oom in ReadMCMCState()");
   fread(&ir, sizeof(int), 1, f);
   fread(tmp, (s-1)* sizeof(double), 1, f);
   for (i = 0; i < s - 1; i++) stree.nodes[s + i].age = tmp[i];
   if (com.clock > 1)  /* rates for clock2 & clock3 */
      fread(stree.nodes[0].rates, (s * 2 - 1)*g * sizeof(double), 1, f);
   for (i = 0; i < g; i++) {
      fread(tmp, (data.ns[i] - 1) * sizeof(double), 1, f);
      for (j = 0; j < data.ns[i] - 1; j++)
         gnodes[i][data.ns[i] + j].age = tmp[j];
   }

   fread(data.rgene, g1 * sizeof(double), 1, f);
   if (com.clock > 1)
      fread(data.sigma2, g1 * sizeof(double), 1, f);
   if (mcmc.usedata == 1 && !com.fix_kappa)
      fread(data.kappa, g * sizeof(double), 1, f);
   if (mcmc.usedata == 1 && !com.fix_alpha)
      fread(data.alpha, g * sizeof(double), 1, f);

   fread(mcmc.steplength, mcmc.nsteplength*sizeof(double), 1, f);
   fread(&lnpR, sizeof(double), 1, f);
   fread(&lnL, sizeof(double), 1, f);
   fclose(f);
   free(tmp);
   printf("lnpR = %.6f lnL = %.6f\n", lnpR, lnL);

   return(0);
}


int UseLocus(int locus, int copyconP, int setModel, int setSeqName)
{
   /* MCMCtree:
      This point nodes to the gene tree at locus gnodes[locus] and set com.z[]
      etc. for likelihood calculation for the locus.  Note that the gene tree
      topology (gnodes[]) is never copied, but nodes[].conP are repositioned in the
      algorithm.  The pointer for root gnodes[][com.ns].conP is assumed to be the
      start of the whole block for the locus.
      If (copyconP && mcmc.useData), the conP for internal nodes point
      to a fixed place (indicated by data.conP_offset[locus]) in the alternative
      space com.conPin[!com.curconP].  Note that the conP for each locus uses the
      correct space so that this routine can be used by all the proposal steps,
      some of which operates on one locus and some change all loci.

      Try to replace this with UseLocus() for B&C.
   */
   int i, s1 = com.ncode*data.npatt[locus], sG = (com.conPSiteClass ? com.ncatG : 1);
   double *conPt = com.conPin[!com.curconP] + data.conP_offset[locus];

   com.ns = data.ns[locus];
   com.ls = data.ls[locus];
   tree.root = data.root[locus];
   tree.nnode = 2 * com.ns - 1;
   nodes = gnodes[locus];
   if (copyconP && mcmc.usedata == 1) { /* this preserves the old conP. */
      memcpy(conPt, gnodes[locus][com.ns].conP, sG*s1*(com.ns - 1) * sizeof(double));
      for (i = com.ns; i < tree.nnode; i++)
         nodes[i].conP = conPt + (i - com.ns)*s1;
   }

   if (setModel && mcmc.usedata == 1) {
      com.cleandata = data.cleandata[locus];
      for (i = 0; i < com.ns; i++)
         com.z[i] = data.z[locus][i];
      com.npatt = com.posG[1] = data.npatt[locus];
      com.posG[0] = 0;
      com.fpatt = data.fpatt[locus];

      /* The following is model-dependent */
      if (data.datatype[locus] == BASE) {
         com.kappa = data.kappa[locus];
         com.alpha = data.alpha[locus];

         xtoy(data.pi[locus], com.pi, com.ncode);
         if (com.model <= TN93)
            eigenTN93(com.model, com.kappa, com.kappa, com.pi, &nR, Root, Cijk);

         if (com.alpha)
            DiscreteGamma(com.freqK, com.rK, com.alpha, com.alpha, com.ncatG, DGammaUseMedian);
      }
      /*
            com.NnodeScale = data.NnodeScale[locus];
            com.nodeScale = data.nodeScale[locus];
            nS = com.NnodeScale*com.npatt * (com.conPSiteClass?com.ncatG:1);
            for(i=0; i<nS; i++)
               com.nodeScaleF[i]=0;
      */
   }
   if (setSeqName)
      for (i = 0; i < com.ns; i++) com.spname[i] = stree.nodes[nodes[i].ipop].name;
   return(0);
}



int AcceptRejectLocus(int locus, int accept)
{
   /* This accepts or rejects the proposal at one locus.
      This works for proposals that change one locus only.
      After UseLocus(), gnodes[locus][ns].conP points to the alternative
      conP space.  If the proposal is accepted, this copies the newly calculated
      conP into gnodes[locus][ns].conP.  In either case, gnodes[].conP is
      repositioned.
      Proposals that change all loci use switchconP() to accept the proposal.
   */
   int i, ns = data.ns[locus], s1 = com.ncode*data.npatt[locus], sG = 1;
   double *conP = com.conPin[com.curconP] + data.conP_offset[locus];

   if (mcmc.usedata == 1) {
      if (com.conPSiteClass) sG = com.ncatG;
      if (accept)
         memcpy(conP, gnodes[locus][ns].conP, sG*s1*(ns - 1) * sizeof(double));
      for (i = ns; i < ns * 2 - 1; i++)
         gnodes[locus][i].conP = conP + (i - ns)*s1;
   }
   return(0);
}

void switchconPin(void)
{
   /* This reposition pointers gnodes[locus].conP to the alternative com.conPin, to avoid
      recalculation of conditional probabilities, when a proposal is accepted in steps that
      change all loci in one go, such as UpdateTimes() and UpdateParameters().
      Note that for site-class models (com.conPSiteClass), gnodes[].conP points to the space
      for class 0, and the space for class 1 starts (ns-1)*ncode*npatt later.  Such
      repositioning for site classes is achieved in fx_r().
   */
   int i, locus;
   double *conP;

   com.curconP = !com.curconP;

   for (locus = 0; locus < data.ngene; locus++) {
      conP = com.conPin[com.curconP] + data.conP_offset[locus];
      for (i = data.ns[locus]; i < data.ns[locus] * 2 - 1; i++, conP += com.ncode*data.npatt[locus])
         gnodes[locus][i].conP = conP;
   }
}

int SetPGene(int igene, int _pi, int _UVRoot, int _alpha, double xcom[])
{
   /* This is not used. */
   if (com.ngene != 1) error2("com.ngene==1?");
   return (0);
}


int SetParameters(double x[])
{
   /* This is not used. */
   return(0);
}

int GetPMatBranch2(double PMat[], double t)
{
   /* This calculates the transition probability matrix.
   */
   double Qrates[2], T, C, A, G, Y, R, mr;

   NPMat++;
   Qrates[0] = Qrates[1] = com.kappa;
   if (com.seqtype == 0) {
      if (com.model <= K80)
         PMatK80(PMat, t, com.kappa);
      else if (com.model <= TN93) {
         T = com.pi[0]; C = com.pi[1]; A = com.pi[2]; G = com.pi[3]; Y = T + C; R = A + G;
         if (com.model == F84) {
            Qrates[0] = 1 + com.kappa / Y;   /* kappa1 */
            Qrates[1] = 1 + com.kappa / R;   /* kappa2 */
         }
         else if (com.model <= HKY85) Qrates[1] = Qrates[0];
         mr = 1 / (2 * T*C*Qrates[0] + 2 * A*G*Qrates[1] + 2 * Y*R);
         PMatTN93(PMat, t*mr*Qrates[0], t*mr*Qrates[1], t*mr, com.pi);
      }
   }
   return(0);
}


int ConditionalPNode(int inode, int igene, double x[])
{
   int n = com.ncode, i, j, k, h, ison;
   double t;

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      if (nodes[ison].nson > 0 && !com.oldconP[ison])
         ConditionalPNode(ison, igene, x);
   }
   for (i = 0; i < com.npatt*n; i++) nodes[inode].conP[i] = 1;

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];

      t = nodes[ison].branch*_rateSite;
      if (t < 0) {
         printf("\nt =%12.6f ratesite = %.6f", nodes[ison].branch, _rateSite);
         error2("blength < 0");
      }

      GetPMatBranch2(PMat, t);

      if (nodes[ison].nson < 1 && com.cleandata) {        /* tip && clean */
         for (h = 0; h < com.npatt; h++)
            for (j = 0; j < n; j++)
               nodes[inode].conP[h*n + j] *= PMat[j*n + com.z[ison][h]];
      }
      else if (nodes[ison].nson < 1 && !com.cleandata) {  /* tip & unclean */
         for (h = 0; h < com.npatt; h++)
            for (j = 0; j < n; j++) {
               for (k = 0, t = 0; k < nChara[com.z[ison][h]]; k++)
                  t += PMat[j*n + CharaMap[com.z[ison][h]][k]];
               nodes[inode].conP[h*n + j] *= t;
            }
      }
      else {
         for (h = 0; h < com.npatt; h++)
            for (j = 0; j < n; j++) {
               for (k = 0, t = 0; k < n; k++)
                  t += PMat[j*n + k] * nodes[ison].conP[h*n + k];
               nodes[inode].conP[h*n + j] *= t;
            }
      }

   }  /*  for (ison)  */

   /* node scaling.  Is this coded?  Please check.  */
   if (com.NnodeScale && com.nodeScale[inode])
      NodeScale(inode, 0, com.npatt);
   return (0);
}

double lnLmorphF73(int inode, int locus, double *vp)
{
   /* This calculates the lnL for a morphological locus.
      data.zmorph[locus][] has the morphological measurements.
      data.Rmorph[locus][i*com.ls+j] is correlation between characters i and j.
      data.rgene[locus] is the locus rate or mu_locus for sequence data.
   */
   int s = com.ns, j, h, *sons = nodes[inode].sons;
   double v[2] = { -1,-1 }, x[2], zz, vv, y, t, nu, lnL = 0;
   int debug = 0;

   for (j = 0; j < 2; j++) {
      if (nodes[sons[j]].nson)
         lnL += lnLmorphF73(sons[j], locus, &v[j]);
      else {
         t = nodes[inode].age - nodes[sons[j]].age;
         nu = (com.clock == 1 ? data.rgene[locus] : stree.nodes[sons[j]].rates[locus]);
         v[j] = t*nu;
         /* Correct branch length for population variance. MdR. */
         if (nodes[sons[j]].nson == 0)  v[j] += data.zpopvar[locus];
      }
   }
   vv = v[0] + v[1];
   *vp = nodes[inode].branch + v[0] * v[1] / vv;
   for (h = 0, zz = 0; h < com.ls; h++) {
      x[0] = data.zmorph[locus][sons[0]][h];
      x[1] = data.zmorph[locus][sons[1]][h];
      data.zmorph[locus][inode][h] = (v[0] * x[1] + v[1] * x[0]) / vv;
      zz += square(x[0] - x[1]);
   }
   lnL += y = -0.5*com.ls*log(2 * Pi*vv) - zz / (2 * vv) - data.ldetRm[locus] / 2; /* MdR */

   if (debug) {
      printf("\nnode %2d sons %2d %2d v %6.4f %7.4f vp %7.4f, x %7.4f lnL = %7.4f %8.4f",
         inode + 1, sons[0] + 1, sons[1] + 1, v[0], v[1], *vp, data.zmorph[locus][inode][0], y, lnL);

      if (inode == data.root[locus]) exit(0);
   }
   return(lnL);
}

double lnpD_locus(int locus)
{
   /* This calculates ln{Di|Gi, Bi} using times in stree.nodes[].age and rates.
      Rates are in data.rgene[] if (clock==1) and in stree.nodes[].rates if
      (clock==0).  Branch lengths in the gene tree, nodes[].branch, are calculated
      in this routine but they may not be up-to-date before calling this routine.
      UseLocus() is called before this routine.
   */
   int  i, j, dad;
   double lnL = 0, b, t;

   if (mcmc.usedata == 0)  return(0);
   for (i = 0; i < tree.nnode; i++)   /* age in gene tree */
      nodes[i].age = stree.nodes[nodes[i].ipop].age;
   if (com.clock == 1) {  /* global clock */
      for (i = 0; i < tree.nnode; i++) {
         if (i == tree.root) continue;
         nodes[i].branch = (nodes[nodes[i].father].age - nodes[i].age) * data.rgene[locus];
         if (nodes[i].branch < 0)
            error2("blength < 0");
      }
   }
   else {              /* independent & correlated rates */
      for (i = 0; i < tree.nnode; i++) {
         if (i == tree.root) continue;
         for (j = nodes[i].ipop, b = 0; j != nodes[nodes[i].father].ipop; j = dad) {
            dad = stree.nodes[j].father;
            t = stree.nodes[dad].age - stree.nodes[j].age;
            b += t * stree.nodes[j].rates[locus];
         }
         nodes[i].branch = b;
      }
   }
   if (mcmc.usedata == 1 && data.datatype[locus] == MORPHC)
      lnL = lnLmorphF73(data.root[locus], locus, &t);
   else if (mcmc.usedata == 1)
      lnL = -com.plfun(NULL, -1);
   else if (mcmc.usedata == 2)
      lnL = lnpD_locus_Approx(locus);

   return (lnL*BFbeta);
}

double lnpData(double lnpDi[])
{
   /* This calculates the log likelihood, the log of the probability of the data
      given gtree[] for each locus.
      This updates gnodes[locus][].conP for every node.
   */
   int j, locus;
   double lnL = 0, y;

   if (mcmc.saveconP)
      for (j = 0; j < stree.nspecies * 2 - 1; j++)
         com.oldconP[j] = 0;
   for (locus = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 0, mcmc.usedata, 1);
      y = lnpD_locus(locus);

      if (testlnL && fabs(lnpDi[locus] - y) > 1e-5)
         printf("\tlnLi %.6f != %.6f at locus %d\n", lnpDi[locus], y, locus + 1);
      lnpDi[locus] = y;
      lnL += y;
   }
   return(lnL);
}

double lnLTaylor(int nb, double z[], double g[], double H[])
{
   int i, j;
   double lnL = 0;

   for (i = 0; i < nb; i++)
      lnL += z[i] * g[i];
   for (i = 0; i < nb; i++) {
      lnL += z[i] * z[i] * H[i*nb + i] / 2;
      for (j = 0; j < i; j++)
         lnL += z[i] * z[j] * H[i*nb + j];
   }
   return(lnL);
}

double lnpD_locus_Approx(int locus)
{
   /* This calculates the likelihood using the normal approxiamtion (Thorne et al.
      1998).  The branch lengths on the unrooted tree estimated without clock have
      been read into nodes[][].label and the gradient and Hessian are in data.Gradient[]
      & data.Hessian[].  The branch lengths predicted from the rate-evolution
      model (that is, products of rates and times) are in nodes[].branch.
      The tree is rooted, and the two branch lengths around the root are summed up
      and treated as one branch length, compared with the MLE, which is stored as
      the branch length to the first son of the root.
      This is different from funSS_AHRS(), which has the data branch lengths in
      nodes[].branch and calculate the predicted likelihood by multiplying nodes[].age
      and nodes[].label (which holds the rates).  Think about a redesign to avoid confusion.
   */
   int debug = 0, j, ib, nb = tree.nnode - 1 - 1;  /* don't use tree.nbranch, which is not up to date. */
   int son1 = nodes[tree.root].sons[0], son2 = nodes[tree.root].sons[1];
   double lnL = 0, z[NS * 2 - 1], *g = data.Gradient[locus], *H = data.Hessian[locus], cJC = (com.ncode - 1.0) / com.ncode;
   double bTlog = 1e-5, elog = 0.1, e;   /* change the same line in ReadBlengthGH() as well  */

   /* construct branch lengths */
   for (j = 0; j < tree.nnode; j++) {
      if (j == tree.root || j == son2) continue;
      ib = nodes[j].ibranch;
      if (j == son1)
         z[ib] = nodes[son1].branch + nodes[son2].branch;
      else
         z[ib] = nodes[j].branch;

      if (data.transform == LOG_B) {         /* logarithm transform, MLE not transformed */
         e = (data.blMLE[locus][ib] < bTlog ? elog : 0);
         z[ib] = log((z[ib] + e) / (data.blMLE[locus][ib] + e));
      }
      else {
         if (data.transform == SQRT_B)           /* sqrt transform, MLE transformed */
            z[ib] = sqrt(z[ib]);
         else if (data.transform == ARCSIN_B)    /* arcsin transform, MLE transformed */
            z[ib] = 2 * asin(sqrt(cJC - cJC*exp(-z[ib] / cJC)));
         z[ib] -= data.blMLE[locus][ib];
      }
   }

   if (debug) {
      OutTreeN(F0, 1, 1);  FPN(F0);
      matout(F0, z, 1, nb);
      matout(F0, data.blMLE[locus], 1, nb);
   }

   lnL = lnLTaylor(nb, z, g, H);

   return (lnL);
}


int ReadBlengthGH(char infile[])
{
   /* this reads the MLEs of branch lengths under no clock and their SEs, for
      approximate calculation of sequence likelihood.  The MLEs of branch lengths
      are stored in data.blMLE[locus] as well as nodes[].label, and the gradient
      and Hessian in data.Gradient[locsu][] and data.Hessian[locus][].
      The branch length (sum of 2 branch lengths) around root in rooted tree is
      placed on 1st son of root in rooted tree.
      For the log transform, data.blMLE[] holds the MLEs of branch lengths, while for
      other transforms (sqrt, jc), the transformed values are stored.

      This also frees up memory for sequences.
   */
   FILE* fBGH = gfopen(infile, "r");
   char line[100];
   int locus, i, j, nb, son1, son2, leftsingle;
   double small = 1e-20;
   double dbu[NBRANCH], dbu2[NBRANCH], u, cJC = (com.ncode - 1.0) / com.ncode, sin2u, cos2u;
   double bTlog = 1e-5, elog = 0.1, e;   /* change the line in lnpD_locus_Approx() as well */
   int debug = 0, longbranches = 0;

   if (noisy) printf("\nReading branch lengths, Gradient & Hessian from %s.\n", infile);

   for (locus = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 0, 0, 1);
      printf("Locus %d: %d species\r", locus + 1, com.ns);

      nodes = nodes_t;
      fscanf(fBGH, "%d", &i);
      if (i != com.ns) {
         printf("ns = %d, expecting %d", i, com.ns);
         error2("ns not correct in ReadBlengthGH()");
      }
      if (ReadTreeN(fBGH, &i, 0, 1))
         error2("error when reading gene tree");
      if (i == 0) error2("gene tree should have branch lengths for error checking");
      if ((nb = tree.nbranch) != com.ns * 2 - 3)
         error2("nb = ns * 2 -3 ?");

      if (debug) {
         FPN(F0); OutTreeN(F0, 1, 1);  FPN(F0);
         OutTreeB(F0);
         for (i = 0; i < tree.nnode; i++) {
            printf("\nnode %2d: branch %2d (%9.5f) dad %2d  ", i + 1, nodes[i].ibranch + 1, nodes[i].branch, nodes[i].father + 1);
            for (j = 0; j < nodes[i].nson; j++) printf(" %2d", nodes[i].sons[j] + 1);
         }
      }

      for (i = 0; i < nb; i++) {
         fscanf(fBGH, "%lf", &data.blMLE[locus][i]);
         if (data.blMLE[locus][i] > 10) {
            printf("very large branch length %10.6f\n", data.blMLE[locus][i]);
            longbranches++;
         }
      }
      for (i = 0; i < nb; i++)
         fscanf(fBGH, "%lf", &data.Gradient[locus][i]);

      fscanf(fBGH, "%s", line);
      if (!strstr(line, "Hessian")) error2("expecting Hessian in in.BV");
      for (i = 0; i < nb; i++)
         for (j = 0; j < nb; j++)
            if (fscanf(fBGH, "%lf", &data.Hessian[locus][i*nb + j]) != 1)
               error2("err when reading the hessian matrix in.BV");

      UseLocus(locus, 0, 0, 1);
      NodeToBranch();
      son1 = nodes[tree.root].sons[0];
      son2 = nodes[tree.root].sons[1];
      leftsingle = (nodes[son1].nson == 0);
      if (leftsingle) {  /* Root in unrooted tree is 2nd son of root in rooted tree */
         nodes[son1].ibranch = com.ns * 2 - 3 - 1;  /* last branch in unrooted tree, assuming binary tree */
         for (i = 0; i < tree.nnode; i++) {
            if (i != tree.root && i != son1 && i != son2)
               nodes[i].ibranch -= 2;
         }
      }
      else {  /* Root in unrooted tree is 1st son of root in rooted tree */
         nodes[son1].ibranch = nodes[son2].ibranch - 1;
         for (i = 0; i < tree.nnode; i++) {
            if (i != tree.root && i != son1 && i != son2)
               nodes[i].ibranch--;
         }
      }
      nodes[tree.root].ibranch = nodes[son2].ibranch = -1;
      for (i = 0; i < tree.nnode; i++)
         nodes[i].branch = (nodes[i].ibranch == -1 ? 0 : data.blMLE[locus][nodes[i].ibranch]);

      if (debug) {
         FPN(F0);  FPN(F0);  OutTreeN(F0, 1, 1);
         for (i = 0; i < tree.nnode; i++) {
            printf("\nnode %2d: branch %2d (%9.5f) dad %2d  ", i + 1, nodes[i].ibranch + 1, nodes[i].branch, nodes[i].father + 1);
            for (j = 0; j < nodes[i].nson; j++) printf(" %2d", nodes[i].sons[j] + 1);
         }
         FPN(F0);
      }

      if (data.transform == SQRT_B)       /* sqrt transform on branch lengths */
         for (i = 0; i < nb; i++) {
            dbu[i] = 2 * sqrt(data.blMLE[locus][i]);
            dbu2[i] = 2;
         }
      else if (data.transform == LOG_B)   /* logarithm transform on branch lengths */
         for (i = 0; i < nb; i++) {
            e = (data.blMLE[locus][i] < bTlog ? elog : 0);
            dbu[i] = data.blMLE[locus][i] + e;
            dbu2[i] = dbu[i];
         }
      else if (data.transform == ARCSIN_B) {  /* arcsin transform on branch lengths */
         for (i = 0; i < nb; i++) {
            /*
            if(data.blMLE[locus][i] < 1e-20)
               error2("blength should be > 0 for arcsin transform");
            */
            u = 2 * asin(sqrt(cJC - cJC*exp(-data.blMLE[locus][i] / cJC)));
            sin2u = sin(u / 2);  cos2u = cos(u / 2);
            dbu[i] = sin2u*cos2u / (1 - sin2u*sin2u / cJC);
            dbu2[i] = (cos2u*cos2u - sin2u*sin2u) / 2 / (1 - sin2u*sin2u / cJC) + dbu[i] * dbu[i] / cJC;
         }
      }

      if (data.transform) {          /* this part is common to all transforms */
         for (i = 0; i < nb; i++) {
            for (j = 0; j < i; j++)
               data.Hessian[locus][j*nb + i] = data.Hessian[locus][i*nb + j] *= dbu[i] * dbu[j];
            data.Hessian[locus][i*nb + i] = data.Hessian[locus][i*nb + i] * dbu[i] * dbu[i]
               + data.Gradient[locus][i] * dbu2[i];
         }
         for (i = 0; i < nb; i++)
            data.Gradient[locus][i] *= dbu[i];
      }

      if (data.transform == SQRT_B)       /* sqrt transform on branch lengths */
         for (i = 0; i < nb; i++)
            data.blMLE[locus][i] = sqrt(data.blMLE[locus][i]);
      else if (data.transform == LOG_B) { /* logarithm transform on branch lengths */
         ;
      }
      else if (data.transform == ARCSIN_B)    /* arcsin transform on branch lengths */
         for (i = 0; i < nb; i++)
            data.blMLE[locus][i] = 2 * asin(sqrt(cJC - cJC*exp(-data.blMLE[locus][i] / cJC)));
   }    /* for(locus) */

   fclose(fBGH);
   /* free up memory for sequences */
   for (locus = 0; locus < data.ngene; locus++) {
      free(data.fpatt[locus]);
      for (j = 0; j < data.ns[locus]; j++)
         free(data.z[locus][j]);
   }

   if (longbranches) printf("\a\n%d branch lengths are >10.  Check that you are not using random sequences.", longbranches);
   return(0);
}


int GenerateBlengthGH(char infile[])
{
   /* This generates the sequence alignment and tree files for calculating branch
      lengths, gradient, and hessian.
      This mimics Jeff Thorne's estbranches program.
   */
   FILE *fseq, *ftree, *fctl;
   FILE *fBGH = gfopen(infile, "w");
   char tmpseqf[32], tmptreef[32], ctlf[32], outf[32], line[10000];
   int i, locus;

   mcmc.usedata = 1;
   for (locus = 0; locus < data.ngene; locus++) {
      printf("\n\n*** Locus %d ***\n", locus + 1);
      sprintf(tmpseqf, "tmp%04d.txt", locus + 1);
      sprintf(tmptreef, "tmp%04d.trees", locus + 1);
      sprintf(ctlf, "tmp%04d.ctl", locus + 1);
      sprintf(outf, "tmp%04d.out", locus + 1);
      fseq = gfopen(tmpseqf, "w");
      ftree = gfopen(tmptreef, "w");
      fctl = gfopen(ctlf, "w");

      UseLocus(locus, 0, 0, 1);
      com.cleandata = data.cleandata[locus];
      for (i = 0; i < com.ns; i++)
         com.z[i] = data.z[locus][i];
      com.npatt = com.posG[1] = data.npatt[locus];
      com.posG[0] = 0;
      com.fpatt = data.fpatt[locus];

      DeRoot();

      printPatterns(fseq);
      fprintf(ftree, "  1\n\n");
      OutTreeN(ftree, 1, 0);   FPN(ftree);

      fprintf(fctl, "seqfile = %s\n", tmpseqf);
      fprintf(fctl, "treefile = %s\n", tmptreef);
      fprintf(fctl, "outfile = %s\nnoisy = 3\n", outf);
      if (com.seqtype) {
         fprintf(fctl, "seqtype = %d\n", com.seqtype);
      }
      fprintf(fctl, "model = %d\n", com.model);
      if (com.seqtype == 1) {
         fprintf(fctl, "icode = %d\n", com.icode);
         fprintf(fctl, "fix_kappa = 0\n kappa = 2\n");
      }
      if (com.seqtype == 2) {
         fprintf(fctl, "aaRatefile = %s\n", com.daafile);
      }
      if (com.alpha)
         fprintf(fctl, "fix_alpha = 0\nalpha = 0.5\nncatG = %d\n", com.ncatG);
      fprintf(fctl, "Small_Diff = 0.1e-6\ngetSE = 2\n");
      fprintf(fctl, "method = %d\n", (com.alpha == 0 || com.ns > 100 ? 1 : 0));

      fclose(fseq);  fclose(ftree);  fclose(fctl);

      if (com.seqtype) sprintf(line, "codeml %s", ctlf);
      else            sprintf(line, "baseml %s", ctlf);
      printf("running %s\n", line);
      system(line);

      appendfile(fBGH, "rst2");
   }
   fclose(fBGH);
   return(0);
}


int GetOptions(char *ctlf)
{
   int  transform0 = ARCSIN_B; /* default transform: SQRT_B, LOG_B, ARCSIN_B */
   int  iopt, i, j, nopt = 32, lline = 4096;
   char line[4096], *pline, opt[33], *comment = "*#";
   char *optstr[] = { "seed", "seqfile","treefile", "outfile", "mcmcfile", "checkpoint", "BayesFactorBeta",
        "seqtype", "aaRatefile", "icode", "noisy", "usedata", "ndata", "duplication", "model", "clock",
        "TipDate", "RootAge", "fossilerror", "alpha", "ncatG", "cleandata",
        "BDparas", "kappa_gamma", "alpha_gamma",
        "rgene_gamma", "sigma2_gamma", "print", "burnin", "sampfreq",
        "nsample", "finetune" };
   double t = 1, *eps = mcmc.steplength;
   FILE  *fctl = gfopen(ctlf, "r");

   data.rgeneprior = 0;  /* default rate prior is gamma-Dirichlet. */
   data.transform = transform0;
   strcpy(com.checkpointf, "mcmctree.ckpt");
   if (fctl) {
      if (noisy) printf("\nReading options from %s..\n", ctlf);
      for (;;) {
         if (fgets(line, lline, fctl) == NULL) break;
         for (i = 0, t = 0, pline = line; i < lline&&line[i]; i++)
            if (isalnum(line[i])) { t = 1; break; }
            else if (strchr(comment, line[i])) break;
            if (t == 0) continue;
            if ((pline = strstr(line, "=")) == NULL)
               error2("option file.");
            *pline = '\0';   sscanf(line, "%s", opt);  *pline = '=';
            sscanf(pline + 1, "%lf", &t);

            for (iopt = 0; iopt < nopt; iopt++) {
               if (strncmp(opt, optstr[iopt], 8) == 0) {
                  if (noisy >= 9)
                     printf("\n%3d %15s | %-20s %6.2f", iopt + 1, optstr[iopt], opt, t);
                  switch (iopt) {
                  case (0): SetSeed((int)t, 1);                 break;
                  case (1): sscanf(pline + 1, "%s", com.seqf);    break;
                  case (2): sscanf(pline + 1, "%s", com.treef);   break;
                  case (3): sscanf(pline + 1, "%s", com.outf);    break;
                  case (4): sscanf(pline + 1, "%s", com.mcmcf);   break;
                  case (5): sscanf(pline + 1, "%d", &com.checkpoint); break;
                  case (6): sscanf(pline + 1, "%lf", &BFbeta);    break; /* beta for marginal likelihood */
                  case (7): com.seqtype = (int)t;    break;
                  case (8): sscanf(pline + 2, "%s", com.daafile);  break;
                  case (9): com.icode = (int)t;      break;
                  case (10): noisy = (int)t;          break;
                  case (11):
                     j = sscanf(pline + 1, "%d %s%d", &mcmc.usedata, com.inBVf, &data.transform);
                     if (mcmc.usedata == 2)
                        if (strchr(com.inBVf, '*')) { strcpy(com.inBVf, "in.BV"); data.transform = transform0; }
                        else if (j == 2)              data.transform = transform0;
                        break;
                  case (12): com.ndata = (int)t;           break;
                  case (13): stree.duplication = (int)t;   break;
                  case (14): com.model = (int)t;           break;
                  case (15): com.clock = (int)t;           break;
                  case (16):
                     sscanf(pline + 2, "%lf%lf", &com.TipDate, &com.TipDate_TimeUnit);
                     if (com.TipDate && com.TipDate_TimeUnit == 0) error2("should set com.TipDate_TimeUnit");
                     data.transform = SQRT_B;  /* SQRT_B, LOG_B, ARCSIN_B */
                     break;
                  case (17):
                     stree.RootAge[2] = stree.RootAge[3] = 0.025;  /* default tail probs */
                     if ((strchr(line, '>') || strchr(line, '<')) && (strstr(line, "U(") || strstr(line, "B(")))
                        error2("don't mix < U B on the RootAge line");

                     if ((pline = strchr(line, '>')))
                        sscanf(pline + 1, "%lf", &stree.RootAge[0]);
                     if ((pline = strchr(line, '<'))) {  /* RootAge[0]=0 */
                        sscanf(pline + 1, "%lf", &stree.RootAge[1]);
                     }
                     if ((pline = strstr(line, "U(")))
                        sscanf(pline + 2, "%lf,%lf", &stree.RootAge[1], &stree.RootAge[2]);
                     else if ((pline = strstr(line, "B(")))
                        sscanf(pline + 2, "%lf,%lf,%lf,%lf", &stree.RootAge[0], &stree.RootAge[1], &stree.RootAge[2], &stree.RootAge[3]);
                     break;
                  case (18):
                     data.pfossilerror[0] = 0.0;
                     data.pfossilerror[2] = 1;  /* default: minimum 2 good fossils */
                     sscanf(pline + 1, "%lf%lf%lf", data.pfossilerror, data.pfossilerror + 1, data.pfossilerror + 2);
                     break;
                  case (19): com.alpha = t;           break;
                  case (20): com.ncatG = (int)t;      break;
                  case (21): com.cleandata = (int)t;  break;
                  case (22):
                     sscanf(pline + 1, "%lf%lf%lf%lf", &data.BDS[0], &data.BDS[1], &data.BDS[2], &data.BDS[3]);
                     break;
                  case (23):
                     sscanf(pline + 1, "%lf%lf", data.kappagamma, data.kappagamma + 1); break;
                  case (24):
                     sscanf(pline + 1, "%lf%lf", data.alphagamma, data.alphagamma + 1); break;
                  case (25):
                     sscanf(pline + 1, "%lf%lf%lf%d", &data.rgenepara[0], &data.rgenepara[1], &data.rgenepara[2], &data.rgeneprior);
                     if (data.rgenepara[2] <= 0)  data.rgenepara[2] = 1;
                     if (data.rgeneprior < 0)     data.rgeneprior = 0;
                     break;
                  case (26):
                     sscanf(pline + 1, "%lf%lf%lf", data.sigma2para, data.sigma2para + 1, data.sigma2para + 2);
                     if (data.sigma2para[2] <= 0) data.sigma2para[2] = 1;
                     break;
                  case (27): mcmc.print = (int)t;     break;
                  case (28): mcmc.burnin = (int)t;    break;
                  case (29): mcmc.sampfreq = (int)t;  break;
                  case (30): mcmc.nsample = (int)t;   break;
                  case (31):
                     puts("finetune is deprecated now.");
                     break;
                     sscanf(pline + 1, "%d:%lf%lf%lf%lf%lf%lf", &j, eps, eps + 1, eps + 2, eps + 3, eps + 4, eps + 5);
                     break;
                  }
                  break;
               }
            }
            if (iopt == nopt) {
               printf("\noption %s in %s\n", opt, ctlf);  exit(-1);
            }
      }
      fclose(fctl);
   }
   else
      if (noisy) error2("\nno ctl file...");

   if (com.ndata > NGENE) error2("raise NGENE?");
   else if (com.ndata <= 0) com.ndata = 1;
   if (mcmc.usedata == 2) {
      com.model = 0;  com.fix_kappa = com.fix_alpha = 1; com.alpha = 0;
      if (com.seqtype == 1) com.ncode = 61;
   }
   if (com.seqtype == 2)
      com.ncode = 20;

   if (com.alpha == 0) { com.fix_alpha = 1; com.nalpha = 0; }
   if (com.clock < 1 || com.clock>3) error2("clock should be 1, 2, 3?");
   if (mcmc.burnin <= 0) puts("burnin=0: no automatic step adjustment?");

   if (BFbeta && mcmc.usedata == 0)
      error2("marginal like for prior with usedata =0?");
   else if (BFbeta == 0)
      BFbeta = 1;

   return(0);
}


double lnPDFInfinitesitesClock(double t1, double FixedDs[])
{
   /* This calculates the ln of the joint pdf, which is proportional to the
      posterior for given root age, assuming infinite sites.  Fixed distances are
      in FixedDs[]: d11, d12, ..., d(1,s-1), d21, d31, ..., d_g1.
      Note that the posterior is one dimensional, and the variable used is root age.
   */
   int s = stree.nspecies, g = data.ngene, i, j;
   double lnp, summu = 0, prodmu = 1, *gD = data.rgenepara;

   stree.nodes[stree.root].age = t1;
   for (j = s; j < stree.nnode; j++)
      if (j != stree.root)
         stree.nodes[j].age = t1*FixedDs[j - s] / FixedDs[0];

   lnp = lnpriorTimes();

   data.rgene[0] = summu = prodmu = FixedDs[0] / t1;
   for (i = 1; i < g; i++) {
      data.rgene[i] = FixedDs[s - 1 + i - 1] / t1;
      summu += data.rgene[j];
      prodmu *= data.rgene[j];
   }
   lnp += (gD[0] - gD[2] * g)*log(summu) - gD[1] / g*summu + (gD[2] - 1)*log(prodmu); /* f(mu_i) */

   lnp += (2 - s)*log(FixedDs[0] / t1) - g*log(t1);                               /* Jacobi */

   return (lnp);
}

double InfinitesitesClock(double *FixedDs)
{
   /* This runs MCMC to calculate the posterior density of the root age
      when there are infinite sites at each locus.  The clock is assumed, so that
      the posterior is one-dimensional.
   */
   int i, j, ir, nround = 0, nsaved = 0;
   int s = stree.nspecies, g = data.ngene;
   double t, tnew, naccept = 0;
   double e = 0.01 + 0.1*rndu(), lnp, lnpnew, lnacceptance, c, *x, Pjump = 0;
   double tmean, t025, t975, tL, tU;
   char timestr[32];

   matout2(F0, FixedDs, 1, s - 1 + g - 1, 9, 5);
   printf("\nRunning MCMC to get %d samples for t0 (root age)\n", mcmc.nsample);
   t = stree.nodes[stree.root].age;
   lnp = lnPDFInfinitesitesClock(t, FixedDs);
   x = (double*)malloc(max2(mcmc.nsample, (s + g) * 3) * sizeof(double));
   if (x == NULL) error2("oom x");

   for (ir = -mcmc.burnin, tmean = 0; ir < mcmc.nsample*mcmc.sampfreq; ir++) {
      if (ir == 0 || (nround >= 100 && ir < 0 && ir % (mcmc.burnin / 4) == 0)) {
         nround = 0; naccept = 0; tmean = 0;
         ResetStepLengths(NULL, &Pjump, &e, 1);
      }
      lnacceptance = e*rndSymmetrical();
      c = exp(lnacceptance);
      tnew = t*c;
      lnpnew = lnPDFInfinitesitesClock(tnew, FixedDs);
      lnacceptance += lnpnew - lnp;

      if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) {
         t = tnew; lnp = lnpnew;  naccept++;
      }
      nround++;
      tmean = (tmean*(nround - 1) + t) / nround;

      if (ir >= 0 && (ir + 1) % mcmc.sampfreq == 0)
         x[nsaved++] = t;

      Pjump = naccept / nround;
      if ((ir + 1) % max2(mcmc.sampfreq, mcmc.sampfreq*mcmc.nsample / 10000) == 0)
         printf("\r%3.0f%%  %7.2f mean t0 = %9.6f", (ir + 1.) / (mcmc.nsample*mcmc.sampfreq) * 100, Pjump, tmean);
      if (mcmc.sampfreq*mcmc.nsample > 20 && (ir + 1) % (mcmc.sampfreq*mcmc.nsample / 10) == 0)
         printf(" %s\n", printtime(timestr));
   }

   qsort(x, (size_t)mcmc.nsample, sizeof(double), comparedouble);
   t025 = x[(int)(mcmc.nsample*.025 + .5)];
   t975 = x[(int)(mcmc.nsample*.975 + .5)];

   /* Below x[] is used to collect the posterior means and 95% CIs */
   for (i = 0; i < 3; i++) {
      t = (i == 0 ? tmean : (i == 1 ? t025 : t975));
      lnPDFInfinitesitesClock(t, FixedDs);  /* calculates t and mu */
      for (j = s; j < stree.nnode; j++)
         x[i*(s + g) + (j - s)] = stree.nodes[j].age;
      for (j = 0; j < g; j++)
         x[i*(s + g) + s - 1 + j] = data.rgene[j];
   }
   printf("\nmean (95%% CI) CI-width for times\n\n");
   for (j = s; j < stree.nnode; j++) {
      tL = x[(s + g) + j - s];
      tU = x[2 * (s + g) + j - s];
      printf("Node %2d: %9.6f (%9.6f, %9.6f) %9.6f\n", j + 1, x[j - s], tL, tU, tU - tL);
   }
   printf("\nmean & 95%% CI for rates\n\n");
   for (j = 0; j < g; j++)
      printf("gene %2d: %9.6f (%9.6f, %9.6f)\n", j + 1, x[s - 1 + j], x[2 * (s + g) + s - 1 + j], x[(s + g) + s - 1 + j]);

   printf("\nNote: the posterior has only one dimension.\n");
   free(x);
   exit(0);
}


double lnPDFInfinitesitesClock23(double FixedDs[])
{
   /* This calculates the log joint pdf, which is proportional to the posterior
      when there are infinite sites at each locus.  The data are fixed branch lengths,
      stored in FixedDs, locus by locus.
      The variables in the posterior are node ages, rate r0 for root son0, and mu & sigma2.
      This PDF includes f(t1,...,t_{s-1})*f(r_ij), but not f(mu_i) or f(sigma2_i) for loci.
      The latter are used to calculate ratios in the MCMC algorithm.
   */
   int s = stree.nspecies, g = data.ngene, locus, j;
   int root = stree.root, *sons = stree.nodes[root].sons;
   double lnJ, lnp, b, r0;  /* r0 is rate for root son0, fixed. */
   double t, t0 = stree.nodes[root].age - stree.nodes[sons[0]].age;
   double t1 = stree.nodes[root].age - stree.nodes[sons[1]].age;
   double summu = 0, prodmu = 1, *gD = data.rgenepara, *gDs = data.sigma2para;

   /* compute branch rates using times and branch lengths */
   for (locus = 0; locus < g; locus++) {
      for (j = 0; j < stree.nnode; j++) {
         if (j == root || j == sons[0]) continue;  /* r0 is in the posterior */
         t = stree.nodes[nodes[j].father].age - stree.nodes[j].age;

         if (t <= 0)
            error2("t<=0");
         if (j == sons[1]) {
            b = FixedDs[locus*stree.nnode + sons[0]];
            r0 = stree.nodes[sons[0]].rates[locus];
            stree.nodes[j].rates[locus] = (b - r0*t0) / t1;
            if (r0 <= 0 || t1 <= 0 || b - r0*t0 <= 0) {
               printf("\nr0 = %.6f t1 = %.6f b-t0*t0 = %.6f", r0, t1, b - r0*t0);
               error2("r<=0 || t1<=0 || b-r0*t0<=0...");
            }
         }
         else
            stree.nodes[j].rates[locus] = FixedDs[locus*stree.nnode + j] / t;
      }
   }
   if (debug == 9) {
      printf("\n   (age tlength)        rates & branchlengths\n");
      for (j = 0; j < stree.nnode; j++, FPN(F0)) {
         t = (j == root ? -1 : stree.nodes[nodes[j].father].age - stree.nodes[j].age);
         printf("%2d (%7.4f%9.4f)  ", j + 1, stree.nodes[j].age, t);
         for (locus = 0; locus < g; locus++) printf(" %9.4f", stree.nodes[j].rates[locus]);
         printf("  ");
         for (locus = 0; locus < g; locus++) printf(" %9.4f", FixedDs[locus*stree.nnode + j]);
      }
      sleep2(2);
   }

   lnp = lnpriorTimes() + lnpriorRates();
   for (j = 0, lnJ = -log(t1); j < stree.nnode; j++)
      if (j != root && j != sons[0] && j != sons[1])
         lnJ -= log(stree.nodes[nodes[j].father].age - stree.nodes[j].age);
   lnp += g*lnJ;

   return (lnp);
}


double Infinitesites(FILE *fout)
{
   /* This reads the fixed branch lengths and runs MCMC to calculate the posterior
      of times and rates when there are infinite sites at each locus.
      This is implemented for the global clock model (clock = 1) and the
      independent-rates model (clock = 2).  I think this works for clock3 as well.
      mcmc.steplength[] is used to propose changes: 0: t, 1:mu, 2:r; 3:mixing, 4:sigma2.

      For clock=2, the MCMC samples the following variables:
          s-1 times, rate r0 for left root branch at each locus,
          mu at each locus & sigma2 at each locus.

       Eevery species has to be present at every locus.

   */
   int root = stree.root, *sons = stree.nodes[root].sons;
   int lline = 10000, locus, nround, ir, i, j, k, ip, s = stree.nspecies, g = data.ngene;
   int nxpr[2] = { 8,3 }, MixingStep = 1, rooted;
   char timestr[36];
   double e[6] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
   double y, ynew, yb[2], sumold, sumnew, *gD, naccept[5] = { 0 }, Pjump[5] = { 0 };
   double *x, *mx, *FixedDs, *b, maxt0, tson[2], lnL, lnLnew, lnacceptance, c, lnc;
   char *FidedDf[2] = { "FixedDsClock1.txt", "FixedDsClock23.txt" };
   FILE *fin = gfopen(FidedDf[com.clock > 1], "r"), *fmcmc = gfopen(com.mcmcf, "w");

   if (BFbeta != 1) error2("BFbeta should not be used for Infinitesites?");
   com.model = 0;  com.alpha = 0;
   mcmc.usedata = 0;
   if (data.rgeneprior == 1) puts("\aInfiniteSites, not working for cond i.i.d. locus rate prior?");

   printStree();
   if (com.clock > 1) GetMem();
   GetInitials();  /* This sets root age.  This works for TipDate models? */

   printf("\nInfiniteSites, reading branch lengths from %s (s = %d  g = %d):\n\n", FidedDf[com.clock > 1], s, g);
   com.np = s - 1 + g + g + g;   /* times, mu, sigma2, rates */
   FixedDs = (double*)malloc((g*stree.nnode + com.np * 2) * sizeof(double));
   if (FixedDs == NULL) error2("oom");
   x = FixedDs + g*stree.nnode;
   mx = x + com.np;
   fscanf(fin, "%d", &i);
   if (i != s) error2("wrong number of species in FixedDs.txt");

   if (data.pfossilerror[0]) {
      puts("model of fossil errors for infinite data not tested yet.");
      getchar();
      SetupPriorTimesFossilErrors();
   }

   if (com.clock == 1) { /* global clock: read FixedDs[] and close file. */
      for (i = 0; i < s - 1 + g - 1; i++)
         fscanf(fin, "%lf", &FixedDs[i]);
      fclose(fin);
      InfinitesitesClock(FixedDs);
      return(0);
   }

   /* print header line in the mcmc.txt sample file */
   fprintf(fmcmc, "Gen");
   for (i = s; i < stree.nnode; i++) fprintf(fmcmc, "\tt_n%d", i + 1);
   for (j = 0; j < g; j++) fprintf(fmcmc, "\tmu_L%d", j + 1);
   for (j = 0; j < g; j++) fprintf(fmcmc, "\tsigma2_L%d", j + 1);
   for (j = 0; j < g; j++) fprintf(fmcmc, "\tr_left_L%d", j + 1);
   fprintf(fmcmc, "\n");

   for (locus = 0, b = FixedDs, nodes = nodes_t; locus < g; locus++, b += stree.nnode) {
      ReadTreeN(fin, &i, 0, 1);
      OutTreeN(F0, 1, 1); FPN(F0);
      rooted = 1;
      if (tree.nnode == 2 * s - 2 && nodes[tree.root].nson == 3)
         rooted = 0;
      else if (tree.nnode != 2 * s - 1)
         error2("use species tree for each locus!");
      b[root] = -1;
      if (rooted) b[sons[0]] = nodes[sons[0]].branch + nodes[sons[1]].branch;
      b[sons[1]] = -1;
      for (j = 0; j < stree.nnode; j++) {
         if (j != root && j != sons[0] && j != sons[1])
            b[j] = nodes[j].branch;
      }
   }
   fclose(fin);

   printf("\nFixed distances at %d loci\n", g);
   for (j = 0; j < stree.nnode; j++, FPN(F0)) {
      printf("node %3d  ", j + 1);
      for (locus = 0; locus < g; locus++)
         printf(" %9.6f", FixedDs[locus*stree.nnode + j]);
   }

   for (i = 0; i < g; i++) { /* Reset r0 so that r0*t0<b0.  */
      y = FixedDs[i*stree.nnode + sons[0]] / (stree.nodes[root].age - stree.nodes[sons[0]].age);
      stree.nodes[sons[0]].rates[i] = y*rndu();
   }

   for (i = 0; i < com.np; i++) mx[i] = 0;
   for (i = 0, k = 0; i < s - 1; i++) x[k++] = stree.nodes[s + i].age;
   for (i = 0; i < g; i++) x[k++] = data.rgene[i];
   for (i = 0; i < g; i++) x[k++] = data.sigma2[i];
   for (i = 0; i < g; i++) x[k++] = stree.nodes[sons[0]].rates[i];

   lnL = lnPDFInfinitesitesClock23(FixedDs);

   printf("\nStarting MCMC (np=%d) lnp = %9.3f\nInitials:            ", com.np, lnL);
   for (i = 0; i < com.np; i++) printf(" %5.3f", x[i]);
   printf("\n\nparas: %d times, %d mu, %d sigma2, %d rates r0 (left daughter of root)", s - 1, g, g, g);
   printf("\nUsing steplength parameters from the control file\n");

   /* MCMC proposals: t (0), mu (1), sigma2 (4), r0 (2), mixing (3) */
   for (ir = -mcmc.burnin, nround = 0; ir < mcmc.sampfreq*mcmc.nsample; ir++) {
      if (ir == 0 || (ir < 0 && nround >= 100 && ir % (mcmc.burnin / 4) == 0)) {
         nround = 0; zero(naccept, 5);  zero(mx, com.np);
         ResetStepLengths(NULL, Pjump, e, 5);
      }
      for (ip = 0; ip < com.np; ip++) {
         lnacceptance = 0;
         if (ip < s - 1) {  /* times */
            y = stree.nodes[s + ip].age;
            for (i = 0; i < 2; i++) tson[i] = stree.nodes[stree.nodes[s + ip].sons[i]].age;
            yb[0] = max2(tson[0], tson[1]);
            yb[1] = (s + ip == root ? OldAge : stree.nodes[stree.nodes[s + ip].father].age);
            if (s + ip == root) {
               for (i = 0; i < g; i++) {
                  maxt0 = FixedDs[i*stree.nnode + sons[0]] / stree.nodes[sons[0]].rates[i];
                  yb[1] = min2(yb[1], tson[0] + maxt0);
               }
            }
            else if (s + ip == sons[0]) {
               for (i = 0; i < g; i++) {
                  maxt0 = FixedDs[i*stree.nnode + sons[0]] / stree.nodes[sons[0]].rates[i];
                  yb[0] = max2(yb[0], stree.nodes[root].age - maxt0);
               }
            }
            ynew = y + e[0] * rndSymmetrical();
            ynew = stree.nodes[s + ip].age = reflect(ynew, yb[0], yb[1]);
         }
         else if (ip - (s - 1) < g) {    /* mu for loci */
            gD = data.rgenepara;
            for (j = 0, sumold = 0; j < g; j++) sumold += data.rgene[j];
            lnacceptance = lnc = e[1] * rndSymmetrical();
            c = exp(lnc);
            y = data.rgene[ip - (s - 1)];
            ynew = data.rgene[ip - (s - 1)] *= c;

            sumnew = sumold + ynew - y;
            lnacceptance += (gD[0] - gD[2] * g)*log(sumnew / sumold) - gD[1] / g*(ynew - y) + (gD[2] - 1)*lnc;
         }
         else if (ip - (s - 1 + g) < g) { /* sigma2 for loci */
            gD = data.sigma2para;
            for (j = 0, sumold = 0; j < g; j++) sumold += data.sigma2[j];
            lnacceptance = lnc = e[4] * rndSymmetrical();
            c = exp(lnc);
            y = data.sigma2[ip - (s - 1 + g)];
            ynew = data.sigma2[ip - (s - 1 + g)] *= c;
            sumnew = sumold + ynew - y;
            lnacceptance += (gD[0] - gD[2] * g)*log(sumnew / sumold) - gD[1] / g*(ynew - y) + (gD[2] - 1)*lnc;
         }
         else {                   /* rate r0 for root son0 for loci (r0*t0<b0) */
            y = stree.nodes[root].age - stree.nodes[sons[0]].age;
            if (y <= 0)
               printf("age error");
            yb[0] = 0;
            yb[1] = FixedDs[(ip - (s - 1 + g * 2))*stree.nnode + sons[0]] / y;
            y = stree.nodes[sons[0]].rates[ip - (s - 1 + g * 2)];
            ynew = y + e[2] * rndSymmetrical();
            ynew = stree.nodes[sons[0]].rates[ip - (s - 1 + g * 2)] = reflect(ynew, yb[0], yb[1]);
         }


         lnLnew = lnPDFInfinitesitesClock23(FixedDs);
         lnacceptance += lnLnew - lnL;

         if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) {
            x[ip] = ynew;  lnL = lnLnew;
            if (ip < s - 1)           naccept[0] += 1.0 / (s - 1);    /* t */
            else if (ip < s - 1 + g)    naccept[1] += 1.0 / g;        /* mu */
            else if (ip < s - 1 + g * 2)  naccept[4] += 1.0 / g;      /* sigma2 */
            else                 naccept[2] += 1.0 / g;               /* r0 for son0 */
         }
         else {
            if (ip < s - 1)                                       /* t */
               stree.nodes[s + ip].age = y;
            else if (ip - (s - 1) < g)                            /* mu */
               data.rgene[ip - (s - 1)] = y;
            else if (ip - (s - 1 + g) < g)                        /* sigma2 */
               data.sigma2[ip - (s - 1 + g)] = y;
            else                                                  /* r0 for son0 */
               stree.nodes[sons[0]].rates[ip - (s - 1 + g * 2)] = y;
         }
      }  /* for(ip, com.np) */

      if (MixingStep) {  /* this multiplies times by c and divides mu and r by c. */
         lnc = e[3] * rndSymmetrical();
         c = exp(lnc);
         lnacceptance = (s - 1 - g - g)*(lnc);

         gD = data.rgenepara;
         for (j = 0, sumold = 0; j < g; j++)    sumold += data.rgene[j];
         sumnew = sumold / c;
         lnacceptance += (gD[0] - gD[2] * g)*log(sumnew / sumold) - gD[1] / g*(sumnew - sumold) - (gD[2] - 1)*g*lnc;

         for (j = s; j < stree.nnode; j++)  stree.nodes[j].age *= c;
         for (i = 0; i < g; i++)             data.rgene[i] /= c;
         for (i = 0; i < g; i++)             stree.nodes[sons[0]].rates[i] /= c;
         lnLnew = lnPDFInfinitesitesClock23(FixedDs);
         lnacceptance += lnLnew - lnL;
         if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) {
            lnL = lnLnew;
            naccept[3]++;
         }
         else {
            for (j = s; j < stree.nnode; j++)  stree.nodes[j].age /= c;
            for (i = 0; i < g; i++)  data.rgene[i] *= c;
            for (i = 0; i < g; i++)  stree.nodes[sons[0]].rates[i] *= c;
         }
      }
      nround++;

      for (j = s, k = 0; j < stree.nnode; j++)  x[k++] = stree.nodes[j].age;
      for (i = 0; i < g; i++)             x[k++] = data.rgene[i];
      for (i = 0; i < g; i++)             x[k++] = data.sigma2[i];
      for (i = 0; i < g; i++)             x[k++] = stree.nodes[sons[0]].rates[i];
      for (i = 0; i < com.np; i++) mx[i] = (mx[i] * (nround - 1) + x[i]) / nround;
      for (i = 0; i < 5; i++) Pjump[i] = naccept[i] / nround;

      k = mcmc.sampfreq*mcmc.nsample;
      if (noisy && (k <= 10000 || (ir + 1) % (k / 2000) == 0)) {
         printf("\r%3.0f%%", (ir + 1.) / (mcmc.nsample*mcmc.sampfreq)*100.);
         for (j = 0; j < 5; j++) printf(" %4.2f", Pjump[j]);  printf(" ");
         if (com.np < nxpr[0] + nxpr[1]) { nxpr[0] = com.np; nxpr[1] = 0; }
         for (j = 0; j < nxpr[0]; j++) printf(" %5.3f", mx[j]);
         if (com.np > nxpr[0] + nxpr[1] && nxpr[1]) printf(" -");
         for (j = 0; j < nxpr[1]; j++) printf(" %5.3f", mx[com.np - nxpr[1] + j]);
         printf(" %5.2f", lnL);
      }
      if (mcmc.sampfreq*mcmc.nsample > 20 && (ir + 1) % (mcmc.sampfreq*mcmc.nsample / 20) == 0)
         printf(" %s\n", printtime(timestr));
      if (ir > 0 && (ir + 1) % mcmc.sampfreq == 0) {
         fprintf(fmcmc, "%d", ir + 1);
         for (i = 0; i < com.np; i++)
            fprintf(fmcmc, "\t%.5f", x[i]);  FPN(fmcmc);
      }
   }
   free(FixedDs);

   DescriptiveStatisticsSimpleMCMCTREE(fout, com.mcmcf, 1);

   exit(0);
}


int GetInitialsTimes(void)
{
   /* This checks the calibrations for consistency as follows: 
      If i is descendent to j, and tmin[i]>tmax[j], an error is printed out.
      This sets ages for nodes with bounds first.  
      It then goes through a loop looking for nodes with the largest min calibration, and 
      another loop looking for the longest path, each case assigning node ages on the path until 
      reaching an ancestral node with either an assigned age or max bound.
      In the case of shared/mirrored node ages due to gene duplication, conflicts may still arise,
      so that a loop of corrections is used.

      stree.duplication:
      stree.nodes[i].label = 0:  usual node, which may and may not have calibration.
      stree.nodes[j].label =-1:  node i is the driver node, whose age is shared by other nodes.
                                 It may and may not have calibration.
      stree.nodes[j].label = i:  node j shares the age of node i, and doesn't have calibration.
      */
   int s = stree.nspecies, s21 = s * 2 - 1, root = stree.root;
   int i, j, k, k1, ir, ncorrections, driver, from, to, depth, maxdepth;
   double *p, a, b, d, tbpath[2], *tmin, *tmax, *r;
   char *pptable, *flag;
   int debug = 0;

   if (com.TipDate && stree.nodes[stree.root].fossil != BOUND_F)
      error2("\nTipDate model requires bounds on the root age..\n");

   /* set up pop-pop table of ancestor-descendent relationships. */
   pptable = (char*)malloc(s21*(s21 + 1) * sizeof(char));
   if (pptable == NULL) error2("oom pptable");
   flag = pptable + s21*s21;
   memset(pptable, 0, s21*(s21 + 1) * sizeof(char));
   tmin = (double*)malloc(s21 * 2 * sizeof(double));
   if (tmin == NULL) error2("oom tmin");
   memset(tmin, 0, s21 * 2 * sizeof(double));
   tmax = tmin + s21;
   r = (double*)malloc(s * sizeof(double));
   if (r == NULL) error2("oom r");

   for (i = 0; i < s21; i++) pptable[i*s21 + i] = pptable[i*s21 + root] = (char)1;
   for (i = 0; i < s21; i++) {
      for (j = i; j != root; j = stree.nodes[j].father)
         pptable[i*s21 + j] = (char)1;
   }
   for (i = s; i < s21; i++) {
      if (stree.nodes[i].label != -1) continue;    /* -1 means driver */
      for (j = s; j < s21; j++) {
         if (stree.nodes[j].label == i) /* if j mirrors i, ancestors of j are also ancestors of i. */
            for (k = j; k != root && pptable[i*s21 + k] == 0; k = stree.nodes[k].father) {
               pptable[i*s21 + k] = (char)1;
            }
      }
      /* if j mirrors i, copy ancestors of i to j */
      for (j = s; j < s21; j++) {
         if (stree.nodes[j].label == i)
            memcpy(pptable + j*s21 + s, pptable + i*s21 + s, (s - 1) * sizeof(char));
      }
   }
   if (debug && s<20) {
      printf("\n\nSetting up initial node ages...\npop-pop table:\n     ");
      for (j = 0; j < s21; j++) printf(" %2d", j);
      for (i = 0; i < s21; i++) {
         printf("\n%2d:  ", i);
         for (j = 0; j < s21; j++) printf(" %2d", pptable[i*s21 + j]);
         if (i >= s) {
            printf(" %3s ", (stree.nodes[i].fossil ? fossils[stree.nodes[i].fossil] : ""));
            printf("  label = %2d, ", stree.nodes[i].label);
            if (stree.nodes[i].label == -1)      printf("driver");
            else if (stree.nodes[i].label > 0)   printf("mirror");
         }
      }
   }

   /* retrieve tmin & rmax from calibrations.  These are guestimates for gamma, SN, ST, S2N */
   for (i = s; i < s21; i++) stree.nodes[i].age = -1;
   for (i = s; i < s21; i++) {
      if (stree.nodes[i].fossil == 0)  continue;
      p = stree.nodes[i].pfossil;

      switch (stree.nodes[i].fossil) {
      case (LOWER_F):  tmin[i] = p[0];                  break;
      case (UPPER_F):  tmax[i] = p[1];                  break;
      case (BOUND_F):  tmin[i] = p[0]; tmax[i] = p[1];  break;
      case (GAMMA_F):
         tmin[i] = QuantileGamma(0.025, p[0], p[1]); tmax[i] = QuantileGamma(0.975, p[0], p[1]);
         break;
      case (SKEWN_F):
         tmin[i] = p[0] - 1.5*p[1];  tmax[i] = p[0] + 1.5*p[1];
         break;
      case (SKEWT_F):
         tmin[i] = p[0] - 1.5*p[1];  tmax[i] = p[0] + 1.5*p[1];
         break;
      case (S2N_F):
         a = p[0] * p[1] + (1 - p[0])*p[4];
         b = p[0] * p[2] + (1 - p[0])*p[5];
         tmin[i] = a - 1.5*b;  tmax[i] = a + 1.5*b;
         break;
      }
      if (tmin[i] < 0)      tmin[i] = 0;
      if (tmax[i] > OldAge) tmax[i] = OldAge/10;
   }
   for (i = s; i < s21; i++) {
      for (j = s; j < s21; j++)
         if (j != i && pptable[i*s21 + j] && tmax[j] > 0 && tmin[i] > tmax[j]) {
            printf("\n\ndaughter node %3d minimun age %9.4f \nmother node %3d max age %9.4f\n", i + 1, tmin[i], j + 1, tmax[j]);
            error2("strange calibration?");
         }
   }
   for (i = s; i < s21; i++) {
      if (tmin[i])   /* remove duplicated min bounds for ancestors of i */
         for (j = s; j < s21; j++)
            if (j != i && pptable[i*s21 + j] && tmin[j] <= tmin[i])
               tmin[j] = 0;
      if (tmax[i])   /* remove duplicated max bounds for descendents of i */
         for (j = s; j < s21; j++) {
            if (j != i && pptable[j*s21 + i] && tmax[j] >= tmax[i])
               tmax[j] = 0;
         }
   }
   if (debug) {
      printf("\n\nbounds on nodes after removing duplicated min and max bounds");
      for (i = s; i < s21; i++) {
         if(tmin[i] || tmax[i]) 
            printf("\n%2d:  (%8.4f, %8.4f)", i+1, tmin[i], tmax[i]);
      }
   }
   /* generate ages for bounded nodes if possible. */
   for (i = s; i < s21; i++) {
      if (tmin[i] && tmax[i])
         stree.nodes[i].age = tmin[i] + (tmax[i] - tmin[i])*rndu();
   }
   if (tmax[root] == 0) {
      if(stree.RootAge[1]>0) tmax[root] = stree.RootAge[1];
      else 
         error2("we need a maximum-age bound on root or RootAge?");
   }
   /* look for largest unused tmin, and assign ages on path from node to ancestor with age or with tmax. */
   for (ir = 0; ir < s; ir++) {
      tbpath[0] = 0;
      for (i = s; i < s21; i++) {  /* find max tmin and store in tbpath[0]. */
         if (stree.nodes[i].age == -1 && tbpath[0] < tmin[i]) {
            tbpath[0] = tmin[i];  from = i;
         }
      }
      if (tbpath[0] == 0) break;
      maxdepth = 1;
      for (j = stree.nodes[from].father; j != -1; j = stree.nodes[j].father, maxdepth++) {
         if (stree.nodes[j].age > 0) 
            break;
         else if (tmax[j] > 0) {
            maxdepth++;
            break;
         }
      }
      to = j;
      tbpath[1] = (stree.nodes[to].age > 0 ? stree.nodes[to].age : tmax[to]);

      for (j = 0; j < maxdepth; j++) r[j] = rndu();
      if (maxdepth>1) qsort(r, (size_t)maxdepth, sizeof(double), comparedouble);
      for (j = from, k = 0; k<maxdepth; j = stree.nodes[j].father, k++) {
         stree.nodes[j].age = tbpath[0] + (tbpath[1] - tbpath[0])*r[k];
         driver = (stree.nodes[j].label >= s ? stree.nodes[j].label : (stree.nodes[j].label == -1 ? j : 0));
         if (driver) {
            for (k1 = s; k1 < s21; k1++) {
               if (k1 != j && (k1 == driver || stree.nodes[k1].label == driver))
                  stree.nodes[k1].age = stree.nodes[j].age;
            }
         }
      }
      printf("\n*** min-age round %2d: from nodes %3d (%9.5f) to %3d, %3d nodes ", ir+1, from+1, stree.nodes[from].age, to+1, maxdepth);
      if(debug && s<200) printStree();
   }

   /* look for longest path, and assign ages on path from node to ancestor with age or with tmax. */
   for (ir = 0; ir < s; ir++) {
      maxdepth = 0;
      for (i = s; i < s21; i++) {
         if (stree.nodes[i].age > 0) continue;
         for (j = stree.nodes[i].father, depth = 1; j != -1; j = stree.nodes[j].father, depth++) {
            if (stree.nodes[j].age > 0) 
               break;
            else if (tmax[j] > 0) {
               depth++;
               break;
            }
         }
         if (maxdepth < depth) {
            maxdepth = depth;  from = i;  to = j;
         }
      }
      if (maxdepth == 0) break;  /* all ages are assigned. */
      int *sons = stree.nodes[from].sons;
      tbpath[0] = max2(stree.nodes[sons[0]].age, stree.nodes[sons[1]].age);
      tbpath[1] = (stree.nodes[to].age > 0 ? stree.nodes[to].age : tmax[to]);

      for (j = 0; j < maxdepth; j++) r[j] = rndu();
      if(maxdepth>1) qsort(r, (size_t)maxdepth, sizeof(double), comparedouble);
      for (j = from, k = 0; k<maxdepth; j = stree.nodes[j].father, k++) {
         stree.nodes[j].age = tbpath[0] + (tbpath[1] - tbpath[0])*r[k];
         driver = (stree.nodes[j].label >= s ? stree.nodes[j].label : (stree.nodes[j].label == -1 ? j : 0));
         if (driver) {
            for (k1 = s; k1 < s21; k1++) {
               if (k1 != j && (k1 == driver || stree.nodes[k1].label == driver))
                  stree.nodes[k1].age = stree.nodes[j].age;
            }
         }
      }
      if (debug && s < 200) {
         printf("\n*** longest-path round %2d: from nodes %3d (%9.5f) to %3d, %3d nodes ", ir + 1, from + 1, stree.nodes[from].age, to + 1, maxdepth);
         printStree();
      }
   }

   /* check for consistency for all nodes, and correct if needed. */
   for (ir = 0; ir < 100; ir++) {
      ncorrections = 0;
      for (i = 0; i < s21; i++) {
         j = stree.nodes[i].father;
         if (j == -1 || stree.nodes[i].age < stree.nodes[j].age) continue;
         if(debug) 
            printf("\n*** correction round %2d nodes %2d & %2d, ages %9.5f %9.5f", ir+1, i+1, j+1, stree.nodes[i].age, stree.nodes[j].age);
         ncorrections++;
         stree.nodes[j].age = stree.nodes[i].age * (1.001 + 0.01*rndu());
         /* copy nodes[j].age to all mirrored nodes */
         driver = (stree.nodes[j].label >= s ? stree.nodes[j].label : (stree.nodes[j].label == -1 ? j : 0));
         if (driver) {
            for (k = s; k < s21; k++) {
               if (k != j && (k == driver || stree.nodes[k].label == driver))
                  stree.nodes[k].age = stree.nodes[j].age;
            }
         }
      }
      if (ncorrections == 0) break;
   }

   free(pptable);  free(tmin);  free(r);
   return (0);
}

int GetInitials(void)
{
/* This sets the initial values for starting the MCMC, and returns np, the
   number of parameters in the MCMC, to be collected in collectx().
*/
   int np = 0, i, j, k, jj, dad, nchanges, g = data.ngene, g1 = g + (data.rgeneprior == 1);
   double a_r = data.rgenepara[0], b_r = data.rgenepara[1], a, b, smallr = 1e-3;

   com.rgene[0] = -1;  /* com.rgene[] is not used.  -1 to force crash */
   puts("\ngetting initial values to start MCMC.");

   GetInitialsTimes();

   /* initial mu (mean rates) for genes */
   np = stree.nspecies - 1 + g1;
   for (i = 0; i < g1; i++)
      data.rgene[i] = smallr + rndgamma(a_r) / b_r;   /* mu_i & mu_0 */

   if (com.clock > 1) {               /* sigma2, rates for nodes or branches */
      np += g1;
      if (mcmc.print >= 2) np += g*(stree.nnode - 1);

      /* sigma2 in lnrates among loci */
      for (i = 0; i < g1; i++)
         data.sigma2[i] = smallr + rndgamma(data.sigma2para[0]) / data.sigma2para[1];
      /* rates at nodes */
      for (j = 0; j < stree.nnode; j++) {
         if (j == stree.root) {
            for (i = 0; i < g; i++)  stree.nodes[j].rates[i] = -99;
         }
         else {
            for (i = 0; i < g; i++) {
               stree.nodes[j].rates[i] = smallr + rndgamma(a_r) / b_r;
            }
         }
      }
   }

   /* set up substitution parameters */
   if (mcmc.usedata == 1) {
      for (i = 0; i < g; i++)
         if (data.datatype[i] == BASE && com.model >= K80 && !com.fix_kappa) {
            data.kappa[i] = rndgamma(data.kappagamma[0]) / data.kappagamma[1] + 0.5;
            np++;
         }
      for (i = 0; i < g; i++)
         if (data.datatype[i] == BASE && !com.fix_alpha) {
            data.alpha[i] = rndgamma(data.alphagamma[0]) / data.alphagamma[1] + 0.1;
            np++;
         }
   }

   if (data.pfossilerror[0]) {
      a = data.pfossilerror[0];
      b = data.pfossilerror[1];
      data.Pfossilerr = a / (a + b)*(0.4 + 0.6*rndu());
      np++;
   }

   return(np);
}

int collectx(FILE* fout, double x[])
{
   /* this collects parameters into x[] for printing and summarizing.
      It returns the number of parameters.

        clock=1: times, rates for genes, kappa, alpha
        clock=0: times, rates or rates by node by gene, sigma2, rho_ij, kappa, alpha
   */
   int i, j, np = 0, g = data.ngene, g1 = g + (data.rgeneprior == 1);
   static int firsttime = 1;

   if (firsttime && fout)  fprintf(fout, "Gen");
   for (i = stree.nspecies; i < stree.nspecies * 2 - 1; i++) {
      if (firsttime && fout)  fprintf(fout, "\tt_n%d", i + 1);
      x[np++] = stree.nodes[i].age;
   }
   for (i = 0; i < g1; i++) {
      if (firsttime && fout) {
         if (i == g)      fprintf(fout, "\tmu_bar");
         else if (g > 1)  fprintf(fout, "\tmu%d", i + 1);
         else          fprintf(fout, "\tmu");
      }
      x[np++] = data.rgene[i];
   }
   if (com.clock > 1) {
      for (i = 0; i < g1; i++) {
         if (firsttime && fout) {
            if (i == g)      fprintf(fout, "\tsigma2_bar");
            else if (g > 1)  fprintf(fout, "\tsigma2_%d", i + 1);
            else          fprintf(fout, "\tsigma2");
         }
         x[np++] = data.sigma2[i];
      }
      if (mcmc.print >= 2)
         for (i = 0; i < g; i++) {
            for (j = 0; j < stree.nnode; j++) {
               if (j == stree.root) continue;
               if (firsttime && fout) {
                  if (g > 1) fprintf(fout, "\tr_g%d_n%d", i + 1, j + 1);
                  else    fprintf(fout, "\tr_n%d", j + 1);
               }
               x[np++] = stree.nodes[j].rates[i];
            }
         }
   }
   if (mcmc.usedata == 1) {
      if (!com.fix_kappa)
         for (i = 0; i < g; i++) {
            if (data.datatype[i] != 0) continue;  /* kappa in nucleotide alignments */
            if (firsttime && fout) {
               if (g > 1) fprintf(fout, "\tkappa_%d", i + 1);
               else    fprintf(fout, "\tkappa");
            }
            x[np++] = data.kappa[i];
         }

      if (!com.fix_alpha)
         for (i = 0; i < g; i++) {
            if (data.datatype[i] != 0) continue;  /* alpha in nucleotide alignments */
            if (firsttime && fout) {
               if (g > 1) fprintf(fout, "\talpha_%d", i + 1);
               else    fprintf(fout, "\talpha");
            }
            x[np++] = data.alpha[i];
         }
   }
   if (data.pfossilerror[0]) {
      if (firsttime && fout)  fprintf(fout, "\tpFossilErr");
      x[np++] = data.Pfossilerr;
   }

   if (np != com.np) {
      printf("np in collectx is %d != %d\n", np, com.np);
      if (!mcmc.usedata && (com.model || com.alpha)) printf("\nUse JC69 for no data");
      error2("");
   }
   if (firsttime && fout && mcmc.usedata) fprintf(fout, "\tlnL");
   if (firsttime && fout) fprintf(fout, "\n");

   firsttime = 0;
   return(0);
}


double lnpriorTimesBDS_Approach1(void)
{
   /* Approach 1 or Theorem #.# in Stadler & Yang (2013 Syst Biol).  Based on Tanja's notes of 27 July 2011.
      This calculates g(z*), which is constant throughout the MCMC, and thus wastes time.
   */
   int j, k;
   double lambda = data.BDS[0], mu = data.BDS[1], rho = data.BDS[2], psi = data.BDS[3];
   double lnp = 0, t, t1 = stree.nodes[stree.root].age, gt, gt1, a, e, c1, c2;
   double zstar, z0, z1;

   if (lambda <= 0 || mu < 0 || (rho <= 0 && psi <= 0))
      error2("B-D-S parameters:\n lambda > 0, mu >= 0, either rho > 0 or psi > 0");

   /* if(psi==0) the density is the same as in YR1997 */
   if (psi == 0 && fabs(lambda - mu) < 1e-20) {
      c1 = 1 / t1 + rho*lambda;
      for (j = stree.nspecies; j < stree.nnode; j++)
         if (j != stree.root) {
            c2 = 1 + rho*lambda*stree.nodes[j].age;
            lnp += log(c1 / (c2*c2));
         }
   }
   else if (psi == 0 && fabs(lambda - mu) > 1e-20) {
      a = lambda - rho*lambda - mu;
      e = exp((mu - lambda)*t1);
      c1 = (rho*lambda + a*e) / (1 - e);
      if (fabs(1 - e) < 1E-100)
         printf("e too close to 1..");
      /* a & c1 are constant over loop */
      for (j = stree.nspecies; j < stree.nnode; j++)
         if (j != stree.root) {
            e = exp((mu - lambda)*stree.nodes[j].age);
            c2 = (lambda - mu) / (rho*lambda + a*e);
            c2 *= c2*e*c1;
            if (c2 < 1E-300 || c2>1E300)
               printf("c2 not numeric.");
            lnp += log(c2);
         }
   }
   else {
      c1 = sqrt(square(lambda - mu - psi) + 4 * lambda*psi);
      c2 = -(lambda - mu - 2 * lambda*rho - psi) / c1;
      gt1 = 1 / (exp(-c1*t1)*(1 - c2) + (1 + c2));
      if (gt1<-1E-300 || gt1>1E300)
         printf("gt1 not numeric.");

      for (j = stree.nspecies; j < stree.nnode; j++) {
         if (j == stree.root) continue;
         for (k = stree.nodes[j].sons[0]; k >= stree.nspecies; k = stree.nodes[k].sons[1]);
         z0 = stree.nodes[k].age;
         /* printf("node %2d z[%d] = %7.4f ", j+1, k+1, z0); */

         for (k = stree.nodes[j].sons[1]; k >= stree.nspecies; k = stree.nodes[k].sons[0]);
         z1 = stree.nodes[k].age;
         zstar = max2(z0, z1);
         /* printf("z[%d] = %7.4f -> %7.4f\n", k+1, z1, zstar); */

         zstar = 1 / (exp(-c1*zstar)*(1 - c2) + (1 + c2));

         t = stree.nodes[j].age;
         gt = exp(-c1*t)*(1 - c2) + (1 + c2);
         lnp += -c1*t + log(c1*(1 - c2) / (gt*gt*(gt1 - zstar)));
      }
   }

   /* prior on t1 */
   if (stree.nodes[stree.nspecies].fossil != BOUND_F)
      error2("Only bounds for the root age are implemented.");
   lnp += lnptCalibrationDensity(t1, stree.nodes[stree.nspecies].fossil, stree.nodes[stree.nspecies].pfossil);

   return(lnp);
}

double logqtp0_BDS_Approach2(double t, double lambda, double mu, double rho, double psi, double *p0t)
{
   /* This calculates p0t and log(qt) for Approach 2 of Stadler & Yang (2013 Syst Biol).
   */
   double c1, c2, e = 0, r, logqt;

   if (psi == 0 && fabs(lambda - mu) <= 1e-20) {
      r = 1 + rho*lambda*t;
      *p0t = 1 - rho / (1 + rho*lambda*t);
      logqt = log(4 * r*r);
   }
   else if (psi == 0 && fabs(lambda - mu) > 1e-20) {
      e = exp((mu - lambda)*t);
      r = (rho*lambda + (lambda - rho*lambda - mu)*e) / (lambda - mu);
      *p0t = 1 - rho*(lambda - mu) / (rho*lambda + (lambda - rho*lambda - mu)*e);
      logqt = log(4 * r*r) - (mu - lambda)*t;
   }
   else {
      c1 = sqrt(square(lambda - mu - psi) + 4 * lambda*psi);
      if (c1 == 0) error2("c1==0");
      c2 = -(lambda - mu - 2 * lambda*rho - psi) / c1;
      e = exp(-c1*t);
      r = (e*(1 - c2) - (1 + c2)) / (e*(1 - c2) + (1 + c2));

      *p0t = (lambda + mu + psi + c1*r) / (2 * lambda);
      logqt = c1*t;
      logqt += log(e*e*square(1 - c2) + e * 2 * (1 - c2*c2) + square(1 + c2));
   }
   /*
   if(e<=1e-300)
      printf("lmrp %8.4f %8.4f %7.4f %7.4f t=%7.4f p0 logqt= %7.4f %7.4e\n", lambda, mu, rho, psi, t, *p0t, logqt);
   */
   return(logqt);
}

double lnpriorTimesTipDate_Approach2(void)
{
   /* Approach 2 of Stadler & Yang (2013 Syst Biol).  The BDS parameters are assumed to be fixed.
      This ignores terms involving y, which are constants when the BDSS parameters are fixed.
      t1 is root age.
   */
   int i, m = 0, n = stree.nspecies, k = 0;
   double lambda = data.BDS[0], mu = data.BDS[1], rho = data.BDS[2], psi = data.BDS[3];
   double lnp = 0, t1 = stree.nodes[stree.root].age, p0, p1, logqt, logqt1, p0t1, z, e;

   if (rho == 0) {
      m = stree.nspecies;  n = 0;
   }
   else if (com.TipDate) {
      for (i = 0, m = 0; i < stree.nspecies; i++)
         m += (stree.nodes[i].age > 0);
      n = stree.nspecies - m;
   }

   if (lambda <= 0 || mu < 0 || (rho <= 0 && psi <= 0))
      error2("wrong B-D-S parameters");
   if (rho > 0 && n < 2)
      error2("we want more than 2 extant sampled species");

   /* the numerator f(x, z | L(tmrca) = 2).   */
   logqt1 = logqtp0_BDS_Approach2(t1, lambda, mu, rho, psi, &p0t1);
   lnp -= 2 * logqt1;
   for (i = stree.nspecies; i < stree.nnode; i++) { /* loop over n+m-1 internal node ages except t1  */
      if (i == stree.root)  continue;
      logqt = logqtp0_BDS_Approach2(stree.nodes[i].age, lambda, mu, rho, psi, &p0);
      lnp -= logqt;
   }

   /* the denominator, h(tmrca, n) */
   if (rho) {
      if (fabs(lambda - mu - psi) > 1e-20) {
         e = exp(-(lambda - mu - psi)*t1);   /* this causes overflow for large t */
         z = lambda*rho + (lambda - lambda*rho - mu - psi)*e;
         p1 = rho*square(lambda - mu - psi)*e / (z*z);
         if (p1 <= 0)
            printf("lnpriorTimesTipDate: p1 = %.6f <=0 (t1 = %9.5g)\n", p1, t1);
         lnp -= log(p1*p1);
         if (n > 2) {
            z = rho*lambda*(1 - e) / z;
            if (z <= 0)
               printf("z = %.4f < 0\n", z);
            lnp -= (n - 2)*log(z);
         }
      }
      else {
         if (n <= 1) error2("we don't like n <= 1");
         if (n > 2) lnp -= (n - 2)*log(t1);
         lnp += (n + 2)*log(1 + rho*lambda*t1);
      }
   }
   else {  /* rho = 0 for viruses */
      lnp -= 2 * log(1 - p0t1);
   }

   /* prior on t1 */
   if (stree.nodes[stree.nspecies].fossil != BOUND_F)
      error2("Only bounds for the root age are implemented.");
   lnp += lnptCalibrationDensity(t1, stree.nodes[stree.nspecies].fossil, stree.nodes[stree.nspecies].pfossil);

   return(lnp);
}

double p01t_BDSEquation6Stadler2010(double t, double lambda, double mu, double rho, double psi, double *p0t)
{
   /* This calculates p0t and p1t, defined in equations 1 & 2 in Stadler (2010).
   */
   double c1, c2, e, r, p1t;

   if (fabs(lambda - mu) < 1e-20 && psi == 0) {
      r = 1 / (1 / rho + lambda*t);
      *p0t = 1 - r;
      p1t = r*r / rho;
   }
   else if (fabs(lambda - mu) > 1e-20 && psi == 0) {
      e = exp((mu - lambda)*t);
      r = rho*(lambda - mu) / (rho*lambda + (lambda - mu - rho*lambda)*e);
      *p0t = 1 - r;
      p1t = r*r / rho*e;
   }
   else {
      c1 = sqrt(square(lambda - mu - psi) + 4 * lambda*psi);
      if (c1 == 0) error2("c1==0");
      c2 = -(lambda - mu - 2 * lambda*rho - psi) / c1;
      e = exp(-c1*t);
      r = (e*(1 - c2) - (1 + c2)) / (e*(1 - c2) + (1 + c2));

      *p0t = (lambda + mu + psi + c1*r) / (2 * lambda);
      p1t = 4 * rho / (2 * (1 - c2*c2) + e*square(1 - c2) + square(1 + c2) / e);
   }

   /* printf("lmrp %8.6f %8.6f %7.4f %7.4f t=%7.4f p01= %7.4f %7.4f\n", lambda, mu, rho, psi, t, *p0t, p1t); */

   return(p1t);
}

double lnpriorTimesTipDateEquation6Stadler2010(void)
{
   /* Equation 6 in Stadler T. 2010. Sampling-through-time in birth-death trees.
      J Theor Biol 267:396-404.
      t1 is root age.
   */
   int i, m, n, k = 0;
   double lambda = data.BDS[0], mu = data.BDS[1], rho = data.BDS[2], psi = data.BDS[3];
   double lnp = 0, t1 = stree.nodes[stree.root].age, p0, p1, z, e;

   if (lambda <= 0 || mu <= 0 || (rho <= 0 && psi <= 0))
      error2("wrong B-D-S parameters..");
   for (i = 0, m = 0; i < stree.nspecies; i++)
      m += (stree.nodes[i].age > 0);
   n = stree.nspecies - m;

   if (n > 2) {
      if (fabs(lambda - mu) < 1e-20)
         z = 1 + 1 / (rho*lambda*t1);
      else {
         e = exp((mu - lambda)*t1);
         z = (lambda*rho + (lambda - mu - lambda*rho)*e) / (rho*lambda*(1 - e));
      }
      lnp = (n - 2)*log(z*z);
   }
   if (n > 1) {
      p1 = p01t_BDSEquation6Stadler2010(t1, lambda, mu, rho, 0, &p0);
      lnp -= log((n - 1)*p1*p1);
   }
   /*  this term is constant when the BDS parameters are fixed.
   lnp += (n+m-2)*log(lambda) + (k+m)*log(psi);
   */

   p1 = p01t_BDSEquation6Stadler2010(t1, lambda, mu, rho, psi, &p0);
   lnp += 2 * log(p1);
   for (i = stree.nspecies; i < stree.nnode; i++)  /* loop over n+m-1 internal node ages except t1  */
      if (i != stree.root)
         lnp += log(p01t_BDSEquation6Stadler2010(stree.nodes[i].age, lambda, mu, rho, psi, &p0));

   /* the y terms do not change when the BDS parameters are fixed.
   for(i=0; i<stree.nspecies; i++) {
      if( (t=stree.nodes[i].age) ) {
         p1 = p01t_BDSEquation6Stadler2010(t, lambda, mu, rho, psi, &p0);
         lnp += log(p0/p1);
      }
   }
   */

   /* prior on t1 */
   if (stree.nodes[stree.nspecies].fossil != BOUND_F)
      error2("Only bounds for the root age are implemented.");
   lnp += lnptCalibrationDensity(t1, stree.nodes[stree.nspecies].fossil, stree.nodes[stree.nspecies].pfossil);
   return(lnp);
}


#define P0t_YR07(expmlt) (rho*(lambda-mu)/(rho*lambda +(lambda*(1-rho)-mu)*expmlt))

double lnPDFkernelBDS_YR07(double t, double t1, double vt1, double lambda, double mu, double rho)
{
   /* this calculates the log of the pdf of the BDS kernel.
   */
   double lnp, mlt = (mu - lambda)*t, expmlt, small = 1e-20;

   if (fabs(mu - lambda) < small)
      lnp = log((1 + rho*lambda*t1) / (t1*square(1 + rho*lambda*t)));
   else {
      expmlt = exp(mlt);
      lnp = P0t_YR07(expmlt);
      if (lnp < 0) puts("this should not be negative, in lnPDFkernelBDS_YR07");
      lnp = log(lnp*lnp * lambda / (vt1*rho) * expmlt);
   }
   return(lnp);
}

double CDFkernelBDS_YR07(double t, double t1, double vt1, double lambda, double mu, double rho)
{
   /* this calculates the CDF for the kernel density from the BDS model
   */
   double cdf, expmlt, small = 1e-20;

   if (fabs(lambda - mu) < small)
      cdf = (1 + rho*lambda*t1)*t / (t1*(1 + rho*lambda*t));
   else {
      expmlt = exp((mu - lambda)*t);
      if (expmlt < 1e10)
         cdf = rho*lambda / vt1 * (1 - expmlt) / (rho*lambda + (lambda*(1 - rho) - mu)*expmlt);
      else {
         expmlt = 1 / expmlt;
         cdf = rho*lambda / vt1 * (expmlt - 1) / (rho*lambda*expmlt + (lambda*(1 - rho) - mu));
      }
   }
   return(cdf);
}


double lnPDFkernelBeta(double t, double t1, double p, double q)
{
   /* This calculates the PDF for the beta kernel.
      The normalizing constant is calculated outside this routine.
   */
   double lnp, x = t / t1;

   if (x < 0 || x>1) error2("t outside of range (0, t1)");
   lnp = (p - 1)*log(x) + (q - 1)*log(1 - x);
   lnp -= log(t1);
   return(lnp);
}


double lnptNCgiventC(void)
{
   /* This calculates f_BDS(tNC|tC), the conditional probability of no-calibration
      node ages given calibration node ages, that is, the first term in equation 3
      in Yang & Rannala (2006).  This is called by lnpriorTimes().

      The beta kernel is added on 5 October 2009, specified by data.priortime=1.

      The routine uses stree.nodes[].pfossil, stree.nodes[].fossil,
      and lambda, mu, rho from the birth-death process with species sampling
      (data.BDS[0-3]).
      It sorts the node ages in the species tree and then uses the
      birth-death prior conditional on the calibration points.
      t[0] is t1 in Yang & Rannala (2006).

      rankt[3]=1 means that the 3rd yougest age is node number ns+1.  This is set up
      for the calibration nodes only, = -1 for non-calibration nodes.  First we collect
      all times into t[] and sort them. Next we collect all calibration times into tc[]
      and sort them.  Third, we find the ranks of calibration times in tc[], that is,
      i1, i2, ..., ic in YB06.

      The term for root is not in this routine.  The root is excluded from the ranking.
   */
   int  i, j, k, n1 = stree.nspecies - 1, rankprev, nfossil = 0;
   int  ranktc[MaxNFossils]; /* ranks of calibration nodes */
   double t[NS - 1], tc[MaxNFossils], t1 = stree.nodes[stree.root].age;
   double lnpt = 0, expmlt = 1, vt1, P0t1, cdf, cdfprev, small = 1e-20;
   double lambda = data.BDS[0], mu = data.BDS[1], rho = data.BDS[2];
   double p = data.BDS[0], q = data.BDS[1], lnbeta = 0;
   int debug = 0;

   /* (A) Calculate f(tNC), joint of the non-calibration node ages.
          The root age is excluded, so the loop starts from j = 1.
          In the calculation of (eq.9)/(eq.11), (s - 2)! cancels out, as does
          the densities g(t_i_1) etc. in eq.11.  After cancellation, only the
          densities of the non-calibration node ages remain in eq.9.
   */
   if (data.priortime == 0) {  /* BDS kernel */
      if (rho < 0) error2("rho < 0");
      /* vt1 is needed only if (lambda != mu) */
      if (fabs(lambda - mu) > small) {
         expmlt = exp((mu - lambda)*t1);
         P0t1 = P0t_YR07(expmlt);
         vt1 = 1 - P0t1 / rho*expmlt;
      }
      else {
         P0t1 = rho / (1 + rho*mu*t1);
         vt1 = mu*t1*P0t1;
      }

      for (j = 1, lnpt = 0; j < n1; j++) {
         k = stree.nspecies + j;
         if (!stree.nodes[k].usefossil)   /* non-calibration node age */
            lnpt += lnPDFkernelBDS_YR07(stree.nodes[k].age, t1, vt1, lambda, mu, rho);
      }
   }
   else {                   /* beta kernel */
      lnbeta = LnGamma(p) + LnGamma(q) - LnGamma(p + q);
      for (j = 1, lnpt = 0; j < n1; j++) {
         k = stree.nspecies + j;
         if (!stree.nodes[k].usefossil)   /* non-calibration node age */
            lnpt += lnPDFkernelBeta(stree.nodes[k].age, t1, p, q) + lnbeta;
      }
   }

   /* (B) Divide by f(tC), marginal of calibration node ages (eq.9/eq.11).
      This goes through the calibration nodes in the order of their ages,
      with sorted node ages stored in tc[].
   */
   for (j = stree.nspecies; j < stree.nnode; j++) {
      t[j - stree.nspecies] = stree.nodes[j].age;
      if (j != stree.root && stree.nodes[j].usefossil)
         tc[nfossil++] = stree.nodes[j].age;
   }
   if (nfossil > MaxNFossils) error2("raise MaxNFossils?");

   if (nfossil) {
      /* The only reason for sorting t[] is to construct ranktc[].  */
      qsort(t, (size_t)n1, sizeof(double), comparedouble);
      qsort(tc, (size_t)nfossil, sizeof(double), comparedouble);

      for (i = j = 0; i < nfossil; i++) {
         if (i) j = ranktc[i - 1] + 1;
         for (; j < n1; j++)
            if (tc[i] < t[j]) break;
         ranktc[i] = j;
      }
      if (debug == 1) {
         matout2(F0, t, 1, n1, 9, 5);
         matout2(F0, tc, 1, nfossil, 9, 5);
         for (i = 0; i < nfossil; i++)
            printf("%9d", ranktc[i]);
         FPN(F0);
      }
   }

   for (i = 0, rankprev = 0, cdfprev = 0; i < nfossil + 1; i++) {
      if (i < nfossil) {
         if (data.priortime == 0)   /* BDS kernel */
            cdf = CDFkernelBDS_YR07(tc[i], t1, vt1, lambda, mu, rho);
         else                    /* beta kernel */
            cdf = CDFBeta(tc[i] / t1, p, q, lnbeta);

         k = ranktc[i] - rankprev - 1;
      }
      else {
         cdf = 1;
         k = n1 - rankprev - 1;
      }
      if (k > 0) {
         if (cdf <= cdfprev) error2("cdf diff<0 in lnptNCgiventC");
         lnpt += LnGamma(k + 1.0) - k*log(cdf - cdfprev);
      }
      rankprev = ranktc[i];
      cdfprev = cdf;
   }
   if (debug == 1) printf("\npdf = %.12f\n", exp(lnpt));


   /*
   {
   double t1=stree.nodes[4].age, t2=stree.nodes[5].age, t3=stree.nodes[6].age, Gt2, Gt3;
   Gt2 = CDFkernelBDS_YR07(t2, t1, vt1, lambda, mu, rho);
   Gt3 = CDFkernelBDS_YR07(t3, t1, vt1, lambda, mu, rho);
   if(stree.nodes[5].usefossil==0 && stree.nodes[6].usefossil==0)
      lnpt += log(1.0/2);
   if(stree.nodes[5].usefossil==1 && stree.nodes[6].usefossil==0)
      if(t3<t2) lnpt += log(Gt2);
      else      lnpt += log(1-Gt2);
   if(stree.nodes[5].usefossil==0 && stree.nodes[6].usefossil==1)
      if(t2<t3) lnpt += log(Gt3);
      else      lnpt += log(1-Gt3);

   }
   */
   return (lnpt);
}


double lnptCalibrationDensity(double t, int fossil, double p[7])
{
   /* This calculates the log of calibration density
   */
   double lnpt, a = p[0], b = p[1], thetaL, thetaR, tailL = 99, tailR = 99, s, z, t0, P, c, A;

   switch (fossil) {
   case (LOWER_F):  /* truncated Cauchy C(a(1+p), s^2). tL = a. */
      P = p[1];
      c = p[2];
      tailL = p[3];
      t0 = a*(1 + P);
      s = a*c;
      A = 0.5 + 1 / Pi*atan(P / c);
      if (t > a) {
         z = (t - t0) / s;
         lnpt = log((1 - tailL) / (Pi*A*s*(1 + z*z)));
      }
      else {
         z = P / c;
         thetaL = (1 / tailL - 1) / (Pi*A*c*(1 + z*z));
         lnpt = log(tailL*thetaL / a) + (thetaL - 1)*log(t / a);
      }
      break;
   case (UPPER_F):
      /* equation (16) in Yang and Rannala (2006) */
      tailR = p[2];
      if (t < b)
         lnpt = log((1 - tailR) / b);
      else {
         thetaR = (1 - tailR) / (tailR*b);
         lnpt = log(tailR*thetaR) - thetaR*(t - b);
      }
      break;
   case (BOUND_F):
      tailL = p[2];
      tailR = p[3];
      if (t > a && t < b)
         lnpt = log((1 - tailL - tailR) / (b - a));
      else if (t < a) {
         thetaL = (1 - tailL - tailR)*a / (tailL*(b - a));
         lnpt = log(tailL*thetaL / a) + (thetaL - 1)*log(t / a);
      }
      else {
         thetaR = (1 - tailL - tailR) / (tailR*(b - a));
         lnpt = log(tailR*thetaR) - thetaR*(t - b);
      }
      break;
   case (GAMMA_F):
      lnpt = a*log(b) - b*t + (a - 1)*log(t) - LnGamma(a);
      break;
   case (SKEWN_F):
      lnpt = logPDFSkewN(t, p[0], p[1], p[2]);
      break;
   case (SKEWT_F):
      lnpt = log(PDFSkewT(t, p[0], p[1], p[2], p[3]));
      break;
   case (S2N_F):
      a = PDFSkewN(t, p[1], p[2], p[3]);
      b = PDFSkewN(t, p[4], p[5], p[6]);
      lnpt = log(p[0] * a + (1 - p[0])*b);
      break;
   }
   return(lnpt);
}


double lnptC(void)
{
/* This calculates the prior density of times at calibration nodes as specified
   by the fossil calibration information, the second term in equation 3 in
   Yang & Rannala (2006).
   a=tL, b=tU.

   The term for root is always calculated in this routine.
   If there is a fossil at root and if it is a lower bound, it is re-set to a
   pair of bounds.
*/
   int i, j, nfossil = 0, fossil;
   double p[7], t, lnpt = 0;
   int debug = 0;

   /* RootAge[] will cause stree.nodes[stree.root].fossil to be set before this routine. */
   if (stree.nfossil == 0 && stree.nodes[stree.root].fossil == 0)
      return(0);

   for (j = stree.nspecies; j < stree.nnode; j++) {
      if (j != stree.root && !stree.nodes[j].usefossil)
         continue;
      nfossil++;
      t = stree.nodes[j].age;
      fossil = stree.nodes[j].fossil;
      for (i = 0; i < 7; i++) p[i] = stree.nodes[j].pfossil[i];

      if (j != stree.root || (stree.nodes[j].usefossil && fossil != LOWER_F)) {
         if (fossil == LOWER_F)
            p[3] = stree.nodes[j].pfossil[3];
         else if (fossil == UPPER_F)
            p[2] = stree.nodes[j].pfossil[2];
         else if (fossil == BOUND_F) {
            p[2] = stree.nodes[j].pfossil[2];
            p[3] = stree.nodes[j].pfossil[3];
         }
      }
      else if (!stree.nodes[j].usefossil) {  /* root fossil absent or deleted, using RootAge */
         p[0] = stree.RootAge[0];
         p[1] = stree.RootAge[1];
         if (stree.RootAge[0] > 0) {
            fossil = BOUND_F;
            p[2] = stree.RootAge[2];
            p[3] = stree.RootAge[3];
         }
         else {
            fossil = UPPER_F;
            p[2] = stree.RootAge[2];
         }
      }
      else if (fossil == LOWER_F) {
         /* root fossil is L.  B(a,b,tailL, tailR) is used.  p & c ignored. */
         fossil = BOUND_F;
         p[1] = stree.RootAge[1];
         p[3] = stree.nodes[j].pfossil[3];
         p[2] = (stree.RootAge[0] > 0 ? stree.RootAge[3] : stree.RootAge[2]);
      }

      lnpt += lnptCalibrationDensity(t, fossil, p);
   }

   return(lnpt);
}


double getScaleFossilCombination(void);

double getScaleFossilCombination(void)
{
   /* This uses Monte Carlo integration to calculate the normalizing constant for the
      joint prior on fossil times, when some fossils are assumed to be in error and not
      used.  The constant is the integral of the density over the feasible region, where
      the times satisfy the age constraints on the tree.
      It is assumed that the ancestral nodes have smaller node numbers than descendent
      nodes, as is the case if the node numbers are assigned by ReadTreeN().
      nfs = nfossilused if a fossil is used for the root or nfs = nfossilused + 1
      if otherwise.
      CLower[] contains constants for lower bounds, to avoid duplicated computation.
   */
   int nrepl = 5000000, ir, i, j, k, nfs, ifs[MaxNFossils] = { 0 }, feasible;
   int root = stree.root, rootfossil = stree.nodes[root].fossil, fossil;
   signed char ConstraintTab[MaxNFossils][MaxNFossils] = { {0} };  /* 1: row>col; 0: irrelevant */
   double C, sC, accept, mt[MaxNFossils] = { 0 }, t[MaxNFossils], *pfossil[MaxNFossils], pfossilroot[4];
   double importance, CLower[MaxNFossils][3];
   double sNormal = 1, r, thetaL, thetaR, a, b, c, p, s, A, zc, zn, tailL, tailR;
   int debug = 1;

   /* Set up constraint table */
   for (i = stree.nspecies, nfs = 0; i < stree.nnode; i++) {
      if (i == root || stree.nodes[i].usefossil)
         ifs[nfs++] = i;
   }
   for (i = 1; i < nfs; i++) {
      for (j = 0; j < i; j++) {
         for (k = ifs[i]; k != -1; k = stree.nodes[k].father) {
            if (k == ifs[j]) { ConstraintTab[i][j] = 1; break; }
         }
      }
   }
   if (debug) {
      for (i = 0; i < nfs; i++) {
         printf("\n%d (%2d %s): ", i + 1, ifs[i] + 1, fossils[stree.nodes[ifs[i]].fossil]);
         for (j = 0; j < i; j++)
            printf("%4d", ConstraintTab[i][j]);
      }
      FPN(F0);
   }

   /* Set up calibration info.  Save constants for lower bounds */
   for (i = 0; i < 4; i++)
      pfossilroot[i] = stree.nodes[root].pfossil[i];
   if (ifs[0] == root) {                    /* copy info for root into pfossilroot[] */
      if (!stree.nodes[root].usefossil) {  /* root fossil absent or excluded */
         for (i = 0; i < 4; i++)
            pfossilroot[i] = stree.RootAge[i];
         rootfossil = (stree.RootAge[0] > 0 ? BOUND_F : UPPER_F);
      }
      else if (stree.nodes[root].fossil == LOWER_F) {   /* root fossil is lower bound */
         pfossilroot[1] = stree.RootAge[1];
         rootfossil = BOUND_F;
      }
   }
   for (i = 0; i < nfs; i++) {
      j = ifs[i];
      pfossil[i] = (j == root ? pfossilroot : stree.nodes[j].pfossil);

      if (j != root && stree.nodes[j].fossil == LOWER_F) {
         a = pfossil[i][0];
         p = pfossil[i][1];
         c = pfossil[i][2];
         tailL = pfossil[i][3];
         s = a*c*sNormal;
         A = 0.5 + 1 / Pi*atan(p / c);
         CLower[i][0] = (1 / tailL - 1) / (Pi*A*c*(1 + (p / c)*(p / c)));  /* thetaCauchy */
         CLower[i][1] = (1 / tailL - 1) * a / (s*sqrt(Pi / 2));            /* thetaNormal  */
         CLower[i][2] = s / (A*c*a*sqrt(2 * Pi));                      /* s/Act*sqrt(2Pi) */
      }
   }

   /* Take samples */
   for (ir = 0, C = sC = accept = 0; ir < nrepl; ir++) {
      for (i = 0, importance = 1; i < nfs; i++) {
         j = ifs[i];
         fossil = (j == root ? rootfossil : stree.nodes[j].fossil);
         a = pfossil[i][0];
         b = pfossil[i][1];
         r = rndu();
         switch (fossil) {
         case(LOWER_F):  /* simulating from folded normal, the importance density */
            tailL = pfossil[i][3];
            if (r < tailL) {   /* left tail, CLower[i][1] has thetaNormal */
               thetaL = CLower[i][1];
               t[i] = a * pow(rndu(), 1 / thetaL);
               importance *= CLower[i][0] / thetaL * pow(t[i] / a, CLower[i][0] - thetaL);
            }
            else {
               c = pfossil[i][2];
               s = a*c*sNormal;
               t[i] = a + fabs(rndNormal()) * s;
               zc = (t[i] - a*(1 + b)) / (c*a);
               zn = (t[i] - a) / s;
               importance *= CLower[i][2] * exp(zn*zn / 2) / (1 + zc*zc);
            }
            break;
         case(UPPER_F):
            tailR = pfossil[i][2];
            if (r > tailR) { /* flat part */
               t[i] = b*rndu();
            }
            else {  /* right tail */
               thetaR = (1 - tailR) / (tailR*b);
               t[i] = b - log(rndu()) / thetaR;
            }
            break;
         case (BOUND_F):
            tailL = pfossil[i][2];
            tailR = pfossil[i][3];
            if (r > tailL + tailR)  /* flat part */
               t[i] = a + (b - a) * rndu();
            else if (r < tailL) {  /* left tail */
               thetaL = (1 - tailL - tailR)*a / (tailL*(b - a));
               t[i] = a * pow(rndu(), 1 / thetaL);
            }
            else {                 /* right tail */
               thetaR = (1 - tailL - tailR) / (tailR*(b - a));
               t[i] = b - log(rndu()) / thetaR;
            }
            break;
         case(GAMMA_F):
            t[i] = rndgamma(a) / b;
            break;
         default:
            printf("\nfossil = %d (%s) not implemented.\n", fossil, fossils[fossil]);
            exit(-1);
         }
      }

      feasible = 1;
      for (i = 1; i < nfs; i++) {
         for (j = 0; j < i; j++)
            if (ConstraintTab[i][j] && t[i] > t[j])
            {
               feasible = 0; break;
            }
         if (feasible == 0) break;
      }
      if (feasible) {
         accept++;
         C += importance;
         sC += importance*importance;
         for (i = 0; i < nfs; i++)
            mt[i] += t[i] * importance;
      }
      if ((ir + 1) % 100000 == 0) {
         a = ir + 1;
         s = (sC / a - C / a*C / a) / a;
         printf("\r%7d %.2f%% C = %.6f +- %.6f mt", ir + 1, accept / (ir + 1) * 100, C / (ir + 1), sqrt(s));
         for (i = 0; i < nfs; i++)
            printf(" %6.4f", mt[i] / C);
      }
   }   /* for(ir) */
   return(C / nrepl);
}


int SetupPriorTimesFossilErrors(void)
{
   /* This prepares for lnpriorTimes() under models of fossil errors.  It calculates
      the scaling factor for the probability density of times for a given combination
      of the indicator variables, indicating which fossil is in error and not used.
      nMinCorrect is the minimum number of correct fossils.
   */
   int nMinCorrect = (int)data.pfossilerror[2], ncomFerr, is, i, icom, nused = 0, it;
   double t;
   int debug = 1;

   if (data.pfossilerror[0] == 0 || stree.nfossil > MaxNFossils || stree.nfossil < nMinCorrect)
      error2("Fossil-error model is inappropriate.");

   for (i = nMinCorrect, ncomFerr = 0; i <= stree.nfossil; i++) {
      ncomFerr += (int)(Binomial((double)stree.nfossil, i, &t) + .1);
      if (t != 0) error2("too many fossils to deal with? ");
   }

   data.CcomFossilErr = (double*)malloc(ncomFerr * sizeof(double));
   if (data.CcomFossilErr == NULL) error2("oom for CcomFossilErr");

   /* cycle through combinations of indicators. */
   for (i = 0, icom = 0; i < (1 << stree.nfossil); i++) {
      it = i;
      for (is = stree.nspecies, nused = 0; is < stree.nnode; is++) {
         if (stree.nodes[is].fossil) {
            stree.nodes[is].usefossil = 1 - it % 2;
            nused += stree.nodes[is].usefossil;
            it /= 2;
            if (debug) printf("%2d (%2d)", stree.nodes[is].usefossil, is + 1);
         }
      }
      if (nused < nMinCorrect) continue;
      if (debug) printf("\n\n********* Combination %2d/%2d:  %2d fossils used", icom + 1, ncomFerr, nused);
      data.CcomFossilErr[icom++] = log(getScaleFossilCombination());
   }
   return(0);
}



double lnpriorTimes(void)
{
   /* This calculates the prior density of node times in the master species tree.
         Node ages are in stree.nodes[].age.
   */
   int nMinCorrect = (int)data.pfossilerror[2];
   int  is, i, icom, nused, it;
   double pE = data.Pfossilerr, lnpE, ln1pE, pEconst = 0, lnC, lnNC = 0;
   double lnpt = 0, scaleF = -1e300, y;
   int debug = 0;

   if (stree.nfossil == 0 || com.TipDate) {
      if (1)        lnpt = lnpriorTimesBDS_Approach1();
      else if (0)   lnpt = lnpriorTimesTipDate_Approach2();
      else if (0)   lnpt = lnpriorTimesTipDateEquation6Stadler2010();
   }
   else if (data.pfossilerror[0] == 0) { /* no fossil errors in model */
      lnC = lnptC();
      if (stree.duplication == 0) lnNC = lnptNCgiventC();
      lnpt = lnC + lnNC;

      if (debug) printf("\nftC ftNC ft = %9.5f%9.5f%9.5f", exp(lnC), exp(lnNC), exp(lnC)*exp(lnNC));
   }
   else {   /* fossil errors: cycle through combinations of indicators, using nMinCorrect. */
      if (nMinCorrect) {
         pEconst = y = pow(pE, (double)stree.nfossil);
         if (y < 1e-300) error2("many fossils & small pE.  Rewrite the code?");
         for (i = 1; i < nMinCorrect; i++) /* i is the number of correct fossils used */
            pEconst += y *= (stree.nfossil - i + 1)*(1 - pE) / (i*pE);
         pEconst = log(1 - pEconst);
      }
      lnpE = log(pE);
      ln1pE = log(1 - pE);
      for (i = 0, icom = 0; i < (1 << stree.nfossil); i++) { /* sum over U */
         it = i;
         for (is = stree.nspecies, nused = 0; is < stree.nnode; is++) {
            if (stree.nodes[is].fossil) {
               stree.nodes[is].usefossil = 1 - it % 2;
               nused += stree.nodes[is].usefossil;
               it /= 2;
            }
         }
         if (nused < nMinCorrect) continue;
         lnC = lnptC();
         if (stree.duplication == 0) lnNC = lnptNCgiventC();

         y = nused*ln1pE + (stree.nfossil - nused)*lnpE - pEconst - data.CcomFossilErr[icom++]
            + lnC + lnNC;

         if (debug)
            printf("\nU%d%d ftC ftNC ft = %9.5f%9.5f%9.5f", stree.nodes[5].usefossil, stree.nodes[6].usefossil,
               exp(lnC), exp(lnNC), exp(lnC)*exp(lnNC));

         if (y < scaleF + 100)
            lnpt += exp(y - scaleF);
         else {
            lnpt = 1;
            scaleF = y;
         }
         lnpt += exp(y - scaleF);
      }
      lnpt = scaleF + log(lnpt);
   }
   if (debug) exit(0);

   return (lnpt);
}

int getPfossilerr(double postEFossil[], double nround)
{
   /* This is modified from lnpriorTimes(), and prints out the conditonal Perror
      for each fossil given the times and pE.
   */
   int nMinCorrect = (int)data.pfossilerror[2];
   int  is, i, icom, k, nf = stree.nfossil, nused, it;
   double pE = data.Pfossilerr, lnpE, ln1pE, pEconst = 0;
   double pt, pU[MaxNFossils] = { 0 }, scaleF = -1e300, y;
   char U[MaxNFossils];  /* indicators of fossil errors */

   if (data.priortime == 1)
      error2("beta kernel for prior time not yet implemented for model of fossil errors.");
   if (nMinCorrect) {
      pEconst = y = pow(pE, (double)stree.nfossil);
      for (i = 1; i < nMinCorrect; i++) /* i is the number of correct fossils used */
         pEconst += y *= (stree.nfossil - i + 1)*(1 - pE) / (i*pE);
      pEconst = log(1 - pEconst);
   }

   lnpE = log(pE); ln1pE = log(1 - pE);
   for (i = 0, icom = 0, pt = 0; i < (1 << stree.nfossil); i++) { /* sum over U */
      it = i;
      for (is = stree.nspecies, k = 0, nused = 0; is < stree.nnode; is++) {
         if (stree.nodes[is].fossil) {
            stree.nodes[is].usefossil = 1 - it % 2;
            nused += stree.nodes[is].usefossil;
            U[k++] = !stree.nodes[is].usefossil;
            it /= 2;
         }
      }
      if (nused < nMinCorrect) continue;
      if (k != nf) error2("k != nf in getPfossilerr()?");

      y = nused*ln1pE + (nf - nused)*lnpE - pEconst - data.CcomFossilErr[icom++] + lnptC() + lnptNCgiventC();

      if (y < scaleF + 100) {
         pt += y = exp(y - scaleF);
         for (k = 0; k < nf; k++)
            if (U[k]) pU[k] += y;
      }
      else {         /* change scaleF */
         scaleF = y;
         pt = 1;
         for (k = 0; k < nf; k++) pU[k] = U[k];  /*  1 if U[k]=1 or 0 if U[k]=0 */
      }
   }
   for (k = 0; k < nf; k++) pU[k] /= pt;
   for (k = 0; k < nf; k++)
      postEFossil[k] = (postEFossil[k] * (nround - 1) + pU[k]) / nround;

   return (0);
}



int LabelOldCondP(int spnode)
{
   /* This sets com.oldconP[j]=0 if node j in the gene tree needs updating, after
      either rates or times have changed for spnode in the species tree.  This is to
      avoid duplicated computation of conditional probabilities.  The routine workes
      on the current gene tree and accounts for the fact that some species may be
      missing at some loci.
      The routine first finds spnode or its first ancestor that is present in the
      gene tree, then identifies that node in the genetree.  This reveals the
      oldest node j in the gene tree that is descendent of spnode.  Node j and all
      its ancestors in the gene tree need updating.

      Before calling this routine, set com.oldconP[]=1.  This routine changes some
      com.oldconP[] into 0 but do not change any 0 to 1.

      The gene tree is in nodes[], as UseLocus has been called prior to this.
      This is called by UpdateTimes and UpdateRates.
   */
   int i, j = spnode;

   if (j >= tree.nnode || j != nodes[j].ipop) {

      /* From among spnode and its ancestors in the species tree, find the
         first node, i, that is in genetree.
      */
      for (i = spnode; i != -1; i = stree.nodes[i].father) {

         /* Find that node in genetree that is node i in species tree.
            Its descendent, node j, is the oldest node in gene tree that is
            descendent of spnode.
         */
         for (j = 0; j < tree.nnode && nodes[j].ipop != i; j++);

         if (j < tree.nnode) break;
      }
   }

   if (j < tree.nnode)
      for (; ; j = nodes[j].father) {
         com.oldconP[j] = 0;
         if (j == tree.root) break;
      }

   return(0);
}


int UpdateTimes(double *lnL, double steplength[], char accept[])
{
   /* This updates the node ages in the master species tree stree.nodes[].age,
      stree.duplication:
      stree.nodes[i].label = 0:  usual node, which may and may not have calibration.
                            -1:  node i is the driver node, whose age is shared by other nodes.
                                 It may and may not have calibration.
      stree.nodes[j].label = i:  node j shares the age of node i, and doesn't have calibration.
   */
   int  locus, is, i, j, *sons, dad, nmirror, *mirrors;
   double t, tnew, tmin, tmax, yb[2], y, ynew, e;
   double lnacceptance, lnLd, lnpDinew[NGENE], lnpTnew, lnpRnew = -1;

   if (stree.duplication) {
      mirrors = (int*)malloc(stree.nspecies * sizeof(int));
      if (mirrors == NULL) error2("oom mirrors");
   }

   for (is = stree.nspecies; is < stree.nnode; is++) {
      if (debug == 2) printf("\nnode %2d age %9.6f", is + 1, stree.nodes[is].age);
      if (stree.nodes[is].label > 0) continue;  /* mirror node, age specified by driver. */
      nmirror = 0;
      if (stree.nodes[is].label == -1) {        /* driver node, with age shared by other nodes */
         for (j = is + 1; j < stree.nnode; j++) {
            if (stree.nodes[j].label == is)
               mirrors[nmirror++] = j;
         }
         if (debug == 2) {
            printf(" mirrors for node %2d: ", is + 1);
            for (j = 0; j < nmirror; j++)
               printf("%2d", mirrors[j] + 1);
            printf("\n");
         }
      }

      t = stree.nodes[is].age;
      sons = stree.nodes[is].sons;
      dad = stree.nodes[is].father;
      tmin = max2(stree.nodes[sons[0]].age, stree.nodes[sons[1]].age);
      tmax = (is == stree.root ? 1e9 : stree.nodes[dad].age);
      for (j = 0; j < nmirror; j++) { /* is is driver node */
         sons = stree.nodes[mirrors[j]].sons;
         dad = stree.nodes[mirrors[j]].father;
         tmin = max2(tmin, stree.nodes[sons[0]].age);
         tmin = max2(tmin, stree.nodes[sons[1]].age);
         tmax = min2(tmax, stree.nodes[dad].age);
      }
      yb[0] = (tmin > 0 ? log(tmin) : -99);
      yb[1] = log(tmax);

      y = log(t);
      e = steplength[is - stree.nspecies];
      if (e <= 0) error2("steplength = 0 in UpdateTimes");

      ynew = y + e*rndSymmetrical();
      ynew = reflect(ynew, yb[0], yb[1]);
      stree.nodes[is].age = tnew = exp(ynew);

      if (debug == 2) printf("\nnode %2d age %9.6f %9.6f", is + 1, t, tnew);

      for (j = 0; j < nmirror; j++)  /* copy new age to mirror nodes */
         stree.nodes[mirrors[j]].age = tnew;

      lnacceptance = ynew - y;
      lnpTnew = lnpriorTimes();
      lnacceptance += lnpTnew - data.lnpT;
      if (com.clock == 3) {
         lnpRnew = lnpriorRates();
         lnacceptance += lnpRnew - data.lnpR;
      }

      for (locus = 0, lnLd = 0; locus < data.ngene; locus++) {
         UseLocus(locus, 1, mcmc.usedata, 0);

         if (mcmc.saveconP) {
            for (i = 0; i < stree.nnode; i++) com.oldconP[i] = 1;
            LabelOldCondP(is);
            for (j = 0; j < nmirror; j++)   /* update all mirror nodes. */
               LabelOldCondP(mirrors[j]);
         }
         if (com.oldconP[tree.root])
            lnpDinew[locus] = data.lnpDi[locus];
         else
            lnpDinew[locus] = lnpD_locus(locus);
         lnLd += lnpDinew[locus] - data.lnpDi[locus];
      }
      lnacceptance += lnLd;

      if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) {
         accept[is - stree.nspecies] = (char)1;
         data.lnpT = lnpTnew;
         if (com.clock == 3) data.lnpR = lnpRnew;
         for (locus = 0; locus < data.ngene; locus++)
            data.lnpDi[locus] = lnpDinew[locus];
         *lnL += lnLd;
         if (mcmc.usedata == 1) switchconPin();
         if (debug == 2) printf(" Y (lnL = %.4f)\n", *lnL);
      }
      else {
         stree.nodes[is].age = t;
         for (j = 0; j < nmirror; j++)  /* driver node, with age shared by mirror nodes */
            stree.nodes[mirrors[j]].age = t;

         if (debug == 2) printf(" N (lnL = %.4f)\n", *lnL);
         for (locus = 0; locus < data.ngene; locus++)
            AcceptRejectLocus(locus, 0);  /* reposition conP */
      }

      if (debug) printStree();
   }  /* for(is) */

   if (stree.duplication)
      free(mirrors);
   return(0);
}


#if (0)  /*  this is not used now. */

int UpdateTimesClock23(double *lnL, double steplength[], char accept[])
{
   /* This updates the node ages in the master species tree stree.nodes[].age,
      one by one.  It simultaneously changes the rates at the three branches
      around the node so that the branch lengths remain the same, and there is
      no need to update the lnL.  Sliding window on the logarithm of ages.
   */
   int  is, i, *sons, dad;
   double tb[2], t, tnew, tson[2], yb[2], y, ynew, e, rateratio[3];
   double lnacceptance, lnpTnew, lnpRnew;

   if (debug == 2) puts("\nUpdateTimesClock23 ");

   for (is = stree.nspecies; is < stree.nnode; is++) {
      t = stree.nodes[is].age;
      sons = stree.nodes[is].sons;
      tson[0] = stree.nodes[sons[0]].age;
      tson[1] = stree.nodes[sons[1]].age;
      dad = stree.nodes[is].father;
      tb[0] = max2(tson[0], tson[1]);
      yb[0] = (tb[0] > 0 ? log(tb[0]) : -99);
      if (is != stree.root)  yb[1] = log(tb[1] = stree.nodes[dad].age);
      else                   yb[1] = 99;

      y = log(t);
      e = steplength[is - stree.nspecies];
      if (e <= 0) error2("steplength = 0 in UpdateTimesClock23");
      ynew = y + e*rndSymmetrical();
      ynew = reflect(ynew, yb[0], yb[1]);
      stree.nodes[is].age = tnew = exp(ynew);
      lnacceptance = ynew - y;

      /* Thorne et al. (1998) equation 9. */
      rateratio[0] = (t - tson[0]) / (tnew - tson[0]);
      rateratio[1] = (t - tson[1]) / (tnew - tson[1]);
      for (i = 0; i < data.ngene; i++) {
         stree.nodes[sons[0]].rates[i] *= rateratio[0];
         stree.nodes[sons[1]].rates[i] *= rateratio[1];
      }
      lnacceptance += data.ngene*log(rateratio[0] * rateratio[1]);
      if (is != stree.root) {
         rateratio[2] = (t - tb[1]) / (tnew - tb[1]);
         for (i = 0; i < data.ngene; i++)
            stree.nodes[is].rates[i] *= rateratio[2];
         lnacceptance += data.ngene*log(rateratio[2]);
      }

      lnpTnew = lnpriorTimes();
      lnacceptance += lnpTnew - data.lnpT;
      lnpRnew = lnpriorRates();
      lnacceptance += lnpRnew - data.lnpR;

      if (debug == 2) printf("species %2d t %8.4f %8.4f", is, t, tnew);

      if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) {
         accept[is - stree.nspecies] = (char)1;
         data.lnpT = lnpTnew;
         data.lnpR = lnpRnew;
         if (debug == 2) printf(" Y (%4d)\n", NPMat);
      }
      else {
         stree.nodes[is].age = t;
         for (i = 0; i < data.ngene; i++) {
            stree.nodes[sons[0]].rates[i] /= rateratio[0];
            stree.nodes[sons[1]].rates[i] /= rateratio[1];
            if (is != stree.root)  stree.nodes[is].rates[i] /= rateratio[2];
         }
         if (debug == 2) printf(" N (%4d)\n", NPMat);
      }
   }  /* for(is) */
   return(0);
}

#endif

int UpdateParaRates(double *lnL, double steplength[], char accept[], double space[])
{
   /* This updates mu (mean rate) and sigma2 (variance of lnrates) for each locus.
   If (data.rgeneprior==0), they have gamma-dirichlet (dos reis et al. 2014: eq. 5).
   If (data.rgeneprior==1), locus rates mu_i and sigma2_i have conditional iid prior
         (Zhu et al. 2015 Sys Biol, p.279 equation 8).
      For the clock (clock1), this changes locus rate (mu_i) and changes the likelihood
      but there is no prior for branch rates.
      For relaxed clocks (clock2 and clock3), this changes mu and sigma2 for loci;
      it changes the rate prior but not the likelihood.
   */
   int g = data.ngene, g1 = g + (data.rgeneprior == 1);
   int locus, ip, np = 1 + (com.clock > 1), j;
   char *parastr[2] = { "mu", "sigma2" };
   double lnacceptance, lnpDinew = -1e99, lnpRnew = 0, lnLd = 0, naccept = 0;
   double yb[2] = { -99,99 }, y, ynew, pold, pnew, sumold, sumnew, e;
   double *para, *gD, a, b, bold, bnew;   /* gamma-Dirichlet */

   if (debug == 5) puts("\nUpdateParaRates (rgene & sigma2)");
   if (mcmc.saveconP) FOR(j, stree.nspecies * 2 - 1) com.oldconP[j] = 0;
   for (ip = 0; ip < np; ip++) {  /* mu (rgene) and sigma2 for each locus and for overall */
      if (ip == 0) { para = data.rgene;  gD = data.rgenepara; } /* rgene (mu) */
      else         { para = data.sigma2; gD = data.sigma2para; } /* sigma2 */

      for (locus = 0; locus < g1; locus++) {
         e = steplength[ip*g1 + locus];
         if (e <= 0) error2("steplength = 0 in UpdateParaRates");

         for (j = 0, sumold = 0; j < g; j++) sumold += para[j];
         pold = para[locus];
         y = log(pold);
         ynew = y + e*rndSymmetrical();
         ynew = reflect(ynew, yb[0], yb[1]);
         para[locus] = pnew = exp(ynew);
         sumnew = sumold + pnew - pold;
         lnacceptance = ynew - y;

         if (ip == 0 && com.clock == 1 && locus < g) { /* rgene (mu) */
            UseLocus(locus, 1, mcmc.usedata, 0);
            lnpDinew = lnpD_locus(locus);
            lnLd = lnpDinew - data.lnpDi[locus];  /* likelihood ratio */
            lnacceptance += lnLd;
         }

         if (data.rgeneprior == 0)  /* gamma-dirichlet prior (dos reis et al. 2014, eq 5) */
            lnacceptance += (gD[0] - gD[2] * g)*log(sumnew / sumold) - gD[1] / g*(sumnew - sumold) + (gD[2] - 1)*(ynew - y);
         else {                     /* conditional iid prior (Zhu et al. 2015 SB, p.279 eq. 8) */
            if (locus < g) {          /* mu_i & sigma2_i */
               a = gD[2];  b = gD[2] / para[g];
               lnacceptance += logPriorRatioGamma(pnew, pold, a, b); /* f(mu_i) */
            }
            else {                 /* mu_bar & sigma2_bar */
               a = gD[0];  b = gD[1];
               lnacceptance += logPriorRatioGamma(pnew, pold, a, b); /* f(mu_bar) */
               a = gD[2];  bnew = gD[2] / pnew; bold = gD[2] / pold;
               lnacceptance += g*a*log(bnew / bold);
               for (j = 0; j < g; j++)
                  lnacceptance -= (bnew - bold)*para[j];               /* f(mu_i) */
            }
         }

         if (debug == 5) printf("%-7s locus %2d %9.5f -> %9.5f ", parastr[ip], locus, pold, pnew);

         if (com.clock > 1) {
            lnpRnew = lnpriorRates();  /* this allows correlation of branch rates among loci */
            lnacceptance += lnpRnew - data.lnpR;
         }
         if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) {
            accept[ip*g1 + locus] = (char)1;
            if (com.clock == 1 && locus < g) {    /* mu_i but not sigma2_i exist in model. */
               *lnL += lnLd;
               data.lnpDi[locus] = lnpDinew;
               if (mcmc.usedata == 1) AcceptRejectLocus(locus, 1);
            }
            else
               data.lnpR = lnpRnew;
         }
         else {
            para[locus] = pold;
            if (ip == 0 && com.clock == 1 && mcmc.usedata == 1)
               AcceptRejectLocus(locus, 0);  /* reposition conP */
         }
      }  /* for(locus) */
   }     /* for(ip) */
   return(0);
}


#if (0)
void getSinvDetS(double space[])
{
   /* This uses the variance-correlation matrix data.correl[g*g] to constructs
      Sinv (inverse of S) and detS (|S|).  It also copies the variances into
      data.sigma2[g].  This is called every time data.correl is updated.

      What restrictions should be placed on data.correl[]???

      space[data.ngene]
   */
   int i, j, g = data.ngene;
   int debug = 0;

   for (i = 0; i < g; i++)
      data.Sinv[i*g + i] = data.sigma2[i] = data.correl[i*g + i];
   for (i = 0; i < g; i++) {
      for (j = 0; j < i; j++)
         data.Sinv[i*g + j] = data.Sinv[j*g + i]
         = data.correl[i*g + j] * sqrt(data.sigma2[i] * data.sigma2[j]);
   }

   if (debug) {
      printf("\ncorrel & S & Sinv ");
      matout(F0, data.correl, g, g);
      matout(F0, data.Sinv, g, g);
   }

   matinv(data.Sinv, g, g, space);
   data.detS = fabs(space[0]);

   if (debug) {
      matout(F0, data.Sinv, g, g);
      printf("|S| = %.6g\n", data.detS);
   }
   if (data.detS < 0)
      error2("detS < 0");
}
#endif


double lnpriorRates(void)
{
   /* This calculates the log of the prior of branch rates under the two rate-drift models:
      the independent rates (clock=2) and the geometric Brownian motion model (clock=3).

      clock=2: the algorithm cycles through the branches, and add up the log
      probabilities.
      clock=3: the root rate is mu or data.rgene[].  The algorithm cycles through
      the ancestral nodes and deals with the two daughter branches.
   */
   int i, inode, locus, dad = -1, g = data.ngene, s = stree.nspecies, sons[2];
   double lnpR = -log(2 * Pi) / 2.0*(2 * s - 2)*g, t, tA, t1, t2, Tinv[4], detT;
   double zz, r = -1, rA, r1, r2, y1, y2;
   double a, b;

   if (com.clock == 3 && data.priorrate == 1)
      error2("gamma prior for rates for clock3 not implemented yet.");
   else if (com.clock == 2 && data.priorrate == 1) {   /* clock2, gamma rate prior */
      lnpR = 0;
      for (locus = 0; locus < g; locus++) {
         a = data.rgene[locus] * data.rgene[locus] / data.sigma2[locus];
         b = data.rgene[locus] / data.sigma2[locus];
         lnpR += (a*log(b) - LnGamma(a)) * (2 * s - 2);
         for (inode = 0; inode < stree.nnode; inode++) {
            if (inode == stree.root) continue;
            r = stree.nodes[inode].rates[locus];
            lnpR += -b*r + (a - 1)*log(r);
         }
      }
   }
   else if (com.clock == 2 && data.priorrate == 0) {  /* clock2, LN rate prior */
      for (locus = 0; locus < g; locus++)
         lnpR -= log(data.sigma2[locus]) / 2.*(2 * s - 2);
      for (inode = 0; inode < stree.nnode; inode++) {
         if (inode == stree.root) continue;
         for (locus = 0; locus < g; locus++) {
            r = stree.nodes[inode].rates[locus];
            zz = log(r / data.rgene[locus]) + data.sigma2[locus] / 2;
            lnpR += -zz*zz / (2 * data.sigma2[locus]) - log(r);
         }
      }
   }
   else if (com.clock == 3 && data.priorrate == 0) {  /* clock3, LN rate prior */
      for (inode = 0; inode < stree.nnode; inode++) {
         if (stree.nodes[inode].nson == 0) continue; /* skip the tips */
         dad = stree.nodes[inode].father;
         for (i = 0; i < 2; i++) sons[i] = stree.nodes[inode].sons[i];
         t = stree.nodes[inode].age;
         if (inode == stree.root) tA = 0;
         else                   tA = (stree.nodes[dad].age - t) / 2;
         t1 = (t - stree.nodes[sons[0]].age) / 2;
         t2 = (t - stree.nodes[sons[1]].age) / 2;
         detT = t1*t2 + tA*(t1 + t2);
         Tinv[0] = (tA + t2) / detT;
         Tinv[1] = Tinv[2] = -tA / detT;
         Tinv[3] = (tA + t1) / detT;
         for (locus = 0; locus < g; locus++) {
            rA = (inode == stree.root ? data.rgene[locus] : stree.nodes[inode].rates[locus]);
            r1 = stree.nodes[sons[0]].rates[locus];
            r2 = stree.nodes[sons[1]].rates[locus];
            y1 = log(r1 / rA) + (tA + t1)*data.sigma2[locus] / 2;
            y2 = log(r2 / rA) + (tA + t2)*data.sigma2[locus] / 2;
            zz = (y1*y1*Tinv[0] + 2 * y1*y2*Tinv[1] + y2*y2*Tinv[3]);
            lnpR -= zz / (2 * data.sigma2[locus]) + log(detT*square(data.sigma2[locus])) / 2 + log(r1*r2);
         }
      }
   }
   return lnpR;
}

double lnpriorRatioRates(int locus, int inodeChanged, double rold)
{
   /* This calculates the lnpriorRatio when one rate is changed.
      If inodeChanged is tip, we sum over 1 term (equation 7 in RY2007).
      If inodeChanged is not tip, we sum over 2 terms.
   */
   double rnew = stree.nodes[inodeChanged].rates[locus], lnpRd = 0, a, b, z, znew;
   double zz, r = -1, rA, r1, r2, y1, y2, t, tA, t1, t2, Tinv[4], detT;
   int i, inode, ir, dad = -1, sons[2], OldNew;

   if (com.clock == 2 && data.priorrate == 0) {         /* clock2, LN rate prior */
      z = log(rold / data.rgene[locus]) + data.sigma2[locus] / 2;;
      znew = log(rnew / data.rgene[locus]) + data.sigma2[locus] / 2;
      lnpRd = -log(rnew / rold) - (znew*znew - z*z) / (2 * data.sigma2[locus]);
   }
   else if (com.clock == 2 && data.priorrate == 1) {   /* clock2, gamma rate prior */
      a = data.rgene[locus] * data.rgene[locus] / data.sigma2[locus];
      b = data.rgene[locus] / data.sigma2[locus];
      lnpRd = -b*(rnew - rold) + (a - 1)*log(rnew / rold);
   }
   else {                                          /* clock3 */
      for (ir = 0; ir < (stree.nodes[inodeChanged].nson == 0 ? 1 : 2); ir++) {
         inode = (ir == 0 ? stree.nodes[inodeChanged].father : inodeChanged);
         dad = stree.nodes[inode].father;
         for (i = 0; i < 2; i++) sons[i] = stree.nodes[inode].sons[i];
         t = stree.nodes[inode].age;
         if (inode == stree.root) tA = 0;
         else                   tA = (stree.nodes[dad].age - t) / 2;
         t1 = (t - stree.nodes[sons[0]].age) / 2;
         t2 = (t - stree.nodes[sons[1]].age) / 2;
         detT = t1*t2 + tA*(t1 + t2);
         Tinv[0] = (tA + t2) / detT;
         Tinv[1] = Tinv[2] = -tA / detT;
         Tinv[3] = (tA + t1) / detT;
         for (OldNew = 0; OldNew < 2; OldNew++) {  /* old rate & new rate */
            stree.nodes[inodeChanged].rates[locus] = (OldNew == 0 ? rold : rnew);
            rA = (inode == stree.root ? data.rgene[locus] : stree.nodes[inode].rates[locus]);
            r1 = stree.nodes[sons[0]].rates[locus];
            r2 = stree.nodes[sons[1]].rates[locus];
            y1 = log(r1 / rA) + (tA + t1)*data.sigma2[locus] / 2;
            y2 = log(r2 / rA) + (tA + t2)*data.sigma2[locus] / 2;
            zz = (y1*y1*Tinv[0] + 2 * y1*y2*Tinv[1] + y2*y2*Tinv[3]);
            zz = zz / (2 * data.sigma2[locus]) + log(r1*r2);
            lnpRd -= (OldNew == 0 ? -1 : 1) * zz;
         }
      }
   }
   return(lnpRd);
}


int UpdateRates(double *lnL, double steplength[], char accept[])
{
   /* This updates rates under the rate-drift models (clock=2 or 3).
      It cycles through nodes in the species tree at each locus.

      Slight waste of computation: For every proposal to change rates, lnpriorRates()
      is called to calculate the prior for all rates at all loci on the whole
      tree.  This wasting computation if rates are not correlated across loci.
   */
   int locus, inode, j, g = data.ngene;
   double naccept = 0, lnpDinew, lnacceptance, lnLd, r, rnew, lnpRd, e;
   double yb[2] = { -99,99 }, y, ynew;

   for (locus = 0; locus < g; locus++) {
      for (inode = 0; inode < stree.nnode; inode++) {
         if (inode == stree.root) continue;
         r = stree.nodes[inode].rates[locus];
         y = log(r);
         e = steplength[locus*(2 * stree.nspecies - 2) + inode - (inode >= stree.root)];
         if (e <= 0) error2("steplength = 0 in UpdateRates");
         ynew = y + e*rndSymmetrical();
         ynew = reflect(ynew, yb[0], yb[1]);
         stree.nodes[inode].rates[locus] = rnew = exp(ynew);
         lnacceptance = ynew - y;              /* proposal ratio */

         UseLocus(locus, 1, mcmc.usedata, 0);  /* copyconP=1 */

         if (mcmc.saveconP) {
            FOR(j, stree.nspecies * 2 - 1) com.oldconP[j] = 1;
            LabelOldCondP(inode);
         }
         lnpRd = lnpriorRatioRates(locus, inode, r);
         lnacceptance += lnpRd;
         lnpDinew = lnpD_locus(locus);
         lnLd = lnpDinew - data.lnpDi[locus];
         lnacceptance += lnLd;
         if (lnacceptance > 0 || rndu() < exp(lnacceptance)) {
            accept[locus*(2 * stree.nspecies - 2) + inode - (inode >= stree.root)] = (char)1;
            if (mcmc.usedata == 1) AcceptRejectLocus(locus, 1);
            data.lnpR += lnpRd;
            data.lnpDi[locus] = lnpDinew;
            *lnL += lnLd;
         }
         else {
            if (mcmc.usedata == 1) AcceptRejectLocus(locus, 0);
            stree.nodes[inode].rates[locus] = r;
         }
      }
   }
   return(0);
}


int UpdateParameters(double *lnL, double steplength[], char accept[])
{
   /* This updates parameters in the substitution model such as kappa and alpha for loci.
      Should we update the birth-death process parameters here as well?
   */
   int locus, j, ip, np = !com.fix_kappa + !com.fix_alpha;
   double naccept = 0, lnLd, lnpDinew, lnacceptance, c = 1;
   double yb[2] = { -99,99 }, y, ynew, pold, pnew, e, *gammaprior;

   if (np == 0) return(0);
   if (debug == 4) puts("\nUpdateParameters");
   if (mcmc.saveconP) FOR(j, stree.nspecies * 2 - 1) com.oldconP[j] = 0;
   for (locus = 0; locus < data.ngene; locus++) {
      for (ip = 0; ip < np; ip++) {
         e = steplength[locus*np + ip];
         if (e <= 0) error2("steplength = 0 in UpdateParameters");

         if (ip == 0 && !com.fix_kappa) {  /* kappa */
            pold = data.kappa[locus];
            y = log(pold);
            ynew = y + e*rndSymmetrical();
            ynew = reflect(ynew, yb[0], yb[1]);
            data.kappa[locus] = pnew = exp(ynew);
            gammaprior = data.kappagamma;
         }
         else {  /* alpha */
            pold = data.alpha[locus];
            y = log(pold);
            ynew = y + e*rndSymmetrical();
            ynew = reflect(ynew, yb[0], yb[1]);
            data.alpha[locus] = pnew = exp(ynew);
            gammaprior = data.alphagamma;
         }
         lnacceptance = ynew - y;

         UseLocus(locus, 1, mcmc.usedata, 0); /* this copies parameter from data.[] to com. */

         lnpDinew = lnpD_locus(locus);
         lnLd = lnpDinew - data.lnpDi[locus];
         lnacceptance += lnLd;
         lnacceptance += logPriorRatioGamma(pnew, pold, gammaprior[0], gammaprior[1]);

         if (debug == 4)
            printf("\nLocus %2d para%d %9.4f%9.4f %10.5f", locus + 1, ip + 1, pold, pnew, lnLd);

         if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) {
            accept[locus*np + ip] = (char)1;
            *lnL += lnLd;
            data.lnpDi[locus] = lnpDinew;
            if (mcmc.usedata == 1) AcceptRejectLocus(locus, 1);
            if (debug == 4) printf(" Y\n");
         }
         else {
            if (ip == 0 && !com.fix_kappa)
               data.kappa[locus] = pold;
            else
               data.alpha[locus] = pold;

            if (mcmc.usedata == 1) AcceptRejectLocus(locus, 0);

            if (debug == 4) printf(" N\n");
         }
      }
   }
   return(0);
}



int mixingTipDate(double *lnL, double steplength, char *accept)
{
   /* This increases or decreases all node ages in a 1-D move, applied to TipDated data.
      1 September 2011, Aguas de Lindoia
      Znote 2014-11-21:  Please check that changemu=0 is correct and gives the same result as changemu=1.
   */
   int s=stree.nspecies, changemu = 0, g = data.ngene, i, j, k, jj, locus, ndivide;
   double c, lnc, lnacceptance = 0, lnpTnew, lnpRnew = -9999, lnpDinew[NGENE], lnLd;
   double *gD = data.rgenepara, summuold = 0, summunew;
   double *minages, *x, tz;

   minages = (double*)malloc(s * 2 * sizeof(double));
   if (minages == NULL) error2("oom");
   memset(minages, 0, s * 2 * sizeof(double));
   x = minages + s;
   x[0] = 1;

   if (com.clock == 1) changemu = 1;
   ndivide = (changemu ? g : 0);

   /* find minages[] and x[] for the interior nodes */
   if (debug == 6)
      printStree();

   for (j = 0; j < s; j++) {
      tz = stree.nodes[j].age;
      for (k = stree.nodes[j].father; k != -1; k = stree.nodes[k].father)
         if (tz < minages[k - s])
            break;
         else
            minages[k - s] = tz;
   }
   if (debug == 6) {
      for (j = s; j < s*s - 1; j++)
         printf("node %2d age %9.5f minages %9.5f\n", j + 1, stree.nodes[j].age, minages[j - s]);
   }
   for (j = s; j < s*2 - 1; j++) {
      if (j == stree.root) continue;
      jj = j - s;
      x[jj] = (stree.nodes[j].age - minages[jj]) / (stree.nodes[stree.nodes[j].father].age - minages[jj]);
   }
   lnacceptance = lnc = steplength*rndSymmetrical();
   c = exp(lnc);
   stree.nodes[stree.root].age = minages[0] + (stree.nodes[stree.root].age - minages[0]) * c;
   for (j = s; j < s*2-1; j++) {
      if (j == stree.root) continue;
      jj = j - s;
      tz = stree.nodes[j].age;
      stree.nodes[j].age = minages[jj] + x[jj] * (stree.nodes[stree.nodes[j].father].age - minages[jj]);
      lnacceptance += log((stree.nodes[j].age - minages[jj]) / (tz - minages[jj]));
   }

   lnpTnew = lnpriorTimes();
   lnacceptance += lnpTnew - data.lnpT;

   /* dividing locus rates by c.  */
   if (changemu) {
      for (j = 0; j < g; j++) {  /* mu */
         summuold += data.rgene[j];
         data.rgene[j] /= c;
      }
      summunew = summuold / c;
      lnacceptance += (gD[0] - gD[2] * g)*log(summunew / summuold) - gD[1] / g*(summunew - summuold) + (gD[2] - 1)*g*(-lnc);
   }
   if (com.clock > 1) {              /* rate-drift models */
      ndivide += g*(s * 2 - 2);
      for (i = 0; i < stree.nnode; i++)
         if (i != stree.root)
            for (j = 0; j < g; j++)
               stree.nodes[i].rates[j] /= c;
      lnpRnew = lnpriorRates();
      lnacceptance += lnpRnew - data.lnpR;
   }
   lnacceptance -= ndivide*lnc;

   if (mcmc.saveconP)
      for (j = 0; j < s * 2 - 1; j++)
         com.oldconP[j] = 0;
   for (locus = 0, lnLd = 0; locus < g; locus++) {
      UseLocus(locus, 1, mcmc.usedata, 0);
      lnpDinew[locus] = lnpD_locus(locus);
      lnLd += lnpDinew[locus] - data.lnpDi[locus];
   }
   lnacceptance += lnLd;

   if (lnacceptance > 0 || rndu() < exp(lnacceptance)) { /* accept */
      *accept = (char)1;
      data.lnpT = lnpTnew;
      data.lnpR = lnpRnew;
      for (locus = 0; locus < g; locus++)
         data.lnpDi[locus] = lnpDinew[locus];
      if (mcmc.usedata == 1) switchconPin();
      *lnL += lnLd;
   }
   else {   /* reject */
      /* restore all node ages */
      stree.nodes[stree.root].age = minages[0] + (stree.nodes[stree.root].age - minages[0]) / c;
      for (j = s; j < s * 2 - 1; j++) {
         if (j == stree.root) continue;
         jj = j - s;
         stree.nodes[j].age = minages[jj] + x[jj] * (stree.nodes[stree.nodes[j].father].age - minages[jj]);
      }
      if (debug == 6)
         printStree();
      /* restore all rates */
      if (changemu) for (j = 0; j < g; j++) data.rgene[j] *= c;
      if (com.clock > 1) {
         for (i = 0; i < stree.nnode; i++)
            if (i != stree.root)
               for (j = 0; j < g; j++)
                  stree.nodes[i].rates[j] *= c;
      }

      for (locus = 0; locus < g; locus++)
         AcceptRejectLocus(locus, 0);  /* reposition conP */

   }
   free(minages);
   return(0);
}


int mixing(double *lnL, double steplength, char *accept)
{
   /* If com.clock==1, this multiplies all times by c and divides all locus rates
      (mu or data.rgene[]) by c, so that the likelihood is unchanged.

      If com.clock=2 or 3, this multiplies all times by c and divide all branch rates
      by c.  If (changemu), the locus rates mu (data.rgene[]) as well as the branch rates
      stree.nodes[].rates at all loci are divided by c.
      If(changemu==0), the locus rates mu_i are not changed and only the branch rates
      are divided by c.  In both cases, the likelihood is unchanged.
      In the case of cond iid prior (data.rgeneprior==1), mu_bar is not changed.
   */
   int  changemu = 1, g = data.ngene, i, j, ndivide, s = stree.nspecies;
   double *gD = data.rgenepara, a, b, c, lnc, summuold = -1, summunew = -1;
   double lnacceptance = 0, lnpTnew, lnpRnew = -1;

   if (com.clock == 1) changemu = 1;
   ndivide = (changemu ? g : 0);

   if (steplength <= 0) error2("steplength = 0 in mixing");
   if (com.TipDate)
      return mixingTipDate(lnL, steplength, accept);
   /*
      else if(com.clock==2 || com.clock==3)
         return mixingCladeStretch(lnL, steplength, accept);
   */
   lnc = steplength*rndSymmetrical();
   c = exp(lnc);
   for (j = s; j < stree.nnode; j++)
      stree.nodes[j].age *= c;

   if (changemu) {
      for (j = 0; j < g; j++) data.rgene[j] /= c;
      if (data.rgeneprior == 0) {    /* gamma-dirichlet: eq 5 in dos reis et al. (2014) */
         for (j = 0, summunew = 0; j < g; j++) summunew += data.rgene[j];
         summuold *= c;
         lnacceptance += (gD[0] - gD[2] * g)*log(summunew / summuold) - gD[1] / g*(summunew - summuold) + (gD[2] - 1)*g*(-lnc);
      }
      else {  /* conditional iid prior */
         a = gD[2];  b = gD[2] / data.rgene[g];  /* data.rgene[g] is mu_bar. */
         for (j = 0; j < g; j++)     /* f(mu_i) */
            lnacceptance += -(a - 1)*lnc - b*data.rgene[j] * (1 - c);
      }
   }
   if (com.clock > 1) {            /* rate-drift models */
      ndivide += g*(s * 2 - 2);
      for (i = 0; i < stree.nnode; i++)
         if (i != stree.root)
            for (j = 0; j < g; j++)
               stree.nodes[i].rates[j] /= c;
      lnpRnew = lnpriorRates();
      lnacceptance += lnpRnew - data.lnpR;
   }

   lnpTnew = lnpriorTimes();
   lnacceptance += lnpTnew - data.lnpT + (s - 1 - stree.duplication - ndivide)*lnc;

   if (lnacceptance > 0 || rndu() < exp(lnacceptance)) { /* accept */
      *accept = (char)1;
      data.lnpT = lnpTnew;
      data.lnpR = lnpRnew;
   }
   else {   /* reject */
      for (j = s; j < stree.nnode; j++)
         stree.nodes[j].age /= c;
      if (changemu)
         for (j = 0; j < g; j++) data.rgene[j] *= c;
      if (com.clock > 1) {
         for (i = 0; i < stree.nnode; i++)
            if (i != stree.root)
               for (j = 0; j < g; j++)
                  stree.nodes[i].rates[j] *= c;
      }
   }
   return(0);
}

int mixingCladeStretch(double *lnL, double steplength, char *accept)
{
   /* If com.clock==1, this multiplies all times by c and divides all rates (mu or
      data.rgene[]) by c, so that the likelihood does not change.

      If com.clock=2 or 3, this multiplies all times by c and divide by c all rates:
      stree.nodes[].rates and mu (data.rgene[]) at all loci.  The likelihood is
      not changed.
   */
   int  s = stree.nspecies, g = data.ngene, is, locus, i, j, k, nsnode, desc[NS], *sons;
   int  ndivide = 0;
   double *gD = data.rgenepara, c, lnc;
   double lnacceptance = 0, lnpTnew, lnpRnew = -1, lnpDinew[NGENE], lnLd;

   printf("mixingCladeStretch(): this move is incorrect.  Think and fix it?");

   *accept = (char)0;
   is = s + (int)((s - 1)*rndu());
   lnc = steplength*rndSymmetrical();
   c = exp(lnc);
   for (j = s, nsnode = 0; j < s * 2 - 1; j++) { /* collect descendent of is into desc[] */
      for (k = j; k != -1; k = stree.nodes[k].father)
         if (k == is) break;
      if (k != -1) desc[nsnode++] = j;  /* j is descendent of is */
   }
   for (j = 0; j < nsnode; j++) {
      stree.nodes[desc[j]].age *= c;
      sons = stree.nodes[desc[j]].sons;
      for (locus = 0; locus < g; locus++) {
         stree.nodes[sons[0]].rates[locus] /= c;
         stree.nodes[sons[1]].rates[locus] /= c;
      }
   }
   for (locus = 0, lnLd = 0; locus < g; locus++) {
      UseLocus(locus, 1, mcmc.usedata, 0);
      lnpDinew[locus] = lnpD_locus(locus);
      lnLd += lnpDinew[locus] - data.lnpDi[locus];
      lnpRnew = lnpriorRates();
      lnacceptance += lnpRnew - data.lnpR;
   }
   lnacceptance += lnLd;
   lnpTnew = lnpriorTimes();
   lnacceptance += lnpTnew - data.lnpT + (nsnode - (g*nsnode * 2))*lnc;

   if (lnacceptance > 0 || rndu() < exp(lnacceptance)) { /* accept */
      *accept = (char)1;
      data.lnpT = lnpTnew;
      data.lnpR = lnpRnew;
      for (locus = 0; locus < g; locus++)
         data.lnpDi[locus] = lnpDinew[locus];
      if (mcmc.usedata == 1) switchconPin();
      *lnL += lnLd;
   }
   else {   /* reject */
      for (j = 0; j < nsnode; j++) {
         stree.nodes[desc[j]].age /= c;
         sons = stree.nodes[desc[j]].sons;
         for (locus = 0; locus < g; locus++) {
            stree.nodes[sons[0]].rates[locus] *= c;
            stree.nodes[sons[1]].rates[locus] *= c;
         }
      }
   }

   if (debug)
      printf("\nmixingCladeStretch %d (%d) %9.2f %4c ", is + 1, nsnode, lnLd, (*accept ? 'Y' : 'N'));

   return(0);
}

int UpdatePFossilErrors(double steplength, char *accept)
{
   /* This updates the probability of fossil errors stree.Pfossilerr.
      The proposal do not change the likelihood.
   */
   double lnacceptance, lnpTnew, pold, pnew;
   double p = data.pfossilerror[0], q = data.pfossilerror[1];

   if (steplength <= 0) error2("steplength = 0 in UpdatePFossilErrors");
   pold = data.Pfossilerr;
   pnew = pold + steplength*rndSymmetrical();
   data.Pfossilerr = pnew = reflect(pnew, 0, 1);
   lnacceptance = (p - 1)*log(pnew / pold) + (q - 1)*log((1 - pnew) / (1 - pold));
   lnpTnew = lnpriorTimes();
   lnacceptance += lnpTnew - data.lnpT;

   if (lnacceptance >= 0 || rndu() < exp(lnacceptance)) {
      *accept = (char)1;
      data.lnpT = lnpTnew;
   }
   else {
      data.Pfossilerr = pold;
   }
   return(0);
}


int DescriptiveStatisticsSimpleMCMCTREE(FILE *fout, char infile[], int nbin)
{
   FILE *fin = gfopen(infile, "r"), *fFigTree;
   char timestr[32], FigTreef[96] = "FigTree.tre";
   double *dat, *x, *mean, *median, *minx, *maxx, *x005, *x995, *x025, *x975, *xHPD025, *xHPD975, *var;
   double *y, *Tint, tmp[2], rho1;
   int  n, p, i, j, jj;
   size_t sizedata = 0;
   static int lline = 10000000, ifields[MAXNFIELDS], SkipC1 = 1, HasHeader;
   char *line;
   static char varstr[MAXNFIELDS][32] = { "" };

   puts("\nSummarizing MCMC samples . ..");
   if ((line = (char*)malloc(lline * sizeof(char))) == NULL) error2("oom ds");
   scanfile(fin, &n, &p, &HasHeader, line, ifields);
   printf("%d records, %d variables\n", n, p);

   sizedata = (size_t)p*n * sizeof(double);
   if (fabs((double)p*n * sizeof(double) - sizedata) > 1)
      error2("data matrix too big to fit in memory");
   dat = (double*)malloc(sizedata);
   mean = (double*)malloc((p * 13 + n) * sizeof(double));
   if (dat == NULL || mean == NULL) error2("oom DescriptiveStatistics.");
   memset(dat, 0, sizedata);
   memset(mean, 0, (p * 12 + n) * sizeof(double));
   median = mean + p; minx = median + p; maxx = minx + p;
   x005 = maxx + p; x995 = x005 + p; x025 = x995 + p; x975 = x025 + p; xHPD025 = x975 + p; xHPD975 = xHPD025 + p;
   var = xHPD975 + p;   Tint = var + p;  y = Tint + p;

   if (HasHeader) { /* line has the first line of variable names */
      for (i = 0; i < p; i++) {
         if (ifields[i]<0 || ifields[i]>lline - 2) error2("line not long enough?");
         sscanf(line + ifields[i], "%s", varstr[i]);
      }
   }
   for (i = 0; i < n; i++)
      for (j = 0; j < p; j++) {
         fscanf(fin, "%lf", &dat[j*n + i]);
      }
   fclose(fin);

   printf("Collecting mean, median, min, max, percentiles, etc.\n");
   for (j = SkipC1, x = dat + j*n; j < p; j++, x += n) {
      memmove(y, x, n * sizeof(double));
      Tint[j] = 1 / Eff_IntegratedCorrelationTime(y, n, &mean[j], &var[j], &rho1);
      qsort(x, (size_t)n, sizeof(double), comparedouble);
      minx[j] = x[0];  maxx[j] = x[n - 1];
      median[j] = (n % 2 == 0 ? (x[n / 2 - 1] + x[n / 2]) / 2 : x[n / 2]);
      x005[j] = x[(int)(n*.005)];    x995[j] = x[(int)(n*.995)];
      x025[j] = x[(int)(n*.025)];    x975[j] = x[(int)(n*.975)];

      HPDinterval(x, n, tmp, 0.05);
      xHPD025[j] = tmp[0];
      xHPD975[j] = tmp[1];
      if ((j + 1) % 100 == 0 || j == p - 1)
         printf("\r\t\t\t%6d/%6d done  %s", j + 1, p, printtime(timestr));
   }
   FPN(F0);

   if (nodes[stree.root].age == 0) {  /* the ages need be reset if mcmc.print<0 */
      for (j = stree.nspecies; j < stree.nnode; j++)
         nodes[j].age = mean[SkipC1 + j - stree.nspecies];
      for (j = 0; j < stree.nnode; j++) {
         nodes[j].label = -1;  /* Label strings have the 95% CI, so labels are disabled. */
         if (j != stree.root) {
            nodes[j].branch = nodes[nodes[j].father].age - nodes[j].age;
            if (nodes[j].branch < 0)
               puts("blength<0");
         }
      }
   }
   for (j = stree.nspecies; j < stree.nnode; j++) {
      nodes[j].annotation = (char*)malloc(32 * sizeof(char));
      jj = SkipC1 + j - stree.nspecies;
      sprintf(nodes[j].annotation, "[&95%%HPD={%.6g, %.6g}]", xHPD025[jj], xHPD975[jj]);
   }
   FPN(fout);
   OutTreeN(fout, 1, PrBranch | PrLabel | PrNodeStr); FPN(fout); /* prints ages & CIs */

   /* sprintf(FigTreef, "%s.FigTree.tre", com.seqf); */
   if ((fFigTree = (FILE*)fopen(FigTreef, "w")) == NULL)
      error2("FigTree.tre file creation error");
   fprintf(fFigTree, "#NEXUS\nBEGIN TREES;\n\n\tUTREE 1 = ");
   OutTreeN(fFigTree, 1, PrBranch | PrLabel | PrNodeStr);
   fprintf(fFigTree, "\n\nEND;\n");

   if (com.TipDate) {
      printf("\nThe FigTree tree is in FigTree.tre");
      fprintf(fFigTree, "\n[Note for FigTree: Under Time Scale, set Offset = %.1f, Scale factor = -%.1f\n",
         com.TipDate, com.TipDate_TimeUnit);
      fprintf(fFigTree, "Untick Scale Bar, & tick Tip Labels, Node Bars, Scale Axis, Reverse Axis, Show Grid.]\n");
   }
   fclose(fFigTree);

   /* rategrams */
   if (com.clock >= 2 && mcmc.print >= 2) {
      jj = SkipC1 + stree.nspecies - 1 + data.ngene * 2;
      for (i = 0; i < data.ngene; i++) {
         fprintf(fout, "\nrategram locus %d:\n", i + 1);
         for (j = 0; j < stree.nnode; j++) {
            if (j == stree.root) continue;
            nodes[j].branch = mean[jj++];
         }
         OutTreeN(fout, 1, PrBranch); FPN(fout);
      }
   }

   printf("\n\nPosterior means (95%% Equal-tail CI) (95%% HPD CI) HPD-CI-width\n\n");
   for (j = SkipC1; j < p; j++) {
      printf("%-10s ", varstr[j]);
      printf("%10.4f (%7.4f, %7.4f) (%7.4f, %7.4f) %7.4f", mean[j], x025[j], x975[j], xHPD025[j], xHPD975[j], xHPD975[j] - xHPD025[j]);
      if (j < SkipC1 + stree.nspecies - 1) {
         printf("  (Jnode %2d)", 2 * stree.nspecies - 1 - 1 - j + SkipC1);
         if (com.TipDate)
            printf(" time: %7.3f (%6.3f, %6.3f)", com.TipDate - mean[j] * com.TipDate_TimeUnit,
               com.TipDate - x975[j] * com.TipDate_TimeUnit, com.TipDate - x025[j] * com.TipDate_TimeUnit);
      }
      printf("\n");
   }

   /* int correl (double x[], double y[], int n, double *mx, double *my, double *vxx, double *vxy, double *vyy, double *r) */

   if (fout) {
      fprintf(fout, "\n\nPosterior mean (95%% Equal-tail CI) (95%% HPD CI) HPD-CI-width\n\n");
      for (j = SkipC1; j < p; j++) {
         fprintf(fout, "%-10s ", varstr[j]);
         fprintf(fout, "%10.4f (%7.4f, %7.4f) (%7.4f, %7.4f) %7.4f", mean[j], x025[j], x975[j], xHPD025[j], xHPD975[j], xHPD975[j] - xHPD025[j]);
         if (j < SkipC1 + stree.nspecies - 1) {
            fprintf(fout, " (Jnode %2d)", 2 * stree.nspecies - 1 - 1 - j + SkipC1);
            if (com.TipDate)
               fprintf(fout, " time: %7.3f (%6.3f, %6.3f)", com.TipDate - mean[j] * com.TipDate_TimeUnit,
                  com.TipDate - x975[j] * com.TipDate_TimeUnit, com.TipDate - x025[j] * com.TipDate_TimeUnit);
         }
         fprintf(fout, "\n");
      }
#if(1)
      printf("\nmean"); for (j = SkipC1; j < p; j++) printf("\t%.4f", mean[j]);
      printf("\nEff "); for (j = SkipC1; j < p; j++) printf("\t%.4f", 1 / Tint[j]);
#endif
   }

   free(dat); free(mean); free(line);
   for (j = stree.nspecies; j < stree.nnode; j++)
      free(nodes[j].annotation);
   return(0);
}


int MCMC(FILE* fout)
{
   FILE *fmcmc = NULL;
   int nxpr[2] = { 6, 2 }, i, j, k, ir, g = data.ngene;
   double lnL = 0, mlnL = 0, nround = 0, *x, *mx, postEFossil[MaxNFossils] = { 0 };
   double au = data.rgenepara[0], bu = data.rgenepara[1], a = data.rgenepara[2];
   char timestr[36];

   mcmc.saveconP = 1;
   if (mcmc.usedata != 1) mcmc.saveconP = 0;
   for (j = 0; j < stree.nspecies * 2 - 1; j++)  com.oldconP[j] = 0;
   GetMem();
   if (mcmc.usedata == 2)
      ReadBlengthGH(com.inBVf);

   printf("\n%d burnin, sampled every %d, %d samples\n",
      mcmc.burnin, mcmc.sampfreq, mcmc.nsample);
   if (mcmc.usedata) puts("Approximating posterior");
   else              puts("Approximating prior");

   printf("(Settings: cleandata=%d  print=%d  saveconP=%d)\n",
      com.cleandata, mcmc.print, mcmc.saveconP);

   com.np = GetInitials();

   for (j = 0; j < mcmc.nsteplength; j++)
      mcmc.steplength[j] = 0.001 + 0.1*rndu();

   /*
   printf("\nsteplength values: ");
   matout2(F0, mcmc.steplength, 1, mcmc.nsteplength, 9, 5);
   */
   printStree();

   x = (double*)malloc(com.np * 2 * sizeof(double));
   if (x == NULL) error2("oom in mcmc");
   mx = x + com.np;

   if (mcmc.print > 0)
      fmcmc = gfopen(com.mcmcf, "w");
   collectx(fmcmc, x);
   if (!com.fix_kappa && !com.fix_alpha && data.ngene == 2) { nxpr[0] = 6; nxpr[1] = 4; }

   puts("\npriors: ");
   if (mcmc.usedata == 1) {
      if (!com.fix_kappa) printf("\tG(%7.4f, %7.4f) for kappa\n", data.kappagamma[0], data.kappagamma[1]);
      if (!com.fix_alpha) printf("\tG(%7.4f, %7.4f) for alpha\n", data.alphagamma[0], data.alphagamma[1]);
   }
   if (data.rgeneprior == 0) {
      printf("\tmu_i ~ gammaDir(%.3f, %.3f, %.3f), ", au, bu, a);
      printf("SD(mu_i) = %8.5f, corr(mu_i,mu_j) = %7.5f\n", sqrt(au / (bu*bu)*((g - 1) / (g*a + 1)*(au + 1) + 1)), (g*a - au) / (g*(a + 1) + (g - 1)*au));
      if (com.clock > 1)
         printf("\tsigma2 ~ gammaDir(%.4f, %.4f, %.4f)\n", data.sigma2para[0], data.sigma2para[1], data.sigma2para[2]);
   }
   else {
      printf("\tmu_i ~ Cond iid(%.3f, %.3f, %.3f), ", au, bu, a);
      printf("SD(mu_i) = %8.5f, corr(mu_i,mu_j) = %7.5f\n", sqrt(au / (bu*bu)*(1 + (1 + au) / a)), a / (a + 1 + au));
      if (com.clock > 1)
         printf("\tsigma2 ~ Cond iid(%.4f, %.4f, %.4f)\n", data.sigma2para[0], data.sigma2para[1], data.sigma2para[2]);
   }

   if (com.checkpoint == 2) {
      ReadMCMCstate(com.checkpointf);
      mcmc.burnin = 0;
      printf("\nInitial parameters, np = %d.\nSptree & Genetrees read from %s.\n", com.np, com.checkpointf);
   }

   /* calculates prior for times and likelihood for each locus */
   if (data.pfossilerror[0])
      SetupPriorTimesFossilErrors();
   data.lnpT = lnpriorTimes();

   for (j = 0; j < data.ngene; j++) com.rgene[j] = -1;  /* com.rgene[] is not used. */
   if (com.clock > 1)
      data.lnpR = lnpriorRates();
   printf("\nInitial parameters (np = %d):\n", com.np);
   for (j = 0; j < com.np; j++) printf(" %9.6f", x[j]); FPN(F0);
   lnL = lnpData(data.lnpDi);
   printf("\nlnL0 = %9.2f\n", lnL);

   printf("\nStarting MCMC (np = %d) . . .\n", com.np);
   if (com.clock == 1)
      printf("paras: %d times, %d mu, (& kappa, alpha)\n", stree.nspecies - 1, data.ngene);
   else
      printf("paras: %d times, %d mu, %d sigma2 (& rates, kappa, alpha)\n", stree.nspecies - 1, data.ngene, data.ngene);

   zero(mx, com.np);
   if (com.np < nxpr[0] + nxpr[1]) { nxpr[0] = com.np; nxpr[1] = 0; }
   for (ir = -mcmc.burnin; ir < mcmc.sampfreq*mcmc.nsample; ir++) {
      if (ir == 0 || (nround >= 100 && mcmc.burnin >= 200 && ir < 0 && ir % (mcmc.burnin / 5) == 0)) {
         if (mcmc.burnin >= 200)
            ResetStepLengths(NULL, mcmc.Pjump, mcmc.steplength, mcmc.nsteplength);

         nround = 0;
         zero(mcmc.Pjump, mcmc.nsteplength);
         zero(mx, com.np);
         mlnL = 0;
         testlnL = 1;
         if (fabs(lnL - lnpData(data.lnpDi)) > 0.001) {
            printf("\n%12.6f = %12.6f?  Resetting lnL\n", lnL, lnpData(data.lnpDi));
            lnL = lnpData(data.lnpDi);
            exit(-1);
         }
         testlnL = 0;
      }
      nround++;
      memset(mcmc.accept, 0, mcmc.nsteplength * sizeof(char));

      UpdateTimes(&lnL, mcmc.steplength, mcmc.accept);

      UpdateParaRates(&lnL, mcmc.steplength + mcmc.steplengthOffset[1], mcmc.accept + mcmc.steplengthOffset[1], com.space);

      if (com.clock > 1)
         UpdateRates(&lnL, mcmc.steplength + mcmc.steplengthOffset[2], mcmc.accept + mcmc.steplengthOffset[2]);

      if (mcmc.usedata == 1)
         UpdateParameters(&lnL, mcmc.steplength + mcmc.steplengthOffset[3], mcmc.accept + mcmc.steplengthOffset[3]);

      mixing(&lnL, mcmc.steplength[mcmc.steplengthOffset[4]], mcmc.accept + mcmc.steplengthOffset[4]);

      if (data.pfossilerror[0])
         UpdatePFossilErrors(mcmc.steplength[mcmc.steplengthOffset[5]], mcmc.accept + mcmc.steplengthOffset[5]);

      /* printf("root age = %.4f\n", stree.nodes[stree.root].age); */

      for (j = 0; j < mcmc.nsteplength; j++)
         mcmc.Pjump[j] = (mcmc.Pjump[j] * (nround - 1) + mcmc.accept[j]) / nround;
      if (mcmc.print) collectx(fmcmc, x);
      for (j = 0; j < com.np; j++) mx[j] = (mx[j] * (nround - 1) + x[j]) / nround;
      mlnL = (mlnL*(nround - 1) + lnL / BFbeta) / nround;

      if (data.pfossilerror[0])
         getPfossilerr(postEFossil, nround);

      if (mcmc.print && ir >= 0 && (ir == 0 || (ir + 1) % mcmc.sampfreq == 0)) {
         fprintf(fmcmc, "%d", ir + 1);
         for (j = 0; j < com.np; j++) fprintf(fmcmc, "\t%.7f", x[j]);
         if (mcmc.usedata) fprintf(fmcmc, "\t%.3f", lnL / BFbeta);
         FPN(fmcmc);
      }
      if ((ir + 1) % max2(mcmc.sampfreq, mcmc.sampfreq*mcmc.nsample / 100) == 0) {
         printf("\r%3.0f%%", (ir + 1.) / (mcmc.nsample*mcmc.sampfreq)*100.);

         for (j = 0; j < 5; j++) printf(" %4.2f", mcmc.Pjump[j]);
         printf(" ");

         FOR(j, nxpr[0]) printf(" %5.3f", mx[j]);
         if (com.np > nxpr[0] + nxpr[1] && nxpr[1]) printf(" -");
         FOR(j, nxpr[1]) printf(" %5.3f", mx[com.np - nxpr[1] + j]);
         if (mcmc.usedata) printf(" %4.1f", mlnL);
      }

      k = mcmc.sampfreq*mcmc.nsample;
      if (k > 20 && (ir+1) % (k/20) == 0) {
         testlnL = 1;
         if (fabs(lnL - lnpData(data.lnpDi)) > 0.001) {
            printf("\n%12.6f = %12.6f?  Resetting lnL\n", lnL, lnpData(data.lnpDi));
            lnL = lnpData(data.lnpDi);
            exit(-1);
         }
         testlnL = 0;
         if (mcmc.print) fflush(fmcmc);
         printf(" %s\n", printtime(timestr));
         if (com.checkpoint == 1 && ir >= 0 && (ir + 1) % (k / 10) == 0) /* save MCMC state */
            SaveMCMCstate(com.checkpointf, ir, data.lnpR, lnL);
      }
   }  /* for(ir) */

   if (mcmc.print) fclose(fmcmc);

   if (BFbeta != 1) printf("\nBFbeta = %8.6f  E_b(lnf(X)) = %9.4f\n", BFbeta, mlnL);
   printf("\nTime used: %s", printtime(timestr));
   fprintf(fout, "\nTime used: %s", printtime(timestr));

   fprintf(fout, "\n\nmean of parameters using all iterations\n");
   for (i = 0; i < com.np; i++)  fprintf(fout, " %9.5f", mx[i]);
   FPN(fout);

   copySptree();
   for (j = 0; j < stree.nspecies - 1; j++)
      nodes[stree.nspecies + j].age = mx[j];
   for (j = 0; j < stree.nnode; j++)
      if (j != stree.root)
         nodes[j].branch = nodes[nodes[j].father].age - nodes[j].age;
   fputs("\nSpecies tree for FigTree.  Branch lengths = posterior mean times; 95% CIs = labels", fout);
   FPN(fout); OutTreeN(fout, 1, PrNodeNum); FPN(fout);
   FPN(fout); OutTreeN(fout, 1, 1); FPN(fout);

   if (mcmc.print)
      DescriptiveStatisticsSimpleMCMCTREE(fout, com.mcmcf, 1);
   if (data.pfossilerror[0]) {
      free(data.CcomFossilErr);
      printf("\nPosterior probabilities that each fossil is in error.\n");
      for (i = k = 0; i < stree.nspecies * 2 - 1; i++) {
         if ((j = stree.nodes[i].fossil))
            printf("Node %2d: %s (%9.4f, %9.4f) %8.4f\n", i + 1, fossils[j],
               stree.nodes[i].pfossil[0], stree.nodes[i].pfossil[1], postEFossil[k++]);
      }
      fprintf(fout, "\nPosterior probabilities that each fossil is in error.\n");
      for (i = k = 0; i < stree.nspecies * 2 - 1; i++) {
         if ((j = stree.nodes[i].fossil))
            fprintf(fout, "Node %2d: %s (%9.4f, %9.4f) %8.4f\n", i + 1, fossils[j],
               stree.nodes[i].pfossil[0], stree.nodes[i].pfossil[1], postEFossil[k++]);
      }
   }

   free(x);
   FreeMem();
   printf("\ntime prior: %s\n", (data.priortime ? "beta" : "Birth-Death-Sampling"));
   printf("rate prior: %s\n", (data.priorrate ? "gamma" : "Log-Normal"));
   if (mcmc.usedata == 2 && data.transform)
      printf("%s transform is used in approx like calculation.\n", Btransform[data.transform]);
   printf("\nTime used: %s\n", printtime(timestr));
   return(0);
}
