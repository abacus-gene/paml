/* baseml.c
   Maximum likelihood parameter estimation e aligned DNA (RNA) sequences,
   combined with phylogenetic tree estimation.
   Copyright, Ziheng YANG, July 1992 onwards

   cc -o baseml -fast baseml.c tools.c -lm
   cl -O2 baseml.c tools.c
   baseml <ControlFileName>
*/

#include "paml.h"

#define NS            7000
#define NBRANCH       (NS*2-2)
#define NNODE         (NS*2-1)
#define MAXNSONS      64
#define NGENE         500
#define LSPNAME       96
#define NCODE         5
#define NCATG         100

#define spaceming2(n) ((n)*((n)*2+9+2)*sizeof(double))

#define NP            (NBRANCH+NGENE+11)
/*
#define NP            (NBRANCH+3*NNODE+NGENE+11+NCATG*NCATG-1)
*/
extern int noisy, NFunCall, NEigenQ, NPMatUVRoot, *ancestor, GeneticCode[][64];
extern double Small_Diff, *SeqDistance;

int Forestry(FILE *fout, FILE* ftree);
void DetailOutput(FILE *fout, double x[], double var[]);
int GetOptions(char *ctlf);
int GetInitials(double x[], int *fromfile);
int GetInitialsTimes(double x[]);
int SetxInitials(int np, double x[], double xb[][2]);
int SetParameters(double x[]);
int SetPGene(int igene, int _pi, int _UVRoot, int _alpha, double x[]);
int SetPSiteClass(int iclass, double x[]);
int testx(double x[], int np);
double GetBranchRate(int igene, int ibrate, double x[], int *ix);
int GetPMatBranch(double Pt[], double x[], double t, int inode);
int ConditionalPNode(int inode, int igene, double x[]);

int TransformxBack(double x[]);
int AdHocRateSmoothing(FILE*fout, double x[NS * 3], double xb[NS * 3][2], double space[]);
void DatingHeteroData(FILE* fout);
void InitializeNodeScale(void);

int TestModel(FILE *fout, double x[], int nsep, double space[]);
int OldDistributions(int inode, double AncientFreqs[]);
int SubData(int argc, char *argv[]);
int GroupDistances();


struct CommonInfo {
   char *z[NS];
   char *spname[NS], seqf[2048], outf[2048], treef[2048], cleandata;
   char oldconP[NNODE];  /* update conP for nodes to save computation (0 yes; 1 no) */
   int seqtype, ns, ls, ngene, posG[NGENE + 1], lgene[NGENE], *pose, npatt, readpattern;
   int np, ntime, nrgene, nrate, nalpha, npi, nhomo, ncatG, ncode, Mgene;
   size_t sspace, sconP;
   int fix_kappa, fix_rgene, fix_alpha, fix_rho, nparK, fix_blength;
   int clock, model, getSE, runmode, print, verbose, ndata, idata, bootstrap;
   int icode, coding, method, nbtype;
   int ndata_trees_opt;  /* 0 1 2 3 */
   double *fpatt, kappa, alpha, rho, rgene[NGENE], pi[4], piG[NGENE][4];
   double freqK[NCATG], rK[NCATG], MK[NCATG*NCATG];
   double *blengths0;
   double(*plfun)(double x[], int np), *conP, *fhK, *space;
   int conPSiteClass;        /* is conP memory allocated for each class? */
   int NnodeScale;
   char *nodeScale;          /* nScale[ns-1] for interior nodes */
   double *nodeScaleF;       /* nScaleF[npatt] for scale factors */
   int fix_omega;
   double omega;
}  com;

struct TREEB {
   int  nbranch, nnode, root, branches[NBRANCH][2];
   double lnL;
}  tree;

struct TREEN {
   int father, nson, sons[MAXNSONS], ibranch, ipop;
   double branch, age, *pkappa, pi[4], *conP, label, label2;
   char fossil, *name, *annotation;
}  *nodes, **gnodes, nodes_t[2 * NS - 1];


/* for stree.nodes[].fossil: lower, upper, bounds, gamma, inverse-gamma */
enum { LOWER_F = 1, UPPER_F, BOUND_F } FOSSIL_FLAGS;
char *fossils[] = { " ", "L", "U", "B" };


struct TIPDATE {
   int flag, ymd;
   double timeunit, youngest;
}  tipdate;

struct SPECIESTREE {
   int nbranch, nnode, root, nspecies, nfossil;
   struct TREESPN {
      char name[LSPNAME + 1], fossil, usefossil;  /* fossil: 0, 1, 2, 3 */
      int father, nson, sons[2];
      double label, age, pfossil[7];   /* lower and upper bounds or alpha & beta */
      double *lnrates;          /* log rates for loci */
   } nodes[2 * NS - 1];
}  stree;
/* all trees are binary & rooted, with ancestors unknown. */

struct DATA { /* locus-specific data and tree information */
   int ns[NGENE], ls[NGENE], npatt[NGENE], ngene, lgene[NGENE];
   int root[NGENE + 1], BlengthMethod, fix_nu, nbrate[NGENE];
   char *z[NGENE][NS], cleandata[NGENE];
   double *fpatt[NGENE], lnpT, lnpR, lnpDi[NGENE];
   double Qfactor[NGENE], pi[NGENE][NCODE], nu[NGENE];
   double rgene[NGENE], kappa[NGENE], alpha[NGENE], omega[1];
   int NnodeScale[NGENE];
   char *nodeScale[NGENE];    /* nScale[data.ns[locus]-1] for interior nodes */
}  data;

int nR = NCODE, LASTROUND, BayesEB;
double PMat[NCODE*NCODE], Cijk[NCODE*NCODE*NCODE], Root[NCODE];
int StepMatrix[NCODE*NCODE];


FILE *frub, *flnf, *frst, *frst1, *frst2 = NULL, *finitials = NULL;
char *ratef = "rates";
char *models[] = { "JC69","K80","F81","F84","HKY85","T92","TN93","REV","UNREST", "REVu","UNRESTu" };
enum { JC69, K80, F81, F84, HKY85, T92, TN93, REV, UNREST, REVu, UNRESTu } MODELS;
char *clockstr[] = { "", "Global clock", "Local clock", "ClockCombined" };
enum { GlobalClock = 1, LocalClock, ClockCombined } ClockModels;

double _rateSite = 1; /* rate for a site */
int N_PMatUVRoot = 0;


#define BASEML 1
#include "treesub.c"
#include "treespace.c"

int main (int argc, char *argv[])
{
   FILE* fout, * fseq = NULL, * ftree = NULL, * fpair[6];
   char pairfs[1][32] = { "2base.t" };
   char rstf[2048] = "rst", ctlf[4096] = "baseml.ctl", timestr[64];
   char *Mgenestr[] = { "diff. rate", "separate data", "diff. rate & pi",
                     "diff. rate & kappa", "diff. rate & pi & kappa" };
   int getdistance = 1, i, k;
   size_t s2 = 0;
   
   if (argc > 2 && !strcmp(argv[argc - 1], "--stdout-no-buf"))
      setvbuf(stdout, NULL, _IONBF, 0);

   starttimer();
   SetSeed(-1, 0);

   com.ndata = 1;
   com.cleandata = 0;  noisy = 0;  com.runmode = 0;

   com.clock = 0;
   com.fix_rgene = 0;  /* 0: estimate rate factors for genes */
   com.nhomo = 0;
   com.getSE = 0;      /* 1: want S.E.s of estimates;  0: don't want them */

   com.seqtype = 0;   com.model = 0;
   com.fix_kappa = 0; com.kappa = 5;
   com.fix_alpha = 1; com.alpha = 0.;  com.ncatG = 4;    /* alpha=0 := inf */
   com.fix_rho = 1;   com.rho = 0;     com.nparK = 0;
   com.ncode = 4;     com.icode = 0;
   com.print = 0;     com.verbose = 1;  com.fix_blength = 0;
   com.method = 0;    com.space = NULL;

   printf("BASEML in %s\n", pamlVerStr);
   frub = gfopen("rub", "w");  frst = gfopen(rstf, "w"); frst1 = gfopen("rst1", "w");

   if (argc > 1)  strncpy(ctlf, argv[1], 4095);
   GetOptions(ctlf);
   if (com.runmode != -2) finitials = fopen("in.baseml", "r");
   if (com.getSE == 2)    frst2 = fopen("rst2", "w");

   fprintf(frst, "Supplemental results for BASEML\n\nseqf:  %s\ntreef: %s\n",
      com.seqf, com.treef);
   fout = gfopen(com.outf, "w");
   fpair[0] = (FILE*)gfopen(pairfs[0], "w");

   /* for stepwise addition, com.sspace should be calculated using com.ns. */
   com.sspace = 1000000 * sizeof(double);
   if ((com.space = (double*)malloc(com.sspace)) == NULL)
      error2("oom space");

   if (com.clock == 5 || com.clock == 6)
      DatingHeteroData(fout);

   if (com.seqf[0] == '\0' || (fseq = fopen(com.seqf, "r")) == NULL) {
      printf("\n\nSequence file %s not found!\n", com.seqf);
      exit(-1);
   }
   if ((ftree = fopen(com.treef, "r")) == NULL) {
      printf("\ntree file %s not found.\n", com.treef);
      exit(-1);
   }
   if (com.ndata_trees_opt >= 2) {
      ReadMainTree(ftree, &stree);
      fclose(ftree);
      if (com.ndata_trees_opt == 3) {
         printf("\nPrinting trees for datasets into the file genetrees.trees\n\n");
         GenerateGtrees(fout, fseq, "genetrees.trees");
         fclose(fseq);   fclose(fout);
         exit(0);
      }
   }

   gnodes = malloc(sizeof(struct TREEN*));
   for (com.idata = 0; com.idata < com.ndata; com.idata++) {
      if (com.ndata > 1) {
         printf("\nData set %4d\n", com.idata + 1);
         fprintf(fout, "\n\nData set %4d\n", com.idata + 1);

         fprintf(frst1, "%4d", com.idata + 1);
      }
      if (com.idata)  GetOptions(ctlf); /* com.cleandata is read from here again. */
      ReadSeq((com.verbose ? fout : NULL), fseq, com.cleandata, 0, 0);
      if (com.ngene > 1 && (com.fix_blength == 2 || com.fix_blength == 3))
         error2("fix_blength = 2 or 3 does not work for partitioned data or Mgene models");
      if (com.fix_blength == 3) {
         printf("\nRelative branch lengths in tree fixed; estimating a scale factor.\n");
         if ((com.blengths0 = (double*)malloc((com.ns * 2 - 2) * sizeof(double))) == NULL)
            error2("oom blengths0");
      }

      SetMapAmbiguity(com.seqtype, 0);

      /* AllPatterns(fout); */

      if (com.rho && com.readpattern) error2("rho doesn't work with readpattern.");
      if (com.ndata == 1) fclose(fseq);
      k = (com.ns * 2 - 1) * sizeof(struct TREEN);
      if ((nodes = (struct TREEN*)malloc(k)) == NULL) error2("oom");
      gnodes[0] = nodes;

      if (com.coding) {
         if (com.ls % 3 != 0 || (com.ngene != 1 && com.ngene != 3))
            error2("this is not a coding sequence.  Remove icode?");
      }

      if (com.Mgene && com.ngene == 1)  error2("option Mgene for 1 gene?");
      if (com.ngene > 1 && com.nhomo) error2("nhomo for mutliple genes?");
      if (com.nalpha && (!com.alpha || com.ngene == 1 || com.fix_alpha))
         error2("Malpha");
      if (com.nalpha > 1 && com.rho) error2("Malpha or rho");
      if (com.alpha == 0)  com.nalpha = 0;
      else                 com.nalpha = (com.nalpha ? com.ngene : !com.fix_alpha);
      if (com.Mgene == 1)  com.nalpha = !com.fix_alpha;

      if (com.ngene == 1) com.Mgene = 0;
      if ((com.nhomo == 1 && com.ngene > 1) || (com.Mgene > 1 && com.nhomo >= 1))
         error2("nhomo does not work with Mgene options");

      if ((com.Mgene >= 2 && com.model == JC69) || (com.Mgene >= 3 && com.model == F81)
         || ((com.Mgene == 2 || com.Mgene == 4) && com.model == K80)
         || (com.Mgene > 1 && com.nhomo > 1)
         || (com.Mgene >= 2 && (com.model == UNREST || com.model == UNRESTu)))
         error2("model || Mgene");
      fprintf(fout, "\nBASEML (in %s)  %s  %s ", pamlVerStr, com.seqf, models[com.model]);
      if (com.clock) fprintf(fout, " %s ", clockstr[com.clock]);
      if (!com.nparK && com.alpha && com.rho) fprintf(fout, "  Auto-");
      if (com.alpha) fprintf(fout, "dGamma (ncatG=%d)", com.ncatG);
      if (com.nalpha > 1) fprintf(fout, "(%d gamma)", com.nalpha);
      if (com.ngene > 1)
         fprintf(fout, " (%d genes: %s)  ", com.ngene, Mgenestr[com.Mgene]);
      if (com.nhomo > 1)
         fprintf(fout, "\nNonhomo:%2d  fix_kappa%2d\n", com.nhomo, com.fix_kappa);
      if (com.nparK && com.ncatG > 6)
         printf("\a\n%d rate categories for nonparametric model?\n", com.ncatG);
      if (com.nparK) fprintf(fout, "\nnparK:%4d  K:%4d\n", com.nparK, com.ncatG);

      if (getdistance) {
         i = com.ns*(com.ns - 1) / 2;
         SeqDistance = (double*)realloc(SeqDistance, i * sizeof(double));
         ancestor = (int*)realloc(ancestor, i * sizeof(int));
         if (SeqDistance == NULL || ancestor == NULL) error2("oom distance&ancestor");
      }
      InitializeBaseAA(fout);

      if (com.Mgene == 3)
         for (i = 0; i < com.ngene; i++)
            xtoy(com.pi, com.piG[i], com.ncode);
      if (com.model == JC69 && com.ngene <= 1 && !com.readpattern && !com.print) {
         PatternWeightJC69like();
         if (fout) {
            fprintf(fout, "\n\nPrinting out site pattern counts\n");
            printPatterns(fout);
         }
         if (com.verbose >= 2) {
            fprintf(fout, "\nSite-to-pattern map: ");
            for (i = 0; i < com.ls; i++)
               fprintf(fout, " %2d", com.pose[i] + 1);
            fprintf(fout, "\n");

            fprintf(fout, "\n**** Alignment mutated for JC69-like models ****\n");
            fprintf(fout, "\n%6d %6d\n", com.ns, com.ls);
            for (i = 0; i < com.ns; i++) {
               fprintf(fout, "\n%-30s  ", com.spname[i]);
               print1seq(fout, com.z[i], com.ls, com.pose);
            }
            fprintf(fout, "\n");
         }
      }
      com.sconP = (com.ns - 1)*com.ncode*(size_t)com.npatt * sizeof(double);
      if ((com.conP = (double*)realloc(com.conP, com.sconP)) == NULL)
         error2("oom conP");
      if (com.alpha || com.nparK) {
         s2 = com.npatt*com.ncatG * sizeof(double);
         if ((com.fhK = (double*)realloc(com.fhK, s2)) == NULL) error2("oom");
      }

      printf("\n%9zu bytes for distance ", com.ns*(com.ns - 1) / 2 * (sizeof(double) + sizeof(int)));
      printf("\n%9zu bytes for conP\n", com.sconP);
      printf("%9zu bytes for fhK\n%9zu bytes for space\n", s2, com.sspace);

      /* for(i=0; i<com.ns*2-1; i++) xtoy(com.pi,nodes[i].pi, 4);  check this? */

      if (getdistance)
         DistanceMatNuc(fout, fpair[0], com.model, com.alpha);


#if(0)
      /* Selecting the most divergent sequences using the distance matrix
      */
      {
         char keep[NS];
         FILE *ftmp = gfopen("newdata.txt", "w");
         int is, js, i, j, nskeep, isbest = 0, isk, chosen[NS];
         double dmax, dminmax, d, dbig = 9;

         if (com.model == 0 && com.print == 0)
            error2("choose RateAncestor = 1 to print the seqs correctly.");
         nskeep = 5;
         printf("\nPicking up the most different sequences.\nHow many do you want? ");
         scanf("%d", &nskeep);

         if (nskeep >= com.ns) error2("nothing to do");

         for (is = 0; is < com.ns; is++) {
            keep[is] = 0;
            chosen[is] = -1;
         }

         for (is = 0, dmax = 0; is < com.ns; is++) {
            for (js = 0; js < is; js++) {
               d = SeqDistance[is*(is - 1) / 2 + js];
               if (dmax < d) { dmax = d; chosen[1] = is; chosen[0] = js; }
            }
         }
         keep[chosen[0]] = keep[chosen[1]] = 1;
         printf("selected seq %3d %s\n", chosen[0] + 1, com.spname[chosen[0]]);
         printf("selected seq %3d %s\n", chosen[1] + 1, com.spname[chosen[1]]);
         for (isk = 2; isk < nskeep; isk++) {
            for (is = 0, dminmax = 0; is < com.ns; is++) {
               if (keep[is]) continue;
               /* d is the smallest distance to those chosen */
               for (js = 0, d = dbig; chosen[js] != -1; js++) {
                  i = max2(is, chosen[js]);
                  j = min2(is, chosen[js]);
                  if (d > SeqDistance[i*(i - 1) / 2 + j]) d = SeqDistance[i*(i - 1) / 2 + j];
               }
               if (dminmax < d) {
                  dminmax = d;  isbest = is;
               }
            }
            keep[isbest] = 1;
            chosen[isk] = isbest;
            printf("selected seq %5d (dmin = %8.4f): %s\n", isbest + 1, dminmax, com.spname[isbest]);
         }

         fprintf(ftmp, "%6d %6d\n", nskeep, com.ls);
         for (j = 0; j < com.ns; j++) {
            if (keep[j] == 0) continue;
            fprintf(ftmp, "%-40s  ", com.spname[j]);
            print1seq(ftmp, com.z[j], com.ls, com.pose);
            fprintf(ftmp, "\n");
         }
         fclose(ftmp);
         return(0);

      }

#endif

      if (com.Mgene == 1)        MultipleGenes(fout, ftree, fpair, com.space);
      else if (com.runmode == 0) Forestry(fout, ftree);
      else if (com.runmode == 2) fprintf(frst, "%2d", StarDecomposition(fout, com.space));
      else if (com.runmode == 3) StepwiseAddition(fout, com.space);
      else if (com.runmode >= 4) Perturbation(fout, (com.runmode == 4), com.space);

      fprintf(frst, "\n"); 
      if ((com.idata + 1) % 10 == 0) fflush(frst);
      if (com.ndata > 1 && com.runmode) {
         fprintf(frst1, "\t");
         OutTreeN(frst1, 1, 0);
      }
      fprintf(frst1, "\n"); 
      fflush(frst1);
      free(gnodes[0]);
      memset(com.spname, 0, NS * sizeof(char));
      if (com.fix_blength == 3)  free(com.blengths0);

      printf("\nTime used: %s\n", printtime(timestr));
   }   /* for(com.idata) */
   free(gnodes);
   if (com.ndata > 1 && fseq) fclose(fseq);
   if (com.ndata_trees_opt <= 1 && ftree)
      fclose(ftree);
   free(com.space);
   fclose(fout);  fclose(frst);   fclose(fpair[0]);
   if (finitials) { fclose(finitials);  finitials = NULL; }

   return (0);
}


/* x[]: t[ntime], rgene[ngene-1], kappa[nbranch], pi[nnode*3],
        { alpha, rho || rK[], fK[] || rK[], MK[] }
*/

void PrintBestTree(FILE *fout, FILE *ftree, int btree);

int Forestry(FILE *fout, FILE *ftree)
{
   static int ages = 0;
   int status = 0, inbasemlg = 0, i, j = 0, itree = 0, ntree=1, np, iteration = 1;
   int pauptree = 0, btree = 0, haslength;
   double x[NP], xcom[NP - NBRANCH], lnL, lnL0 = 0, lnLbest = 0, e = 1e-7, nchange = -1;
   double xb[NP][2], tl = 0, *g = NULL, *H = NULL;
   FILE *finbasemlg = NULL, *frate = NULL;

   if (com.ndata_trees_opt <= 1) {
      if (com.ndata_trees_opt == 0) rewind(ftree);
      j = GetTreeFileType(ftree, &ntree, &pauptree, 0);
   }
   if (com.alpha)
      frate = gfopen(ratef, "w");
   if (com.alpha && com.rho == 0 && com.nhomo == 0 && com.nparK == 0 && com.ns < 15) {
      inbasemlg = 1;  finbasemlg = gfopen("in.basemlg", "w");
   }
   flnf = gfopen("lnf", "w+");
   fprintf(flnf, "%6d %6d %6d\n", ntree, com.ls, com.npatt);

   for (itree = 0; ntree == -1 || itree < ntree; itree++, iteration = 1) {
      if (com.ndata_trees_opt <= 1 && ReadTreeN(ftree, &haslength, 0, 1))
         error2("end of tree file.");
      if (com.ndata_trees_opt>=2) {
         printf("\nExtracting gene tree for dataset #%2d from main tree...\n", com.idata + 1);
         GenerateGtree_locus(com.idata, com.ns, 0);
      }
      ProcessNodeAnnotation(&i);

      if (noisy) printf("\nTREE # %2d: ", itree + 1);
      fprintf(fout, "\nTREE # %2d:  ", itree + 1);
      fprintf(frub, "\n\nTREE # %2d\n", itree + 1);
      if (com.print) fprintf(frst, "\n\nTREE # %2d\n", itree + 1);
      fprintf(flnf, "\n\n%2d\n", itree + 1);

      LASTROUND = 0;
      if ((com.fix_blength == 2 || com.fix_blength == 3) && !haslength)
         error2("We need branch lengths in tree");
      if (com.fix_blength > 0 && !haslength) com.fix_blength = 0;
      if (ages++ == 0 && com.fix_blength > 0 && haslength) {
         if (com.clock) puts("\nBranch lengths in tree are ignored");
         else {
            if (com.fix_blength == 2) puts("\nBranch lengths in tree are fixed.");
            else if (com.fix_blength == 1)
               puts("\nBranch lengths in tree used as initials.");
            if (com.fix_blength == 1) {
               for (i = 0; i < tree.nnode; i++)
                  if ((x[nodes[i].ibranch] = nodes[i].branch) < 0)
                     x[nodes[i].ibranch] = 1e-5;
            }
         }
      }

      if (com.cleandata)
         nchange = MPScore(com.space);
      if (noisy&&com.ns < 99)  {
         OutTreeN(F0, 0, 0); printf(" MP score: %.2f\n", nchange);
      }
      OutTreeN(fout, 0, 0);  fprintf(fout, "  MP score: %.2f", nchange);
      if (!com.clock && com.model <= REV && com.nhomo <= 2
         && nodes[tree.root].nson <= 2 && com.ns > 2) {
         puts("\nThis is a rooted tree, without clock.  Check.");
         if (com.verbose) fputs("\nThis is a rooted tree.  Please check!", fout);
      }

      fflush(fout);  fflush(flnf);

      GetInitials(x, &i);

      if (i == -1) iteration = 0;

      if ((np = com.np) > NP) error2("raise NP and recompile");

      if (com.sspace < spaceming2(np)) {
         com.sspace = spaceming2(np);
         if ((com.space = (double*)realloc(com.space, com.sspace)) == NULL)
            error2("oom space");
      }

      if (itree && !finitials && (com.nhomo == 0 || com.nhomo == 2))
         for (i = 0; i < np - com.ntime; i++) x[com.ntime + i] = max2(xcom[i], 0.001);

      PointconPnodes();

      lnL = com.plfun(x, np);
      if (noisy) {
         printf("\nntime & nrate & np:%6d%6d%6d\n", com.ntime, com.nrate, com.np);
         if (noisy > 2 && com.ns < 50) { 
            OutTreeB(F0); printf("\n");
            matout(F0, x, 1, np); 
         }
         printf("\nlnL0 = %12.6f\n", -lnL);
      }

      if (iteration && np) {
         SetxBound(np, xb);
         SetxInitials(np, x, xb); /* start within the feasible region */
         if (com.method == 1)
            j = minB(noisy > 2 ? frub : NULL, &lnL, x, xb, e, com.space);
         else if (com.method == 3)
            j = minB2(noisy > 2 ? frub : NULL, &lnL, x, xb, e, com.space);
         else
            j = ming2(noisy > 2 ? frub : NULL, &lnL, com.plfun, NULL, x, xb, com.space, e, np);

         if (j == -1 || lnL <= 0 || lnL > 1e7) status = -1;
         else                                  status = 0;
      }

      if (itree == 0) {
         lnL0 = lnLbest = lnL; btree = 0;
      }
      else if (lnL < lnLbest) {
         lnLbest = lnL;  btree = itree;
      }
      if (noisy) printf("Out...\nlnL  = %12.6f\n", -lnL);
      fprintf(fout, "\nlnL(ntime:%3d  np:%3d): %13.6f %+12.6f\n",
         com.ntime, np, -lnL, -lnL + lnL0);
      if (com.fix_blength < 2) { 
         OutTreeB(fout);  fprintf(fout, "\n");
      }
      if (LASTROUND == 0) {
         LASTROUND = 1;
         if ((com.npi && com.model != T92) || com.nparK >= 2) TransformxBack(x);
      }
      if (com.clock) { /* move ages into x[] */
         for (i = 0, j = !nodes[tree.root].fossil; i < tree.nnode; i++)
            if (i != tree.root && nodes[i].nson && !nodes[i].fossil)
               x[j++] = nodes[i].age;
      }
      for (i = 0; i < np; i++)
         fprintf(fout, " %8.6f", x[i]);
      fprintf(fout, "\n"); 
      fflush(fout);
      if (inbasemlg) matout(finbasemlg, x, 1, np);

      /*
      for(i=0;i<np;i++) fprintf(frst1,"\t%.6f", x[i]);
      fprintf(frst1,"\t%.4f", -lnL);
      */
      if (com.getSE) {
         puts("Calculating SE's");
         if (com.sspace < np*(np + 1) * sizeof(double)) {
            com.sspace = np*(np + 1) * sizeof(double);
            if ((com.space = (double*)realloc(com.space, com.sspace)) == NULL)
               error2("oom space for SE");
         }

         g = com.space;
         H = g + com.np;
         HessianSKT2004(x, lnL, g, H);
         if (com.getSE >= 2 && com.clock == 0 && nodes[tree.root].nson == 3) {  /* g & H */
            fprintf(frst2, "\n %d\n\n", com.ns);
            OutTreeN(frst2, 1, 1);  fprintf(frst2, "\n\n");
            for (i = 0; i < com.ntime; i++)
               if (x[i] > 0.0004 && fabs(g[i]) < 0.005) g[i] = 0;
            for (i = 0; i < com.ntime; i++) fprintf(frst2, " %9.6f", x[i]);
            fprintf(frst2, "\n\n");
            for (i = 0; i < com.ntime; i++) fprintf(frst2, " %9.6f", g[i]);
            fprintf(frst2, "\n\n");
            fprintf(frst2, "\nHessian\n\n");
            for (i = 0; i < com.ntime; i++) {
               for (j = 0; j < com.ntime; j++)
                  fprintf(frst2, " %10.4g", H[i * np + j]);
               fprintf(frst2, "\n");
            }
            fflush(frst2);
         }
         for (i = 0; i < np*np; i++)  H[i] *= -1;
         matinv(H, np, np, H + np*np);
         fprintf(fout, "SEs for parameters:\n");
         for (i = 0; i < np; i++)
            fprintf(fout, " %8.6f", (H[i*np + i] > 0. ? sqrt(H[i*np + i]) : -1));
         fprintf(fout, "\n");
      }

      /* if(com.clock) SetBranch(x); */
      /* GroupDistances(); */

      if (com.nbtype > 1 && com.clock>1)
         fprintf(fout, "\nWarning: branch rates are not yet applied in tree length and branch lengths");
      if (AbsoluteRate) {
         fprintf(fout, "\nNote: mutation rate is not applied to tree length.  Tree has ages, for TreeView & FigTree\n");
         for (i = 0; i < tree.nnode; i++)
            nodes[i].branch *= tipdate.timeunit;
      }

      if (!AbsoluteRate) {
         for (i = 0, tl = 0; i < tree.nnode; i++)
            if (i != tree.root) tl += nodes[i].branch;
         fprintf(fout, "\ntree length = %9.5f%s\n", tl, (com.ngene > 1 ? " (1st gene)" : ""));
      }
      fprintf(fout, "\n"); 
      OutTreeN(fout, 1, 0); fprintf(fout, "\n\n");
      OutTreeN(fout, 1, 1); fprintf(fout, "\n");
      if (com.clock) {
         fprintf(fout, "\n"); 
         OutTreeN(fout, 1, PrNodeNum); 
         fprintf(fout, "\n");
      }

      if (tipdate.flag) {  /* scale back the ages in nodes[].branch */
         for (i = 0; i < tree.nnode; i++)  nodes[i].branch /= tipdate.timeunit;
      }
      if (com.np - com.ntime || com.clock)
         DetailOutput(fout, x, H);
      if (status) {
         printf("convergence?\n");
         fprintf(fout, "check convergence..\n");
      }
      if (itree == 0)
         for (i = 0; i < np - com.ntime; i++) xcom[i] = x[com.ntime + i];
      else if (!j)
         for (i = 0; i < np - com.ntime; i++) xcom[i] = xcom[i] * .8 + x[com.ntime + i] * 0.2;
      /*
            TestModel(fout,x,1,com.space);
            fprintf(fout,"\n\n# of lfun calls:%10d\n", NFunCall);
      */

      if (com.coding && com.Mgene != 1 && !AbsoluteRate) {
         fputs("\nTree with branch lengths for codon models:\n", fout);
         for (i = 0; i < tree.nnode; i++)
            nodes[i].branch *= (1 + com.rgene[1] + com.rgene[2]);
         fprintf(fout, "\n"); OutTreeN(fout, 1, 1); fprintf(fout, "\n");
         for (i = 0; i < tree.nnode; i++)
            nodes[i].branch /= (1 + com.rgene[1] + com.rgene[2]);
      }

      com.print -= 9;  com.plfun(x, np);  com.print += 9;
      if (com.print) {
         if (com.plfun != lfun)
            lfunRates(frate, x, np);

         /** think about more-general models.  Check that this may not be working for the clock models clock>=2 **/
         if ((com.nhomo <= 2 && com.nparK == 0 && com.model <= REV && com.rho == 0)
            || (com.nhomo > 2 && com.alpha == 0 && com.Mgene == 0))
            AncestralSeqs(frst, x);
      }
   }         /* for (itree) */

   if (ntree > 1) {
      fprintf(frst1, "\t%d\t", btree + 1);
      PrintBestTree(frst1, ftree, btree);
   }
   if (finbasemlg) fclose(finbasemlg);
   if (frate) fclose(frate);

   if (ntree == -1)  ntree = itree;
   if (ntree > 1) { rewind(flnf);  rell(flnf, fout, ntree); }
   fclose(flnf);

   return(0);
}


void PrintBestTree(FILE *fout, FILE *ftree, int btree)
{
   int itree, ntree, i;

   rewind(ftree);
   GetTreeFileType(ftree, &ntree, &i, 0);
   for (itree = 0; ntree == -1 || itree < ntree; itree++) {
      if (ReadTreeN(ftree, &i, 0, 1)) {
         puts("\nend of tree file."); break;
      }
      if (itree == btree)
         OutTreeN(fout, 1, 0);
   }
}


int TransformxBack(double x[])
{
   /* transform variables x[] back to their original definition after iteration,
      for output and for calculating SEs.
   */
   int i, k, K = com.ncatG;

   k = com.ntime + com.nrgene + com.nrate;
   for (i = 0; i < com.npi; i++)
      f_and_x(x + k + 3 * i, x + k + 3 * i, 4, 0, 0);

   k += com.npi * 3 + K - 1;        /* K-1 for rK */
   if (com.nparK == 2)          /* rK & fK */
      f_and_x(x + k, x + k, K, 0, 0);
   else if (com.nparK == 3)     /* rK & MK (double stochastic matrix) */
      for (i = 0; i < K - 1; k += K - 1, i++)  f_and_x(x + k, x + k, K, 0, 0);
   else if (com.nparK == 4)     /* rK & MK */
      for (i = 0; i < K; k += K - 1, i++)  f_and_x(x + k, x + k, K, 0, 0);
   return(0);
}


void DetailOutput(FILE *fout, double x[], double var[])
{
   int n = com.ncode, i, j, k = com.ntime, nkappa[] = { 0, 1, 0, 1, 1, 1, 2, 5, 11 };
   int n31pi = (com.model == T92 ? 1 : 3);
   double *p = com.pi, t = 0, k1, k2, S, V, Y = p[0] + p[1], R = p[2] + p[3];
   int inode, a;
   double *Qrate = x + com.ntime + com.nrgene, *AncientFreqs = NULL;
   double EXij[NCODE * NCODE] = { 0 }, * Q = PMat, p0[NCODE], c[NCODE], ya;
   double tnew, treelength;

   fprintf(fout, "\nDetailed output identifying parameters\n");
   if (com.clock) OutputTimesRates(fout, x, var);
   k = com.ntime;

   if (com.nrgene) { /* this used to be:  if (com.nrgene && !com.clock) */
      fprintf(fout, "\nrates for %d genes:%6.0f", com.ngene, 1.);
      for (i = 0; i < com.nrgene; i++)
         fprintf(fout, " %8.5f", x[k++]);
      fprintf(fout, "\n");
   }

   if (com.nhomo == 1) {
      if (com.nrate) fprintf(fout, "rate (kappa or abcde) under %s:", models[com.model]);
      for (i = 0; i < com.nrate; i++) fprintf(fout, " %8.5f", x[k++]);
      fprintf(fout, "\nBase frequencies:\n");
      for (j = 0; j < n; j++) fprintf(fout, " %8.5f", com.pi[j]);
      k += n31pi;
   }
   else if (com.nhomo >= 3) {
      if (com.model >= F84 && com.model <= REV )
         fprintf(fout, "kappa or abcde or Qrates under %s:", models[com.model]);
      for (i = 0; i < com.nrate; i++) {
         if (i == 0 || i % nkappa[com.model] == 0) fprintf(fout, "\n");
         fprintf(fout, " %8.5f", x[k++]);
      }
      SetParameters(x);
      if (com.alpha == 0) {
         if ((AncientFreqs = (double*)malloc(tree.nnode * 4 * sizeof(double))) == NULL)
            error2("out of memory for OldDistributions()");
         OldDistributions(tree.root, AncientFreqs);
      }
      fputs("\n\n(frequency parameters for branches)  [frequencies at nodes] (see Yang & Roberts 1995 fig 1)\n", fout);
      for (i = 0; i < tree.nnode; i++) {
         fprintf(fout, "Node #%2d  (", i + 1);
         for (j = 0; j < 4; j++) fprintf(fout, " %8.6f ", nodes[i].pi[j]);
         fprintf(fout, ")");
         if (com.alpha == 0) {
            fprintf(fout, "  [");
            for (j = 0; j < 4; j++) fprintf(fout, " %8.6f", AncientFreqs[i * 4 + j]);
            fprintf(fout, " GC = %5.3f ]", AncientFreqs[i * 4 + 1] + AncientFreqs[i * 4 + 3]);
         }
         fprintf(fout, "\n");
      }
      fprintf(fout, "\nNote: node %d is root.\n", tree.root + 1);
      k += com.npi*n31pi;

      /* print expected numbers of changes along branches */
      if (com.alpha == 0) {
         fprintf(fout, "\n\nExpected numbers of nucleotide changes on branches\n\n");
         for (inode = 0, treelength = 0; inode < tree.nnode; inode++) {
            if (inode == tree.root) continue;
            t = nodes[inode].branch;
            xtoy(AncientFreqs + nodes[inode].father*n, p0, n);

            if (com.nhomo > 2 && com.model <= TN93) {
               eigenTN93(com.model, *nodes[inode].pkappa, *(nodes[inode].pkappa + 1), nodes[inode].pi, &nR, Root, Cijk);
               QTN93(com.model, Q, *nodes[inode].pkappa, *(nodes[inode].pkappa + 1), nodes[inode].pi);
            }
            else if (com.nhomo > 2 && com.model == REV)
               eigenQREVbase(NULL, Q, nodes[inode].pkappa, nodes[inode].pi, &nR, Root, Cijk);

            for (i = 0; i < n; i++)  { /* calculate correction vector c[4] */
               c[i] = 0;   Q[i*n + i] = 0;
            }

            for (i = 0; i < n; i++) {  /* calculate correction vector c[4] */
               for (a = 1; a < n; a++) {
                  ya = -(1 - exp(Root[a] * t)) / Root[a];
                  for (k = 0; k < n; k++)
                     c[i] += p0[k] * Cijk[k*n*n + i*n + a] * ya;
               }
               for (j = 0; j < n; j++)
                  EXij[i*n + j] = Q[i*n + j] * (nodes[inode].pi[i] * t + c[i]);
            }
            fprintf(fout, "Node #%2d  ", inode + 1);
            if (inode < com.ns) fprintf(fout, "%-10s  ", com.spname[inode]);
            tnew = sum(EXij, n * n);
            treelength += tnew;
            fprintf(fout, "(blength = %9.6f, %9.6f)", t, tnew);
            fprintf(fout, "\n%s:     ", (com.model <= TN93 ? "kappa" : "abcde"));
            for (j = 0; j < nkappa[com.model]; j++)  fprintf(fout, " %9.6f ", *(nodes[inode].pkappa + j));
            fprintf(fout, "\npi source: "); for (j = 0; j < n; j++) fprintf(fout, " %9.6f ", p0[j]);
            fprintf(fout, "\npi model : "); for (j = 0; j < n; j++) fprintf(fout, " %9.6f ", nodes[inode].pi[j]);
            fprintf(fout, "\nExij\n");
            for (i = 0; i < n; i++) {
               for (j = 0; j < n; j++)  fprintf(fout, " %9.6f", EXij[i*n + j]);
               fprintf(fout, "  ");
               for (j = 0; j < n; j++)  fprintf(fout, " %8.1f", EXij[i*n + j] * com.ls);
               fprintf(fout, "\n");
            }
            /* matout(fout, c, 1, n);  */
            fprintf(fout, "\n");
         }  /* for (inode) */
         fprintf(fout, "\ntree length (eq. 6 in Matsumoto et al. 2015) = %9.6f\n", treelength);
      }

      if (com.alpha == 0) free(AncientFreqs);
   }
   else if (!com.fix_kappa) {
      fprintf(fout, "\nParameters %s in the rate matrix", (com.model <= TN93 ? "(kappa)" : ""));
      fprintf(fout, " (%s) (Yang 1994 J Mol Evol 39:105-111):\n", models[com.model]);

      if (com.nhomo == 2) {
         fprintf(fout, "\nbranch         t    kappa      TS     TV\n");
         for (i = 0; i < tree.nbranch; i++) {
            if (com.model == F84) { k1 = 1 + x[k + i] / R;  k2 = 1 + x[k + i] / R; }
            else                   k1 = k2 = x[k + i];
            S = 2 * p[0] * p[1] * k1 + 2 * p[2] * p[3] * k2;
            V = 2 * Y*R;
            /* Qfactor = 1 / (S + V); */
            /* t=(com.clock ? nodes[tree.branches[i][1]].branch : x[i]); */
            t = nodes[tree.branches[i][1]].branch;
            fprintf(fout, "%2d..%-2d %9.5f %8.5f %9.5f %8.5f\n",
               tree.branches[i][0] + 1, tree.branches[i][1] + 1, t, x[k + i],
               t*S / (S + V), t*V / (S + V));
         }
      }
      /* Mgene = 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff*/
      else if (com.Mgene >= 2) {
         for (i = 0; i < com.ngene; i++) {
            fprintf(fout, "\nGene #%d: ", i + 1);
            p = (com.Mgene == 3 ? com.pi : com.piG[i]);
            Qrate = (com.Mgene == 2 ? x + k : x + k + i*nkappa[com.model]);
            if (com.model <= TN93)
               for (j = 0; j < nkappa[com.model]; j++)
                  fprintf(fout, " %8.5f", Qrate[j]);
            else if (com.model == REV || com.model == REVu)
               /* output Q matrix, no eigen calculation */
               eigenQREVbase(fout, Q, Qrate, p, &nR, Root, Cijk);
         }
         if (com.Mgene >= 3) k += com.ngene*nkappa[com.model];
         else              k += nkappa[com.model];
      }
      else {
         if (com.model < REV)
            for (i = 0; i < com.nrate; i++)
               fprintf(fout, " %8.5f", x[k++]);
         else
            k += com.nrate;
      }
      fprintf(fout, "\n");
   }

   if (com.Mgene < 2) {
      if ((com.model == REV || com.model == REVu) && com.nhomo <= 2) /* output Q, no eigen calculation */
         eigenQREVbase(fout, Q, Qrate, com.pi, &nR, Root, Cijk);
      else if (com.model == UNREST || com.model == UNRESTu)
         QUNREST(fout, PMat, x + com.ntime + com.nrgene, com.pi);
   }

   for (j = 0; j < com.nalpha; j++) {
      if (!com.fix_alpha)
         fprintf(fout, "\nalpha (gamma, K=%d) = %8.5f", com.ncatG, (com.alpha = x[k++]));
      if (com.nalpha > 1)
         DiscreteGamma(com.freqK, com.rK, com.alpha, com.alpha, com.ncatG, DGammaUseMedian);
      fprintf(fout, "\nrate: "); for (i = 0; i < com.ncatG; i++) fprintf(fout, " %8.5f", com.rK[i]);
      fprintf(fout, "\nfreq: "); for (i = 0; i < com.ncatG; i++) fprintf(fout, " %8.5f", com.freqK[i]);
      fprintf(fout, "\n");
   }
   if (!com.fix_rho) {
      fprintf(fout, "rho for the auto-discrete-gamma model: %9.5f", x[k]);
      fprintf(fout, "\n");
   }
   if (com.nparK >= 1 && com.nalpha <= 1) {
      fprintf(fout, "\nrate:");  for (i = 0; i<com.ncatG; i++) fprintf(fout, " %8.5f", com.rK[i]);
      fprintf(fout, "\nfreq:");  for (i = 0; i<com.ncatG; i++) fprintf(fout, " %8.5f", com.freqK[i]);
      fprintf(fout, "\n");
   }
   if (com.rho || (com.nparK >= 3 && com.nalpha <= 1)) {
      fprintf(fout, "transition probabilities between rate categories:\n");
      for (i = 0; i < com.ncatG; i++) {
         for (j = 0; j < com.ncatG; j++)
            fprintf(fout, " %8.5f", com.MK[i * com.ncatG + j]);
         fprintf(fout, "\n");
      }
   }
   fprintf(fout, "\n");
}

int GetStepMatrix(char*line)
{
   /* read user definitions of the REV and UNREST models
      StepMatrix[4*4]:
         -1 at diagonals, 0 for default rate, and positive values for rates
   */
   char *p, *errstr = "StepMatrix specification in the control file";
   int n = com.ncode, i, j, k, b1 = -1, b2;

   p = strchr(line, '[');
   if (p == NULL) error2("model specification.  Expecting '['.");
   sscanf(++p, "%d", &com.nrate);
   if (com.nrate < 0 || (com.model == REVu && com.nrate > 5) || (com.model == UNRESTu && com.nrate > 11))
      error2(errstr);
   for (i = 0; i < n; i++) for (j = 0; j < n; j++)
      StepMatrix[i*n + j] = (i == j ? -1 : 0);
   for (i = 0; i < com.nrate; i++) {
      while (*p && *p != '(') p++;
      if (*p++ != '(') error2("expecting (");
      for (k = 0; k < 12; k++) {
         while (isspace(*p)) p++;
         if (*p == ')') break;
         b1 = CodeChara(*p++, 0);
         b2 = CodeChara(*p++, 0);
         if (b1 < 0 || b1>3 || b2 < 0 || b2>3) error2("bases out of range.");
         if (b1 == b2 || StepMatrix[b1*n + b2] > 0) {
            printf("pair %c%c already specified.\n", BASEs[b1], BASEs[b2]);
         }
         if (com.model == REVu) StepMatrix[b1*n + b2] = StepMatrix[b2*n + b1] = i + 1;
         else                StepMatrix[b1*n + b2] = i + 1;
      }
      printf("rate %d: %d pairs\n", i + 1, k);
   }
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++)
         printf("%3d", StepMatrix[i * n + j]);
      printf("\n");
   }

   return(0);
}

int GetOptions(char* ctlf)
{
   int iopt, i, j, nopt = 31, lline = 4096;
   char line[4096], * pline, opt[32], * comment = "*#", str[1024];
   char* optstr[] = { "seqfile","outfile","treefile","noisy", "ndata", "cleandata",
        "verbose","runmode", "method", "clock", "TipDate", "fix_rgene","Mgene", "nhomo",
        "getSE","RateAncestor", "model","fix_kappa","kappa",
        "fix_alpha","alpha","Malpha","ncatG", "fix_rho","rho",
        "nparK", "bootstrap", "Small_Diff","icode", "fix_blength", "seqtype" };
   double t;
   FILE* fctl;
   int ng = -1, markgenes[NGENE], it=1;

   com.ndata_trees_opt = 0;
   com.nalpha = 0;
   fctl = gfopen(ctlf, "r");
   if (noisy) printf("Reading options from %s..\n", ctlf);
   for (; ; ) {
      if (fgets(line, lline, fctl) == NULL) break;
      for (i = 0, t = 0, pline = line; i < lline && line[i]; i++) {
         if (isalnum(line[i])) {
            t = 1; break;
         }
         else if (strchr(comment, line[i]))
            break;
      }
      if (t == 0) continue;
      sscanf(line, "%s%*s%lf", opt, &t);
      if ((pline = strstr(line, "=")) == NULL)
         error2("option file.");

      for (iopt = 0; iopt < nopt; iopt++) {
         if (strncmp(opt, optstr[iopt], 8) == 0) {
            if (noisy >= 9)
               printf("\n%3d %15s | %-20s %6.2f", iopt + 1, optstr[iopt], opt, t);
            switch (iopt) {
            case (0): sscanf(pline + 1, "%s", com.seqf);    break;
            case (1): sscanf(pline + 1, "%s", com.outf);    break;
            case (2): sscanf(pline + 1, "%s", com.treef);    break;
            case (3): noisy = (int)t;           break;
            case (4): 
               sscanf(pline + 1, "%d%s%d", &com.ndata, str, &it);
               if (strstr(str, "separate_trees"))
                  com.ndata_trees_opt = 1;               /* separate trees for ndata */
               else if (strstr(str, "maintree")) {
                  com.ndata_trees_opt = 2;               /* generate genetree.trees & run ML */
                  if (it == 0)  com.ndata_trees_opt = 3; /* dry run, no ML */
               }
               break;
            case (5): com.cleandata = (char)t;  break;
            case (6): com.verbose = (int)t;     break;
            case (7): com.runmode = (int)t;     break;
            case (8): com.method = (int)t;      break;
            case (9): com.clock = (int)t;       break;
            case (10):
               sscanf(pline + 1, "%d%lf", &tipdate.flag, &tipdate.timeunit);
               break;
            case (11): com.fix_rgene = (int)t;   break;
            case (12): com.Mgene = (int)t;       break;
            case (13): com.nhomo = (int)t;       break;
            case (14): com.getSE = (int)t;       break;
            case (15): com.print = (int)t;       break;
            case (16): com.model = (int)t;
               if (com.model > UNREST) GetStepMatrix(line);
               break;
            case (17): com.fix_kappa = (int)t;   break;
            case (18):
               com.kappa = t;
               if (com.fix_kappa && (com.clock == 5 || com.clock == 6)
                  && com.model != 0 && com.model != 2) {
                  ng = splitline(++pline, min2(ng, com.ndata), markgenes);
                  for (j = 0; j < min2(ng, com.ndata); j++)
                     if (!sscanf(pline + markgenes[j], "%lf", &data.kappa[j])) break;
                  /*
                  matout(F0, data.kappa, 1, min2(ng,com.ndata));
                  */
               }
               break;
            case (19): com.fix_alpha = (int)t;   break;
            case (20):
               com.alpha = t;
               if (com.fix_alpha && t && (com.clock == 5 || com.clock == 6)) {
                  ng = splitline(++pline, min2(ng, com.ndata), markgenes);
                  for (j = 0; j < min2(ng, com.ndata); j++)
                     if (!sscanf(pline + markgenes[j], "%lf", &data.alpha[j])) break;
                  /*
                  matout(F0, data.alpha, 1, min2(ng,com.ndata));
                  */
               }
               break;
            case (21): com.nalpha = (int)t;      break;
            case (22): com.ncatG = (int)t;       break;
            case (23): com.fix_rho = (int)t;     break;
            case (24): com.rho = t;              break;
            case (25): com.nparK = (int)t;       break;
            case (26): com.bootstrap = (int)t;   break;
            case (27): Small_Diff = t;           break;
            case (28): com.icode = (int)t; com.coding = 1; break;
            case (29): com.fix_blength = (int)t; break;
            case (30):
               com.seqtype = (int)t;
               if (com.seqtype == 2)      com.ncode = 2;
               else if (com.seqtype == 5) com.ncode = 5;
               break;
            }
            break;
         }
      }
      if (iopt == nopt) {
         printf("\noption %s in %s not recognised\n", opt, ctlf); 
         exit(-1);
      }
   }
   fclose(fctl);

   if ((com.fix_kappa || (com.fix_alpha&&com.alpha)) && (com.clock == 5 || com.clock == 6))
      printf("Using parameters from the control file.");


   if (com.seqtype != 0 && com.seqtype != 4 && com.seqtype != 5) error2("seqtype?");
   if (com.fix_alpha == 1 && com.alpha == 0) {
      if (!com.fix_rho || com.rho) error2("fix rho to 0 if alpha=0.");
   }
   if (com.nparK >= 1) {
      com.fix_alpha = com.fix_rho = 1;
      if (com.alpha == 0) com.alpha = 0.5;  /* used to generate rK as initial values */
      if (com.nparK <= 2) com.rho = 0;
      else             com.rho = 0.4;
      if (com.nhomo >= 1)
         error2("nhomo & nparK");
   }
   if (com.model != F84 && com.kappa <= 0)  error2("init kappa..");
   if (!com.fix_alpha && com.alpha <= 0)  error2("init alpha..");
   if (!com.fix_rho && com.rho == 0) { com.rho = 0.001;  puts("init rho reset"); }

   if (com.alpha)
   {
      if (com.ncatG<2 || com.ncatG>NCATG) error2("ncatG");
   }
   else if (com.ncatG > 1) com.ncatG = 1;

   if (com.method && (com.clock || com.rho))
   {
      com.method = 0; puts("Iteration method reset: method = 0");
   }
   if (com.method && (com.model == UNREST || com.model == UNRESTu))
      error2("I did not implemented method = 1 for UNREST.  Use method = 0");

   if (com.model > UNREST && com.Mgene > 2)
      error2("u models don't work with Mgene");

   if (com.nhomo > 5)
      error2("nhomo");
   else if (com.nhomo == 2) {
      if (com.model != K80 && com.model != F84 && com.model != HKY85) 
         error2("nhomo = 2 works with F84 or HKY85");
   }
   else if (com.nhomo > 2 && !(com.model >= F84 && com.model <= REV)) error2("nhomo & model");
   else if (com.nhomo > 2 && com.method) error2("nhomo & method.");
   else {
      if (com.nhomo == 1 && !(com.model >= F81 && com.model <= REV) && com.model != REVu)
         error2("nhomo=1 and model");
   }

   if (com.nhomo >= 2 && com.clock == 2) error2("clock=2 & nhomo imcompatible");

   if (com.clock >= 5 && com.model >= TN93) error2("model & clock imcompatible");

   if (com.nhomo > 1 && com.runmode > 0)  error2("nhomo incompatible with runmode");
   if (com.runmode == -2) error2("runmode = -2 not implemented in baseml.");
   if (com.clock && com.runmode > 2)  error2("runmode & clock?");
   if (com.runmode)  com.fix_blength = 0;
   if (com.runmode == 3 && (com.npi || com.nparK))
      error2("runmode incompatible with nparK or nhomo.");

   if (com.model == JC69 || com.model == K80 || com.model == UNREST)
      if (com.nhomo != 2)  com.nhomo = 0;
   if (com.model == JC69 || com.model == F81) { com.fix_kappa = 1; com.kappa = 1; }
   if ((com.model == TN93 || com.model == REV || com.model == REVu) && com.nhomo <= 2)
      com.fix_kappa = 0;
   if (com.seqtype == 5 && com.model != REVu)
      error2("RNA-editing model requires REVu");

   if (com.nparK == 3) {
      puts("\n\nnparK==3, double stochastic, may not work.  Use nparK=4?\n");
      getchar();
   }
   com.fix_omega = 1; com.omega = 0;
   if (com.ndata <= 0) com.ndata = 1;
   if (com.bootstrap && com.ndata != 1) error2("ndata=1 for bootstrap.");

   return (0);
}


int GetInitials(double x[], int *fromfile)
{
   /* This calculates com.np, com.ntime, com.npi, com.nrate etc., and gets
      initial values for ML iteration.  It also calculates some of the
      statistics (such as eigen values and vectors) that won't change during
      the likelihood iteration.  However, think about whether this is worthwhile.
      The readfile option defeats this effort right now as the x[] is read in
      after such calculations.  Some of the effort is duplicated in
      SetParameters().  Needs more careful thinking.
   */
   int i, j, k, K = com.ncatG, n31pi = (com.model == T92 ? 1 : 3);
   int nkappa[] = { 0, 1, 0, 1, 1, 1, 2, 5, 11, -1, -1 }, nkappasets;
   size_t sconP_new = (size_t)(tree.nnode - com.ns)*com.ncode*com.npatt*com.ncatG * sizeof(double);
   double t = -1;

   NFunCall = NPMatUVRoot = NEigenQ = 0;
   if (com.clock == ClockCombined && com.ngene <= 1)
      error2("Combined clock model requires mutliple genes.");
   GetInitialsTimes(x);

   com.plfun = lfunAdG;
   if (com.alpha == 0 && com.nparK == 0)
      com.plfun = lfun;
   else if ((com.alpha && com.rho == 0) || com.nparK == 1 || com.nparK == 2)
      com.plfun = lfundG;

   if (com.clock && com.fix_blength == -1)
      com.fix_blength = 0;

   if (com.method && com.fix_blength != 2 && com.plfun == lfundG) {
      com.conPSiteClass = 1;
      if (com.sconP < sconP_new) {
         com.sconP = sconP_new;
         printf("\n%9zu bytes for conP, adjusted\n", com.sconP);
         if ((com.conP = (double*)realloc(com.conP, com.sconP)) == NULL)
            error2("oom conP");
      }
   }
   InitializeNodeScale();

   com.nrgene = (!com.fix_rgene)*(com.ngene - 1);
   for (j = 0; j < com.nrgene; j++) x[com.ntime + j] = 1;
   if (com.fix_kappa && (com.Mgene == 3 || com.Mgene == 4)) error2("Mgene options");

   if (com.model <= UNREST) {
      com.nrate = 0;
      if (!com.fix_kappa) {
         com.nrate = nkappa[com.model];
         if (com.Mgene >= 3)
            com.nrate *= com.ngene;
      }
   }
   switch (com.nhomo) {
   case (0): com.npi = 0;           break;   /* given 1 pi */
   case (1): com.npi = 1;           break;   /* estimate 1 pi */
   case (2): com.npi = 0;  com.nrate = tree.nbranch;  break;  /* b kappa's */
   case (4):
      com.npi = tree.nnode;     /* nnode pi   */
      com.nrate = nkappa[com.model] * (com.fix_kappa ? 1 : tree.nbranch);
      break;
   case (3):   /* ns+2 pi: tips, internal & root */
   case (5):   /* pi: user-specified sets of pi for branches */
      if (com.nhomo == 3) {
         for (i = 0; i < tree.nnode; i++)
            nodes[i].label = (i < com.ns ? i : com.ns);
         com.npi = com.nbtype = com.ns + (tree.root >= com.ns);
         if (tree.root >= com.ns)
            nodes[tree.root].label = com.npi++;
      }
      else {  /* nbtype is for Qrates & npi for pi. */
         com.npi = com.nbtype;
         if (tree.root >= com.ns) {   /* root may have a new set of pi */
            j = (int)nodes[tree.root].label;
            if (j == com.nbtype)
               nodes[tree.root].label = com.npi++;
            else if (j <0 || j > com.nbtype)
               error2("nhomo = 5: label for root strange?");
         }
      }
      
      printf("%d sets of frequency parameters\n", com.npi);
      /* number of sets of Qrates */
      j = (com.fix_kappa == 0 ? tree.nbranch : (com.fix_kappa == 1 ? 1 : com.nbtype));
      com.nrate = nkappa[com.model] * j;
      break;
   }

   if (com.model <= TN93 && com.Mgene <= 1 && com.nhomo == 0)
      eigenTN93(com.model, com.kappa, com.kappa, com.pi, &nR, Root, Cijk);
   if (com.model == REV || com.model == UNREST)
      for (j = 0; j < (com.Mgene >= 3 ? com.ngene : 1); j++) {
         k = com.ntime + com.nrgene + j*(com.model == REV ? 5 : 11);
         for (i = 0; i < com.nrate; i++)
            x[k + i] = 0.4 + 0.2*rndu();
         if (com.model == UNREST)
            x[k] = x[k + 3] = x[k + 8] = 0.8 + 0.2*rndu();
         if (com.model == REV) {
            nkappasets = com.nrate / 5;
            for (j = 0; j < nkappasets; j++)
               x[k + j * 5] = 1 + 0.1*(j + 1); /*  0.8+0.2*rndu();  */
         }
      }
   else
      for (i = 0; i < com.nrate; i++)
         x[com.ntime + com.nrgene + i] = 0.02 + com.kappa*(.8 + 0.4*rndu());

   for (i = 0; i < com.npi*n31pi; i++)  /* initials for transformed pi's */
      x[com.ntime + com.nrgene + com.nrate + i] = rndu()*.2;
   com.np = k = com.ntime + com.nrgene + com.nrate + com.npi*n31pi;

   if (com.alpha || com.nparK) {
      for (i = 0; i < com.nalpha; i++)
         x[k++] = com.alpha = 0.01 + com.alpha*(.9 + 0.2*rndu());

      if (!com.fix_rho)   x[k++] = com.rho = 0.01 + com.rho*(.8 + 0.4*rndu());

      if (com.rho)
         AutodGamma(com.MK, com.freqK, com.rK, &t, com.alpha, com.rho, K);
      else
         DiscreteGamma(com.freqK, com.rK, com.alpha, com.alpha, K, DGammaUseMedian);

      if (com.nparK) {
         xtoy(com.rK, x + k, K - 1);
         k += K - 1;
      }
      switch (com.nparK) {
      case (2):                            /* rK & fK */
         /* zero(x+k, K-1); */
         for (i = 0; i < K - 1; i++) x[k + i] = 0.1*(0.5 - rndu());
         k += K - 1;
         break;
      case (3):                            /* rK & MK (double stochastic) */
         /* zero(x+k, (K-1)*(K-1));  */
         for (i = 0; i < (K - 1)*(K - 1); i++) x[k + i] = 0.1*(0.5 - rndu());
         k += (K - 1)*(K - 1);
         break;
      case (4):                            /* rK & MK */
         /* zero(x+k, K*(K-1)); */
         for (i = 0; i < K*(K - 1); i++) x[k + i] = 0.1*(0.5 - rndu());
         k += K*(K - 1);
         break;
      }
      com.np = k;
   }

   if (com.fix_blength == -1)
      for (i = 0; i < com.np; i++)  x[i] = (i < com.ntime ? .1 + rndu() : .5 + rndu());

   if (finitials) readx(x, fromfile);
   else    *fromfile = 0;

   return (0);
}




int SetParameters(double x[])
{
   /* This sets parameters in com., nodes[], etc and is called before lfun() etc.
      Iinitialize U, V, Root etc, if necessary.
      For nhomo models (nhomo=1,3,4)
         x[] has frequencies if (LASTROUND==1) or exp(pi)/(1+SUM(exp[pi])) if otherwise
   */
   int i, j, k, K = com.ncatG, status = 0, n31pi = (com.model == T92 ? 1 : 3);
   int nkappa[] = { 0, 1, 0, 1, 1, 1, 2, 5, 11, -1, -1 };
   double k1 = com.kappa, k2 = com.kappa, t, space[NCATG*(NCATG + 1)];

   if (com.clock >= 5) return(0);
   if (com.fix_blength != 2) SetBranch(x);

   for (i = 0; i < com.nrgene; i++) com.rgene[i + 1] = x[com.ntime + i];
   if (com.clock && com.clock < 5 && AbsoluteRate)
      com.rgene[0] = x[0]; /* so that rgene are absolute rates */

   if (!com.fix_kappa && com.model <= TN93 && com.clock < 5) {
      com.kappa = k1 = k2 = x[com.ntime + com.nrgene];
      if (com.model == TN93) k2 = x[com.ntime + com.nrgene + 1];
   }
   if (com.nhomo == 1) {
      k = com.ntime + com.nrgene + com.nrate;
      if (com.model == T92) {  /* T=A, C=G */
         com.pi[0] = com.pi[2] = (1 - x[k]) / 2;  com.pi[1] = com.pi[3] = x[k] / 2;
      }
      else {
         if (!LASTROUND) f_and_x(x + k, com.pi, 4, 0, 0);
         else            xtoy(x + k, com.pi, 3);
         com.pi[3] = 1 - sum(com.pi, 3);
      }
      if (com.model <= TN93)
         eigenTN93(com.model, k1, k2, com.pi, &nR, Root, Cijk);
   }
   else if (com.nhomo == 2)
      for (i = 0, k = com.ntime + com.nrgene; i < tree.nbranch; i++)
         nodes[tree.branches[i][1]].pkappa = x + k + i;

   if (com.model <= TN93 && com.nhomo == 0 && com.Mgene <= 1)
      RootTN93(com.model, k1, k2, com.pi, &t, Root);
   else if (com.nhomo >= 3) {
      for (i = 0, k = com.ntime + com.nrgene; i < tree.nbranch; i++) {  /* Qrate for branches */
         if (com.fix_kappa == 0)       /* each branch has its Qrate */
            nodes[tree.branches[i][1]].pkappa = x + k + i*nkappa[com.model]; 
         else if (com.fix_kappa == 1)  /* same Qrate for all branches */
            nodes[tree.branches[i][1]].pkappa = x + k;
         else if (com.fix_kappa == 2)  /* sets of Qrate */
            nodes[tree.branches[i][1]].pkappa = x + k + (int)nodes[tree.branches[i][1]].label*nkappa[com.model];
      }
      k += com.nrate;

      for (i = 0; i < tree.nnode; i++) {          /* pi for nodes */
         j = (com.nhomo == 4 ? i : (int)nodes[i].label);
         if (com.model == T92) {
            nodes[i].pi[0] = nodes[i].pi[2] = (1 - x[k + j]) / 2;  /* TA */
            nodes[i].pi[1] = nodes[i].pi[3] = x[k + j] / 2;        /* CG */
         }
         else {
            if (!LASTROUND) f_and_x(x + k + j * 3, nodes[i].pi, 4, 0, 0);
            else            xtoy(x + k + j * 3, nodes[i].pi, 3);
            nodes[i].pi[3] = 1 - sum(nodes[i].pi, 3);
         }
      }
      xtoy(nodes[tree.root].pi, com.pi, 4);
   }
   else if ((com.model == REV || com.model == REVu) && com.Mgene <= 1)
      eigenQREVbase(NULL, PMat, x + com.ntime + com.nrgene, com.pi, &nR, Root, Cijk);
   /*
   else if ((com.model==UNREST || com.model==UNRESTu) && com.Mgene<=1)
      eigenQunrest (NULL, x+com.ntime+com.nrgene,com.pi,&nR,cRoot,cU,cV);
   */

   if (com.nparK == 0 && (com.alpha == 0 || com.fix_alpha*com.fix_rho == 1))
      return(status);
   if (com.nalpha > 1) return (status);
   k = com.ntime + com.nrate + com.nrgene + com.npi*n31pi;
   if (!com.fix_alpha) {
      com.alpha = x[k++];
      if (com.fix_rho)
         DiscreteGamma(com.freqK, com.rK, com.alpha, com.alpha, K, DGammaUseMedian);
   }
   if (!com.fix_rho) {
      com.rho = x[k++];
      AutodGamma(com.MK, com.freqK, com.rK, &t, com.alpha, com.rho, K);
   }
   if (com.nparK == 0) return(status);

   /* nparK models */
   xtoy(x + k, com.rK, K - 1);

   if (com.nparK == 2) {
      if (!LASTROUND)  f_and_x(x + k + K - 1, com.freqK, K, 0, 0);
      else             xtoy(x + k + K - 1, com.freqK, K - 1);
      com.freqK[K - 1] = 1 - sum(com.freqK, K - 1);
   }
   else if (com.nparK == 3) {   /* rK & MK (double stochastic matrix) */
      for (i = 0, k += K - 1; i < K - 1; k += K - 1, i++) {
         if (!LASTROUND) f_and_x(x + k, com.MK + i*K, K, 0, 0);
         else            xtoy(x + k, com.MK + i*K, K - 1);
         com.MK[i*K + K - 1] = 1 - sum(com.MK + i*K, K - 1);
      }
      for (j = 0; j < K; j++) {
         for (i = 0, com.MK[(K - 1)*K + j] = 1; i < K - 1; i++)
            com.MK[(K - 1)*K + j] -= com.MK[i*K + j];
         if (com.MK[(K - 1)*K + j] < 0)
            printf("SetPar: MK[K-1][j]=%.5f<0\n", com.MK[(K - 1)*K + j]);
      }
   }
   else if (com.nparK == 4) { /* rK & MK */
      for (i = 0, k += K - 1; i < K; k += K - 1, i++) {
         if (!LASTROUND) f_and_x(x + k, com.MK + i*K, K, 0, 0);
         else            xtoy(x + k, com.MK + i*K, K - 1);
         com.MK[i*K + K - 1] = 1 - sum(com.MK + i*K, K - 1);
      }
      PtoPi(com.MK, com.freqK, K, space);
   }
   com.rK[K - 1] = (1 - innerp(com.freqK, com.rK, K - 1)) / com.freqK[K - 1];
   return (status);
}


int SetPGene(int igene, int _pi, int _UVRoot, int _alpha, double x[])
{
   /* xcom[] does not contain time parameters
      Note that com.piG[][] have been homogeneized if (com.Mgene==3)
   */
   int nr[] = { 0, 1, 0, 1, 1, 1, 2, 5, 11 };
   int k = com.nrgene + (com.Mgene >= 3)*igene*nr[com.model];
   double *xcom = x + com.ntime;
   double ka1 = xcom[k], ka2 = (com.model == TN93 ? xcom[k + 1] : -1);

   if (com.Mgene == 2 && com.fix_kappa) ka1 = ka2 = com.kappa;

   if (_pi) {
      xtoy(com.piG[igene], com.pi, 4);
   }
   if (_UVRoot) {
      if (com.model == K80) com.kappa = ka1;
      else if (com.model <= TN93)
         eigenTN93(com.model, ka1, ka2, com.pi, &nR, Root, Cijk);
      else if (com.model == REV || com.model == REVu)
         eigenQREVbase(NULL, PMat, xcom + k, com.pi, &nR, Root, Cijk);
   }
   if (_alpha) {
      com.alpha = xcom[com.nrgene + com.nrate + com.npi + igene]; /* check?? */
      DiscreteGamma(com.freqK, com.rK, com.alpha, com.alpha, com.ncatG, DGammaUseMedian);
   }
   return(0);
}


int SetxBound(int np, double xb[][2])
{
   /* sets lower and upper bounds for variables during iteration
   */
   int i, j, k = 0, nf = 0, n31pi = (com.model == T92 ? 1 : 3);
   double rateb[] = { 1e-5, 999 }, rgeneb[] = { 1e-4, 999 };
   double alphab[] = { .005, 999 }, rhob[] = { -0.2, 0.99 }, pb[] = { .00001, .99999 };
   double fb[] = { -19,9 }; /* transformed freqs.*/

   SetxBoundTimes(xb);
   for (i = com.ntime; i < np; i++)
      for (j = 0; j < 2; j++)  xb[i][j] = rateb[j];

   for (i = 0; i < com.nrgene; i++)
      for (j = 0; j < 2; j++)  xb[com.ntime + i][j] = rgeneb[j];
   for (i = 0; i < com.nrate; i++)
      for (j = 0; j < 2; j++)  xb[com.ntime + com.nrgene + i][j] = rateb[j];
   k = com.ntime + com.nrgene + com.nrate;
   for (i = 0; i < com.npi*n31pi; i++) {
      xb[k][0] = (com.model == T92 ? pb[0] : fb[0]);
      xb[k++][1] = (com.model == T92 ? pb[1] : fb[1]);
   }
   for (i = 0; i < com.nalpha; i++, k++)  
      for (j = 0; j < 2; j++) 
         xb[k][j] = alphab[j];
   if (!com.fix_rho)   
      for (j = 0; j < 2; j++) xb[np - 1][j] = rhob[j];
   if (com.nparK) {
      for (i = 0; i < com.ncatG - 1; i++) { xb[k][0] = rateb[0]; xb[k++][1] = rateb[1]; }
      if (com.nparK == 2) nf = com.ncatG - 1;
      else if (com.nparK == 3) nf = (com.ncatG - 1)*(com.ncatG - 1);
      else if (com.nparK == 4) nf = (com.ncatG - 1)*com.ncatG;
      for (i = 0; i < nf; i++) { xb[k][0] = fb[0]; xb[k++][1] = fb[1]; }
   }
   if (noisy > 2 && np < 50) {
      printf("\nBounds (np=%d):\n", np);
      for (i = 0; i < np; i++) printf(" %10.6f", xb[i][0]);
      printf("\n");
      for (i = 0; i < np; i++) printf(" %10.6f", xb[i][1]);
      printf("\n");
   }
   return(0);
}

int testx(double x[], int np)
{
   /* This is used for LS branch length estimation by nls2, called only if(clock==0)
   */
   int i;
   double tb[] = { .4e-6, 99 };

   for (i = 0; i < com.ntime; i++)
      if (x[i]<tb[0] || x[i]>tb[1])
         return (-1);
   return (0);
}



int ConditionalPNode(int inode, int igene, double x[])
{
   int n = com.ncode, i, j, k, h, ison, pos0 = com.posG[igene], pos1 = com.posG[igene + 1];
   double t;

   for (i = 0; i < nodes[inode].nson; i++)
      if (nodes[nodes[inode].sons[i]].nson > 0 && !com.oldconP[nodes[inode].sons[i]])
         ConditionalPNode(nodes[inode].sons[i], igene, x);
   if (inode < com.ns) {  /* young ancestor */
      for (h = pos0*n; h < pos1*n; h++)
         nodes[inode].conP[h] = 0;
   }
   else
      for (h = pos0*n; h < pos1*n; h++)
         nodes[inode].conP[h] = 1;
   if (com.cleandata && inode < com.ns) /* young ancestor */
      for (h = pos0; h < pos1; h++)
         nodes[inode].conP[h*n + com.z[inode][h]] = 1;

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      t = nodes[ison].branch*_rateSite;

      if (com.clock < 5) {
         if (com.clock)  t *= GetBranchRate(igene, (int)nodes[ison].label, x, NULL);
         else            t *= com.rgene[igene];
      }
      GetPMatBranch(PMat, x, t, ison);

      if (nodes[ison].nson < 1 && com.cleandata) {        /* tip && clean */
         for (h = pos0; h < pos1; h++)
            for (j = 0; j < n; j++)
               nodes[inode].conP[h*n + j] *= PMat[j*n + com.z[ison][h]];
      }
      else if (nodes[ison].nson < 1 && !com.cleandata) {  /* tip & unclean */
         for (h = pos0; h < pos1; h++)
            for (j = 0; j < n; j++) {
               for (k = 0, t = 0; k < nChara[(int)com.z[ison][h]]; k++)
                  t += PMat[j*n + CharaMap[(int)com.z[ison][h]][k]];
               nodes[inode].conP[h*n + j] *= t;
            }
      }
      else {
         for (h = pos0; h < pos1; h++)
            for (j = 0; j < n; j++) {
               for (k = 0, t = 0; k < n; k++)
                  t += PMat[j*n + k] * nodes[ison].conP[h*n + k];
               nodes[inode].conP[h*n + j] *= t;
            }
      }
   }        /*  for (ison)  */
   if (com.NnodeScale && com.nodeScale[inode])  NodeScale(inode, pos0, pos1);
   return (0);
}

int PMatCijk(double P[], double t)
{
/* P(t)ij = SUM Cijk * exp{Root*t}
*/
   int i, j, k, n = com.ncode, nr = nR;
   double exptm1[5] = {0};

   memset(P, 0, n*n * sizeof(double));
   for (k = 1; k < nr; k++) exptm1[k] = expm1(t*Root[k]);
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         for (k = 0; k < nr; k++)
            P[i*n + j] += Cijk[i*n*nr + j*nr + k] * exptm1[k];
      }
      P[i*n + i] ++;
   }
   return (0);
}



#ifdef UNDEFINED

int CollapsSite(FILE *fout, int nsep, int ns, int *ncat, int SiteCat[])
{
   int j, k, it, h, b[NS], ndiff, n1, bit1;
   /* n1: # of 1's   ...  bit1: the 1st nonzero bit */

   *ncat = 5 + (1 << (ns - 1)) - 1 + nsep * 11;
   if (fout) fprintf(fout, "\n# cat:%5d  # sep:%5d\n\n", *ncat, nsep);

   for (h = 0; h < (1 << 2 * ns); h++) {
      for (j = 0, it = h; j < ns; b[ns - 1 - j] = it % 4, it /= 4, j++);
      for (j = 1, ndiff = 0; j < ns; j++) {
         for (k = 0; k < j; k++) if (b[j] == b[k]) break;
         if (k == j) ndiff++;
      }
      switch (ndiff) {
      default: SiteCat[h] = 0;      break;
      case (0): SiteCat[h] = b[0] + 1; break;
      case (1):
         for (j = 1, it = 0, n1 = 0, bit1 = 0; j < ns; j++) {
            k = (b[j] != b[0]);
            it = it * 2 + k;
            n1 += k;
            if (bit1 == 0 && k) bit1 = ns - 1 - j;
         }
         it = 5 + it - 1;
         if (nsep == 0) { SiteCat[h] = it; break; }

         SiteCat[h] = it + min2(bit1 + 1, nsep) * 11;
         if (n1 == 1 && bit1 < nsep) {
            SiteCat[h] -= 11;
            SiteCat[h] += (b[0] * 4 + b[ns - 1 - bit1] - b[0] - (b[0] <= b[ns - 1 - bit1]));
         }
         break;
      }
      if (fout) {
         for (j = 0; j < ns; k++) fprintf(fout, "%1c", BASEs[b[j]]);
         fprintf(fout, "%5d    ", SiteCat[h]);
         if (h % 4 == 3) fprintf(fout, "\n");
      }
   }
   return (0);
}

int GetPexpML(double x[], int ncat, int SiteCat[], double pexp[])
{
   int  j, it, h, nodeb[NNODE];
   int  isum, nsum = 1 << (2 * (tree.nbranch - com.ns + 1));
   double fh, y, Pt[NBRANCH][16];

   if (com.ngene > 1 || com.nhomo || com.alpha || com.nparK || com.model > REV)
      error2("Pexp()");
   SetParameters(x);
   for(j=0; j<tree.nbranch; j++)
      PMatCijk(Pt[j], nodes[tree.branches[j][1]].branch);

   for (h = 0; h < (1 << 2 * com.ns); h++) {
      if (SiteCat[h] == 0) continue;
      for (j = 0, it = h; j < com.ns; nodeb[com.ns - 1 - j] = it & 3, it >>= 2, j++);
      for (isum = 0, fh = 0; isum < nsum; isum++) {
         for (j = 0, it = isum; j < tree.nbranch - com.ns + 1; j++)
         {
            nodeb[com.ns + j] = it % 4; it /= 4;
         }
         for (j = 0, y = com.pi[nodeb[tree.root]]; j < tree.nbranch; j++)
            y *= Pt[j][nodeb[tree.branches[j][0]] * 4 + nodeb[tree.branches[j][1]]];
         fh += y;
      }
      pexp[SiteCat[h]] += fh;
   }
   pexp[0] = 1 - sum(pexp + 1, ncat - 1);
   return (0);
}


int TestModel(FILE *fout, double x[], int nsep, double space[])
{
   /* test of models, using com.
   */
   int j, h, it, ls = com.ls, ncat, *SiteCat;
   double *pexp = space, *nobs, lmax0, lnL0, X2, ef, de;

   SiteCat = (int*)malloc((1 << 2 * com.ns) * sizeof(int));
   if (SiteCat == NULL)  error2("oom");
   CollapsSite(F0, nsep, com.ns, &ncat, SiteCat);
   fprintf(fout, "\n\nAppr. test of model.. ncat%6d  nsep%6d\n", ncat, nsep);

   nobs = pexp + ncat;
   zero(pexp, 2 * ncat);
   /* nobs */
   for (h = 0; h < com.npatt; h++) {
      for (j = 0, it = 0; j < com.ns; j++) it = it * 4 + (com.z[j][h] - 1);
      nobs[SiteCat[it]] += com.fpatt[h];
   }
   GetPexpML(x, ncat, SiteCat, pexp);

   for (h = 0, lnL0 = 0, X2 = 0, lmax0 = -(double)ls*log((double)ls); h < ncat; h++) {
      if (nobs[h] > 1) {
         lmax0 += nobs[h] * log(nobs[h]);
         lnL0 += nobs[h] * log(pexp[h]);
      }
      ef = com.ls*pexp[h];
      de = square(nobs[h] - ef) / ef;
      X2 += de;
      fprintf(fout, "\nCat #%3d%9.0f%9.2f%9.2f", h, nobs[h], ef, de);
   }
   fprintf(fout, "\n\nlmax0:%12.4f  D2:%12.4f   X2:%12.4f\n",
      lmax0, 2 * (lmax0 - lnL0), X2);
   free(SiteCat);

   return (0);
}


#endif


int OldDistributions(int inode, double AncientFreqs[])
{
   /* reconstruct nucleotide frequencies at and down inode for nonhomogeneous models com.nhomo==3 or 4.
      AncientFreqs[tree.nnode*4]
   */
   int i, n = 4;

   if (com.alpha || com.model > REV) {
      puts("OldDistributions() does not run when alpha > 0 or model >= TN93");
      return(-1);
   }
   if (inode == tree.root)
      xtoy(nodes[inode].pi, AncientFreqs + inode*n, n);
   else {
      if (com.nhomo > 2 && com.model <= TN93)
         eigenTN93(com.model, *nodes[inode].pkappa, *(nodes[inode].pkappa + 1), nodes[inode].pi, &nR, Root, Cijk);
      else if (com.nhomo > 2 && com.model == REV)
         eigenQREVbase(NULL, PMat, nodes[inode].pkappa, nodes[inode].pi, &nR, Root, Cijk);

      PMatCijk(PMat, nodes[inode].branch);
      matby(AncientFreqs + nodes[inode].father*n, PMat, AncientFreqs + inode*n, 1, n, n);
   }
   for (i = 0; i < nodes[inode].nson; i++)
      OldDistributions(nodes[inode].sons[i], AncientFreqs);
   return (0);
}





/* problems and notes

   (1) AdG model: generation of com.MK[K*K] is not independent
       of com.alpha or DiscreteGamma().

non-homogeneous process models:
  nhomo            fix_kappa         models

  0 (1 pi given)   0 (solve 1 kappa)   JC69, K80, F81, F84, HKY85,
                                       REV(5), UNREST(11)
                   1 (1 kappa given)   K80, F84, HKY85

  1 (solve 1 pi)   0 (as above)        F84, HKY85, REV(5)
                   1                   F81(0), F84, HKY85

  2 (b kappa's)    ?                   K80, F84, HKY85

  3 (ns+2 pi)      0 (solve 1 kappa)   F84 & HKY85
                   1 (nbranch kappa)   F84 & HKY85

  4 (nnode pi)     0,1  (as above)


space-time process models:
  nparK     fix_alpha       fix_rho        parameters
  0         (0,1)          (0,1)           alpha & rho for AdG

  1         set to 1       set to 1        rK[]        (freqK=1/K) (K-1)
  2         set to 1       set to 1        rK[] & freqK[]          2(K-1)
  3         set to 1       set to 1        rK[] & MK[] (freqK=1/K) K(K-1)
  4         set to 1       set to 1        rK[] & MK[]             K*K-1


Local clock models
   parameters under local clock (com.clock=2)
      com.ntime = (#ancestral nodes) - 1 + (#rates) - 1
   Parameters include (ns-1) node times t[] and rates for branches.
      x[0]: t0*r0
      x[1..(nid-1)]: ti/t0 ancestral node times expressed as ratios
      x[ns-1 .. ntime-1]: rates for branches

*/
