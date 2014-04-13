/* mcmctree.c 

   Markov chain Monte Carlo on trees (Bayes phylogenetic analysis)
   
                   Ziheng YANG, December 2002

   cc -o mcmctree -march=pentiumpro -mcpu=pentiumpro -O4 -funroll-loops -fomit-frame-pointer -finline-functions mcmctree.c tools.c -lm
   cc -o mcmctree -march=athlon -mcpu=athlon -O4 -funroll-loops -fomit-frame-pointer -finline-functions mcmctree.c tools.c -lm
           cc -o mcmctree -fast mcmctree.c tools.c -lm

               cl -O2 -W4 mcmctree.c tools.c
                  mcmctree <ControlFileName>
*/

#include "paml.h"

#define NS            40
#define NBRANCH      (NS*2-2)
#define NNODE        (NS*2-1)
#define NGENE         10          /* used for gnodes[NGENE] */
#define LSPNAME       30
#define NCODE         4
#define NCATG         8

extern int noisy, NFunCall;
extern char BASEs[];

int GetOptions(char *ctlf);
int ReadTreeSeqs(FILE*fout);
int GetMem(void);
void FreeMem(void);
int UseLocus(int locus, int copyconP, int setSeqName);
int AcceptRejectLocus(int locus, int accept);
void switchconPin(void);
int SetPGene(int igene, int _pi, int _UVRoot, int _alpha, double xcom[]);
int DownSptreeSetTime(int inode);
int GetRgenes (void);
void getSinvDetS (double space[]);
int GetInitials(double space[]);
int GetGtree(int locus);
void printGtree(int printBlength);
int SetParameters(double x[]);
int getQfactor(int locus, double kappa);
int ConditionalPNode(int inode, int igene, double x[]);
double BranchLength (double t, double y0, double yt, double sigma2);
double lnpData(double lnpDi[]);
double lnpD_locus(int locus);
double lnpriorTimes(void);
double lnpriorRates(void);
void copySptree(void);
void printSptree(void);
int MCMC(FILE* fout);
int SetOldCondP (int spnode);
double UpdateTimes(double *lnL, double finetune);
double UpdateRates(double *lnL, double finetune);
double UpdateRatesClock(double *lnL, double finetune);
double UpdateParameters(double *lnL, double finetune);
double mixing(double finetune);

struct CommonInfo {
   char *z[NS], *spname[NS], seqf[96],outf[96],treef[96];
   char oldconP[NNODE];       /* update conP for node? (0 yes; 1 no) */
   int seqtype, ns, ls, ngene, posG[2],lgene[1], *pose, npatt;
   int np, ncode, ntime, nrate, nrgene, nalpha, npi, ncatG, print, seed;
   int cleandata;
   int model, clock, fix_alpha, fix_kappa, fix_rgene, Mgene;
   double *fpatt, kappa, alpha;
   double rgene[NGENE],piG[NGENE][NCODE];  /* not used */
   double (*plfun)(double x[],int np), freqK[NCATG], rK[NCATG], *conP, *fhK;
   double pi[NCODE], pi_sqrt[NCODE], Qfactor;
   int    sconP, curconP;                /* curconP = 0 or 1 */
   double *conPin[2], *conP0, space[10000];  /* space used for S^-1 and |S| */
   int npi0;
   int conPSiteClass, readpattf;         /* not used */
   int nnodeScale;
   char nodeScale[NNODE];    /* nScale[ns-1] for interior nodes */
   double *nodeScaleF;       /* nScaleF[npatt] for scale factors */
}  com;
struct TREEB {
   int  nbranch, nnode, root, branches[NBRANCH][2];
}  tree;
struct TREEN { /* ipop is node number in species tree.  label is rate. */
   int father, nson, sons[2], ibranch, ipop;
   double branch, age, label, *conP;
}  *nodes, *gnodes[NGENE], nodes_t[2*NS-1];

/* nodes_t[] is working space.  nodes is a pointer and moves around.  
   gnodes[] holds the gene trees, subtrees constructed from the master species 
   tree.  Each locus has a full set of rates (lnrates) for all branches on the 
   master tree, held in sptree.nodes[].lnrates.  Branch lengths in the gene 
   tree are calculated by using those rates and the divergence dates.
*/


/* for sptree.nodes[].fossil: lower, upper, bounds, gamma, inverse-gamma */
enum {FREE_F=0, LOWER_F, UPPER_F, BOUND_F, GAMMA_F, IGAMMA_F} FOSSIL_FLAGS;
char *fossils[]={" ", "L", "U", "B", "G", "IG"}; 
enum {GAMMA, IGAMMA} DISTRIBUTIONS;


struct SPECIESTREE {
   int nbranch, nnode, root, nspecies, nfossil;
   struct TREESPN {
      char name[LSPNAME+1];  /* fossil: 0, 1, 2, 3 */
      int father, nson, sons[2], fossil;
      double age, tfossil[2];   /* lower and upper bounds or alpha & beta */
      double *lnrates;              /* log rates for loci */
   } nodes[2*NS-1];
}  sptree;
/* all trees are binary & rooted, with ancestors unknown. */


struct DATA { /* locus-specific data and tree information */
   int ns[NGENE], ls[NGENE], npatt[NGENE], ngene, lgene[NGENE];
   int root[NGENE+1], conP_offset[NGENE], BlengthMethod;
   char *z[NGENE][NS], cleandata[NGENE];
   double *fpatt[NGENE], lnpT, lnpR, lnpDi[NGENE];
   double Qfactor[NGENE], pi[NGENE][NCODE];
   double rgene[NGENE], kappa[NGENE];
   double birth, death, sampling;
   double rgenegamma[2], kappagamma[2], sigma2gamma[2];
   /* correl[g*g] for rho and sigma2, S^-1[g][g], |S| */
   double *correl, *Sinv, detS;
}  data;

struct MCMCPARAMETERS {
   int burnin, nsample, sampfreq, usedata, saveconP, print;
   double finetune[5];
}  mcmc; /* control parameters */


char *models[]={"JC69","K80","F81","F84","HKY85","T92","TN93","REV"};
enum {JC69, K80, F81, F84, HKY85, T92, TN93, REV} MODELS;

int nR=4;
double PMat[16], Cijk[64], Root[4];
double _rateSite=1, OLDAGE=99;
int LASTROUND=0, debug=0, testlnL=0, NPMat=0; /* no use for this */

/* char *ratef="rates";  */ /* not used in this program */

#define MCMCTREE  1
#include "treesub.c"

int main (int argc, char *argv[])
{
   char ctlf[]="mcmctree.ctl";
   FILE  *fout;

   error2("mcmctree is disabled right now");

   noisy=3;
   com.alpha=0.;     com.ncatG=1;
   com.ncode=4;      com.clock=1;

   printf("MCMCTREE in %s\n", VerStr);
   if (argc>1) 
      { strcpy(ctlf, argv[1]); printf ("\nctlfile reset to %s.\n", ctlf); }

   data.birth=2; data.death=1; data.sampling=0.05; 
   com.cleandata=0; mcmc.usedata=1;

   GetOptions (ctlf);
   fout=gfopen(com.outf,"w");

   ReadTreeSeqs(fout);

   GetMem();
   MCMC(fout);
   FreeMem();
   return (0);
}



extern double gamma_pexp[];

int ReadTreeSeqs (FILE*fout)
{
/* This reads the combined species tree, the fossil calibration information, 
   and sequence data at each locus.  sptree.nodes[].tfossil[] has tL, tU for 
   bounds or alpha and beta for the gamma prior.  
   This also constructs the gene tree at each locus, by pruning the master 
   species tree..
*/
   FILE *fseq, *ftree;
   int i,j, locus, lspname=LSPNAME, clean0=com.cleandata;
   double *pe=gamma_pexp, modetiles[3]; /* mode, 2.5%, 97.5% */

   ftree=gfopen(com.treef,"r");
  
   /* read master species tree and process fossil calibration info */
   fscanf(ftree, "%d%d", &sptree.nspecies, &i);
   com.ns=sptree.nspecies;
   /* to read master species names into sptree.nodes[].name */
   if(noisy) puts("Reading master tree.");
   FOR(j,sptree.nspecies) com.spname[j]=sptree.nodes[j].name;
   nodes=nodes_t;
   ReadaTreeN(ftree,&i,&i,1,1);
   OutaTreeN(F0,0,0); FPN(F0);
   OutaTreeN(F0,1,0); FPN(F0);
   /* copy master tree into sptree */
   sptree.nnode=tree.nnode;  sptree.nbranch=tree.nbranch; 
   sptree.root=tree.root;    sptree.nfossil=0;
   for(i=0; i<sptree.nspecies*2-1; i++) {
      sptree.nodes[i].father=nodes[i].father;
      sptree.nodes[i].nson=nodes[i].nson;
      if(nodes[i].nson!=0 && nodes[i].nson!=2) 
         error2("master tree has to be binary.");
      for(j=0;j<sptree.nodes[i].nson;j++) sptree.nodes[i].sons[j]=nodes[i].sons[j];

      /* collect fossil information */
      sptree.nodes[i].fossil=0;
      sptree.nodes[i].tfossil[0]=nodes[i].branch; /* ">": Lower bound */
      sptree.nodes[i].tfossil[1]=nodes[i].label;  /* "<": Upper bound */

      if(nodes[i].branch && nodes[i].label) {
         modetiles[0]=nodes[i].age;   /* mode */
         modetiles[1]=nodes[i].branch;    /* ">": Lower bound */
         modetiles[2]=nodes[i].label;     /* "<": Upper bound */
         if(modetiles[0]==0)
            sptree.nodes[i].fossil = BOUND_F;
         else {
            getab_gamma(&sptree.nodes[i].tfossil[0], &sptree.nodes[i].tfossil[1], modetiles);
            printf("\na = %8.6f  b = %8.6f  CDFs (%8.6f, %8.6f) %8.6f\n", sptree.nodes[i].tfossil[0],sptree.nodes[i].tfossil[1], pe[0],pe[1], pe[1]-pe[0]);
            sptree.nodes[i].fossil = GAMMA_F; 
         }
         sptree.nfossil++;
      }
      else if(nodes[i].branch) 
         { sptree.nodes[i].fossil = LOWER_F; sptree.nfossil++; }
      else if(nodes[i].label) 
         { sptree.nodes[i].fossil = UPPER_F; sptree.nfossil++; }

      nodes[i].branch=nodes[i].label=0;
   }

   /* read sequences at each locus, construct gene tree by pruning sptree */
   fseq=gfopen(com.seqf,"r");
   fscanf(fseq, "%d", &data.ngene);
   printf("\nReading sequence data..  %d loci\n", data.ngene);
   for(locus=0; locus<data.ngene; locus++) {
      fprintf(fout, "\n\n*** Locus %d ***\n", locus+1);
      printf("\n\n*** Locus %d ***\n", locus+1);

      com.cleandata=clean0;
      FOR(j,sptree.nspecies) com.spname[j]=NULL; /* points to nowhere */
      ReadSeq (NULL, fseq);                      /* allocates com.spname[] */
      Initialize (fout);
      if(com.model==JC69) PatternJC69like(fout);

      xtoy(com.pi, data.pi[locus], com.ncode);

      data.cleandata[locus] = (char)com.cleandata;
      data.ns[locus]=com.ns;  
      data.ls[locus]=com.ls;  
      data.npatt[locus]=com.npatt;
      data.fpatt[locus]=com.fpatt; com.fpatt=NULL;
      for(i=0; i<com.ns; i++) { data.z[locus][i]=com.z[i]; com.z[i]=NULL; }

      printf("%3d patterns, %s\n", com.npatt,(com.cleandata?"clean":"messy"));
      GetGtree(locus);      /* free com.spname[] */
   }
   for(i=0,com.cleandata=1; i<data.ngene; i++) 
      if(!data.cleandata[i]) com.cleandata=0;

   fclose(ftree); fclose(fseq);
   return(0);
}


int getQfactor (int locus, double kappa)
{
   double T, C, A, G, Y, R, *pi=data.pi[locus];
   double Qrates[2];

   T=pi[0]; C=pi[1]; A=pi[2]; G=pi[3]; Y=T+C; R=A+G;
   Qrates[0]=Qrates[1]=com.kappa;
   if (com.model==F84) { 
      Qrates[0]=1+com.kappa/Y;   /* kappa1 */
      Qrates[1]=1+com.kappa/R;   /* kappa2 */
   }
   data.Qfactor[locus]=1/(2*T*C*Qrates[0] + 2*A*G*Qrates[1] + 2*Y*R);

   return(0);
}



int GetGtree (int locus)
{
/* construct the gene tree at locus by pruning tips in the master species 
   tree.  com.spname[] have names of species at the current locus and 
   the routine use them to compare with sptree.nodes[].name to decide which 
   species to keep for the locus.  See GetSubTreeN() for more info.
*/
   int ns=data.ns[locus], i,j, ipop[NS], keep[NS], newnodeNO[2*NS-1];

   FOR(j,sptree.nspecies) keep[j]=0;
   FOR(i,ns) {
      FOR(j,sptree.nspecies) 
         if(!strcmp(com.spname[i], sptree.nodes[j].name)) break;
      if(j==sptree.nspecies) {
         printf("species %s not found in master tree\n", com.spname[i]);
         exit(-1);
      }
      keep[j]=i+1; ipop[i]=j;
      free(com.spname[i]);
   }
   /* copy master species tree and then prune it. */
   copySptree();
   GetSubTreeN(keep, newnodeNO);
   com.ns=ns;

   FOR(i,sptree.nnode)  
      if(newnodeNO[i]!=-1) nodes[newnodeNO[i]].ipop=i;
   /* if(debug) */ 
      printGtree(0);

   gnodes[locus]=(struct TREEN*)malloc((ns*2-1)*sizeof(struct TREEN));
   if(gnodes[locus]==NULL) error2("oom gtree");
   memcpy(gnodes[locus], nodes, (ns*2-1)*sizeof(struct TREEN));
   data.root[locus]=tree.root;

   return(0);
}


void printGtree (int printBlength)
{
   int i,j;

   FOR(i,com.ns) com.spname[i]=sptree.nodes[nodes[i].ipop].name;
   FOR (i,tree.nnode) if(i!=tree.root) 
      nodes[i].branch=nodes[nodes[i].father].age-nodes[i].age;
   printf("\nns = %d  nnode = %d", com.ns, tree.nnode);
   printf("\n%7s%7s %8s %7s%7s","father","node","(ipop)","nson:","sons");
   FOR (i, tree.nnode) {
      printf ("\n%7d%7d   (%2d) %7d  ",
         nodes[i].father, i, nodes[i].ipop, nodes[i].nson);
      FOR(j, nodes[i].nson) printf (" %2d", nodes[i].sons[j]);
   }
   FPN(F0); OutaTreeN(F0,0,0); FPN(F0); OutaTreeN(F0,1,0); FPN(F0); 
   if(printBlength) { OutaTreeN(F0,1,1); FPN(F0); }
}


void copySptree (void)
{
/* This copies sptree into nodes = nodes_t, for printing or editing
*/
   int i,j;

   nodes=nodes_t;
   com.ns=sptree.nspecies;   tree.root=sptree.root;
   tree.nnode=sptree.nnode;  tree.nbranch=sptree.nbranch; 
   for(i=0; i<com.ns*2-1; i++) {
      if(i<com.ns) com.spname[i]=sptree.nodes[i].name;
      nodes[i].father=sptree.nodes[i].father;
      nodes[i].nson=sptree.nodes[i].nson;
      for(j=0;j<nodes[i].nson;j++) nodes[i].sons[j]=sptree.nodes[i].sons[j];
      nodes[i].age=sptree.nodes[i].age;
      if(i!=tree.root) 
         nodes[i].branch=sptree.nodes[nodes[i].father].age-sptree.nodes[i].age;
   }
}

void printSptree (void)
{
   int i;

   printf("\n************\nSpecies tree\nns = %d  nnode = %d", sptree.nspecies, sptree.nnode);
   printf("\n%7s%7s  %-8s %12s %12s%16s\n","father","node","name","time","fossil","sons");
   FOR (i, sptree.nnode) {
      printf("%7d%7d  %-14s %9.5f ", 
         sptree.nodes[i].father, i, sptree.nodes[i].name, sptree.nodes[i].age);
      printf("  %2s %6.2f %6.2f", 
         fossils[sptree.nodes[i].fossil], sptree.nodes[i].tfossil[0], sptree.nodes[i].tfossil[1]);
      if(sptree.nodes[i].nson)
         printf("  (%2d %2d)", sptree.nodes[i].sons[0], sptree.nodes[i].sons[1]);
      FPN(F0);
   }
   copySptree();
   FPN(F0); OutaTreeN(F0,0,0); FPN(F0); OutaTreeN(F0,1,0);  FPN(F0); 
   OutaTreeN(F0,1,1); FPN(F0);
}



int GetMem (void)
{
/* This allocates memory for conditional probabilities (conP).  gnodes[locus] is 
   not allocated here but in GetGtree().

   conP0 is the space for the tips.  Only one copy is needed as the conditional 
   probabilities for tips are fixed according to the data and never change.  
   Space for conP0 is allocated for a locus only if the locus has unclean 
   (ambiguity) data, that is, if(data.cleandata[locus]==0).
   Conditional probabilities for internal nodes are com.conPin[2], allocated 
   according to data.ns[locus] and data.npatt[locus] at all loci.  Two copies 
   of the space are allocated, hence the [2].  The copy used for the current 
   gene trees is com.conPin[com.curconP] while the copy for proposed gene trees 
   is com.conPin[!com.curconP].  data.conP_offset[locus] marks the starting 
   position in conPin[] for each locus.
*/
   /* sconP0: tips; sconP: internal nodes */
   int locus,j, sconP0, s1, g=data.ngene;
   double *conP0, *conP, *rates;

   /* get mem for conP (internal nodes), and conP0 (tips) */
   if(mcmc.usedata) {
      data.conP_offset[0]=0;
      for(locus=0,com.sconP=sconP0=0; locus<g; locus++) {
         s1= com.ncode * data.npatt[locus];
         com.sconP += s1*(data.ns[locus]-1)*sizeof(double);
         if(!data.cleandata[locus])  
            sconP0 += s1*data.ns[locus]*sizeof(double);
         if(locus<g-1)
            data.conP_offset[locus+1]=data.conP_offset[locus]+(data.ns[locus]-1)*s1;
      }
      if((com.conPin[0]=(double*)malloc(2*com.sconP+sconP0))==NULL) 
         error2("oom conP");
      com.conPin[1]      =com.conPin[0]+com.sconP/sizeof(double);
      if(sconP0) com.conP0=com.conPin[1]+com.sconP/sizeof(double);
      printf("\n%d bytes for conP, %d bytes for conP0\n", 2*com.sconP,sconP0);

      /* set gnodes[locus][].conP for tips and internal nodes */
      com.curconP=0; conP=com.conPin[0]; conP0=com.conP0;
      for(locus=0; locus<g; locus++) {
         for(j=data.ns[locus]; j<data.ns[locus]*2-1; j++,conP+=com.ncode*data.npatt[locus])
            gnodes[locus][j].conP = conP;
         if(!data.cleandata[locus]) {
            for(j=0; j<data.ns[locus]; j++,conP0+=com.ncode*data.npatt[locus])
               gnodes[locus][j].conP = conP0;
            UseLocus(locus, 0, 0);
            InitConditionalPNode ();
         }
      }
   }
   if(com.clock==0) {
      s1=(sptree.nspecies*2-1)*g*sizeof(double);
      if(noisy) printf("%d bytes for rates.\n", s1);
      if((rates=(double*)malloc(s1))==NULL) error2("oom for rates");
      for(j=0; j<sptree.nspecies*2-1; j++) 
         sptree.nodes[j].lnrates = rates+g*j;
      if((data.correl=(double*)malloc(2*g*g*sizeof(double)))==NULL)
         error2("oom when getting variance-covariance matrix.");
      data.Sinv = data.correl+g*g;
   }
   return(0);
}


void FreeMem (void)
{
   int locus, j;

   FOR(locus,data.ngene) free(gnodes[locus]);
   if(mcmc.usedata) free(com.conPin[0]);
   for(locus=0; locus<data.ngene; locus++) {
      free(data.fpatt[locus]);
      for(j=0;j<data.ns[locus]; j++)
         free(data.z[locus][j]);
   }
   if(!com.clock) {
      free(sptree.nodes[0].lnrates);
      free(data.correl);
   }
}



int UseLocus (int locus, int copyconP, int setSeqName)
{
/* This point nodes to the gene tree at locus gnodes[locus] and set com.z[] 
   etc. for likelihood calculation for the locus.  Note that the gene tree 
   topology (gnodes[]) is never copied, but nodes[].conP are repositioned in the 
   algorithm.  If (copyconP && mcmc.useData), the conP for internal nodes point 
   to a fixed place (indicated by data.conP_offset[locus]) in the alternative 
   space com.conPin[!com.curconP].  Note that the conP for each locus uses the 
   correct space so that this routine can be used by all the proposal steps, 
   some of which operates on one locus and some change all loci.
   The conP for tips always point to com.conP0 (when com.cleandata=0).
*/
   int i, s1=com.ncode*data.npatt[locus];
   double *conPt=com.conPin[!com.curconP]+data.conP_offset[locus];

   com.ns=data.ns[locus]; com.ls=data.ls[locus];
   tree.root=data.root[locus]; tree.nnode=2*com.ns-1;
   nodes=gnodes[locus];
   if(copyconP && mcmc.usedata) { /* this preserves the old conP. */
      memcpy(conPt, gnodes[locus][com.ns].conP, s1*(com.ns-1)*sizeof(double));
      for(i=com.ns; i<com.ns*2-1; i++)
         nodes[i].conP = conPt+(i-com.ns)*s1;
   }

   if(mcmc.usedata) {
      com.cleandata=data.cleandata[locus];
      com.npatt=com.posG[1]=data.npatt[locus];  com.posG[0]=0;
      com.fpatt=data.fpatt[locus];
      for(i=0; i<com.ns; i++) com.z[i]=data.z[locus][i];

      /* The following is model-dependent */
      if(!com.fix_kappa) com.kappa=data.kappa[locus];
      xtoy(data.pi[locus], com.pi, com.ncode);
      if(com.model==REV)
         getpi_sqrt (com.pi, com.pi_sqrt, com.ncode, &com.npi0);
      getQfactor (locus, com.kappa);
      com.Qfactor=data.Qfactor[locus];
      /* Should call EigenQ () here, if necessary. */

   }
   if(setSeqName)
      FOR(i,com.ns) com.spname[i]=sptree.nodes[nodes[i].ipop].name;
   return(0);
}


int AcceptRejectLocus (int locus, int accept)
{
/* This accepts or rejects the proposal at one locus.  
   This works for proposals that change one locus only.  
   After UseLocus(), gnodes[locus][ns].conP points to the alternative 
   conP space.  If the proposal is accepted, this copies the newly calculated 
   conP into gnodes[locus][ns].conP.  In either case, gnodes[].conP is 
   repositioned.
   Proposals that change all loci use switchconP() to accept the proposal.
*/
   int i, ns=data.ns[locus], s1=com.ncode*data.npatt[locus];
   double *conP=com.conPin[com.curconP]+data.conP_offset[locus];

   if(mcmc.usedata) {
      if(accept)
         memcpy(conP, gnodes[locus][ns].conP, s1*(ns-1)*sizeof(double));
      for(i=ns; i<ns*2-1; i++)
         gnodes[locus][i].conP = conP+(i-ns)*s1;
   }
   return(0);
}

void switchconPin(void)
{
/* This resets gnodes[locus].conP to the alternative com.conPin, to avoid 
   recalculation of conP, when a proposal is accepted in steps that change all 
   loci in one go, such as UpdateTimes() and UpdateParameters().
*/
   int i,locus;
   double *conP=com.conPin[com.curconP=!com.curconP];
   
   for(locus=0; locus<data.ngene; locus++)
      for(i=data.ns[locus]; i<data.ns[locus]*2-1; i++,conP+=com.ncode*data.npatt[locus])
         gnodes[locus][i].conP = conP;
}



int SetPGene(int igene, int _pi, int _UVRoot, int _alpha, double xcom[])
{
/* This is not used. */
   if(com.ngene!=1) error2("com.ngene==1?");
   return (0);
}


int SetParameters (double x[])
{
/* This is not used. */
   return(0);
}

int GetPMatBranch2(double PMat[], double t)
{
/* This gets the transition probability matrix.
   Should remove duplicated calculation of Qfactor, for all branches for the 
   same locus.
*/
   double Qrates[2], T,C,A,G,Y,R;

   NPMat++;
   Qrates[0]=Qrates[1]=com.kappa;
   if(com.seqtype==0) {
      if (com.model<=K80)
         PMatK80(PMat, t, com.kappa);
      else if(com.model<REV) {
         T=com.pi[0]; C=com.pi[1]; A=com.pi[2]; G=com.pi[3]; Y=T+C; R=A+G;
         if (com.model==F84) { 
            Qrates[0]=1+com.kappa/Y;   /* kappa1 */
            Qrates[1]=1+com.kappa/R;   /* kappa2 */
         }
         else if (com.model<=HKY85) Qrates[1]=Qrates[0];
         com.Qfactor=1/(2*T*C*Qrates[0] + 2*A*G*Qrates[1] + 2*Y*R);
         PMatTN93(PMat, t*com.Qfactor*Qrates[0], t*com.Qfactor*Qrates[1], t*com.Qfactor, com.pi);
      }
   }
   return(0);
}


int ConditionalPNode (int inode, int igene, double x[])
{
   int n=com.ncode, i,j,k,h, ison;
   double t;

   FOR(i,nodes[inode].nson) {
      ison=nodes[inode].sons[i];
      if (nodes[ison].nson>0  && !com.oldconP[ison])
         ConditionalPNode (nodes[inode].sons[i], igene, x);
   }
   FOR(i,com.npatt*n) nodes[inode].conP[i]=1;

   FOR(i, nodes[inode].nson) {
      ison=nodes[inode].sons[i];

      t=nodes[ison].branch*_rateSite;
      if(t<0) error2("blength<0");

      GetPMatBranch2(PMat, t);
      if(com.cleandata && nodes[ison].nson<1) /* tip */
         for(h=0; h<com.npatt; h++) 
            FOR(j,n) 
               nodes[inode].conP[h*n+j]*=PMat[j*n+com.z[ison][h]];
      else {
         for(h=0; h<com.npatt; h++) 
            FOR(j,n) {
               for(k=0,t=0; k<n; k++)    /* t is used as temp */
                  t+=PMat[j*n+k]*nodes[ison].conP[h*n+k];  /* expensive line! */
               nodes[inode].conP[h*n+j]*=t;
            }
      }
   }  /*  for (ison)  */

   /* node scaling, not done yet */
   if(com.nnodeScale && com.nodeScale[inode])  NodeScale(inode, 0, com.npatt);
   return (0);
}


double lnpData (double lnpDi[])
{
/* This calculates the log likelihood, the log of the probability of the data 
   given gtree[] for each locus.
   This updates gnodes[locus][].conP.
*/
   int j,locus;
   double lnL=0, y;

   if(mcmc.saveconP) FOR(j,sptree.nspecies*2-1) com.oldconP[j]=0;
   for(locus=0; locus<data.ngene; locus++) {
      UseLocus(locus,0, 1);
      /* UseLocus(locus,0, 0); */
      y=lnpD_locus(locus);

      if(testlnL && fabs(lnpDi[locus]-y)>1e-5)
         printf("\tlnLi %.6f != %.6f at locus %d\n", lnpDi[locus],y,locus+1);
      lnpDi[locus]=y;
      lnL+=y;
   }

   return(lnL);
}

double lnpD_locus (int locus)
{
/* This calculates ln{Di|Gi, Bi} using times in sptree.nodes[].age and rates.
   Rates are in data.rgene[] if (clock==1) and in sptree.nodes[].lnrates if 
   (clock==0).  Branch lengths in the gene tree, nodes[].branch, are calculated 
   in this routine but they may not be up-to-date before calling this routine.
   UseLocus() is called before this routine.
*/
   int  i,j, dad;
   double lnL,b,t;

   if(!mcmc.usedata) return(0);
   if(com.clock) {
      FOR(i,tree.nnode)  /* age in gene tree */
         nodes[i].age=sptree.nodes[nodes[i].ipop].age;
      FOR(i,tree.nnode) {
         if(i==tree.root) continue;
         dad = nodes[i].father;
         nodes[i].branch = (nodes[dad].age-nodes[i].age) * data.rgene[locus];
      }
   }
   else {
      FOR(i,tree.nnode) {
         if(i==tree.root) continue;
         for(j=nodes[i].ipop,b=0; j!=nodes[nodes[i].father].ipop; j=dad) {
            dad = sptree.nodes[j].father;
            t = sptree.nodes[dad].age-sptree.nodes[j].age;

            b += BranchLength (t, sptree.nodes[dad].lnrates[locus], 
                                  sptree.nodes[j].lnrates[locus], 
                                  data.correl[locus*data.ngene+locus]);
         }
         nodes[i].branch=b;
      }
   }
   lnL = -com.plfun(NULL, -1);
   return (lnL);
}



int GetOptions (char *ctlf)
{
   int i, nopt=20, lline=255, seed;
   char line[255],*pline, opt[20], *comment="*#";
   char *optstr[] = {"seed", "seqfile","outfile","treefile", "model", "clock",
        "alpha", "ncatG", "usedata", "cleandata", "BlengthMethod", "BDparas", 
        "kappa_gamma", "rgene_gamma", "sigma2_gamma", 
        "print", "burnin", "sampfreq", "nsample", "finetune"};
   double t=1, *eps=mcmc.finetune;
   FILE  *fctl=gfopen (ctlf, "r");

   if (fctl) {
      if (noisy) printf ("\nReading options from %s..\n", ctlf);
      for (;;) {
         if (fgets (line, lline, fctl) == NULL) break;
         for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
            if (isalnum(line[i]))  { t=1; break; }
            else if (strchr(comment,line[i])) break;
         if (t==0) continue;
         sscanf (line, "%s%*s%lf", opt, &t);
         if ((pline=strstr(line, "="))==NULL) error2 ("option file.");

         for (i=0; i<nopt; i++) {
            if (strncmp(opt, optstr[i], 8)==0)  {
               if (noisy>=9)
                  printf ("\n%3d %15s | %-20s %6.2f", i+1,optstr[i],opt,t);
               switch (i) {
                  case ( 0): 
                     seed=(int)t;
                     if(seed<=0) seed=(unsigned int)(2*(int)time(NULL)+1);
                     SetSeed(seed);
                     break;
                  case ( 1): sscanf(pline+1, "%s", com.seqf);    break;
                  case ( 2): sscanf(pline+1, "%s", com.outf);    break;
                  case ( 3): sscanf(pline+1, "%s", com.treef);   break;
                  case ( 4): com.model=(int)t;      break;
                  case ( 5): com.clock=(int)t;      break;
                  case ( 6): com.alpha=t;           break;
                  case ( 7): com.ncatG=(int)t;      break;
                  case ( 8): mcmc.usedata=(int)t;   break;
                  case ( 9): com.cleandata=(int)t;  break;
                  case (10): data.BlengthMethod=(int)t;  break;
                  case (11): 
                     sscanf(pline+1,"%lf%lf%lf", &data.birth,&data.death,&data.sampling);
                     break;
                  case (12): 
                     sscanf(pline+1,"%lf%lf", data.kappagamma, data.kappagamma+1); break;
                  case (13): 
                     sscanf(pline+1,"%lf%lf", data.rgenegamma, data.rgenegamma+1); break;
                  case (14): 
                     sscanf(pline+1,"%lf%lf", data.sigma2gamma,data.sigma2gamma+1); break;
                  case (15): mcmc.print=(int)t;     break;
                  case (16): mcmc.burnin=(int)t;    break;
                  case (17): mcmc.sampfreq=(int)t;  break;
                  case (18): mcmc.nsample=(int)t;   break;
                  case (19): 
                     sscanf(pline+1,"%lf%lf%lf%lf%lf", eps,eps+1,eps+2,eps+3,eps+4);
                     break;
               }
               break;
            }
         }
         if (i==nopt)
            { printf ("\noption %s in %s\n", opt, ctlf);  exit (-1); }
      }
      fclose (fctl);
   }
   else
      if (noisy) error2("\nno ctl file..");


   if (com.model>HKY85)  error2("model.");
   if (com.alpha==0)  {
      com.nalpha=0;    com.plfun=lfun;
   }
   else {
      if (com.ncatG<2 || com.ncatG>NCATG) error2 ("ncatG");
      com.plfun=lfundG;
   }
   if (com.model==JC69 || com.model==F81) { com.fix_kappa=1; com.kappa=1; }

   return (0);
}



int DownSptreeSetTime (int inode)
{
/* This goes down the species tree (to the tips) to specify the node ages, 
   starting from the root age.
*/
   int j,ison;

   for (j=0; j<sptree.nodes[inode].nson; j++) {
      ison=sptree.nodes[inode].sons[j];
      if(sptree.nodes[ison].nson) {
         sptree.nodes[ison].age = sptree.nodes[inode].age*(.5+.5*rndu());
         DownSptreeSetTime(ison);
      }
   }
   return(0);
}

int GetRgenes (void)
{
/* This routine is for ad hoc determinatino of the mean lograte at each locus, 
   \mu, which is the mean for the normal lograte at the root.  The routine uses 
   the calibration nodes with bounds B or gamma G to calulate rough rates.  
   The log rate at the root is drawn from N(\mu, t_root*sigma^2).

   PROBLEM: this uses com.alpha to calculate pairwise distances, while the 
   parameter is being estimated.
*/
   int inode, i,j, ipop, locus, nleft,nright,marks[NS], sons[2];
   int nf, nf_low, nf_up;  /* fossil bounds */
   double d,md, kappa, bigD=5, mu,mu_low,mu_up, *tf, rootage;

   if(mcmc.usedata==0) {
      for(i=0; i<data.ngene; i++) data.rgene[i]=0;
      return(0);
   }

   puts("Estimating rough mean rates at each locus...");

   for(locus=0; locus<data.ngene; locus++) {
      UseLocus(locus, 0, 0);
      nf=nf_low=nf_up=0; mu=mu_low=mu_up=0;
      for(inode=com.ns; inode<tree.nnode; inode++) {
         ipop=nodes[inode].ipop;  tf=sptree.nodes[ipop].tfossil;
         if(sptree.nodes[ipop].fossil==0) continue;
         sons[0]=nodes[inode].sons[0];
         sons[1]=nodes[inode].sons[1];
         for(i=0,nleft=nright=0; i<com.ns; i++) {
            for(j=i,marks[i]=0; j!=tree.root; j=nodes[j].father) {
               if(j==sons[0])       { marks[i]=1; nleft++;  break; }
               else if (j==sons[1]) { marks[i]=2; nright++; break; }
            }
         }
         if(nleft==0 || nright==0) error2("we should not be here.");

printf("\nLocus %2d node %2d (ipop %2d) left %2d right %2d",
       locus+1,inode, ipop, nleft,nright);

         for(i=0,md=0; i<com.ns; i++) {
            for(j=0; j<com.ns; j++) {
               if(marks[i]==1 && marks[j]==2) {
                  d = DistanceIJ(i, j, com.model, com.alpha, &kappa);
                  if(d<0) { puts("big distance"); d=bigD; }
                  md+=d;
               }
            }
         }
         md/=(nleft*nright);  /* average distance for the node */
         switch(sptree.nodes[ipop].fossil) {
         case (LOWER_F): { mu_low+=md/tf[0];         nf_low++; break; }
         case (UPPER_F): { mu_up+=md/tf[1];          nf_up++;  break; }
         case (BOUND_F): { mu+=md/((tf[0]+tf[1])/2);  nf++;    break; }
         case (GAMMA_F): { mu+=md/((tf[0]-1)/tf[1]);  nf++;    break; }
         }

printf(" %s d = %9.5f u=%9.5f", fossils[sptree.nodes[ipop].fossil], md, mu);
matIout(F0, marks, 1, com.ns);

      }
      if(nf==0) { printf("need fossils for locus %d\n",locus+1); exit(-1); }
      mu/=nf; 
      if(nf_low) {
         mu_low/=nf_low; 
         /* if(mu>mu_low) mu=mu_low;  */
      }
      if(nf_up) {
         mu_up/=nf_up;
         /* if(mu<mu_up) mu=mu_up; */
      }
      data.rgene[locus]=log(mu);
      /* initial log rates at the root of the species tree at each locus */
      rootage=sptree.nodes[sptree.root].age;
      sptree.nodes[sptree.root].lnrates[locus] 
         = data.rgene[locus] 
         + sqrt(data.correl[locus*data.ngene+locus]*rootage)*rndnormal();
   }

   for(i=0;i<data.ngene;i++) 
      printf("mu for locus %d = %9.5f\n", i+1, data.rgene[i]);

printf("Enter rates (not logrates) for root.. ");
for(i=0;i<data.ngene;i++) 
scanf("%lf", &data.rgene[i]);
for(i=0;i<data.ngene;i++) data.rgene[i]=log(data.rgene[i]);
for(i=0;i<data.ngene;i++) 
printf("%12.6f", data.rgene[i]); FPN(F0);
getchar();

   return(0);
}



int GetInitials(double space[])
{
/* This sets the initial values for starting the MCMC, and returns np, the 
   number of parameters in the MCMC, to be collected in collectx().
   The routine guarantees that each node age is younger than its ancestor's age.
   It does not check the consistency of divergence times against the fossil 
   constraints.  As the model assumes soft bounds, any divergence times are 
   possible, even though this means that the chain might start from a poor 
   place.
   space[data.ngene]
*/
   int np=0, locus,i,j, g=data.ngene, nroughtime=0;
   double roughtime=0; /* rought age for root */


   com.rgene[0]=-1;  /* com.rgene[] is not used. */

   puts("\ngetting initial values to start MCMC.");

   /* set up rough time unit by looking at the fossil info */
   for(j=sptree.nspecies; j<sptree.nnode; j++) {
      if(!sptree.nodes[j].fossil) continue;
      nroughtime++;
      if(sptree.nodes[j].fossil==LOWER_F)
         roughtime+=sptree.nodes[j].tfossil[0]*2;
      else if(sptree.nodes[j].fossil==UPPER_F)
         roughtime+=sptree.nodes[j].tfossil[1]/2;
      else if(sptree.nodes[j].fossil==BOUND_F)
         roughtime+=(sptree.nodes[j].tfossil[0]+sptree.nodes[j].tfossil[1])/2;
      else if(sptree.nodes[j].fossil==GAMMA_F)
         roughtime+=(sptree.nodes[j].tfossil[0]-1)/sptree.nodes[j].tfossil[1];
   }
   roughtime/=nroughtime;
   sptree.nodes[sptree.root].age = roughtime*(0.5+2*rndu());
   DownSptreeSetTime(sptree.root);

   /* set up initial rates */
   if(com.clock) {
      np = sptree.nspecies-1 + g;
      for(i=0; i<g; i++) data.rgene[i]=rndu();
   }
   else {
      np = sptree.nspecies-1 + g*(sptree.nspecies*2-1) + g*(g+1)/2;

      for(i=0; i<g; i++)
         data.correl[i*g+i] = 0.1+rndu();
      for(i=0; i<g; i++)
         for(j=0; j<i; j++)
            data.correl[i*g+j]=data.correl[j*g+i] = .8-rndu();
      getSinvDetS(space);

      /* initializes \mu from data, and resets lnrates at root */
      GetRgenes();
      for(j=0; j<sptree.nnode; j++)
         for(locus=0; locus<data.ngene; locus++) 
            sptree.nodes[j].lnrates[locus] = rndu()-.5;
   }

   /* set up substitution parameters */
   for(i=0; i<g; i++) {
      if(com.model>=K80 && !com.fix_kappa) {
         data.kappa[i] = 1+2*rndu();
         np++;
      }
   }
   return(np);
}

int collectx (double x[])
{
/* this collects parameters into x[] for printing and summarizing.
   It returns the number of parameters.
*/
   int uselograte=0, i,j, np=0, g=data.ngene;

   for(i=sptree.nspecies; i<sptree.nspecies*2-1; i++) 
      x[np++] = sptree.nodes[i].age;
   if(com.clock)
      for(i=0; i<g; i++) x[np++] = data.rgene[i];
   else {
      if(uselograte) 
         for(i=0; i<g; i++)  for(j=0; j<sptree.nnode; j++)
            x[np++] = sptree.nodes[j].lnrates[i];
      else
         for(i=0; i<g; i++)  for(j=0; j<sptree.nnode; j++)
            x[np++] = exp(sptree.nodes[j].lnrates[i]);

      for(i=0; i<g; i++)  /* variance sigma^2 */
         x[np++] = data.correl[i*g+i];

      for(i=0; i<g; i++)  /* correlation rho_ij */
         for(j=0; j<i; j++)
            x[np++] = data.correl[i*g+j];
   }
   if(!com.fix_kappa)
      for(i=0; i<g; i++)
         x[np++]=data.kappa[i];

   if(np!=com.np) error2("np error in collectx.");
   return(0);
}


#define P0t(expmlt) (rho*(lambda-mu)/(rho*lambda +(lambda*(1-rho)-mu)*expmlt))

double PDFkernelBD(double t, double t1, double vt1, double lambda, double mu, double rho)
{
/* this calculates the kernel density from the birth-death process with 
   species sampling.
*/
   double pdf, expmlt, small=1e-20;

   if(fabs(lambda-mu)<small)
      pdf = (1+rho*lambda*t1)/(t1*square(1+rho*lambda*t));
   else {
      expmlt=exp((mu-lambda)*t);
      pdf = P0t(expmlt);
      pdf = pdf*pdf * lambda/(vt1*rho) * expmlt;
   }
   return(pdf);
}

double CDFkernelBD(double t, double t1, double vt1, double lambda, double mu, double rho)
{
/* this calculates the CDF for the kernel density from the birth-death process with 
   species sampling.
*/
   double cdf, expmlt, small=1e-20;

   if(fabs(lambda-mu)<small)
      cdf = (1+rho*lambda*t1)*t/(t1*(1+rho*lambda*t));
   else {
      expmlt=exp((mu-lambda)*t);
      cdf = rho*lambda/vt1 * (1-expmlt)/(rho*lambda +(lambda*(1-rho)-mu)*expmlt);
   }
   return(cdf);
}

double lnptC_Fossil(void)
{
/* This calculates the prior density of times at calibration nodes as specified 
   by the fossil calibration information.  a=tL, b=tU.
*/
   int j;
   double thetaL=log(0.1)/log(0.9), tail=0.025, t, lnpt=0, a,b,theta;

   if(sptree.nfossil==0) return(0);
   for(j=sptree.nspecies; j<sptree.nnode; j++) {
      if(!sptree.nodes[j].fossil) continue;
      t=sptree.nodes[j].age;
      a=sptree.nodes[j].tfossil[0]; b=sptree.nodes[j].tfossil[1]; 
      if(sptree.nodes[j].fossil==LOWER_F) {
         if(t<a) lnpt += log(tail*thetaL/a) + (thetaL-1)*log(t/a);
         else    lnpt += log(tail*thetaL/a);
      }
      else if(sptree.nodes[j].fossil==UPPER_F) {
         theta = (1-tail)/(tail*b);
         lnpt += (t<b ? log((1-tail)/b) : log(tail*theta)-theta*(t-b));
      }
      else if(sptree.nodes[j].fossil==BOUND_F) {
         if(t<a) {
            theta = (1-tail*2)*a/(tail*(b-a));
            lnpt += log(tail*theta/a) + (theta-1)*log(t/a);
         }
         else if(t<b)
            lnpt += log((1-tail*2)/(b-a));
         else {
            theta = (1-tail*2)/(tail*(b-a));
            lnpt += log(tail*theta) - theta*(t-b);
         }
      }
      else if(sptree.nodes[j].fossil==GAMMA_F) {
         lnpt += a*log(b)-b*t+(a-1)*log(t)-LnGamma(a);
      }
   }
   return(lnpt);
}

double lnpriorTimes (void)
{
/* This calculates the prior density of node times in the master species tree:
   sptree.nodes[].age.  It uses sptree.nodes[].tfossil, sptree.nodes[].fossil[],
   and lamdba, mu, rho from the birth-death process with species sampling 
   (data.birth, data.death, data.sampling).     
   The routine sorts the node ages in the species tree and then uses the 
   birth-death prior conditional on the calibration points.  t[0] is t1 in the 
   paper.

   rank[1]=3: means that age of node [ns+1] is the 3rd youngest.
   nodesort[3]=1 means that the 3rd yougest age is node number ns+1.
   Root (node ns) is excluded in the ranking.
*/
   int  i,j,k, n1=sptree.nspecies-1, rank[NS-1], rankprev, nfossil;
   int  nodesort[NS-1]; /* nodes with times sorted */
   double t[NS-1], lnpt, Scale, vt1=0, expmlt, cdf, cdfprev, small=1e-20;

   double lambda=data.birth, mu=data.death, rho=data.sampling;

   if(sptree.root!=sptree.nspecies) error2("node number for root fixed.");
   for(j=sptree.nspecies; j<sptree.nnode; j++) 
      { t[j-sptree.nspecies]=sptree.nodes[j].age; }

   /* ranking the (n-2) node ages */
   for(i=1,rank[0]=-1,nodesort[0]=-1;i<n1;i++) {
      for(j=1,k=1;j<n1;j++)  if(j!=i && t[j]<=t[i]) k++;
      rank[i]=k;
      nodesort[k]=i;
   }

   if(debug==1) {
      matout2(F0, t, 1, n1, 9, 5);
      FOR(j,n1) printf("%9d", rank[j]);  FPN(F0);
      FOR(j,n1) printf("%9d", nodesort[j]);  FPN(F0);
   }

   /* calculate vt1, needed only if (lambda!=mu) */
   if(fabs(lambda-mu)>small) {
      expmlt=exp((mu-lambda)*t[0]);
      vt1 = 1-P0t(expmlt)/rho*expmlt;
   }
   /* calculate f_BD(t_{_C)}, joint of the remaining times */
   for(j=1,lnpt=1,Scale=0; j<n1; j++) {
      if(!sptree.nodes[sptree.nspecies+j].fossil)
         lnpt *= PDFkernelBD(t[j], t[0], vt1, lambda, mu, rho);
      if(j%50==0) { Scale+=log(lnpt); lnpt=1; }
   }
   lnpt=Scale+log(lnpt);

   /* Now calculate f_BD(t_C), marginal for calibration nodes.
      This goes through the nodes in the order of their ages, so that node j 
      is the k-th youngest.
   */
   for(k=1,nfossil=0,rankprev=0,cdfprev=0; k<n1; k++) { 
      if(!sptree.nodes[j=sptree.nspecies+nodesort[k]].fossil)
         continue;
      cdf = CDFkernelBD(t[nodesort[k]], t[0], vt1, lambda, mu, rho);
      if(k-rankprev-1>0) {
         lnpt -= (k-rankprev-1)*log(cdf-cdfprev);
         lnpt += LnGamma((double)k-rankprev-1+1);
      }
      if(debug==1)
         printf("Fossil at node %d age %9.4f  rank diff %d - %d cdf %9.5f\n", j,t[nodesort[k]], k, rankprev, cdf);
      rankprev=k;  cdfprev=cdf;
      nfossil++;
   }
   if(nfossil && n1-1-rankprev>0) {
      lnpt -= (n1-1.-rankprev)*log(1-cdfprev);
      lnpt += LnGamma(n1-1.-rankprev+1);
   }

   lnpt += lnptC_Fossil();
   if(debug==1) printf("\npdf = %.12f\n", exp(lnpt));

   return (lnpt);
}


int SetOldCondP (int spnode)
{
/* This sets com.oldconP[j]=0 if node j in the gene tree needs updating, after 
   either rate or time has changed for spnode in the species tree.  This is to 
   avoid duplicated computation of conditional probabilities.
   The routine first finds spnode or its first ancestor that is present in the 
   gene tree, then identifies that node in the genetree.  This reveals the 
   oldest node j in the gene tree that is descendent of spnode.  Node j and all 
   its ancestors in the gene tree need updating.

   All nodes were set to have com.oldconP[i]=1, before calling this routine.
   This routine changes some com.oldconP[] into 0 but do not change any 0 to 1.

   The gene tree is in nodes[], as UseLocus has been called prior to this.
   This is called by UpdateTimes and UpdateRates.
*/
   int i, j=spnode, nblength=0;

   if(j>=tree.nnode || j!=nodes[j].ipop) {

      /* From among spnode and its ancestors in the species tree, find the 
         first node, i, that is in genetree.
      */
      for(i=spnode; i!=-1; i=sptree.nodes[i].father) {

         /* Find that node in genetree that is node i in species tree.  
            Its descendent, node j, is the oldest node in gene tree that is 
            descendent of spnode.
         */
         for(j=0; j<tree.nnode && nodes[j].ipop!=i; j++) ;

         if(j<tree.nnode) break;
      }
   }

   if(j<tree.nnode)
      for( ; com.oldconP[j]=0,j!=tree.root; j=nodes[j].father)
         nblength++;
   return(nblength);
}

double UpdateTimes (double *lnL, double finetune)
{
/* This updates the node ages in the master species tree: sptree.nodes[].age.
   It changes one node age at a time, within the bounds of mother and daughter 
   node ages.
*/
   int  locus, is, i;
   double accepted=0, t, tnew, tb[2];
   double lnacceptance=0, lnLd=0, lnpDinew[NGENE], lnpTnew;

   if(debug==2) puts("\nUpdateTimes ");
   for(is=sptree.nspecies; is<sptree.nnode; is++) {
      t=sptree.nodes[is].age;
      lnLd=0;   tb[0]=0;  tb[1]=OLDAGE;
      tb[0]=max2(sptree.nodes[sptree.nodes[is].sons[0]].age,
                 sptree.nodes[sptree.nodes[is].sons[1]].age);
      if(is!=sptree.root) 
         tb[1]=sptree.nodes[sptree.nodes[is].father].age;
      tnew=rnduab(t-finetune/2,t+finetune/2);
      tnew=reflect(tnew,tb[0],tb[1]);

      sptree.nodes[is].age=tnew;  /* the only thing to change in proposal */

      lnpTnew=lnpriorTimes();
      lnacceptance = lnpTnew-data.lnpT;

      for(locus=0,lnLd=0; locus<data.ngene; locus++) {
         UseLocus(locus, 1, 0);

         if(mcmc.saveconP) {
            FOR(i,sptree.nspecies*2-1) com.oldconP[i]=1;
            SetOldCondP(is);
         }
         lnpDinew[locus] = data.lnpDi[locus];
         if(!com.oldconP[tree.root]) {
            lnpDinew[locus] = lnpD_locus(locus);
            lnLd += lnpDinew[locus]-data.lnpDi[locus];
         }
      }
      lnacceptance += lnLd;

      if(debug==2) printf("species %2d tb: %8.5f %8.5f t: %8.5f%8.5f %9.2f", is,tb[0],tb[1],t,tnew,lnLd);

      if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
         accepted++;
         data.lnpT=lnpTnew;
         FOR(locus,data.ngene) data.lnpDi[locus]=lnpDinew[locus];

         *lnL+=lnLd;
         if(mcmc.usedata) switchconPin();
         if(debug==2) printf(" Y (%4d)\n", NPMat);
      }
      else {
         sptree.nodes[is].age = t;
         if(debug==2) printf(" N (%4d)\n", NPMat);
         for(locus=0; locus<data.ngene; locus++)
            AcceptRejectLocus(locus, 0);  /* reposition conP */
      }
   }
   return(accepted/(sptree.nspecies-1.));
}


void getSinvDetS (double space[])
{
/* This constructs Sinv (inverse of S) and detS (|S|) using data.correl[], 
   which has rho_ij off-diagonal and sigma^2 on diagonal.
   This is called every time \sigma^2, \rho_ij are updated.
   space[data.ngene]
*/
   int i,j, g=data.ngene;
   int debug=0;

   for(i=0; i<g; i++) {
      data.Sinv[i*g+i] = data.correl[i*g+i];
      for(j=0; j<i; j++)
         data.Sinv[i*g+j] = data.Sinv[j*g+i]
            = data.correl[i*g+j]*sqrt(data.correl[i*g+i]*data.correl[j*g+j]);
   }

   if(debug) {
      printf("\nCorrel & S & Sinv ");
      matout(F0, data.correl, g, g);
      matout(F0, data.Sinv, g, g);
   }

   matinv (data.Sinv, g, g, space);
   /* check this? */
   printf("\a");
   data.detS = space[0];

   if(debug) {
      matout(F0, data.Sinv, g, g);
      printf("|S| = %.6f\n", data.detS);
      getchar();
   }

}

double lnpriorRates (void)
{
/* This calculates the log of the prior of logrates under the geometric 
   Brownian motion model of rate evolution.  The algorithm cycles through all 
   nodes in the species tree and multiply the conditional g-variate normal 
   densities across the tree, accounting for correlations in lograte between 
   genes.  There is no need for recursion as the order at which sptree.nodes 
   are visited is unimportant.
   The logrates are stored in sptree.nodes[].lnrates[].
   The variance-covariance matrix, S, is implicitly used through its inverse 
   data.Sinv[] and determinant data.detS.
   The root lograte has mean data.rgene[locus] and variance sigma^2*root_age.
   There is no bias correction for the root lograte.

   Space requirement: yd[ngene] uses com.space, for y - \mu.
*/
   int i,j, inode, g=data.ngene, dad;
   double lnpR=0, t, zz, *y, *yd=com.space;

   for(inode=0; inode<sptree.nnode; inode++) {
      y=sptree.nodes[inode].lnrates;
      t = sptree.nodes[inode].age;
      if(inode==sptree.root) { /* prior mean of lnrate at root is data.rgene[] */
         for(i=0; i<g; i++)
            /* yd[i] = y[i] - data.rgene[i]; */
            yd[i] = y[i] - data.rgene[i] + t*data.correl[i*g+i]/2;
      }
      else {
         dad=sptree.nodes[inode].father;
         t = sptree.nodes[dad].age - t;
         for(i=0; i<g; i++)
            yd[i] = y[i] - sptree.nodes[dad].lnrates[i] + t*data.correl[i*g+i]/2;
      }
      /* zz is the quadratic form in exponent. */
      for(i=0,zz=0; i<g; i++) {
         zz    += data.Sinv[i*g+i]*yd[i]*yd[i];
         for(j=0; j<i; j++)
            zz += data.Sinv[i*g+j]*yd[i]*yd[j]*2;
      }
      lnpR += -zz/(2*t) - g/2*log(PI*2) - 1/2.*g*log(t) - log(data.detS)/2;
   }
   return lnpR;
}

double UpdateRates (double *lnL, double finetune)
{
/* This updates logrates under the Brownian motion rate-evolution model.
   The algorithm cycles through the loci.  For each locus, it cycles through 
   the interior nodes in the species tree, and updates the three logrates 
   at that node and at its two descendent nodes.

   Waste of computation: For every proposal to change logrates, lnpriorRates() 
   is called to calculate the prior for all logrates at all loci on the whole 
   tree, thus wasting computation.
*/
   int locus, inode, j, g=data.ngene, sons[2];
   double accepted=0, lnpRnew, lnpDinew, lnacceptance, lnLd, lnRold[3];
   double e=finetune/2;

   if(com.clock)
      return UpdateRatesClock(lnL, finetune);

   if(debug==3) puts("\nUpdateRates");
   for(locus=0; locus<g; locus++) {
      for(inode=sptree.nspecies; inode<sptree.nnode; inode++) {
         for(j=0; j<2; j++) sons[j]=sptree.nodes[inode].sons[j];

         lnRold[0]=sptree.nodes[inode].lnrates[locus];
         lnRold[1]=sptree.nodes[sons[0]].lnrates[locus];
         lnRold[2]=sptree.nodes[sons[1]].lnrates[locus];

         sptree.nodes[inode].lnrates[locus]   = rnduab(lnRold[0]-e, lnRold[0]+e);
         sptree.nodes[sons[0]].lnrates[locus] = rnduab(lnRold[1]-e, lnRold[1]+e);
         sptree.nodes[sons[1]].lnrates[locus] = rnduab(lnRold[2]-e, lnRold[2]+e);

         UseLocus(locus, 1, 0);  /* copyconP=1 */

         if(mcmc.saveconP) {
            FOR(j,sptree.nspecies*2-1) com.oldconP[j]=1;
            for(j=0; j<2; j++) {
               if(sptree.nodes[sons[j]].nson) {
                  SetOldCondP(sptree.nodes[sons[j]].sons[0]);
                  SetOldCondP(sptree.nodes[sons[j]].sons[1]);
               }
               else
                  SetOldCondP(sons[j]);
            }
         }

         /* This priorRates calculation is wasteful but not expensive. */
         lnpRnew=lnpriorRates();

         lnacceptance = lnpRnew-data.lnpR;
         lnpDinew=lnpD_locus(locus);
         lnLd = lnpDinew - data.lnpDi[locus];
         lnacceptance +=  lnLd;

         if(lnacceptance>0 || rndu()<exp(lnacceptance)) {
            accepted++;
            if(mcmc.usedata) AcceptRejectLocus(locus,1);
   
            data.lnpR=lnpRnew;
            data.lnpDi[locus]=lnpDinew;
            *lnL += lnLd;
         }
         else {
            if(mcmc.usedata) AcceptRejectLocus(locus,0);
            sptree.nodes[inode].lnrates[locus]=lnRold[0];
            sptree.nodes[sons[0]].lnrates[locus]=lnRold[1];
            sptree.nodes[sons[1]].lnrates[locus]=lnRold[2];
         }
      }
   }
   return(accepted/(g*(sptree.nspecies-1.)));
}


double logPriorRatio(double xnew, double xold, double a, double b, int distribution)
{
/* This calculates the log of prior ratio when x is updated from xold to xnew.
   x has distribution with parameters a and b.
*/
   double lnr=-99999;

   if(distribution==GAMMA)
      lnr = (a-1)*log(xnew/xold) - b*(xnew-xold);
   else if(distribution==IGAMMA)
      lnr = -9999;
   return(lnr);
}

double UpdateRatesClock (double *lnL, double finetune)
{
/* This updates rates data.rgene[] under the clock by cycling through the loci.
   The proposal affects all branch lengths, so com.oldconP[]=0.
*/
   int locus, j, g=data.ngene;
   double accepted=0, rgene0, lnLd, lnpDinew, lnacceptance, c=1;

   if(debug==3) puts("\nUpdateRatesClock");
   if(mcmc.saveconP) FOR(j,sptree.nspecies*2-1) com.oldconP[j]=0;
   for(locus=0; locus<g; locus++) {
      rgene0=data.rgene[locus];
      data.rgene[locus] = rgene0 * (c=exp(finetune*(.5-rndu())));
  
      UseLocus(locus, 1, 0);
      lnpDinew=lnpD_locus(locus);
      lnLd = lnpDinew - data.lnpDi[locus];   /* likelihood ratio */
      lnacceptance =  lnLd + log(c);         /* proposal ratio */
      /* prior ratio */
      lnacceptance += logPriorRatio(data.rgene[locus],rgene0,data.rgenegamma[0],data.rgenegamma[1],GAMMA);

      if(debug==3)
         printf("\nLocus %2d rgene %9.4f%9.4f %10.5f", locus+1, rgene0, data.rgene[locus], lnLd);

      if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
         accepted++;
         *lnL += lnLd;
         data.lnpDi[locus]=lnpDinew;
         if(mcmc.usedata) AcceptRejectLocus(locus,1);
         if(debug==3) printf(" Y\n");
      }
      else {
         data.rgene[locus]=rgene0;
         if(mcmc.usedata) AcceptRejectLocus(locus,0); /* reposition conP */

         if(debug==3) printf(" N\n");
      }
   }

   return(accepted/g);
}

double UpdateParameters (double *lnL, double finetune)
{
/* This updates parameters in the substitution model such as the ts/tv 
   rate ratio for each locus.

   Should we update the birth-death process parameters here as well?

*/
   int locus, j;
   double accepted=0, kappa0=-1, lnLd,lnpDinew, lnacceptance, c=1;

   if(debug==4) puts("\nUpdateParameters");
   if(com.fix_kappa) return(0);

   if(mcmc.saveconP) FOR(j,sptree.nspecies*2-1) com.oldconP[j]=0;
   for(locus=0; locus<data.ngene; locus++) {
      if(!com.fix_kappa) {
         kappa0=data.kappa[locus];
         data.kappa[locus] = kappa0 * (c=exp(finetune*(.5-rndu())));
      }

      UseLocus(locus, 1, 0); /* this sets com.kappa too. */

      lnpDinew=lnpD_locus(locus);
      lnLd = lnpDinew - data.lnpDi[locus];
      lnacceptance =  lnLd + log(c);
      lnacceptance += logPriorRatio(data.kappa[locus],kappa0,data.kappagamma[0],data.kappagamma[1],GAMMA);

      if(debug==4)
         printf("\nLocus %2d kappa %9.4f%9.4f %10.5f",
            locus+1, kappa0, data.kappa[locus], lnLd);

      if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
         accepted++;
         *lnL += lnLd;
         data.lnpDi[locus]=lnpDinew;
         if(mcmc.usedata) AcceptRejectLocus(locus,1);
         if(debug==4) printf(" Y\n");
      }
      else {
         data.kappa[locus]=kappa0;
         if(mcmc.usedata) AcceptRejectLocus(locus,0);

         if(debug==4) printf(" N\n");
      }
   }
   return(accepted/data.ngene);
}


double UpdateParaRates (double finetune, double space[])
{
/* This updates data.correl[], which has variance \sigma^2 on diagonal and 
   correlation coeffecient \rho_ij off diagonal.  The routine uses the mean 
   lnrates in data.rgene[], while data.rgene[] are fixed in the program.  
   The proposals here do not change the likelihood.
   A large step length is used for changing the variance than the correlation.
   The hyperpriors are the same across loci.
*/
   int i,j, g=data.ngene, nproposal=0;
   double old, r, lnacceptance, lnpRnew, e[2], naccept[2]={0}; /* for V and R */
   double rb[2]={-1, 1};

   if(com.clock) return(0);
   if(debug==5) puts("\nUpdateParaRates");

   e[0]=finetune*4; e[1]=finetune;  /* for V and R */
   for(i=0; i<g; i++) {
      for(j=0; j<=i; j++) {
         old=data.correl[i*g+j];
         r = rnduab(old-e[i!=j], old+e[i!=j]);
         if(i==j) { /* variance sigma^2 */
            data.correl[i*g+i] = fabs(r);
            lnacceptance = logPriorRatio(data.correl[i*g+i],old,data.sigma2gamma[0],data.sigma2gamma[1],GAMMA);
         }
         else { /* correlation coefficent, in (0,1) */
            r = reflect(r,rb[0],rb[1]);
            if(r<rb[0]||r>rb[1]) error2("r?");

            data.correl[i*g+j]=data.correl[j*g+i] = r;
            lnacceptance = 0; /* sliding window, uniform prior */
         }

         if(debug==5) printf("%s %d%d %9.5f -> %9.5f ", (i==j?"V":"R"), i,j, old,r);

         getSinvDetS(space);
         /* prior ratio */
         lnpRnew = lnpriorRates();
         lnacceptance += lnpRnew-data.lnpR;
         if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
            naccept[i!=j]++;
            data.lnpR = lnpRnew;
            if(debug==5) printf(" Y |S| = %.5f\n", data.detS);
         }
         else {
            data.correl[i*g+j]=data.correl[j*g+i] = old;
            getSinvDetS(space);  /* have to recalculate this */
            if(debug==5) printf(" N |S| = %.5f\n", data.detS);
         }
      }
   }
   if(debug==5) {
      printf(" acceptance %.3f (%d) %.3f (%d)\n", 
         naccept[0]/g,g,naccept[1]/(g*(g-1)/2),g*(g-1)/2);
      getchar();
   }
   return((naccept[0]+naccept[1])/(g*(g+1)/2.));
}



double mixing (double finetune)
{
/* This multiplies all times by a constant c and divides all rates by the same 
   constant c, so that the likelihood does not change.  This is true only if 
   the branch length is calculated as the arithmetic or geometric mean, and 
   not if it is calculated using the Bronian motion integral.
   Rates are in data.rgene[] if clock==1 or sptree.nodes[].lnrates[] 
   if clock==0.
*/
   int nmultiply=sptree.nspecies-1, ndivide=(com.clock?data.ngene:0);
   int i,j;
   double accepted=0, a=data.rgenegamma[0], b=data.rgenegamma[1], c, lnc;
   double lnacceptance=0, lnpTnew, lnpRnew=-1;

   lnc=finetune*(rndu()-.5); c=exp(lnc);

   for(j=sptree.nspecies; j<sptree.nnode; j++) 
      sptree.nodes[j].age*=c;

   if(com.clock) {
      for(j=0; j<data.ngene; j++) {
         lnacceptance = -(a-1)*lnc-b*(data.rgene[j]/c-data.rgene[j]);
         data.rgene[j]/=c;
      }
   }
   else {
      for(j=0; j<data.ngene; j++) {
         for(i=0; i<sptree.nnode; i++) 
            sptree.nodes[i].lnrates[j] -= lnc;
      }
      lnpRnew=lnpriorRates();
      lnacceptance = lnpRnew-data.lnpR;
      /* Now we may need to update the likelihood 
         if(data.blengthtype==BROWNIAN)
      */
   }
   lnpTnew=lnpriorTimes();
   lnacceptance += lnpTnew-data.lnpT;
   lnacceptance += (nmultiply-ndivide)*lnc;

   if(lnacceptance>0 || rndu()<exp(lnacceptance)) { /* accept */
      accepted=1;
      data.lnpT=lnpTnew;
      data.lnpR=lnpRnew;
   }
   else {   /* reject */

      for(j=sptree.nspecies; j<sptree.nnode; j++) sptree.nodes[j].age/=c;

      if(com.clock) 
         for(j=0; j<data.ngene; j++) data.rgene[j]*=c;
      else {
         for(j=0; j<data.ngene; j++) {
            for(i=0; i<sptree.nnode; i++) 
               sptree.nodes[i].lnrates[j] += lnc;
         }
      }
   }
   return(accepted);
}


int MCMC (FILE* fout)
{
   char mcmcout[32]="mcmc.out";
   FILE *fmcmcout=gfopen(mcmcout,"w");
   int nsteps=4+!com.clock, nxpr[2]={10,2};
   int i,j,k, ir, g=data.ngene;
   double Paccept[5]={0}, lnL=0, nround=0, *x, *mx, *vx;
   char timestr[32];

   starttime();
   noisy=3;
   mcmc.saveconP=1;
   if(!mcmc.usedata) mcmc.saveconP=0;
   FOR(j,sptree.nspecies*2-1) com.oldconP[j]=0;

   printf("\n%d burnin, sampled every %d, %d samples\n", 
           mcmc.burnin, mcmc.sampfreq, mcmc.nsample);
   if(mcmc.usedata) puts("Approximating posterior, using sequence data");
   else             puts("Approximating prior, not using sequence data");

   printf("(Settings: cleandata=%d.  print=%d  saveconP=%d.)\n", 
          com.cleandata, mcmc.print, mcmc.saveconP);

   com.np=GetInitials(com.space);
   x=(double*)malloc(com.np*(com.np+2)*sizeof(double));
   if(x==NULL) error2("oom in mcmc");
   mx=x+com.np;  vx=mx+com.np;

/*
if(com.clock==0) {
FOR(i,g) { data.rgene[i]=0; data.kappa[i]=2; }
for(i=sptree.nspecies; i<sptree.nnode; i++) 
   sptree.nodes[i].age = (1-(i-sptree.nspecies)*.1);
for(i=0;i<g;i++) {
   data.correl[i*g+i]=0.4;
   for(j=0;j<i;j++) 
      data.correl[i*g+j]=data.correl[j*g+i]=.5;
}
getSinvDetS(com.space);
FOR(i,g) {
   for(j=0; j<sptree.nnode; j++)
      sptree.nodes[j].lnrates[i] = 2*rndu()-1;
   sptree.nodes[sptree.root].lnrates[i] = rndu();
}
}
*/

   printSptree();
   collectx(x);

   puts("\npriors: ");
   printf("G(%.4f, %.4f) for kappa\n", data.kappagamma[0], data.kappagamma[1]);
   if(com.clock)
      printf("G(%.4f, %.4f) for rgene\n", data.rgenegamma[0], data.rgenegamma[1]);
   else
      printf("G(%.4f, %.4f) for sigma^2\n", data.sigma2gamma[0], data.sigma2gamma[1]);

   printf("finetune parameters: ");
   for(j=0;j<nsteps;j++) printf(" %7.4f",mcmc.finetune[j]); FPN(F0);
   printf("\nStarting MCMC... ");
   printf("Initial parameters (np = %d):\n", com.np);
   for(j=0;j<com.np;j++) printf(" %5.3f",x[j]); FPN(F0);

/*
printf("branch length method (0=arithmetic; 1=geometric; 2=Brownian)? ");
scanf("%d", &data.BlengthMethod);
*/


   /* calculates prior for times and likelihood for each locus */
   data.lnpT=lnpriorTimes();
   for(j=0; j<data.ngene; j++) com.rgene[j]=-1;
   if(!com.clock) {
      data.lnpR=lnpriorRates();
   }
   lnL=lnpData(data.lnpDi);
   printf("\nlnL0 = %9.2f\n", lnL);

   zero(mx,com.np);  zero(vx,com.np*com.np);
   if(com.np<nxpr[0]+nxpr[1]) { nxpr[0]=com.np; nxpr[1]=0; }
   for(ir=-mcmc.burnin; ir<mcmc.sampfreq*mcmc.nsample; ir++) {
      if(ir==0) { /* reset after burnin */
#if 0
         if(mcmc.burnin>10) {  /* reset finetune parameters */
            for(j=0;j<nsteps;j++) 
               if(Paccept[j]<.1 || Paccept[j]>.8) {
                  mcmc.finetune[j] *= Paccept[j]/0.3;
                  printf("finetune #%d reset to %.5f\n", j+1,mcmc.finetune[j]);
               }
         }
#endif
         nround=0; zero(Paccept,nsteps);
         zero(mx,com.np); zero(vx,com.np*com.np); 
      }

      Paccept[0] = (Paccept[0]*nround + UpdateTimes(&lnL, mcmc.finetune[0]))/(nround+1);
      Paccept[1] = (Paccept[1]*nround + UpdateRates(&lnL, mcmc.finetune[1]))/(nround+1);
      Paccept[2] = (Paccept[2]*nround + UpdateParameters(&lnL, mcmc.finetune[2]))/(nround+1);
      Paccept[3] = (Paccept[3]*nround + mixing(mcmc.finetune[3]))/(nround+1);
      Paccept[4] = (Paccept[4]*nround + UpdateParaRates(mcmc.finetune[4],com.space))/(nround+1);

      nround++;
      collectx(x);
      FOR(j,com.np) mx[j]=(mx[j]*(nround-1)+x[j])/nround;
      FOR(i,com.np) FOR(j,com.np) 
         vx[i*com.np+j] += (x[i]-mx[i])*(x[j]-mx[j]) * (nround-1.)/nround;
      if(mcmc.print && ir>=0 && (ir+1)%mcmc.sampfreq==0) {
         for(j=0;j<com.np; j++) fprintf(fmcmcout," %7.5f",x[j]);
         fprintf(fmcmcout," %10.3f\n",lnL);
      }
      if((ir+1)%max2(mcmc.sampfreq, mcmc.sampfreq*mcmc.nsample/10000)==0) {
         printf("\r%3.0f%%", (ir+1.)/(mcmc.nsample*mcmc.sampfreq)*100.);
         FOR(j,nsteps) printf(" %4.2f", Paccept[j]);  printf(" ");


         FOR(j,nxpr[0]) printf(" %5.3f", mx[j]);
         if(com.np>nxpr[0]+nxpr[1]) printf(" -");
         FOR(j,nxpr[1]) printf(" %5.3f",mx[com.np-nxpr[1]+j]);
/*
FOR(j,12) printf(" %5.3f", mx[sptree.nspecies-1+j]);
*/
         if(mcmc.usedata) printf(" %4.1f ",lnL);
      }
      if(mcmc.sampfreq*mcmc.nsample>20 && (ir+1)%(mcmc.sampfreq*mcmc.nsample/20)==0) {
         printf(" %s\n", printtime(timestr));
         testlnL=1;
         if(fabs(lnL-lnpData(data.lnpDi))>1e-5) {
            printf("\n%12.6f = %12.6f?\n", lnL, lnpData(data.lnpDi));
            puts("lnL not right?");
         }
         testlnL=0;
      }
   }
   fclose(fmcmcout);
   puts("\nSpecies tree with branch lengths in time units, for TreeView");
   FOR(j,sptree.nspecies-1) sptree.nodes[sptree.nspecies+j].age=mx[j];
   copySptree();
   FPN(F0); OutaTreeN(F0,1,1); FPN(F0);
   if(com.clock==0) {
      printf("Posterior means of S = {\\sigma_ij}:");
      for(i=0,k=sptree.nspecies-1+sptree.nnode*data.ngene; i<data.ngene; i++)
         data.correl[i*data.ngene+i]=mx[k++];
      for(i=0;i<data.ngene;i++) for(j=0;j<i;j++)
         data.correl[i*data.ngene+j]=data.correl[j*data.ngene+i]=mx[k++];
      matout(F0, data.correl, data.ngene, data.ngene);
   }
   printf("\nTime used: %s", printtime(timestr));

   FOR(i,com.np*com.np)  vx[i] /= nround;
   matout(fout, mx, 1, com.np);
   FOR(i,com.np) fprintf(fout, "%12.6f", (vx[i*com.np+i])); FPN(fout);
   matout(fout, vx, com.np, com.np);

   fprintf(fout,"\nVar-cov for rates\n");
   for(i=sptree.nspecies-1; i<sptree.nspecies*5-3; i++,FPN(fout))
      for(j=sptree.nspecies-1; j<=i; j++)
         fprintf(fout, " %9.3f", vx[i*com.np+j]);
   fprintf(fout, "\ncorrel for rates\n");
   for(i=sptree.nspecies-1; i<sptree.nspecies*5-3; i++,FPN(fout))
      for(j=sptree.nspecies-1; j<=i; j++)
         fprintf(fout, " %9.3f", vx[i*com.np+j]/sqrt(vx[i*com.np+i]*vx[j*com.np+j]));

   free(x); 
   if(mcmc.print) {
      printf("\nSummarizing, time reset.");
      fprintf(fout,"\n\nSummary of MCMC results\n");
      DescriptiveStatistics(fout, mcmcout, 20, 20);
   }
   return(0);
}


double BranchLength (double t, double y0, double yt, double sigma2)
{
   int method=data.BlengthMethod;  /* 0=arithmetic; 1=geometric; 2=Brownian */
   int i,n=100;
   double a=1./2*t*sigma2, b=yt-y0+a, rate=0, x;

   if(method==0)        /* arithmetic mean */
      return( (exp(y0)+exp(yt))/2*t );
   else if (method==1)  /* geometric mean */
      return(exp((y0+yt)/2)*t);
   else {               /* brownian motion integral */
      for(i=0; i<n; i++) {
         x=i/(double)n;
         rate += exp(-a*x*x + b*x + y0);
      }
      return (rate*t/n);
   }
}
