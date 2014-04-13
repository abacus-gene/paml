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
#define NGENE         3          /* used for gnodes[NGENE] */
#define LSPNAME       30
#define NCODE         4
#define NCATG         8

extern int noisy, NFunCall;
extern char BASEs[];

int GetOptions(char *ctlf);
int ReadTreeSeqs(FILE*fout);
int GetMem(void);
void FreeMem(void);
int UseLocus(int locus, int copycondP, int setmodel, int setSeqName);
int AcceptRejectLocus(int locus, int accept);
void switchconPin(void);
int SetPGene(int igene, int _pi, int _UVRoot, int _alpha, double xcom[]);
int DownSptreeSetTime(int inode);
int GetRateGenes(FILE* fout);
void getSinvDetS(double space[]);
int GetInitials(FILE*fout, double space[]);
int GetGtree(int locus);
void printGtree(int printBlength);
int SetParameters(double x[]);
int ConditionalPNode(int inode, int igene, double x[]);
double BranchLength(double t, double y0, double yt, int method, double sigma2);
double lnpData(double lnpDi[]);
double lnpD_locus(int locus);
double lnpriorTimes(void);
double lnpriorRates(void);
void copySptree(void);
void printSptree(void);
int MCMC(FILE* fout);
int LabelOldCondP(int spnode);
double UpdateTimes(double *lnL, double finetune);
double UpdateRates(double *lnL, double finetune);
double UpdateRatesClock(double *lnL, double finetune);
double UpdateParameters(double *lnL, double finetune);
double UpdateParaRates(double finetune, double space[]);
double mixingClock (double finetune);
double mixing(double *lnL, double finetune);

struct CommonInfo {
   char *z[NS], *spname[NS], seqf[96],outf[96],treef[96];
   char oldconP[NNODE];       /* update conP for node? (0 yes; 1 no) */
   int seqtype, ns, ls, ngene, posG[2],lgene[1], *pose, npatt;
   int np, ncode, ntime, nrate, nrgene, nalpha, npi, ncatG, print, seed;
   int cleandata, ndata;
   int model, clock, fix_kappa, fix_omega, fix_alpha, fix_rgene, Mgene;
   double *fpatt, kappa, omega, alpha;
   double rgene[NGENE],piG[NGENE][NCODE];  /* not used */
   double (*plfun)(double x[],int np), freqK[NCATG], rK[NCATG], *conP, *fhK;
   double pi[NCODE];
   int    sconP, curconP;                    /* curconP = 0 or 1 */
   double *conPin[2], *conP0, space[10000];  /* space used for S^-1 and |S| */
   int conPSiteClass, readpattf;         /* not used */

   int NnodeScale;
   char *nodeScale;    /* nScale[ns-1] for interior nodes */
   double *nodeScaleF;       /* nScaleF[npatt] for scale factors */
}  com;

struct TREEB {
   int  nbranch, nnode, root, branches[NBRANCH][2];
}  tree;
struct TREEN { /* ipop is node number in species tree */
   int father, nson, sons[2], ibranch, ipop;
   double branch, age, label, *conP;
   char fossil;
}  *nodes, **gnodes, nodes_t[2*NS-1];

/* nodes_t[] is working space.  nodes is a pointer and moves around.  
   gnodes[] holds the gene trees, subtrees constructed from the master species 
   tree.  Each locus has a full set of rates (lnrates) for all branches on the 
   master tree, held in sptree.nodes[].lnrates.  Branch lengths in the gene 
   tree are calculated by using those rates and the divergence times.

   gnodes[][].label in the gene tree is used to store branch lengths estimated 
   under no clock when mcmc.usedata=2 (normal approximation to likelihood).
*/


struct SPECIESTREE {
   int nbranch, nnode, root, nspecies, nfossil;
   struct TREESPN {
      char name[LSPNAME+1], fossil;  /* fossil: 0, 1, 2, 3 */
      int father, nson, sons[2];
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
   double pi[NGENE][NCODE];
   double rgene[NGENE], kappa[NGENE], omega[NGENE], alpha[NGENE];
   double birth, death, sampling;
   double rgenegamma[2], kappagamma[2], omegagamma[2], alphagamma[2], sigma2gamma[2];
   /* correl[g*g] is the variance-correlation matrix in lnrates among loci.  
   correl is updated in the MCMC.  sigma2[g] are the variances, Sinv the inverse 
   of S, and detS is det|S| */
   double *correl, *Sinv, *sigma2, detS;
   double *varb[NGENE];
}  data;

struct MCMCPARAMETERS {
   int burnin, nsample, sampfreq, usedata, saveconP, print;
   double finetune[5];
}  mcmc; /* control parameters */


char *models[]={"JC69","K80","F81","F84","HKY85","T92","TN93","REV"};
enum {JC69, K80, F81, F84, HKY85, T92, TN93, REV} MODELS;

int nR=4;
double PMat[16], Cijk[64], Root[4];
double _rateSite=1, OldAge=99;
int LASTROUND=0, BayesEB, debug=0, testlnL=0, NPMat=0; /* no use for this */

/* for sptree.nodes[].fossil: lower, upper, bounds, gamma, inverse-gamma */
enum {FREE_F=0, LOWER_F, UPPER_F, BOUND_F, GAMMA_F, IGAMMA_F} FOSSIL_FLAGS;
char *fossils[]={" ", "L", "U", "B", "G", "IG"}; 
enum {GAMMA, IGAMMA} DISTRIBUTIONS;
enum {ARITHMETIC, GEOMETRIC, BROWNIAN} BLengthMethods;


#define MCMCTREE  1
#include "treesub.c"

int main (int argc, char *argv[])
{
   char ctlf[]="mcmctree.ctl";
   FILE  *fout;

   noisy=3;
   com.alpha=0.;     com.ncatG=1;
   com.ncode=4;      com.clock=1;

   printf("MCMCTREE in %s\n", VerStr);
   if (argc>1) 
      { strcpy(ctlf, argv[1]); printf ("\nctlfile reset to %s.\n", ctlf); }

   data.birth=2;    data.death=1; data.sampling=0.05; 
   com.cleandata=0; mcmc.usedata=2;

   GetOptions (ctlf);
   fout=gfopen(com.outf,"w");

   fprintf(fout, "MCMCTREE (%s) %s\n", VerStr, com.seqf);
   fprintf(fout, "\nseed = %d\n", com.seed);

   ReadTreeSeqs(fout);

   MCMC(fout);
   return (0);
}

int GetMem (void)
{
/* This allocates memory for conditional probabilities (conP).  
   gnodes[locus] is not allocated here but in GetGtree().

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

   Memory arrangement if(com.conPSiteClass=1):
   ncode*npatt for each node, by node, by iclass, by locus
*/
   /* sconP0: tips; sconP: internal nodes */
   int locus,j,k, sconP0, s1,sG=1, sfhK=0, g=data.ngene;
   double *conP0, *conP, *lnrates, *varb;

   /* get mem for conP (internal nodes), and conP0 (tips) */
   if(mcmc.usedata==1) {
      if(!com.fix_alpha && mcmc.saveconP) {
         com.conPSiteClass=1;  sG=com.ncatG;
      }
      data.conP_offset[0]=0;
      for(locus=0,com.sconP=sconP0=0; locus<g; locus++) {
         s1= com.ncode * data.npatt[locus];
         com.sconP += sG*s1*(data.ns[locus]-1)*sizeof(double);
         if(!data.cleandata[locus])
            sconP0 += s1*data.ns[locus]*sizeof(double);
         if(locus<g-1)
            data.conP_offset[locus+1]=data.conP_offset[locus]+sG*s1*(data.ns[locus]-1);
      }

      if((com.conPin[0]=(double*)malloc(2*com.sconP+sconP0))==NULL) 
         error2("oom conP");

      com.conPin[1]       =com.conPin[0]+com.sconP/sizeof(double);
      if(sconP0) com.conP0=com.conPin[1]+com.sconP/sizeof(double);
      printf("\n%d bytes for conP, %d bytes for conP0\n", 2*com.sconP,sconP0);

      /* set gnodes[locus][].conP for tips and internal nodes */
      com.curconP=0; conP0=com.conP0;
      for(locus=0; locus<g; locus++) {
         conP=com.conPin[0]+data.conP_offset[locus];
         for(j=data.ns[locus]; j<data.ns[locus]*2-1; j++,conP+=com.ncode*data.npatt[locus])
            gnodes[locus][j].conP = conP;
         if(!data.cleandata[locus]) {
            for(j=0; j<data.ns[locus]; j++,conP0+=com.ncode*data.npatt[locus])
               gnodes[locus][j].conP = conP0;
            UseLocus(locus, 0, mcmc.usedata, 0);
            InitConditionalPNode ();
         }
      }

      if(!com.fix_alpha) {
         for(locus=0; locus<g; locus++)
            sfhK = max2(sfhK, data.npatt[locus]);
         sfhK *= com.ncatG*sizeof(double);
         if((com.fhK=(double*)realloc(com.fhK,sfhK))==NULL) error2("oom");
      }

   }
   else if(mcmc.usedata==2) { /* using com.condP for com.varb */
      for(locus=0,k=0; locus<data.ngene; locus++)  k+=2*data.ns[locus]-1;
      if((varb=com.conPin[0]=(double*)malloc(k*sizeof(double)))==NULL)
         error2("oom varb");
      for(j=0; j<k; j++)  varb[j]=-1;
      for(locus=0,j=0; locus<data.ngene; locus++) {
         data.varb[locus]=varb+j;
         j+=2*data.ns[locus]-1;
      }
   }

   if(com.clock==0) {  /* space for lnrates */
      s1=(sptree.nspecies*2-1)*g*sizeof(double);
      if(noisy) printf("%d bytes for lnrates.\n", s1);
      if((lnrates=(double*)malloc(s1))==NULL) error2("oom for lnrates");
      for(j=0; j<sptree.nspecies*2-1; j++) 
         sptree.nodes[j].lnrates = lnrates+g*j;

      /* space for variance-covariance matrix */
      if((data.correl=(double*)malloc((2*g*g+g)*sizeof(double)))==NULL)
         error2("oom when getting sqrt variance-covariance matrix.");
      data.Sinv = data.correl+g*g;
      data.sigma2 = data.Sinv+g*g;
   }
   return(0);
}

void FreeMem (void)
{
   int locus, j;

   FOR(locus,data.ngene) free(gnodes[locus]);
   free(gnodes);
   if(mcmc.usedata) free(com.conPin[0]);
   if(mcmc.usedata==1) {
      for(locus=0; locus<data.ngene; locus++) {
         free(data.fpatt[locus]);
         for(j=0;j<data.ns[locus]; j++)
            free(data.z[locus][j]);
      }
   }
   if(!com.clock) {
      free(sptree.nodes[0].lnrates);
      free(data.correl);
   }

   if(mcmc.usedata==1 && com.alpha) free(com.fhK);
}


int UseLocus (int locus, int copyconP, int setModel, int setSeqName)
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

   The conP for tips always point to com.conP0 (when com.cleandata=0), and are 
   not changed here.

   Try to replace this with UseLocus() for B&C.
*/
   int i, s1=com.ncode*data.npatt[locus], sG=(com.conPSiteClass?com.ncatG:1);
   double *conPt=com.conPin[!com.curconP]+data.conP_offset[locus];

   com.ns=data.ns[locus]; com.ls=data.ls[locus];
   tree.root=data.root[locus]; tree.nnode=2*com.ns-1;
   nodes=gnodes[locus];
   if(copyconP && mcmc.usedata==1) { /* this preserves the old conP. */
      memcpy(conPt, gnodes[locus][com.ns].conP, sG*s1*(com.ns-1)*sizeof(double));
      for(i=com.ns; i<tree.nnode; i++)
         nodes[i].conP = conPt+(i-com.ns)*s1;
   }

   if(setModel && mcmc.usedata==1) {
      com.cleandata=data.cleandata[locus];
      com.npatt=com.posG[1]=data.npatt[locus];  com.posG[0]=0;
      com.fpatt=data.fpatt[locus];
      for(i=0; i<com.ns; i++) com.z[i]=data.z[locus][i];

      /* The following is model-dependent */
      com.kappa=data.kappa[locus];
      com.omega=data.omega[locus];
      com.alpha=data.alpha[locus];

      xtoy(data.pi[locus], com.pi, com.ncode);
      if(com.model<=TN93)
         EigenTN93(com.model, com.kappa, com.kappa, com.pi, &nR, Root, Cijk);

      if(com.alpha)
         DiscreteGamma (com.freqK,com.rK,com.alpha,com.alpha,com.ncatG,0);
/*
      com.NnodeScale=data.NnodeScale[locus];
      com.nodeScale=data.nodeScale[locus];
      nS = com.NnodeScale*com.npatt * (com.conPSiteClass?com.ncatG:1);
      FOR(i,nS) com.nodeScaleF[i]=0;
*/
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
   int i, ns=data.ns[locus], s1=com.ncode*data.npatt[locus], sG=1;
   double *conP=com.conPin[com.curconP]+data.conP_offset[locus];

   if(mcmc.usedata==1) {
      if(com.conPSiteClass) sG=com.ncatG;
      if(accept)
         memcpy(conP, gnodes[locus][ns].conP, sG*s1*(ns-1)*sizeof(double));
      for(i=ns; i<ns*2-1; i++)
         gnodes[locus][i].conP = conP+(i-ns)*s1;
   }
   return(0);
}

void switchconPin(void)
{
/* This reposition pointers gnodes[locus].conP to the alternative com.conPin, 
   to avoid recalculation of conditional probabilities, when a proposal is 
   accepted in steps that change all loci in one go, such as UpdateTimes() 
   and UpdateParameters().
   Note that for site-class models (com.conPSiteClass), gnodes[].conP points 
   to the space for class 0, and the space for class 1 starts (ns-1)*ncode*npatt
   later.  Such repositioning for site classes is achieved in fx_r().
*/
   int i,locus;
   double *conP;

   com.curconP=!com.curconP;
   
   for(locus=0; locus<data.ngene; locus++) {
      conP=com.conPin[com.curconP]+data.conP_offset[locus];
      for(i=data.ns[locus]; i<data.ns[locus]*2-1; i++,conP+=com.ncode*data.npatt[locus])
         gnodes[locus][i].conP = conP;
   }
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

int GetPMatBranch2 (double PMat[], double t)
{
/* This calculates the transition probability matrix.
*/
   double Qrates[2], T,C,A,G,Y,R, mr;

   NPMat++;
   Qrates[0]=Qrates[1]=com.kappa;
   if(com.seqtype==0) {
      if (com.model<=K80)
         PMatK80(PMat, t, com.kappa);
      else if(com.model<=TN93) {
         T=com.pi[0]; C=com.pi[1]; A=com.pi[2]; G=com.pi[3]; Y=T+C; R=A+G;
         if (com.model==F84) { 
            Qrates[0]=1+com.kappa/Y;   /* kappa1 */
            Qrates[1]=1+com.kappa/R;   /* kappa2 */
         }
         else if (com.model<=HKY85) Qrates[1]=Qrates[0];
         mr=1/(2*T*C*Qrates[0] + 2*A*G*Qrates[1] + 2*Y*R);
         PMatTN93(PMat, t*mr*Qrates[0], t*mr*Qrates[1], t*mr, com.pi);
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
         ConditionalPNode(ison, igene, x);
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
   if(com.NnodeScale && com.nodeScale[inode])  NodeScale(inode, 0, com.npatt);
   return (0);
}


double lnpD_locus_Normal (double varb[])
{
/* This calculates the likelihood using the normal approxiamtion (Thorne et al. 
   1998).  The branch lengths estimated under no clock have been read into 
   nodes[][].label and the variances are in data.varb[].  The branch lengths 
   predicted from the rate-evolution model (that is, products of rates and 
   times) are in nodes[].branch.  
   The tree is rooted, and the two branch lengths around the root are summed up
   and treated as one branch length.  This is different from funSS_AHRS(), which 
   has the data branch lengths in nodes[].branch and calculate the predicted 
   likelihood by multiplying nodes[].age and nodes[].label (which holds the 
   rates).  Think about a redesign to avoid confusion.
*/
   int j, son0, son1;
   double lnLb, b,be, v;

   for(j=0,lnLb=0; j<tree.nnode; j++) {
      if(j==son0 || j==son1) continue;
      if(j==tree.root) {
         son0=nodes[j].sons[0];
         son1=nodes[j].sons[1];
         b  = nodes[son0].label +nodes[son1].label;
         be = nodes[son0].branch+nodes[son1].branch;
         v=varb[nodes[tree.root].age-nodes[son0].age>1e-6?son0:son1];
      }
      else {
         b  = nodes[j].label;
         be = nodes[j].branch;
         v = varb[j];
      }
      lnLb += -square(be-b)/(2*v) -log(v)/2;

/* printf("n%-2d lnL=%12.6f\n", j+1,lnLb); */

   }
   return (lnLb);
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

   if (mcmc.usedata==0)  return(0);
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

            b += BranchLength(t, sptree.nodes[dad].lnrates[locus], 
                                 sptree.nodes[j].lnrates[locus], 
                                 data.BlengthMethod, data.sigma2[locus]);
         }
         nodes[i].branch=b;
      }
   }
   if(mcmc.usedata==1)
      lnL = -com.plfun(NULL, -1);
   else if(mcmc.usedata==2)
      lnL = lnpD_locus_Normal(data.varb[locus]);

   return (lnL);
}

double lnpData (double lnpDi[])
{
/* This calculates the log likelihood, the log of the probability of the data 
   given gtree[] for each locus.
   This updates gnodes[locus][].conP for every node.
*/
   int j,locus;
   double lnL=0, y;

   if(mcmc.saveconP) 
      FOR(j,sptree.nspecies*2-1) com.oldconP[j]=0;
   for(locus=0; locus<data.ngene; locus++) {
      UseLocus(locus,0, mcmc.usedata, 1);
      y=lnpD_locus(locus);

      if(testlnL && fabs(lnpDi[locus]-y)>1e-5)
         printf("\tlnLi %.6f != %.6f at locus %d\n", lnpDi[locus],y,locus+1);
      lnpDi[locus]=y;
      lnL+=y;
   }

   return(lnL);
}


int GetOptions (char *ctlf)
{
   int i, nopt=22, lline=255;
   char line[255],*pline, opt[20], *comment="*#";
   char *optstr[] = {"seed", "seqfile","outfile","treefile", "ndata", 
        "model", "clock", 
        "alpha", "ncatG", "usedata", "cleandata", "BlengthMethod", "BDparas", 
        "kappa_gamma", "alpha_gamma", "rgene_gamma", "sigma2_gamma", 
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
                     com.seed=(int)t;
                     if(com.seed<=0) com.seed=abs(2*(int)time(NULL)+1);
                     SetSeed(com.seed);
                     break;
                  case ( 1): sscanf(pline+1, "%s", com.seqf);    break;
                  case ( 2): sscanf(pline+1, "%s", com.outf);    break;
                  case ( 3): sscanf(pline+1, "%s", com.treef);   break;
                  case ( 4): com.ndata=(int)t;      break;
                  case ( 5): com.model=(int)t;      break;
                  case ( 6): com.clock=(int)t;      break;
                  case ( 7): com.alpha=t;           break;
                  case ( 8): com.ncatG=(int)t;      break;
                  case ( 9): mcmc.usedata=(int)t;   break;
                  case (10): com.cleandata=(int)t;  break;
                  case (11): data.BlengthMethod=(int)t;  break;
                  case (12): 
                     sscanf(pline+1,"%lf%lf%lf", &data.birth,&data.death,&data.sampling);
                     break;
                  case (13): 
                     sscanf(pline+1,"%lf%lf", data.kappagamma, data.kappagamma+1); break;
                  case (14): 
                     sscanf(pline+1,"%lf%lf", data.alphagamma, data.alphagamma+1); break;
                  case (15): 
                     sscanf(pline+1,"%lf%lf", data.rgenegamma, data.rgenegamma+1); break;
                  case (16): 
                     sscanf(pline+1,"%lf%lf", data.sigma2gamma, data.sigma2gamma+1); break;
                  case (17): mcmc.print=(int)t;     break;
                  case (18): mcmc.burnin=(int)t;    break;
                  case (19): mcmc.sampfreq=(int)t;  break;
                  case (20): mcmc.nsample=(int)t;   break;
                  case (21):
                     sscanf(pline+1,"%lf%lf%lf%lf%lf", eps,eps+1,eps+2,eps+3,eps+4);
                     break;
               }
               break;
            }
         }
         if (i==nopt)
            { printf ("\noption %s in %s\n", opt, ctlf);  exit (-1); }
      }
      fclose(fctl);
   }
   else
      if (noisy) error2("\nno ctl file..");


   if(com.ndata>NGENE) error2("raise NGENE?");
   else if(com.ndata<=0) com.ndata=1;
   if(com.seqtype==0 || com.seqtype==2) { com.fix_omega=1; com.omega=0; }
   if(com.alpha==0)  { com.fix_alpha=1; com.nalpha=0; }
   if(mcmc.usedata==1) {
      if(com.alpha==0)
         com.plfun=lfun;
      else {
         if (com.ncatG<2 || com.ncatG>NCATG) error2 ("ncatG");
         com.plfun=lfundG;
      }
      if(com.model>HKY85)  error2("model.");
      if (com.model==JC69 || com.model==F81) { com.fix_kappa=1; com.kappa=1; }
   }

   return(0);
}



int DownSptreeSetTime (int inode)
{
/* This goes down the species tree, from the root to the tips, to specify the 
   initial node ages.  If the age of inode is not set already, it will 
   initialize it.
   This is called by GetInitials().
*/
   int j,ison;

   for (j=0; j<sptree.nodes[inode].nson; j++) {
      ison=sptree.nodes[inode].sons[j];
      if(sptree.nodes[ison].nson) {
         if(sptree.nodes[ison].age == -1)
            sptree.nodes[ison].age = sptree.nodes[inode].age*(.6+.4*rndu());
         DownSptreeSetTime(ison);
      }
   }
   return(0);
}

int GetRateGenes (FILE* fout)
{
/* This routine is for ad hoc determination of the mean lnrate at each locus, 
   \mu, which is the mean for the normal lnrate at the root.  The routine uses 
   the calibration nodes with bounds B or gamma G to calulate rough rates.  
   The lnrate at the root is drawn from N(mu, t0*sigma^2).

   PROBLEM: this uses com.alpha to calculate pairwise distances, while the 
   parameter is being estimated.
*/
   int inode, i,j,k, ipop, locus, nleft,nright,marks[NS], sons[2];
   int nf, nf_low, nf_up;  /* fossil bounds */
   double d,md, kappa, bigD=5, mr,mr_low,mr_up, *tf, t, maxage=0;

   if(mcmc.usedata==0) {
      for(i=0; i<data.ngene; i++) data.rgene[i]=0;
      return(0);
   }

   puts("Estimating rough mean rates at each locus...");
   fputs("\nProcessing fossil calibrations to get mean rates at each locus\n", fout);

   for(locus=0; locus<data.ngene; locus++) {
      UseLocus(locus, 0, mcmc.usedata, 0);
      nf=nf_low=nf_up=0; mr=mr_low=mr_up=0;
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

printf("Locus %2d node %2d (ipop %2d) left %2d right %2d",
       locus+1,inode, ipop, nleft,nright);

         for(i=0,md=0; i<com.ns; i++) {
            for(j=0; j<com.ns; j++) {
               if(marks[i]==1 && marks[j]==2) {
                   if(mcmc.usedata==2) {
                      for(k=i,d=0; k!=inode; k=nodes[k].father)
                         d += nodes[k].label;
                      for(k=j; k!=inode; k=nodes[k].father)
                         d += nodes[k].label;
                   }
                   else {
                      /* problematic, as com.alpha is unknown */
                      d = DistanceIJ(i, j, com.model, com.alpha, &kappa);
                      if(d<0) { puts("big distance"); d=bigD; }
                   }
                   md+=d;
               }
            }
         }
         md/=(nleft*nright);  /* average distance for the node */
         switch(sptree.nodes[ipop].fossil) {
         case (LOWER_F): { mr_low+=md/tf[0];          nf_low++; maxage=max2(maxage,tf[0]); break; }
         case (UPPER_F): { mr_up+=md/tf[1];           nf_up++;  maxage=max2(maxage,tf[1]); break; }
         case (BOUND_F): { mr+=md/((tf[0]+tf[1])/2);  nf++;     maxage=max2(maxage,(tf[0]+tf[1])/2); break; }
         case (GAMMA_F): { mr+=md/((tf[0]-1)/tf[1]);  nf++;     maxage=max2(maxage,(tf[0]-1)/tf[1]); break; }
         }

printf(" %s d = %9.5f u=%9.5f  ", fossils[sptree.nodes[ipop].fossil], md, mr);
FOR(i,com.ns) printf("%d", marks[i]); FPN(F0);

      }
      if(nf==0 && !(nf_low&&nf_up)) 
         { printf("need fossils for locus %d\n",locus+1); exit(-1); }
      if(nf)     mr/=nf; 
      if(nf_low) mr_low/=nf_low; 
      if(nf_up)  mr_up/=nf_up;

      if(nf_low && nf_up)           mr = (mr*nf+(mr_low+mr_up)/2)/(nf+1);
      else if(nf_low && mr>mr_low)  mr = (mr+mr_low)/2;
      else if(nf_low && mr<nf_up)   mr = (mr+mr_up)/2;

      data.rgene[locus]=log(mr);
      /* initial lnrate at the root of the species tree at each locus */
      t=sptree.nodes[sptree.root].age;
      sptree.nodes[sptree.root].lnrates[locus] 
         = data.rgene[locus] - data.sigma2[locus]*t/2 
                        + sqrt(data.sigma2[locus]*t)*rndnormal();
   }

   if(maxage>OldAge*.1) 
      error2("Fossil ages too large. Change time units in tree file?");
   for(i=0;i<data.ngene;i++) 
      fprintf(fout, "mr for locus %d =%9.5f log(mr) =%9.5f\n", i+1, exp(data.rgene[i]),data.rgene[i]);

   return(0);
}



int GetInitials(FILE*fout, double space[])
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
   for(j=sptree.nspecies; j<sptree.nnode; j++)  sptree.nodes[j].age=-1;
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
   if(com.clock) { /* times, rates for genes */
      np = sptree.nspecies-1 + g;
      for(i=0; i<g; i++) data.rgene[i]=rndu();
   }
   else {          /* times, rates, signma2, correl */
      np = sptree.nspecies-1 + g*sptree.nnode + g*(g+1)/2;

      /* sigma2 and correlation matrix in lnrates among loci */
      for(i=0; i<g; i++) {
         data.correl[i*g+i] = data.sigma2gamma[0]/data.sigma2gamma[1]*(.8+.2*rndu());
         for(j=0; j<i; j++)
            data.correl[i*g+j]=data.correl[j*g+i] = .5-rndu();
      }
      getSinvDetS(space);

      /* gets mu from data, and resets lnrates at root */
      GetRateGenes(fout);
      /* initial lnrates */
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
      if(!com.fix_omega) {
         data.omega[i] = .1+.2*rndu();
         np++;
      }
      if(!com.fix_alpha) {
         data.alpha[i] = .5+.2*rndu();
         np++;
      }

   }
   return(np);
}

int collectx (double x[])
{
/* this collects parameters into x[] for printing and summarizing.
   It returns the number of parameters.

     clock=1: times, rates for genes, kappa, alpha
     clock=0: times, rates or lnrates by node by gene, sigma2, rho_ij, kappa, alpha
*/
   int printlnrate=0, i,j, np=0, g=data.ngene;

   for(i=sptree.nspecies; i<sptree.nspecies*2-1; i++) 
      x[np++] = sptree.nodes[i].age;
   if(com.clock)
      for(i=0; i<g; i++) x[np++] = data.rgene[i];
   else {
      if(printlnrate) 
         for(i=0; i<g; i++)  for(j=0; j<sptree.nnode; j++)
            x[np++] = sptree.nodes[j].lnrates[i];
      else
         for(i=0; i<g; i++)  for(j=0; j<sptree.nnode; j++)
            x[np++] = exp(sptree.nodes[j].lnrates[i]);

      /* collect variance sigma^2 and correlation rho_ij into x */
      for(i=0; i<g; i++) {
         x[np++] = data.correl[i*g+i];
         for(j=0; j<i; j++)  x[np++] = data.correl[i*g+j];
      }
   }
   if(mcmc.usedata==1) {
      if(!com.fix_kappa)
         for(i=0; i<g; i++)
            x[np++]=data.kappa[i];

      if(!com.fix_omega)
         for(i=0; i<g; i++)
            x[np++]=data.omega[i];

      if(!com.fix_alpha)
         for(i=0; i<g; i++)
            x[np++]=data.alpha[i];
   }
   if(np!=com.np) 
      error2("np error in collectx.");
   return(0);
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

#define P0t_BD(expmlt) (rho*(lambda-mu)/(rho*lambda +(lambda*(1-rho)-mu)*expmlt))

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
      pdf = P0t_BD(expmlt);
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
   double t[NS-1], lnpt, Scale, expmlt=1, vt1, P0t1, cdf, cdfprev, small=1e-20;

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
      expmlt= exp((mu-lambda)*t[0]);
      P0t1  = P0t_BD(expmlt);
      vt1   = 1-P0t1/rho*expmlt;
   }
   else {
      P0t1 = rho/(1+rho*mu*t[0]);
      vt1  = mu*t[0]*P0t1;
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

   /* Adhockery, added 3 May 2004, in Edmonton */
   if(!sptree.nodes[sptree.root].fossil)
      lnpt += 2*log(P0t1*(1-vt1))+(n1-1)*log(vt1);

   lnpt += lnptC_Fossil();
   if(debug==1) printf("\npdf = %.12f\n", exp(lnpt));

   return (lnpt);
}


int LabelOldCondP (int spnode)
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
   node ages.  The age of sptree.root is capped by OldAge.
*/
   int  locus, is, i;
   double naccept=0, t, tnew, tb[2];
   double lnacceptance=0, lnLd=0, lnpDinew[NGENE], lnpTnew,lnpRnew;

   if(debug==2) puts("\nUpdateTimes ");
   for(is=sptree.nspecies; is<sptree.nnode; is++) {
      t=sptree.nodes[is].age;

      lnLd=0;   tb[0]=0;  tb[1]=OldAge;
      tb[0]=max2(sptree.nodes[sptree.nodes[is].sons[0]].age,
                 sptree.nodes[sptree.nodes[is].sons[1]].age);
      if(is!=sptree.root) 
         tb[1]=sptree.nodes[sptree.nodes[is].father].age;
      tnew=rnduab(t-finetune/2,t+finetune/2);
      tnew=reflect(tnew,tb[0],tb[1]);

      sptree.nodes[is].age=tnew;  /* the only thing to change in proposal */

      lnpTnew=lnpriorTimes();
      lnacceptance = lnpTnew-data.lnpT;
      if(!com.clock) {
         lnpRnew=lnpriorRates();
         lnacceptance+= lnpRnew-data.lnpR;
      }


      for(locus=0,lnLd=0; locus<data.ngene; locus++) {
         UseLocus(locus, 1, mcmc.usedata, 0);

         if(mcmc.saveconP) {
            FOR(i,sptree.nnode) com.oldconP[i]=1;
            LabelOldCondP(is);
         }
         if(com.oldconP[tree.root]) error2("strange: no need to update likelihood?");
         lnpDinew[locus] = lnpD_locus(locus);
         lnLd += lnpDinew[locus]-data.lnpDi[locus];

      }
      lnacceptance += lnLd;

      if(debug==2) printf("species %2d tb: %8.5f %8.5f t: %8.5f%8.5f %9.2f", is,tb[0],tb[1],t,tnew,lnLd);

      if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
         naccept++;
         data.lnpT=lnpTnew;
         if(!com.clock) data.lnpR=lnpRnew;
         FOR(locus,data.ngene) data.lnpDi[locus]=lnpDinew[locus];

         *lnL+=lnLd;
         if(mcmc.usedata==1) switchconPin();
         if(debug==2) printf(" Y (%4d)\n", NPMat);
      }
      else {
         sptree.nodes[is].age = t;
         if(debug==2) printf(" N (%4d)\n", NPMat);
         for(locus=0; locus<data.ngene; locus++)
            AcceptRejectLocus(locus, 0);  /* reposition conP */
      }
   }
   return(naccept/(sptree.nspecies-1.));
}


void getSinvDetS (double space[])
{
/* This uses the variance-correlation matrix data.correl[g*g] to constructs 
   Sinv (inverse of S) and detS (|S|).  It also copies the variances into 
   data.sigma2[g].  This is called every time data.correl is updated.

   What restrictions should be placed on data.correl[]???

   space[data.ngene]
*/
   int i, j, g=data.ngene;
   int debug=0;

   for(i=0; i<g; i++)
      data.Sinv[i*g+i] = data.sigma2[i] = data.correl[i*g+i];
   for(i=0; i<g; i++) {
      for(j=0; j<i; j++)
         data.Sinv[i*g+j] = data.Sinv[j*g+i]
            = data.correl[i*g+j]*sqrt(data.sigma2[i]*data.sigma2[j]);
   }

   if(debug) {
      printf("\ncorrel & S & Sinv ");
      matout(F0, data.correl, g, g);
      matout(F0, data.Sinv, g, g);
   }

   matinv (data.Sinv, g, g, space);
   data.detS = fabs(space[0]);

   if(debug) {
      matout(F0, data.Sinv, g, g);
      printf("|S| = %.6g\n", data.detS);
   }
   if(data.detS<0) 
      error2("detS < 0");

}

double lnpriorRates (void)
{
/* This calculates the log of the prior of lnrates under the geometric 
   Brownian motion model of rate evolution.  The algorithm cycles through all 
   nodes in the species tree and multiply the conditional g-variate normal 
   densities across the tree, accounting for correlations in lnrate between 
   genes.  There is no need for recursion down the species tree as the order 
   at which sptree.nodes are visited is unimportant.
   The lnrates are stored in sptree.nodes[].lnrates[].
   The variance-covariance matrix, S, is implicitly used through its inverse 
   data.Sinv[] and determinant data.detS.
   The root lnrate has mean data.rgene[locus] and variance sigma^2*root_age.
   There is no bias correction for the root lnrate.

   yd[ngene] = y - mu.
*/
   int i,j, inode, g=data.ngene, dad;
   double lnpR=0, t, zz, *y, yd[NGENE];

   for(inode=0; inode<sptree.nnode; inode++) {
      y=sptree.nodes[inode].lnrates;
      t = sptree.nodes[inode].age;
      if(inode==sptree.root) { /* prior mean lnrate at root is data.rgene[] */
         for(i=0; i<g; i++)
/*
            yd[i] = y[i] - data.rgene[i];
*/
            yd[i] = y[i] - data.rgene[i] + t*data.sigma2[i]/2;

      }
      else {
         dad=sptree.nodes[inode].father;
         t = sptree.nodes[dad].age - t;
         for(i=0; i<g; i++)
            yd[i] = y[i] - sptree.nodes[dad].lnrates[i] + t*data.sigma2[i]/2;
      }
      /* zz is the quadratic form in exponent. */
      for(i=0,zz=0; i<g; i++) {
         zz    += data.Sinv[i*g+i]*yd[i]*yd[i];
         for(j=0; j<i; j++)
            zz += data.Sinv[i*g+j]*yd[i]*yd[j]*2;
      }
      lnpR += -zz/(2*t) - g/2.*log(PI*2) - 1/2.*g*log(t) - log(data.detS)/2;
   }
   return lnpR;
}

double UpdateRates (double *lnL, double finetune)
{
/* This updates lnrates under the Brownian motion rate-evolution model.
   The algorithm cycles through the loci.  For each locus, it cycles through 
   the nodes in the species tree.

   Waste of computation: For every proposal to change lnrates, lnpriorRates() 
   is called to calculate the prior for all lnrates at all loci on the whole 
   tree, thus wasting computation.
*/
   int locus, inode, j, g=data.ngene;
   double naccept=0, lnpRnew, lnpDinew, lnacceptance, lnLd, lnRold;
   double e=finetune;

   if(com.clock)
      return UpdateRatesClock(lnL, finetune);

   if(debug==3) puts("\nUpdateRates");
   for(locus=0; locus<g; locus++) {
      for(inode=0; inode<sptree.nnode; inode++) {

         lnRold=sptree.nodes[inode].lnrates[locus];
         sptree.nodes[inode].lnrates[locus]   = rnduab(lnRold-e, lnRold+e);
         UseLocus(locus, 1, mcmc.usedata, 0);  /* copyconP=1 */

         if(mcmc.saveconP) {
            FOR(j,sptree.nspecies*2-1) com.oldconP[j]=1;
            LabelOldCondP(inode);
         }

         /* This priorRates calculation is wasteful but not expensive. */
         lnpRnew=lnpriorRates();
         lnacceptance = lnpRnew-data.lnpR;
         lnpDinew=lnpD_locus(locus);
         lnLd = lnpDinew - data.lnpDi[locus];
         lnacceptance +=  lnLd;

         if(lnacceptance>0 || rndu()<exp(lnacceptance)) {
            naccept++;
            if(mcmc.usedata==1) AcceptRejectLocus(locus,1);
   
            data.lnpR=lnpRnew;
            data.lnpDi[locus]=lnpDinew;
            *lnL += lnLd;
         }
         else {
            if(mcmc.usedata==1) AcceptRejectLocus(locus,0);
            sptree.nodes[inode].lnrates[locus]=lnRold;
         }
      }
   }
   return(naccept/(g*sptree.nnode));
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
   double naccept=0, rgene0, lnLd, lnpDinew, lnacceptance, c=1;

   if(debug==3) puts("\nUpdateRatesClock");
   if(mcmc.saveconP) FOR(j,sptree.nspecies*2-1) com.oldconP[j]=0;
   for(locus=0; locus<g; locus++) {
      rgene0=data.rgene[locus];
      data.rgene[locus] = rgene0 * (c=exp(finetune*(.5-rndu())));
  
      UseLocus(locus, 1, mcmc.usedata, 0);
      lnpDinew=lnpD_locus(locus);
      lnLd = lnpDinew - data.lnpDi[locus];   /* likelihood ratio */
      lnacceptance =  lnLd + log(c);         /* proposal ratio */
      /* prior ratio */
      lnacceptance += logPriorRatio(data.rgene[locus],rgene0,data.rgenegamma[0],data.rgenegamma[1],GAMMA);

      if(debug==3)
         printf("\nLocus %2d rgene %9.4f%9.4f %10.5f", locus+1, rgene0, data.rgene[locus], lnLd);

      if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
         naccept++;
         *lnL += lnLd;
         data.lnpDi[locus]=lnpDinew;
         if(mcmc.usedata==1) AcceptRejectLocus(locus,1);
         if(debug==3) printf(" Y\n");
      }
      else {
         data.rgene[locus]=rgene0;
         if(mcmc.usedata==1) AcceptRejectLocus(locus,0); /* reposition conP */

         if(debug==3) printf(" N\n");
      }
   }

   return(naccept/g);
}

double UpdateParameters (double *lnL, double finetune)
{
/* This updates parameters in the substitution model such as the ts/tv 
   rate ratio for each locus.

   Should we update the birth-death process parameters here as well?

*/
   int locus, j, ip, np=!com.fix_kappa+!com.fix_alpha;
   double naccept=0, lnLd,lnpDinew, lnacceptance, c=1;
   double pold, pnew, *gammaprior;

   if(debug==4) puts("\nUpdateParameters");
   if(np==0) return(0);

   if(mcmc.saveconP) FOR(j,sptree.nspecies*2-1) com.oldconP[j]=0;
   for(locus=0; locus<data.ngene; locus++) {
      for(ip=0; ip<np; ip++) {
         if(ip==0 && !com.fix_kappa) {  /* kappa */
            pold=data.kappa[locus];
            data.kappa[locus] = pnew = pold * (c=exp(finetune*(.5-rndu())));
            gammaprior=data.kappagamma;
         }
         else {  /* alpha */
            pold=data.alpha[locus];
            data.alpha[locus] = pnew = pold * (c=exp(finetune*(.5-rndu())));
            gammaprior=data.alphagamma;
         }

         UseLocus(locus, 1, mcmc.usedata, 0); /* this copies parameter from data.[] to com. */

         lnpDinew=lnpD_locus(locus);
         lnLd = lnpDinew - data.lnpDi[locus];
         lnacceptance =  lnLd + log(c);
         lnacceptance += logPriorRatio(pnew,pold,gammaprior[0],gammaprior[1],GAMMA);

         if(debug==4)
            printf("\nLocus %2d para%d %9.4f%9.4f %10.5f", locus+1,ip+1,pold,pnew,lnLd);

         if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
            naccept++;
            *lnL += lnLd;
            data.lnpDi[locus]=lnpDinew;
            if(mcmc.usedata==1) AcceptRejectLocus(locus,1);
            if(debug==4) printf(" Y\n");
         }
         else {
            if(ip==0 && !com.fix_kappa)
               data.kappa[locus]=pold;
            else 
               data.alpha[locus]=pold;
            
            if(mcmc.usedata==1) AcceptRejectLocus(locus,0);

            if(debug==4) printf(" N\n");
         }
      }
   }
   return(naccept/(data.ngene*np));
}


double UpdateParaRates (double finetune, double space[])
{
/* This updates the variance and correlation matrix data.correl[] in lnrates 
   among loci.  The hyperpriors are assumed to be the same across loci.  
   The routine uses the mean lnrates in data.rgene[], but do not change 
   data.rgene[].  The proposals in this routine do not change the likelihood.

   Should check what restrictions should be placed on data.correl[].

*/
   int i,j, g=data.ngene, nproposal=0;
   double old, r, lnacceptance, lnpRnew, naccept[2]={0}; /* for sqrt(S) */
   double rb[2]={-1, 1}, e[2]={0,.1};  /* e[0] for sigma2, e[1] for rho_ij */

   if(com.clock) return(0);
   if(debug==5) puts("\nUpdateParaRates");

   e[0]=finetune;
   for(i=0; i<g; i++) {
      for(j=0; j<=i; j++) {
         old=data.correl[i*g+j];
         r = rnduab(old-e[i!=j], old+e[i!=j]);
         if(i==j) { /* variance sigma^2 */
            data.correl[i*g+i] = fabs(r);
            lnacceptance = logPriorRatio(data.correl[i*g+i],old,data.sigma2gamma[0],data.sigma2gamma[1],GAMMA);
         }
         else { /* correlation coefficent: U(0,1) */
            r = reflect(r,rb[0],rb[1]);
            if(r<rb[0]||r>rb[1]) error2("r?");

            data.correl[i*g+j]=data.correl[j*g+i] = r;
            lnacceptance = 0; /* sliding window, uniform prior */
         }
         if(debug==5) printf("%s %d%d %9.5f -> %9.5f ", (i==j?"sig":"rho"), i,j, old,r);

         getSinvDetS(space);
         /* ratio of priors for rates */
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
      printf(" acceptance %.3f (out of %d) %.3f (out of %d)\n",
         naccept[0]/g,g,naccept[1]/(g*(g-1)/2.),g*(g-1)/2);
      getchar();
   }
   return((naccept[0]+naccept[1])/(g*(g+1)/2.));
}



double mixingClock (double finetune)
{
/* If com.clock==1, this multiplies all times by c and divides all rates 
   (data.rgene[]) by c, so that the likelihood does not change.  This algorithm
   does not really work when com.clock==0, because of the upper bound OldAge.  
*/
   int j;
   double naccept=0, a=data.rgenegamma[0], b=data.rgenegamma[1], c, lnc;
   double lnacceptance=0, lnpTnew, lnpRnew=-1;

   if(!com.clock) error2("strange");
   lnc=finetune*(rndu()-.5);  c=exp(lnc);

   for(j=sptree.nspecies; j<sptree.nnode; j++) 
      sptree.nodes[j].age*=c;

   for(j=0; j<data.ngene; j++) {
      lnacceptance = -(a-1)*lnc-b*(data.rgene[j]/c-data.rgene[j]);
      data.rgene[j]/=c;
   }
   lnpTnew=lnpriorTimes();
   lnacceptance += lnpTnew-data.lnpT;
   lnacceptance += ((sptree.nspecies-1) - data.ngene)*lnc;

   if(lnacceptance>0 || rndu()<exp(lnacceptance)) { /* accept */
      naccept=1;
      data.lnpT=lnpTnew;
      data.lnpR=lnpRnew;
   }
   else {   /* reject */
      for(j=sptree.nspecies; j<sptree.nnode; j++) sptree.nodes[j].age/=c;
      for(j=0; j<data.ngene; j++) data.rgene[j]*=c;
   }
   return(naccept);
}


double mixing (double *lnL, double finetune)
{
/* If com.clock==1, this calls mixingClock().  
   If com.clock==0, this add a random constant to all rates at a locus; the 
   likelihood has to be recalculated.
*/
   int locus, i;
   double naccept=0, lnc;
   double lnacceptance=0, lnLd, lnpRnew=-1, lnpDinew;

   if(com.clock)
      return mixingClock(finetune);

   if(mcmc.saveconP) FOR(i,sptree.nspecies*2-1) com.oldconP[i]=0;
   for(locus=0; locus<data.ngene; locus++) {
      UseLocus(locus, 1, mcmc.usedata, 0);  /* copycondP */

      lnc=finetune*(rndu()-.5)*2;
      for(i=0; i<sptree.nnode; i++) 
         sptree.nodes[i].lnrates[locus] -= lnc;

      lnpRnew=lnpriorRates();
      lnpDinew=lnpD_locus(locus);
      lnLd=lnpDinew-data.lnpDi[locus];
      lnacceptance = lnpRnew-data.lnpR + lnLd;

      if(lnacceptance>0 || rndu()<exp(lnacceptance)) { /* accept */
         naccept++;
         *lnL += lnLd;
         data.lnpDi[locus]=lnpDinew;
         data.lnpR=lnpRnew;
         if(mcmc.usedata==1) AcceptRejectLocus(locus,1);
      }
      else {   /* reject */
         for(i=0; i<sptree.nnode; i++)
            sptree.nodes[i].lnrates[locus] += lnc;
         if(mcmc.usedata==1) AcceptRejectLocus(locus,0);
      }
   }
   return(naccept/data.ngene);
}



static double a_Brownian, b_Brownian, c_Brownian;

double fun_Brownian(double x)
{
   double xt = -a_Brownian/2*x*x + b_Brownian*x + c_Brownian;
   if(xt>300) error2("xt in fun_Brownian too large?");
   return exp(xt);
}


double BranchLength (double t, double y0, double yt, int method, double sigma2)
{
   int i,n=1000000;
   double blength=0, x, a=t*sigma2, b=yt-y0+a, small=1e-20;

   if(method==ARITHMETIC)        /* arithmetic mean */
      blength = (exp(y0)+exp(yt))/2*t;
   else if (method==GEOMETRIC)   /* geometric mean */
      blength = exp((y0+yt)/2)*t;
   else if (method==BROWNIAN) {  /* brownian motion integral */
      a_Brownian=a;  b_Brownian=b;  c_Brownian=y0;
      blength = t*NIntegrateGaussLegendre(fun_Brownian, 0, 1, 10);
   }
   else {
      for(i=0; i<n; i++) {
         x=i/(double)n;
         blength += exp(-a/2*x*x + b*x + y0);
      }
      blength*=t/n;

/*
      sqrta=sqrt(a);
      blength = t*sqrt(2*PI/a)*exp(b*b/(2*a)+y0) 
              * (CDFNormal(sqrta-b/sqrta)-CDFNormal(-b/sqrta));
*/

   }
   return(blength);
}

int ReadBlengthVar(char infile[])
{
/* this reads the MLEs of branch lengths under no clock and their SEs, for 
   approximate calculation of sequence likelihood.  The branch lengths are 
   stored in nodes[].label, and the variances are in data.varb[locsu][].

   This also frees up memory for sequences.
*/
   FILE* fBV=gfopen(infile,"r");
   char line[100000];
   int lline=100000, locus, i,j,inode;
   double small=1e-20, v;

   if(noisy) printf("\n\nReading branch lengths and variances from %s.\n", infile);
   for(locus=0; locus<data.ngene; locus++) {
      UseLocus(locus, 0, 0, 0);
      NodeToBranch ();

      if(noisy) printf("\t\t%d branches.\n", tree.nbranch);
 
      for(;;) {
         if(fgets(line,lline,fBV)==NULL) error2("EOF?");
         if(strstr(line, "Locus")) { 
            printf("%s\n", line);
            fgets(line,lline,fBV);
            break; 
         }
      }
      for(i=0; i<tree.nbranch; i++) {
         fscanf(fBV, "%d..%d", &j, &inode);
         if(tree.branches[i][1] != inode-1) 
            error2("something wrong with node numbering.");
      }
      for(i=0; i<tree.nbranch; i++)
         fscanf(fBV, "%lf", &nodes[tree.branches[i][1]].label);
      for(i=0; i<tree.nbranch; i++) {
         if(fscanf(fBV, "%lf", &v)!=1) error2("error when reading SEs");
         data.varb[locus][tree.branches[i][1]] = v*v;
      }

      for(i=0; i<tree.nnode; i++)
         if(i!=tree.root) {
            v=data.varb[locus][i];
            printf("n%-2d b = %12.9f se&v = %12.9f %12.9f\n", i+1,nodes[i].label,v,sqrt(v));
         }
   }
   fclose(fBV);
   /* free up memory for sequences */
   for(locus=0; locus<data.ngene; locus++) {
      free(data.fpatt[locus]);
      for(j=0;j<data.ns[locus]; j++)
         free(data.z[locus][j]);
   }

   return(0);
}



int MCMC (FILE* fout)
{
   char mcmcout[32]="mcmc.out";
   FILE *fmcmcout=gfopen(mcmcout,"w");
   int nsteps=4+!com.clock, nxpr[2]={8,3};
   int i,j,k, ir, g=data.ngene;
   double Paccept[5]={0}, lnL=0, nround=0, *x, *mx, *vx;
   char timestr[32];

   starttime();
   noisy=3;

   mcmc.saveconP=1;
   if(mcmc.usedata!=1) mcmc.saveconP=0;
   FOR(j,sptree.nspecies*2-1) com.oldconP[j]=0;
   GetMem();

   if(mcmc.usedata==2)
      ReadBlengthVar("in.BV");

   printf("\n%d burnin, sampled every %d, %d samples\n", 
           mcmc.burnin, mcmc.sampfreq, mcmc.nsample);
   if(mcmc.usedata) puts("Approximating posterior");
   else             puts("Approximating prior");

   printf("(Settings: cleandata=%d.  print=%d  saveconP=%d.)\n", 
          com.cleandata, mcmc.print, mcmc.saveconP);

   com.np=GetInitials(fout, com.space);
   x=(double*)malloc(com.np*(com.np+2)*sizeof(double));
   if(x==NULL) error2("oom in mcmc");
   mx=x+com.np;  vx=mx+com.np;

   printSptree();
   collectx(x);
   if(!com.fix_kappa && !com.fix_alpha && data.ngene==2) { nxpr[0]=6; nxpr[1]=4; }

   puts("\npriors: ");
   if(mcmc.usedata==1) {
      if(!com.fix_kappa) printf("G(%.4f, %.4f) for kappa\n", data.kappagamma[0], data.kappagamma[1]);
      if(!com.fix_omega) printf("G(%.4f, %.4f) for omega\n", data.omegagamma[0], data.omegagamma[1]);
      if(!com.fix_alpha) printf("G(%.4f, %.4f) for alpha\n", data.alphagamma[0], data.alphagamma[1]);
   }
   if(com.clock)      printf("G(%.4f, %.4f) for rgene\n", data.rgenegamma[0], data.rgenegamma[1]);
   else               printf("G(%.4f, %.4f) for sigma2\n", data.sigma2gamma[0], data.sigma2gamma[1]);

   printf("finetune parameters: ");
   for(j=0;j<nsteps;j++) printf(" %7.4f",mcmc.finetune[j]); FPN(F0);
   printf("\nStarting MCMC... ");
   printf("Initial parameters (np = %d):\n", com.np);
   for(j=0;j<com.np;j++) printf(" %5.3f",x[j]); FPN(F0);

   /* calculates prior for times and likelihood for each locus */
   data.lnpT=lnpriorTimes();
   for(j=0; j<data.ngene; j++) com.rgene[j]=-1;
   if(!com.clock) {
      data.lnpR=lnpriorRates();
   }
   lnL=lnpData(data.lnpDi);
   printf("\nlnL0 = %9.2f\n", lnL);
   printf("paras: times, rates, sigma2, correl, kappa, alpha\n");

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
      Paccept[3] = (Paccept[3]*nround + mixing(&lnL, mcmc.finetune[3]))/(nround+1);
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
         FOR(j,nxpr[1]) printf(" %5.3f", mx[com.np-nxpr[1]+j]);
         if(mcmc.usedata) printf(" %4.1f ",lnL);
         printf(" %4.1f %4.1f",data.lnpT, data.lnpR);
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
   puts("\nSpecies tree showing node ages, for TreeView");
   FOR(j,sptree.nspecies-1) sptree.nodes[sptree.nspecies+j].age=mx[j];
   copySptree();
   FPN(F0); OutaTreeN(F0,1,1); FPN(F0);
   if(com.clock==0) {
      printf("Posterior means of variance-correlation matrix:");
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
/*
      DescriptiveStatistics(fout, mcmcout, 20, 20);
*/
   }
   FreeMem();
   return(0);
}
