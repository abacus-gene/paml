/* MCML.c
   consistency analysis and estimation of sampling errors 
   by Monte-Carlo simulation.  HKY85+dG or F84+dG and simpler models
   are assumed, for Yang (1994, Systematic Biology, 43:329-342) and later.

                 gcc -o mcml -O2 mcml.c tools.o
                            mcml
*/

#define NS            6

#define NBRANCH       (NS*2-3)
#define NNODE         (NS*2-2)
#define NCATG         4
#define NCODE         4
#define NP            (NBRANCH+2)

#include "tools.h"
int testx (double x[], int np);
int GetInitials (double x[], double space[]);
int CheckModel (void);
double lfunG (double x[], int np);
double lfundG (double x[], int np);
double lfun (double x[], int np);
int PartialLikelihood (int inode);
int RealPi (void);
int DistanceMatFreq (FILE *fout, int model, double alfa);
int ModelTree4 (void);
int ReadPexp (char* pexpf);
int OutPexp (char* pexpf);
int ReadaTreeBranches (char* treef, double x[]);

struct CommonInfo {
   char spname[NS][10];
   int ns, npatt, *fpatt, ncatG, np, ntime, ncode, *pattM;
   int clock, given_kapa, given_alfa, model, get_fobs, runmode, escape;
   double pi[6], lmax, kapa, alfa, *fobs, *chunk, *fhK;
   double freqK[NCATG], rK[NCATG], *ErSave, (*plfun)(double x[],int np);
}  com;
struct TREEB {
   int  nbranch, nnode, origin;
   int  branches[NBRANCH][2];
   double lnL;
}  tree;
struct TREEN {
   int father, nson, sons[NS], ibranch;
   double branch, divtime, *lkl;
}  nodes[2*NS-1];

extern int noisy, NFunCall;
FILE *frst;
char *models[]={"JC69","K80","F81","F84","HKY85", "TN93"};
enum {JC69, K80, F81, F84, HKY85, TN93} MODELS;
static int nR=4, CijkIs0[64];
static double Cijk[64], Root[4];

#define NODESTRUCTURE
#define LSDISTANCE
#define TREESEARCH
#include "treesub.c"

int main(int argc)
{
   FILE *fout;
   int ns=4, modeltrue=F84, readpexp=0, outpexp=1, readatree=0;
   char *pexpfin="pexp.in", *pexpfout="pexp.out", *treef="trees.mc";
   int rt1=0, i,j, ir=0,nr=50, dG=1, itree=0,ntree=2, npattMAX, ipi=1;
   int nl=2, il, ls[]={0, 250, 500};
   double alfa=0, kapa=5, x[NP], Pt[4], *pexp, *space;
   double pi0[][6]={{0.25, 0.25, 0.25, 0.25, 0.5, 0.5},
                     {0.1,  0.2,  0.3,  0.4,  0.3, 0.7}};   /* T C A G Y R */
   double t0[][NBRANCH] = { {0.1, 0.1, 0.5, 0.2, 1.0},
                            {0.1, 1.0, 0.1, 1.0, 0.1}};

   fout=fopen ("mcml.out", "w");   frst=fopen ("rst","w"); 
   noisy=0;
   com.ns=ns;      npattMAX=1<<(2*com.ns);
   com.clock=0;    com.ncatG=4;
   com.runmode=2;  com.escape=0;
   com.ncode=4;

   /* initilize and allocate memory */
   i= npattMAX*(sizeof(int) + sizeof(double)*(4+(com.ns-1)*com.ncode));
   if ((com.chunk=(double*)malloc(i))==NULL) error ("oom");
   printf ("\n# pattern:%10d\n%10d bytes allocated\n", npattMAX, i);
   com.fhK=com.ErSave=(com.chunk+(com.ns-1)*com.ncode*npattMAX);
   com.fobs=com.fhK+npattMAX;  pexp=com.fobs+npattMAX;  space=pexp+npattMAX;
   com.pattM=com.fpatt=(int*)(space+npattMAX);

   for (i=0; i<com.ns; i++)  sprintf(com.spname[i], "seq.%d", i+1);
   SetSeed (117567);

   for (itree=0; itree<ntree; itree++) {

      /* get the tree structure & branch lengths */

      com.get_fobs=1; com.npatt=npattMAX;
      for (i=0; i<npattMAX; i++) com.pattM[i]=i;
      if (readatree)  ReadaTreeBranches (treef, x);
      else  {
         ModelTree4 ();
         FOR (i, tree.nbranch) x[i]=t0[itree][i];
      }
      printf ("\n\nModel tree and branch lengths: ");  
      OutaTreeN (F0, 0, 0);       FPN (F0);  
      OutaTreeB (F0);             FPN (F0);  
      FOR (i, com.ntime) fprintf (F0,"%9.5f", x[i]);    FPN (F0);
      OutaTreeN (fout, 0, 0);     FPN (fout);
      OutaTreeB (fout);           FPN (fout); 
      FOR (i, com.ntime) fprintf (fout,"%9.5f", x[i]);  FPN (fout);      
      PointLklnodes ();

      if (readpexp) { 
         ReadPexp (pexpfin);  
         fprintf (fout, "\nTrue model unknown");
      }
      else {    /* the true model */
         com.model=modeltrue;    xtoy(pi0[ipi], com.pi, 6);
         com.given_kapa=1;       com.kapa=kapa;
         com.given_alfa=1;       com.alfa=alfa;
         CheckModel ();
         com.plfun=(com.alfa==0?lfun:(dG?lfundG:lfunG));

         fprintf (fout, "\nTrue model: %s", models[com.model]);
         if (alfa) fprintf (fout, "+Gamma  alfa=%7.2f  dG=%4d", alfa, dG);
         fprintf (fout, "\nkapa =%7.4f\npi[]", com.kapa);
         matout (fout, com.pi, 1, 4);

         com.lmax=-(*com.plfun) (x, com.np);
      }
      if (outpexp) OutPexp (pexpfout);

      com.get_fobs=0;
      DistanceMatFreq (fout, com.model, com.alfa);
      xtoy (com.fobs, pexp, npattMAX);

      fprintf (fout, "\nSimulation results\n");
      for (il=0; il<nl; il++) {   /* consistency analysis and simulation */
         printf ("\nls =%5d\n", ls[il]);
         fprintf (fout, "\nls =%6d ", ls[il]);
         if (com.ns>4 && ls[il]>0) continue;
         FOR (j,4) Pt[j]=0;
         for (ir=0; ir<(ls[il]==0?1:nr); ir++) {
            if (ls[il]==0) {
               xtoy (pexp, com.fobs, npattMAX);
               for (i=0,com.npatt=npattMAX; i<npattMAX; i++) com.pattM[i]=i;
	    }
            else {
               MultiNomial (ls[il], npattMAX, pexp, com.fpatt, space);
               for (i=0,com.npatt=0; i<npattMAX; i++) 
                  if(com.fpatt[i]>0) {
                     com.fobs[com.npatt]=(double)com.fpatt[i]/ls[il];
                     com.pattM[com.npatt++]=i; 
		  }
	    }

 /*** the model for analysis, change the following 3 lines ***/
            com.model=F84;
            com.given_alfa=1;  com.alfa=0;
            com.given_kapa=0;  com.kapa=kapa;

            if (com.model>K80) RealPi ();
            com.plfun = (com.alfa==0?lfun:lfundG);
            rt1=TreeSearch (NULL, space);   Pt[rt1]++; 
            printf ("\r%4d%5d%6d  P:", ls[il], ir+1, com.npatt);
            FOR (j,4) printf ("%7.1f", Pt[j]/(ir+1.)*100.);
	 }  /* for (ir) */
         FPN (F0);
         FOR (j,4) fprintf (fout, "%7.1f", Pt[j]/ir*100.);
         fprintf (fout, "%7d\n", ir);
         fflush (fout);
      }       /* for (il)  */
   }   /* for (itree) or (ik) */
   return (0);
}

int testx (double x[], int np)
{
   int i;
   double tb[]={1e-4, 99}, kapab[]={1e-4, 99}, alfab[]={0.02, 5};

   if (com.clock && SetBranch (x, 0))   return (-1);
   FOR (i,com.ntime)  if (x[i]<tb[0] || x[i]>tb[1]) return (-1);
   if (np==com.ntime) return(0);
   if (!com.given_kapa)
      if (x[com.ntime]<kapab[0] || x[com.ntime]>kapab[1])  return (-1);
   if (!com.given_alfa)
      if (x[np-1]<alfab[0] || x[np-1]>alfab[1]) return (-1);
   return (0);
}

int CheckModel (void)
{
   if (com.model==JC69 || com.model==F81)  { com.given_kapa=1; com.kapa=1; }
   if (com.model==JC69 || com.model==K80) 
      { fillxc (com.pi, 0.25, 4);   fillxc (com.pi+4, 0.5, 2); }
   return (0);
}

int GetInitials (double x[], double space[])
{
   int j;
   double t;

   CheckModel ();
   com.np=com.ntime+(!com.given_kapa)+(!com.given_alfa);

   FOR (j,com.ntime) x[j]=0.05;

   LSDistance (&t, x, testx, space);

   FOR (j, com.ntime) x[j]=max(x[j], 0.001);
   if (!com.given_kapa)  x[com.ntime]=com.kapa;
   if (!com.given_alfa)  x[com.np-1]=com.alfa;
   return (0);
}

static int Base[NS];
static double PMat[16], Qfactor, kapa1, kapa2;

double lfunG (double x[], int np)
{
/* likelihood with gamma rates at sites, for all site patterns and 
   without using GetMem() and GetSave()
*/
   int  i,j, h, nodeb[NNODE], M[NBRANCH];
   int  im, nm, isum,nsum=1<<(2*(tree.nnode-com.ns)), it,nR4;
   int  ninter=tree.nbranch-com.ns+1;
   double lnL, fh, y, Bmx, S, alfa, *p=com.pi;

   alfa=(com.given_alfa ? com.alfa : x[com.np-1]);
   if (!com.given_kapa)  com.kapa=x[com.ntime];
   if (com.model==F84) {kapa1=1+com.kapa/p[4]; kapa2=1+com.kapa/p[5]; }
   else                 kapa1=kapa2=com.kapa;

   EigenTN93 (com.model, kapa1, kapa2, com.pi, &nR, Root, Cijk);
   nR4=nR*4;

   if (SetBranch (x, 0)) puts ("\nx[] err..");

   for(j=0,nm=1; j<tree.nbranch; nm*=nR,j++) ;
   for (im=0; im<nm; im++) {
      for (j=0,it=im,S=0; j<tree.nbranch; j++) {
         M[j]=it%nR;    it/=nR;
         if (M[j]) S+=nodes[tree.branches[j][1]].branch*Root[M[j]];
      }
      com.ErSave[im] = pow(alfa/(alfa-S),alfa);
   }
   for (h=0,lnL=0; h<com.npatt; h++) {
      for (j=0,it=com.pattM[h]; j<com.ns; nodeb[com.ns-1-j]=it%4,it/=4,j++) ;
      for (im=0,fh=0; im<nm; im++) {
         for (j=0,it=im; j<tree.nbranch;  M[j]=it%nR,it/=nR,j++) ;
         for (isum=0,Bmx=0; isum<nsum; isum++) {
            for (j=0,it=isum; j<ninter; nodeb[com.ns+j]=it&3,it>>=2,j++) ;
            for (j=0,y=com.pi[nodeb[tree.origin]]; j<tree.nbranch; j++) {
               i = nodeb[tree.branches[j][0]]*nR4
                 + nodeb[tree.branches[j][1]]*nR + M[j];
               if (CijkIs0[i]) goto nextone;
               y *= Cijk[i];
            }
            Bmx += y;
            nextone: ;
         }
         fh += Bmx*com.ErSave[im];
      }  /* for (im) */
      if (fh<=0)  printf ("\nlfunG: h=%4d  fh=%9.4f", h, fh);
      if (com.get_fobs) com.fobs[h]=fh;
      lnL -= log(fh)*fh;
   }    /* for (h) */
   return (lnL);
}


double lfundG (double x[], int np)
{
/* discrete-gamma model of variable rates over sites
   using com.fobs[], rather than com.fpatt[]
*/
   int  h, ir, i;
   double lnL, fh, *p=com.pi;

   NFunCall++;
   if (!com.given_kapa) com.kapa=x[com.ntime];
   if (!com.given_alfa)  com.alfa=x[com.np-1];
   DiscreteGamma (com.freqK,com.rK,com.alfa,com.alfa,com.ncatG,0);
   if (SetBranch(x, 0)) puts ("branch len err");

   if (com.model==F84) {kapa1=1+com.kapa/p[4]; kapa2=1+com.kapa/p[5]; }
   else                 kapa1=kapa2=com.kapa;
   Qfactor=1/(2*p[0]*p[1]*kapa1+2*p[2]*p[3]*kapa2 + 2*p[4]*p[5]);

   FOR (i, com.npatt) com.fhK[i]=0;
   for (ir=0; ir<com.ncatG; ir++) {
      FOR (i, tree.nnode)
         nodes[i].branch *= (ir==0?com.rK[ir]:com.rK[ir]/com.rK[ir-1]);
      PartialLikelihood (tree.origin);
      for (h=0; h<com.npatt; h++) {
         for (i=0,fh=0; i<4; i++)
            fh += com.pi[i]*nodes[tree.origin].lkl[h*4+i];
         com.fhK[h] += com.freqK[ir]*fh;
      }
   }
   FOR (i,tree.nnode)  nodes[i].branch /= com.rK[com.ncatG-1];

   for (h=0,lnL=0; h<com.npatt; h++) {
      if (com.fhK[h]<=0) {
         printf ("\nlfun: h=%4d  fhK=%9.4f", h, com.fhK[h]);
         matout (F0, x, 1, com.np);
      }
      if (com.get_fobs) com.fobs[h]=com.fhK[h];
      lnL -= log(com.fhK[h])*com.fobs[h];
   }
   return (lnL);
}

double lfun (double x[], int np)
{
/* using com.fobs[]
*/
   int  h,i;
   double lnL, fh, *p=com.pi;

   NFunCall++;

   if (SetBranch(x, 0)) puts ("branch len err");
   if (!com.given_kapa) com.kapa=x[com.ntime];
   if (com.model==F84) {kapa1=1+com.kapa/p[4]; kapa2=1+com.kapa/p[5]; }
   else                 kapa1=kapa2=com.kapa;
   Qfactor=1/(2*p[0]*p[1]*kapa1+2*p[2]*p[3]*kapa2 + 2*p[4]*p[5]);

   PartialLikelihood (tree.origin);
   for (h=0,lnL=0; h<com.npatt; h++) {
      for (i=0,fh=0; i<4; i++)
         fh += com.pi[i]*nodes[tree.origin].lkl[h*4+i];
      if (fh<=0) 
         printf ("\nlfun: h=%4d  fh=%9.4f", h, fh);
      if (com.get_fobs) com.fobs[h]=fh;
      lnL -= log(fh)*com.fobs[h];
   }
   return (lnL);
}

int PartialLikelihood (int inode)
{
   int it, i,j,k,h, ison;
   double t;

   FOR (i,nodes[inode].nson)
      if (nodes[nodes[inode].sons[i]].nson>0)
         PartialLikelihood (nodes[inode].sons[i]);
   fillxc (nodes[inode].lkl, (double)(inode>=com.ns), com.npatt*4);
   FOR (i, nodes[inode].nson) {
      ison=nodes[inode].sons[i];
      if (com.model<=K80)
         PMatK80 (PMat, nodes[ison].branch, com.kapa);
      else {
         t=nodes[ison].branch*Qfactor;
         PMatTN93 (PMat, t*kapa1, t*kapa2, t, com.pi);
      }
      for (h=0; h<com.npatt; h++) {
         for (j=0,it=com.pattM[h]; j<com.ns; Base[com.ns-1-j]=it%4,it/=4,j++) ;
         if (nodes[ison].nson<1)
            FOR (j,4) nodes[inode].lkl[h*4+j]*=PMat[j*4+Base[ison]];
         else
            FOR (j,4) {
               for (k=0,t=0; k<4; k++)
                  t += PMat[j*4+k] * nodes[ison].lkl[h*4+k];
               nodes[inode].lkl[h*4+j] *= t;
            }
      }    /*  for (h)  */
  }        /*  for (i)  */
  return (0);
}

extern double SeqDistance[];

int DistanceMatFreq (FILE *fout, int model, double alfa)
{
/* using com.fobs[] and com.ns
*/
   int i,j,k, it, h, b[NS];
   double x[16], kapa;

   if (fout) fprintf (fout,"\nDistance matrix, model:%s\n", models[model]);
   FOR (i, com.ns) {
      FOR (j, i) {
         for (h=0,zero(x,16); h<com.npatt; h++) {
            for (k=0,it=h; k<com.ns; b[com.ns-1-k]=it%4,it/=4,k++) ;
            x[b[i]*4+b[j]] += com.fobs[h];
         }
         SeqDistance [i*(i-1)/2+j] = SeqDivergence (x, model, alfa, &kapa);
         if (fout) fprintf (fout, "%10.5f", SeqDistance[i*(i-1)/2+j]);
      }
       if (fout) FPN (fout);
   }
   return (0);
}

int RealPi (void)
{
   int h, j,it;

   FOR (j,6) com.pi[j]=0;
   FOR (h,com.npatt) 
      for (j=0,it=com.pattM[h]; j<com.ns; j++) 
         { com.pi[it%4]+=com.fobs[h]/com.ns; it/=4; }
   com.pi[4]=com.pi[0]+com.pi[1];   com.pi[5]=com.pi[2]+com.pi[3];
   if (fabs(com.pi[4]+com.pi[5]-1)>1e-6) {
      matout (F0, com.pi, 1, 6);
      error ("SUM pi!=1 in RealPi");
   }
   return(0);
}

int ModelTree4 (void)
{
/* get the tree structure ((12)34)
*/
   int i,j, branches[][2]={ {4,5}, {5,0}, {5,1}, {4,2}, {4,3}};

   FOR (i,5) FOR (j,2) tree.branches[i][j]=branches[i][j];
   tree.nnode = (tree.nbranch=5) + 1;
   tree.origin=4;
   com.ntime = com.clock ? tree.nnode-com.ns+(tree.origin<com.ns)
                         : tree.nbranch;
   BranchToNode ();
   com.np=com.ntime;

   return (0);
}

int ReadaTreeBranches (char* treef, double x[])
{
/* get the tree structure and branch lengths
*/ 
   FILE *ftree;
   int i,j;

   if (com.clock) error ("clock not supported in ReadaTreeBranches");
   if ((ftree=fopen(treef,"r"))==NULL) error ("no tree file");

   fscanf (ftree, "%d%d", &i, &tree.nbranch);
   if (com.ns!=i)  error ("change com.ns");
   
   FOR (j, tree.nbranch)  FOR (i,2) {
      if (fscanf(ftree, "%d", &tree.branches[j][i])!=1)  error("err tree");
      tree.branches[j][i]--;
   }
   tree.origin=tree.branches[0][0];
   tree.nnode = tree.nbranch+1;
   com.ntime = com.clock ? tree.nnode-com.ns+(tree.origin<com.ns)
                         : tree.nbranch;
   BranchToNode ();
   com.np=com.ntime;
   for (j=0; j<com.ntime; j++)
      if (fscanf (ftree, "%lf", &x[j])!=1) error ("err: branch lengths");
   return (0);
}

int ReadPexp (char* pexpf)
{
   FILE *fpexp;
   int h, j, ch, it, nc=4, npatt;

   if ((fpexp=fopen(pexpf, "r")) == NULL) error ("ReadPexp");
   fscanf (fpexp, "%d%d", &j, &npatt);
   printf ("\nRead site-pattern probabilities from %s\n", pexpf);
   printf ("\nns:%3d  npattern:%4d", j, npatt);
   if (j!=com.ns) error ("pexp file error");
   if (npatt!=com.npatt)
      { puts ("npattern != 4**ns in ReadPexp..");  zero(com.fobs, com.npatt);}
   for (h=0; h<npatt; h++) {
      for (j=0,it=0; j<com.ns; j++) {
         ch=fgetc(fpexp);
         while (!isalpha(ch)) ch=fgetc(fpexp);
         it = it*nc + changeNUC((char)ch)-1;
      }
      if (fscanf (fpexp, "%lf", &com.fobs[it]) != 1) error ("ReadPexp");
/*
      printf ("\n%6d%20.10f", it+1, com.fobs[it]);
*/
   }
   if (fabs(sum(com.fobs,com.npatt)-1)>1e-6) error("Sum pexp!=1 in ReadPexp");
   fclose (fpexp);
   return (0); 
}

int OutPexp (char* pexpf)
{
   FILE *fpexp;
   int h, it, j, nc=4;

   printf ("\nOutputting site-pattern probabilities into %s\n", pexpf);
   fpexp=fopen(pexpf, "w");
   fprintf (fpexp, "%9d%9d\n\n", com.ns, com.npatt);

   for (h=0; h<com.npatt; h++) {
      for (j=0,it=h; j<com.ns; it/=nc,j++) Base[com.ns-1-j]=it%nc;
      for (j=0; j<com.ns; j++)  fprintf (fpexp, "%1c", NUCs[Base[j]]);
      fprintf (fpexp, "%20.8e\n", com.fobs[h]);
   }
   fclose (fpexp);
   return (0); 
}

