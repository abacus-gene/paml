/* BASEMLG.c 
   ML parameter estimation for models with Gamma-distributed rates
   over sites, combined with tree topology estimation from DNA sequences

                   Copyright, Ziheng YANG, July 1992 onwards
                      gcc -o basemlg -O2 basemlg.c tools.o
                        basemlg <ControlFileName>
*/

#define NS       9
#define NTREE    20    
#define NBRANCH  (NS*2-2)      
#define NNODE    (NS*2-1) 
#define NGENE    4
#define LSPNAME  30
#define NCODE    4
#define NP       (NBRANCH+NGENE+5)

#include "tools.h"
extern char NUCs[];
extern int noisy, NFunCall, *ancestor;
extern double *SeqDistance;

int Forestry (FILE *fout, double space[]);
int GetOptions (char *ctlf);
int testx (double x[], int np);
int GetInitials (double x[], double space[]);
int TestFunction (FILE *fout, double x[], double space[]);
int GetMem (int nbranch, int nR, int nitem);
void GetSave (int nitem, int *nm, int M[], double alfa, double c);
double GetBmx (int ninter, int nsum, int im, int M[], int nodeb[]);
double RhoRate (double x[]);
double lfunG_print (double x[], int np);
double lfunG (double x[], int np);
int lfunG_d (double x[], double *lnL, double dl[], int np);
int lfunG_dd (double x[], double *lnL, double dl[], double ddl[], int np);

struct CommonInfo {
   char *z[NS], spname[NS][LSPNAME+1], seqf[32],outf[32],treef[32];
   int  seqtype, ns, ls, ngene, posG[NGENE+1],lgene[NGENE],*pose,npatt,*fpatt;
   int  clock, given_alfa, given_kapa, given_rgene, method, print, verbose;
   int  model;
   int  runmode, escape, np, ntime, nrate, nrgene, ncode;
   double pi[6], lmax, alfa, kapa, rgene[NGENE], piG[NGENE][6], *chunk;
   double *SSave, *ErSave, *EaSave;
   int  nhomo, nparK, ncatG, given_rho, getSE;  /* unused */
   double rho;                           /* unused */
}  com;
struct TREEB {
   int  nbranch, nnode, origin, branches[NBRANCH][2];
   double lnL;
}  tree;                
struct TREEN {
   int father, nson, sons[NS], ibranch;
   double branch, divtime, *lkl;
}  nodes[2*NS-1];

static int nR=4, CijkIs0[64];
static double Cijk[64], Root[4];

FILE *frub, *flfh;
char *models[]={"JC69", "K80", "F81", "F84", "HKY85", "TN93", "REV", "UNREST"};
enum {JC69, K80, F81, F84, HKY85, TN93, REV, UNREST} MODELS;

/*
#define REV_UNREST
*/
#define BASEMLG
#include "treesub.c"

int main(int argc, char *argv[])
{
   char ctlf[32]="baseml.ctl";
   double  *space, kapat=0;
   FILE   *fout, *fseq;

   noisy=2;
   com.runmode=0;  com.escape=0;
   com.clock=0;    com.given_rgene=0;
   com.method=1;          /* 1: newton, 0: ming1 */

   com.model=F84;
   com.given_kapa=0;      com.kapa=5;
   com.given_alfa=0;      com.alfa=0.2;
   com.print=1;      /* 0: nothing, 1: log(f) & r,  2: r for all sites */
   com.ncode=4;

   frub=fopen ("rub","w");   flfh=fopen ("lfh", "w");

   if ((space=(double*)malloc(50000*sizeof(double)))==NULL) error ("oom");

   if (argc>1) { strcpy(ctlf, argv[1]); printf ("\nctlfile is %s.\n", ctlf); }
   GetOptions (ctlf);
   if (com.clock||com.model==REV)  com.method=0;

   if ((fout=fopen (com.outf, "w"))==NULL) error("outfile creation err.");
   if((fseq=fopen (com.seqf,"r"))==NULL)  error("No sequence file!");
   ReadSeq (NULL, fseq, 0);

   fprintf (fout,"BASEMLG %15s %8s + Gamma", com.seqf, models[com.model]);
   if (com.clock) fprintf (fout, ", Clock");
   if (com.given_kapa)  fprintf (fout,"kapa =%7.3f given\n", com.kapa);
   if (com.ngene>1) fprintf (fout, " (%d genes)  ", com.ngene);

   SeqDistance=(double*)malloc(com.ns*(com.ns-1)/2*sizeof(double));
   ancestor=(int*)malloc(com.ns*(com.ns-1)/2*sizeof(int));
   if (SeqDistance==NULL||ancestor==NULL) error("oom");

   Initialize (fout, space, 0);
   if (com.model==JC69) PatternJC69like (fout);

   DistanceMatNuc (fout, com.model, com.alfa, &kapat);
   if (com.model<=HKY85)
      EigenTN93 (com.model, com.kapa, com.kapa, com.pi, &nR, Root, Cijk);

   if (com.runmode==0)
      Forestry (fout, space);
   else 
      error ("tree search not available with basemlg.");
   fclose(fseq);
   putchar ('\a');
   return (0);
}

/* x[] is arranged in order of t[ntime], rgene[], kapa and alfa (if any) */

int Forestry (FILE *fout, double space[])
{
   int  status=0, i,j, itree, ntree, np;
   double x[NBRANCH+1], lnL[NTREE], e=1e-5;
   double *var=space+NP;
   FILE *ftree, *fInitVal=NULL;

   if ((ftree=fopen(com.treef,"r"))==NULL)   error("treefile err.");
   fscanf (ftree, "%d%d", &i, &ntree);
   if (i!=com.ns || ntree>NTREE) error ("err:ns || ntree..");              
   fprintf (flfh,"%6d%6d%6d\n", ntree, com.ls, com.npatt);

   fInitVal=fopen("InitVal","r");

   for (itree=0; itree<ntree; itree++) {
      printf ("\nTREE # %2d\n", itree+1 );
      fprintf (fout,"\nTREE # %2d:  ", itree+1);
      fprintf (flfh,"\n\n%2d\n", itree+1);
      fprintf (frub,"\n\nTREE # %2d", itree+1);
      com.print*=-1;

      if (ReadaTreeN(ftree, &i, 1)) error ("err tree..");
      OutaTreeN (F0, 0, 0);   
      OutaTreeN (fout, 0, 0);  
      fflush (fout);  fflush (flfh);

      i=(com.clock||com.model>HKY85?1+(com.print!=0):2+!com.given_alfa);
      GetMem (tree.nbranch, nR, i);

      GetInitials (x, space);
      np = com.np;
      printf ("\nnp =%6d", np);
      NFunCall=0;
/*
      TestFunction (fout, x, space);
*/
      if (fInitVal) {
         printf ("\n\nReading initial values from InitVal.. OK?");
         FOR (j, np)
            if (fscanf(fInitVal,"%lf",&x[j])!=1) error ("err InitVal");
      }
      else {
         if (itree==0) printf("\a");
         printf ("\n\nSuggest using BASEML to get initial values..");
         FOR (j, com.ntime) fprintf (fout,"%9.5f", x[j]);
      }

      matout (F0, x, 1, np);
      printf ("\nlnL0=%12.6f", -lfunG(x, np));

      if (com.method==1)
         i = Newton(frub, &lnL[itree], lfunG, lfunG_dd,testx,x,var,e,np);
      else {
         i = ming1 (frub, &lnL[itree], lfunG,
         ((com.clock||com.model>HKY85)?NULL:lfunG_d),testx,x,space,e,np);
         Hessian (np,x,lnL[itree],space,var,lfunG,var+np*np);
      }
      if (noisy) printf ("\nOut...\nlnL:%14.6f\n", -lnL[itree]);

      if (com.method==0) matinv (var, np, np, var+np*np);

      if (i || j<np) { status=-1; fprintf (fout, "\n\ncheck convergence.."); }
      fprintf(fout,"\nlnL(np:%3d):%14.6f%+12.6f\n", com.np,
         -lnL[itree], -lnL[itree]+lnL[0]);
      OutaTreeB (fout);  FPN (fout);

      FOR (i, np) fprintf (fout,"%9.5f", x[i]);   FPN (fout);
      FOR (i, np)
         fprintf (fout,"%9.5f", var[i*np+i]>0.?sqrt(var[i*np+i]):0.);

      fprintf (fout, "\n\n# fun calls:%10d\n", NFunCall);

      com.print*=-1;
      if (com.print)  lfunG_print (x, np);
/*
      RhoRate (x);
*/
   }        /* for (itree) */
   fclose (ftree);

   exit (status);
   return (0);
}

int testx (double x[], int np)
{
   int i,k;
   double tb[]={1e-5,4}, rgeneb[]={0.01,20}, kapab[]={0,80}, alfab[]={0.01,9};

   if (com.clock && SetBranch (x, 0)) return (-1);
   FOR (i,com.ntime)   if (x[i]<tb[0] || x[i]>tb[1])   return (-1);
   if (np==com.ntime)  return (0); 
   for (i=0,k=com.ntime; i<com.nrgene; i++,k++) 
      if (x[k]<rgeneb[0] || x[k]>rgeneb[1])   return (-1);
   for (i=0; i<com.nrate; i++,k++)
      if (x[k]<kapab[0] || x[k]>kapab[1])    return(-1);
   if (!com.given_alfa && (x[np-1]<alfab[0] || x[np-1]>alfab[1]))  return(-1);
   return (0);
}

int GetInitials (double x[], double space[])
{
   int i,j;
   double t;

   com.nrgene = (!com.given_rgene)*(com.ngene-1);
   com.nrate=0;
   if (com.model==REV) {
      com.nrate=5;
      x[com.ntime+com.nrgene] = 1;
      FOR (i,com.nrate-1) x[com.ntime+com.nrgene+i+1]=1/com.kapa;
   }
   else if (!com.given_kapa)
      { com.nrate=1; x[com.ntime+com.nrgene]=com.kapa; }
   if (com.model<=HKY85)
      EigenTN93 (com.model, com.kapa, com.kapa, com.pi, &nR, Root, Cijk);

   com.np = com.ntime+com.nrgene+com.nrate+(!com.given_alfa);
   FOR (j,com.ntime) x[j]=0.05+0.01*(double) (com.ntime-j);
   FOR (j,com.nrgene) x[com.ntime+j]=1;

   if (!com.given_alfa)  x[com.np-1]=com.alfa;
   LSDistance (&t, x, testx);

   FOR (j, com.ntime) if (x[j]<1e-5) x[j]=1e-4;
   return (0);
}

int TestFunction (FILE* fout, double x[], double space[])
{
   int i, np=com.np;
   double lnL;

   printf ("\ntest functions\n");
   SetSeed (23);
   FOR (i,np) x[i]=(double)i*rndu()+0.0001;
   matout (F0, x, 1, np);    matout (fout, x, 1, np);

   lnL=lfunG(x, np);
   printf ("\n\nnp:%6d\nlnL:%12.8f\n", np, lnL);
   fprintf (fout, "\n\nnp:%6d\nlnL:%12.8f\n", np, lnL);

   printf ("\ndl and gradient");     fprintf (fout, "\ndl and gradient");
   lfunG_d (x, &lnL, space, np);
   printf ("\nlnL:%12.8f", lnL);     fprintf (fout, "\nlnL:%12.8f", lnL);
   matout (F0, space, 1, np);        matout (fout, space, 1, np);
   gradient (np, x, lnL, space, lfunG, space+np, 0);
   printf ("\nlnL:%12.8f", lnL);     fprintf (fout, "\nlnL:%12.8f", lnL);
   matout (F0, space, 1, np);        matout (fout, space, 1, np);
   
   printf ("\nddl & Hessian");       fprintf (fout, "\nddl & Hessian");
   lfunG_dd (x, &lnL, space, space+np, np);
   printf ("\nlnL:%12.8f", lnL);     fprintf (fout, "\nlnL:%12.8f", lnL);
   matout (F0, space, 1, np);        matout (fout, space, 1, np);
   matout2 (F0, space+np, np, np, 8, 3);
   matout2 (fout, space+np, np, np, 8, 3);
   fflush (fout);
   Hessian (np, x, lnL, space, space+np, lfunG, space+np+np*np);
   printf ("\nlnL:%12.8f", lnL);     fprintf (fout, "\nlnL:%12.8f", lnL);
   matout2 (F0, space+np, np, np, 8, 3);
   matout2 (fout, space+np, np, np, 8, 3);

   exit (0);
}

int GetMem (int nbranch, int nR, int nitem)
{
/* ns=4: 98KB,  5: 1.6MB,  6: 25MB,  7: 402MB    (for HKY85, 3/3/93)
   nitem=1:ErSave; 2:SSave & ErSave; 3:SSave & ErSave & EaSave 
*/
   static int size=0;
   int nm, j;

   for(j=0,nm=1; j<nbranch; j++)   nm*=nR;
   if (nm*nitem>size && size) free (com.chunk);
   if (nm*nitem>size) {
      size=nm*nitem;
      printf ("\n\nSave %12d bytes\n", size*sizeof(double));
      com.chunk = (double*)malloc(size*sizeof(double));
      if (com.chunk==NULL) error("oom");
   }
   com.ErSave  = com.chunk;
   if (nitem>1) com.SSave=com.ErSave+nm;
   if (nitem>2) com.EaSave=com.SSave+nm;
   return(0);
}

void GetSave (int nitem, int *nm, int M[], double alfa, double c)
{
/* correct for both clock=0 and 1
   nitem=1:ErSave; 2:SSave & ErSave; 3:SSave & ErSave & EaSave 
*/
   int im, j, it;
   double S;

   for(j=0,*nm=1; j<tree.nbranch; j++)   *nm*=nR;
   for (im=0; im< *nm; im++) {
      for (j=0,it=im,S=0; j<tree.nbranch; j++) {
         M[j]=it%nR;    it/=nR;
         if (M[j]) S+=nodes[tree.branches[j][1]].branch*Root[M[j]];
      }
      com.ErSave[im] = pow(alfa/(alfa-c*S),alfa);
      if (nitem>1) com.SSave[im]  = S;
      if (nitem>2) com.EaSave[im] = log(alfa/(alfa-c*S))-c*S/(alfa-c*S);
   }
}

double GetBmx (int ninter, int nsum, int im, int M[], int nodeb[])
{
   int isum, j, it, i, nR4=nR*4;
   double y, Bmx;

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
   return (Bmx);
}

double RhoRate (double x[])
{
/* rate factors for all possible site patterns, and the correlation 
   coefficient (Yang and Wang, in press).
*/
   int  i,j, h, nodeb[NNODE], M[NBRANCH], accurate=(com.ns<8);
   int  im, nm, nsum=1<<(2*(tree.nnode-com.ns)), *fobs;
   int  ninter=tree.nbranch-com.ns+1;
   double lnL, lexp, fh, sumfh=0, rh, Bmx, alfa, kapa;
   double mrh=0, mrh0=0, vrh=0, vrh0=0;

   if (com.ngene>1) error ("ngene>1");
   alfa=(com.given_alfa ? com.alfa : x[com.np-1]);
   kapa=(com.given_kapa ? com.kapa : x[com.ntime]);

   fprintf (flfh, "\n\nmodel:%6d\n  kapa:%9.4f\nalfa:%9.4f\nBranches",
            com.model, kapa, alfa);
   matout (flfh, x, 1, tree.nbranch);
   if (com.model<=HKY85 && !com.given_kapa)  
       RootTN93 (com.model, kapa, kapa, com.pi, &rh, Root);
#ifdef REV_UNREST
   if (com.model==REV)
       EigenREV (NULL, x+com.ntime+com.nrgene, com.pi, &nR, Root, Cijk);
#endif
   if (SetBranch (x, 0)) puts ("\nx[] err..");
   GetSave (2, &nm, M, alfa, 1);
   fobs = (int*) malloc ((1<<(2*com.ns)) * sizeof (int));

   FOR (h,(1<<2*com.ns)) fobs[h]=0;
   for (h=0; h<com.npatt; h++) {
      for (j=0,im=0; j<com.ns; j++) im = im*4+(com.z[j][h]-1);
      fobs[im]=com.fpatt[h];
   }
   for (h=0,lnL=0,lexp=0; h<(1<<2*com.ns); h++) {
      if (accurate==0 && fobs[h]==0) continue;
      for (j=0,im=h; j<com.ns; nodeb[com.ns-1-j]=im%4,im/=4,j++);
      for (im=0,fh=0,rh=0; im<nm; im++) {
         Bmx=GetBmx (ninter, nsum, im, M, nodeb);
         fh += Bmx*com.ErSave[im];
         rh += Bmx*com.ErSave[im]*alfa/(alfa-com.SSave[im]);
      }  /* for (im) */
      if (fh<=0)  printf ("\a\nRhoRate: h=%4d  fh=%9.4f \n", h, fh);
      rh /= fh;      vrh += fh*rh*rh;
      sumfh += fh;   mrh += fh*rh;   mrh0+=rh*(double)fobs[h]/(double)com.ls;
      lexp  -= log(fh)*fh;
      if (fobs[h]) {
         vrh0 += rh*rh*(double)fobs[h]/(double)com.ls;
         lnL  -= log(fh)*(double)fobs[h];
      }
      if (com.print>3)  {
         fprintf (flfh,"\n%6d%9.2f%8.4f  ", fobs[h], fh*com.ls, rh);
         FOR (i,com.ns) fprintf (flfh, "%c", NUCs[nodeb[i]]);
      }
   }    /* for (h) */
   vrh-=1;     vrh0-=mrh0*mrh0;
   fprintf (flfh, "\n%s Vrh", accurate?"accurate":"approximate");
   fprintf (flfh,"\nsumfh = 1 =%12.6f      mrh = 1 =%12.6f\n", sumfh, mrh);
   fprintf (flfh, "\nVr :%12.6f  Vr0:%12.6f   mrh0:%12.6f", vrh, vrh0, mrh0);
   fprintf (flfh, "\nPEV:%12.6f      %12.6f", 1/alfa-vrh, 1/alfa-vrh0);
   fprintf (flfh, "\nRHO:%12.6f      %12.6f", sqrt(vrh*alfa), sqrt(vrh0*alfa));
   fprintf (flfh,"\nLn(L)=%12.6f\nlentropy=%12.6f\n", -lnL, -lexp);
   free (fobs);
   return (accurate ? sqrt(vrh*alfa) : sqrt(vrh0*alfa));
}

double lfunG_print (double x[], int np)
{
   int  i,j, h, igene, lt, nodeb[NNODE], M[NBRANCH], allsites=1;
   int  im, nm, nsum=1<<(2*(tree.nnode-com.ns));
   int  ninter=tree.nbranch-com.ns+1;
   double lnL, fh, rh, mrh0=0,vrh0=0, Bmx, y, alfa,kapa, *rates=NULL;
 
   if (com.print<2) allsites=0;
   if (allsites && (rates=(double*)malloc(com.npatt*sizeof(double)))==NULL)
           error ("oom");

   FOR (i, com.nrgene) com.rgene[i+1]=x[com.ntime+i];
   kapa=(com.given_kapa ? com.kapa : x[com.ntime+com.nrgene]);
   alfa=(com.given_alfa ? com.alfa : x[com.np-1]);

   if (com.model<=HKY85 && !com.given_kapa)  
       RootTN93 (com.model, kapa, kapa, com.pi, &y, Root);
#ifdef REV_UNREST
   if (com.model==REV)
       EigenREV (NULL, x+com.ntime+com.nrgene, com.pi, &nR, Root, Cijk);
#endif
   if (SetBranch (x, 0)) puts ("\nx[] err..");
   for(j=0,nm=1; j<tree.nbranch; j++)   nm*=nR;

   for (h=0,lnL=0,igene=-1,lt=0; h<com.npatt; lt+=com.fpatt[h++]) {
      FOR(j,com.ns) nodeb[j] = com.z[j][h]-1;
      if (h==0 || lt==com.lgene[igene])
         GetSave (2, &nm, M, alfa, com.rgene[++igene]);
      for (im=0,fh=0,rh=0; im<nm; im++) {
         Bmx=GetBmx (ninter, nsum, im, M, nodeb);
         fh += Bmx*com.ErSave[im];
         rh += Bmx*com.ErSave[im]*alfa/(alfa-com.SSave[im]);
      }  /* for (im) */
      if (fh<=0)  printf ("\a\nlfunG_print: h=%4d  fh=%9.4f \n", h, fh);
      rh/=fh;
      vrh0+=rh*rh*(double)com.fpatt[h]/(double)com.ls;
      mrh0+=rh*(double)com.fpatt[h]/(double)com.ls;
      if (allsites) rates[h]=rh;

      fh=log (fh);
      lnL -= fh*(double)com.fpatt[h];
      if (com.print>=1) {
         fprintf (flfh,"\n%4d%6d%14.8f%9.2f%9.4f  ",
            h+1, com.fpatt[h], fh, com.ls*exp(fh), rh);
         FOR (i,com.ns) fprintf (flfh, "%c", NUCs[com.z[i][h]-1]);
      }
   }    /* for (h) */
   if (com.print>=2) 
       fprintf (flfh,"\n\nVrh0:%12.6f\nmrh0:%12.6f\nRHO0:%12.6f\n", 
                (vrh0-=mrh0*mrh0), mrh0, sqrt(vrh0*alfa));
   if (allsites) {
      fprintf (flfh, "\n\nrate factors along sequence.. \n");
      FOR (h, com.ls) {
         FOR (j, com.ns) fprintf (flfh, "%c", NUCs[com.z[j][com.pose[h]]-1]);
         fprintf (flfh,"%8d%6d%9.4f\n",h+1,com.pose[h]+1,rates[com.pose[h]]);
      }
   }
   return (lnL);
}

double lfunG (double x[], int np)
{
/* likelihood with spatial rate variation
*/
   int  i,j, h, igene, lt, nodeb[NNODE], M[NBRANCH];
   int  im, nm, nsum=1<<(2*(tree.nnode-com.ns));
   int  ninter=tree.nbranch-com.ns+1;
   double lnL, fh, rh, Bmx, y, alfa,kapa;

   NFunCall++;
   FOR (i, com.nrgene) com.rgene[i+1]=x[com.ntime+i];
   kapa=(com.given_kapa ? com.kapa : x[com.ntime+com.nrgene]);
   alfa=(com.given_alfa ? com.alfa : x[np-1]);
   if (com.model<=HKY85 && !com.given_kapa)  
       RootTN93 (com.model, kapa, kapa, com.pi, &y, Root);
#ifdef REV_UNREST
   if (com.model==REV)
       EigenREV (NULL, x+com.ntime+com.nrgene, com.pi, &nR, Root, Cijk);
#endif
   if (SetBranch (x, 0)) puts ("\nx[] err..");

   for (h=0,lnL=0,igene=-1,lt=0; h<com.npatt; lt+=com.fpatt[h++]) {
      FOR(j,com.ns) nodeb[j] = com.z[j][h]-1;
      if (h==0 || lt==com.lgene[igene])
         GetSave (1, &nm, M, alfa, com.rgene[++igene]);
      for (im=0,fh=0,rh=0; im<nm; im++) {
         Bmx=GetBmx (ninter, nsum, im, M, nodeb);
         fh += Bmx*com.ErSave[im];
      }  /* for (im) */
      if (fh<=0)  printf ("\a\nlfunG: h=%4d  fh=%9.4f \n", h, fh);
      lnL -= log(fh)*(double)com.fpatt[h];
   }    /* for (h) */
   return (lnL);
}

int lfunG_d (double x[], double *lnL, double dl[], int np)
{
   int  nbranch=tree.nbranch, ninter=nbranch-com.ns+1;
   int  i,j, nodeb[NNODE], M[NBRANCH], h, igene,lt;
   int  im, nm, nsum=1<<(2*(tree.nnode-com.ns));
   double fh, y, Bmx, S, alfa, kapa, Er, dfh[NP];
   double *p=com.pi, drk1[4], c=1;

   if (com.clock||com.model>HKY85) error ("err lfunG_d");
   NFunCall++;
   FOR (i, com.nrgene) com.rgene[i+1]=x[com.ntime+i];
   kapa=(com.given_kapa ? com.kapa : x[com.ntime+com.nrgene]);
   alfa=(com.given_alfa ? com.alfa : x[np-1]);
   if (SetBranch (x, 0)) puts ("\nx[] err..");
   if (!com.given_kapa) {
      RootTN93 (com.model, kapa, kapa, p, &S, Root);
      y=p[0]*p[1]+p[2]*p[3];
      drk1[1]=2*y*S*S;
      drk1[2]=2*(y-p[5]*p[5])*p[4]*S*S;
      drk1[3]=2*(y-p[4]*p[4])*p[5]*S*S;
      if (com.model==F84) {
         y=2*p[0]*p[1]/p[4]+2*p[2]*p[3]/p[5];
         drk1[1]=y*S*S;
         drk1[2]=-S+(1+kapa)*y*S*S;
      }
   }
   *lnL=0; zero(dl,np);
   for (h=0,igene=-1,lt=0; h<com.npatt; lt+=com.fpatt[h++]) {
      FOR(j,com.ns) nodeb[j] = com.z[j][h]-1;
      if (h==0 || lt==com.lgene[igene])
         GetSave (2+!com.given_alfa, &nm, M, alfa, (c=com.rgene[++igene]));
      for (im=0,fh=0,zero(dfh,np); im<nm; im++) {
         Bmx=GetBmx (ninter, nsum, im, M, nodeb);
         S = com.SSave[im];
         Er = com.ErSave[im]*Bmx;
         fh += Er;
         FOR (j, nbranch)
            if (M[j]) dfh[j]+=c*Er*alfa*Root[M[j]]/(alfa-c*S);       /* t */
         if (com.ngene>0 && igene>0)
            dfh[com.ntime+igene-1]+=Er*alfa/(alfa-c*S)*S;            /* c */
         if (!com.given_kapa) {
            for (j=0,y=0; j<nbranch; j++) if (M[j]) y+=x[j]*drk1[M[j]];
            dfh[com.ntime+com.nrgene] += Er*alfa/(alfa-c*S)*y*c;     /* k */
         }
         if (!com.given_alfa)  dfh[np-1] += Er*com.EaSave[im];       /* a */
      }  /* for (im) */
      if (fh<=0)  printf ("\a\nlfunG_d: h=%4d  fh=%9.4f \n", h, fh);
      *lnL -= log(fh)*(double)com.fpatt[h];
      FOR (j, np) dl[j] -= dfh[j]/fh*(double)com.fpatt[h];
   }    /* for (h) */
   return(0);
}

int lfunG_dd (double x[], double *lnL, double dl[], double ddl[], int np)
{
   int  nbranch=tree.nbranch, ninter=nbranch-com.ns+1;
   int  i,j,k, nodeb[NNODE], M[NBRANCH], h, igene,lt;
   int  im, nm, nsum=1<<(2*(tree.nnode-com.ns));
   double fh, y, Bmx,S,alfa,kapa, Er,Ea, y1,y2;
   double dfh[NP], ddfh[NP*NP], *p=com.pi, drk1[4], drk2[4], c=1, s1=0,s2=0;

   if (com.clock||com.model>HKY85) error ("err lfunG_dd");
   NFunCall++;
   FOR (i, com.nrgene) com.rgene[i+1]=x[com.ntime+i];
   kapa=(com.given_kapa ? com.kapa : x[com.ntime+com.nrgene]);
   alfa=(com.given_alfa ? com.alfa : x[np-1]);
   if (SetBranch (x, 0)) puts ("\nx[] err..");

   if (!com.given_kapa) {
      RootTN93 (com.model, kapa, kapa, p, &S, Root);

      y=p[0]*p[1]+p[2]*p[3];
      drk1[1]=2*y*S*S;
      drk1[2]=2*(y-p[5]*p[5])*p[4]*S*S;
      drk1[3]=2*(y-p[4]*p[4])*p[5]*S*S;

      drk2[1]=-8*y*y*S*S*S;
      drk2[2]=-8*y*p[4]*(y-p[5]*p[5])*S*S*S;
      drk2[3]=-8*y*p[5]*(y-p[4]*p[4])*S*S*S;

      if (com.model==F84) {
         y=2*p[0]*p[1]/p[4]+2*p[2]*p[3]/p[5];
         drk1[1]=y*S*S;
         drk1[2]=-S+(1+kapa)*y*S*S;

         drk2[1]=-2*y*y*S*S*S;
         drk2[2]=2*y*S*S*(1-(1+kapa)*y*S);
      }
   }
   *lnL=0, zero(dl,np), zero(ddl,np*np);
   for (h=0,igene=-1,lt=0; h<com.npatt; lt+=com.fpatt[h++]) {
      FOR(j,com.ns) nodeb[j] = com.z[j][h]-1;
      if (h==0 || lt==com.lgene[igene])
         GetSave (2+!com.given_alfa, &nm, M, alfa, (c=com.rgene[++igene]));
      for (im=0,fh=0,zero(dfh,np),zero(ddfh,np*np); im<nm; im++) {
         Bmx=GetBmx (ninter, nsum, im, M, nodeb);
         S = com.SSave[im];
         Er = com.ErSave[im]*Bmx;
         y=1./(alfa-c*S);
         fh += Er;
         y1=Er*alfa*y;
         for (j=0,y2=y1*(1+alfa)*y*c*c; j<nbranch; j++) {
            if (M[j]) {
               dfh[j]+=y1*Root[M[j]]*c;                           /* t  */
               for (i=j; i<nbranch; i++)                          /* tt */
                  if (M[i]) ddfh[i*np+j]+=y2*Root[M[i]]*Root[M[j]];
            }
         }
         if (!com.given_kapa) {
            for (j=0,s1=0,s2=0; j<nbranch; j++)
               if (M[j]) { s1+=x[j]*drk1[M[j]];  s2+=x[j]*drk2[M[j]]; }
            k=com.ntime+com.nrgene;
            dfh[k] += c*y1*s1;                                    /* k  */
            for (j=0,y2=y*s1*(alfa+1); j<nbranch; j++)  if (M[j])
               ddfh[k*np+j] += c*y1*(Root[M[j]]*y2*c+drk1[M[j]]); /* kt */
            ddfh[k*np+k] += y1*c*(c*s1*y2+s2);                    /* kk */
         }
         if (com.ngene>0 && igene>0) {
            k=com.ntime+igene-1;
            dfh[k]+=y1*S;                                         /* c  */
            ddfh[k*np+k]+=y*y1*S*S*(1+alfa);                      /* cc */
            y2=alfa*y*y1*(1+c*S);
            FOR (j,nbranch)  ddfh[k*np+j]+=y2*Root[M[j]];         /* ct */
            if (!com.given_kapa)
               ddfh[(com.ntime+com.nrgene)*np+k]+=y2*s1;          /* ck */
         }
         if (!com.given_alfa) {
            Ea = com.EaSave[im];
            dfh[np-1] += Er*Ea;                                   /* a  */
            ddfh[np*np-1]+=Er*(Ea*Ea+1/alfa-y+c*S*y*y);           /* aa */
            y2=Er*y*(Ea*alfa-c*S*y);
            FOR (j,nbranch)  ddfh[(np-1)*np+j]+=c*Root[M[j]]*y2;  /* at */
            if (com.ngene>0 && igene>0)
               ddfh[(np-1)*np+nbranch+igene-1]+=S*y2;             /* ac */
            if (!com.given_kapa)
               ddfh[(np-1)*np+com.ntime+com.nrgene] += c*y2*s1;   /* ak */
         }
      }  /* for (im) */
      if (fh<=0) printf ("\a\nlfunG_dd: h=%4d  fh=%9.4f \n", h, fh);
      *lnL -= log(fh)*(double)com.fpatt[h];
      FOR (j, np) dl[j] -= dfh[j]/fh*(double)com.fpatt[h];
      FOR (j, np) for (i=j; i<np; i++) {
         ddl[i*np+j] -=
             (fh*ddfh[i*np+j]-dfh[i]*dfh[j])/(fh*fh)*(double)com.fpatt[h];
         ddl[j*np+i] = ddl[i*np+j];
      }
   }    /* for (h) */
   return(0);
}

int GetOptions (char *ctlf)
{
   int i, nopt=23, lline=255;
   char line[255], *pline, opt[20], comment='*';
   char *optstr[] = {"seqfile","outfile","treefile", "noisy", "verbose", "runmode",
        "escape", "clock", "given_rgene", "Mgene", "nhomo", "getSE", "RateAncestor",
        "model", "given_kappa", "kappa", "given_alpha", "alpha", "Malpha", "ncatG", 
        "given_rho", "rho", "nparK"};
   double t;
   FILE  *fctl=fopen (ctlf, "r");

   if (fctl) {
      if (noisy) printf ("\n\nReading options from %s..\n", ctlf);
      for (;;) {
         if (fgets (line, lline, fctl) == NULL) break;
         for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
            if (isalnum(line[i]))  { t=1; break; }
         if (line[0]==comment || t==0) continue;
         sscanf (line, "%s%*s%lf", opt, &t);
         if ((pline=strstr(line, "= "))==NULL) error ("option file.");

         for (i=0; i<nopt; i++) {
            if (strncmp(opt, optstr[i], 8)==0)  {
               if (noisy>2)
                  printf ("\n%3d %15s | %-20s %6.2f", i+1,optstr[i],opt,t);
               switch (i) {
                  case ( 0): sscanf(pline+2, "%s", com.seqf);    break;
                  case ( 1): sscanf(pline+2, "%s", com.outf);    break;
                  case ( 2): sscanf(pline+2, "%s", com.treef);   break;
                  case ( 3): noisy=(int)t;           break;
                  case ( 4): com.verbose=(int)t;     break;
                  case ( 5): com.runmode=(int)t;     break;
                  case ( 6): com.escape=(int)t;      break;
                  case ( 7): com.clock=(int)t;       break;
                  case ( 8):                         break;  /* given_rgene */
                  case ( 9): com.given_rgene=(int)t; break;
                  case (10): com.nhomo=(int)t;       break;
                  case (11): com.getSE=(int)t;       break;
                  case (12): com.print=(int)t;       break;
                  case (13): com.model=(int)t;       break;
                  case (14): com.given_kapa=(int)t;  break;
                  case (15): com.kapa=t;             break;
                  case (16): com.given_alfa=(int)t;  break;
                  case (17): com.alfa=t;             break;
                  case (18):                         break;  /* Malfa */
                  case (19): com.ncatG=(int)t;       break;
                  case (20): com.given_rho=(int)t;   break;
                  case (21): com.rho=t;              break;
                  case (22): com.nparK=(int)t;       break;
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
      if (noisy) printf ("\nno ctl file..");

   if (com.runmode>0)  puts ("tree search using BASEMLG?");
   if (com.alfa==0 || com.rho>0 || com.nhomo>0 || com.nparK>0 
      || com.model>HKY85)
       error ("\noptions in file baseml.ctl inappropriate.. giving up..");

   if (com.model!=F84 && com.kapa<=0)  error("init kapa..");
   if (com.alfa<=0) error("init alfa..");
   if (com.model==JC69 || com.model==F81) { com.given_kapa=1; com.kapa=1; }

   return (0);
}
