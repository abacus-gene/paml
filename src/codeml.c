/* codeml.c  (aaml.c & codonml.c)

   Maximum likelihood parameter estimation for codon sequences (seqtype=1) 
                    or amino-acid sequences (seqtype=2)
                Copyright, Ziheng YANG, 1993-2003

               cc -o codeml -fast codeml.c tools.o -lm
                         codeml <ControlFileName>
*/
/*
#define NSSITESBandits
#define NSSITES_K1_K2_CLASSES

#define DSDN_MC  1
#define DSDN_MC_SITES  1
*/

#include "paml.h"

#define NS            500
#define NBRANCH       (NS*2-2)
#define NNODE         (NS*2-1)
#define NGENE         10000
#define LSPNAME       30
#define NCODE         64
#define NCATG         40

#define NP            (NBRANCH*2+NGENE-1+2)
/*
#define NP            (NBRANCH+NGENE-1+189+2)
*/
extern char BASEs[],AAs[];
extern int noisy, NFunCall, NEigenQ, NPMatUVRoot, *ancestor, GeneticCode[][64];
extern double *SeqDistance;

int Forestry (FILE *fout);
int sortwM3(double x[]);
void DetailOutput(FILE *fout, double x[], double var[]);
int GetOptions (char *ctlf);
int testx (double x[], int np);
int SetxBound (int np, double xb[][2]);
int SetxInitials (int np, double x[], double xb[][2]);
int GetInitials (double x[], int*fromfile);
double *PointKappa (double xcom[], int igene);
double *PointOmega (double xcom[], int igene, int inode, int isiteclass);
int SetParameters (double x[]);
int SetPGene (int igene, int _pi, int _UVRoot, int _alpha, double x[]);
int SetPSiteClass(int iclass, double x[]);
int PMatJC69like (double P[], double t, int n);
int setmark_61_64 (void);
int printfcode (FILE *fout, double fb61[], double space[]);
int InitializeCodon (FILE *fout, double space[]);
int CodonListall(char codon[], int nb[3], int ib[3][4]);
int AA2Codonf (double faa[20], double fcodon[]);
int DistanceMatAA (FILE *fout);
int DistanceMatNG86 (FILE *fout, FILE*fds, FILE*fdn, FILE*dt, double alpha);
int GetDaa(FILE *fout, double daa[]);
void getpcodonClass(double x[], double pcodonClass[]);
int EigenQc(int getstats, double blength, double *S, double *dS, double *dN,
    double Root[], double U[], double V[],
    double kappa[], double omega, double Q[]);
int EigenQaa(FILE *fout, double Root[], double U[], double V[],double rate[]);
int Qcodon2aa(double Qc[], double pic[], double Qaa[], double piaa[]);
int SetAA1STEP(void);
int GetOmegaAA(int OmegaAA[]);
int TestModelQc(FILE *fout, double x[]);
double lfun2dSdN (double x[], int np);
int VariancedSdN(double t, double omega, double vtw[2*2], double vdSdN[2*2]);
int GetCodonFreqs(double pi[]);
int PairwiseCodon(FILE *fout, FILE*fds, FILE*fdn, FILE*dt, double space[]);
int PairwiseAA (FILE *fout, FILE *f2AA);
int lfunNSsites_rate (FILE* fout, double x[], int np);
void marksitesNonsyn(void);
double GetBranchRate(int igene, int ibrate, double x[], int *ix);
int GetPMatBranch (double Pt[], double x[], double t, int inode);
int ConditionalPNode (int inode, int igene, double x[]);
double CDFdN_dS(double x,double par[]);
int DiscreteNSsites(double par[]);
void finishup(void);
int mergeSeqs(FILE*fout);
void Get4foldSites(void);

struct common_info {
   char *z[NS], *spname[NS], seqf[96],outf[96],treef[96],daafile[96], cleandata;
   char oldconP[NNODE];       /* update conP for nodes? to save computation */
   int seqtype, ns, ls, ngene, posG[NGENE+1], lgene[NGENE], npatt,*pose;
   int runmode,clock,verbose,print, codonf,aaDist,model,NSsites;
   int nOmega, nbtype, nOmegaType;  /* branch partition, AA pair (w) partition */
   int method, icode, ncode, Mgene, ndata;
   int fix_rgene,fix_kappa,fix_omega,fix_alpha,fix_rho,nparK,fix_blength,getSE;
   int np, ntime, nrgene, nkappa, nrate, nalpha, ncatG, sspace, sconP1;
   int npi0, hkyREV, readpattf;
   double *fpatt, *space, kappa,omega,alpha,rho,rgene[NGENE];
   double pi[NCODE], pi_sqrt[NCODE], piG[NGENE][64],fb61[64];
   double f3x4MG[NGENE][12],*pf3x4MG;
   double freqK[NCATG], rK[NCATG], MK[NCATG*NCATG],daa[20*20], *conP,*conP0,*fhK;
   double (*plfun)(double x[],int np);
   double omega_fix;  /* fix the last w in the NSbranchB, NSbranch2 models 
          for lineages.  Useful for testing whether w>1 for some lineages. */
   int conPSiteClass; /* conPSiteClass=0 if (method==0) and =1 if (method==1)?? */
   int nnodeScale;
   char nodeScale[NNODE];    /* nScale[ns-1] for interior nodes */
   double *nodeScaleF;  /* nScaleF[npatt] for scale factors */

  /* communication between SetParameters & ConditionalPNode & EigenQc.  
     Remove both? */
   double *pomega, KAPPA[5+3]; /* why not *pkappa? */
}  com;
struct TREEB {
   int  nbranch, nnode, root, branches[NBRANCH][2];
   double lnL;
}  tree;
struct TREEN {
   int father, nson, sons[NS], ibranch;
   double branch, age, omega, *conP, label;
   char fix_age;
}  *nodes;


extern double Small_Diff;
int Nsensecodon, FROM61[64], FROM64[64], FourFold[4][4];
int ChangedInIteration;  /* 1: t changed, update P(t); 2: paras changed, update UVRoot */
double *PMat,*U,*V,*Root, *_UU[6],*_VV[6],*_Root[6];
/* 5 sets for branchsite models (YN2002); 6 sets for Bielawski's branchsite models */

double pcodon0[64],paa0[20], *pcodonClass;  /* for aaDist=FIT1 or FIT2 */

int LASTROUND, UVRootChanged=1;
int IClass=-1;


int OmegaAA[190], AA1STEP[190];

double _rateSite=1;
char *sitesNonsyn=NULL;  /* used to mark sites with nonsyn differences, for NSsites */
double Qfactor_NS, Qfactor_NS_branch[2], freqK_NS;


double AAchem[][20+1]={  /* last element is the max */
{8.1, 10.5, 11.6, 13, 5.5, 10.5, 12.3, 9, 10.4, 5.2, 
 4.9, 11.3,  5.7, 5.2,  8,  9.2,  8.6, 5.4, 6.2, 5.9,    13}, /* p */
{ 31, 124,  56,  54,   55, 85, 83,   3, 96, 111, 
 111, 119, 105, 132, 32.5, 32, 61, 170, 136, 84,        170}, /* v */
{0, 0.65, 1.33, 1.38, 2.75, 0.89, 0.92, 0.74, 0.58,
 0, 0, 0.33, 0, 0, 0.39, 1.42, 0.71, 0.13, 0.2, 0,      -999},/* c */
{-0.11, 0.079, -0.136, -0.285, -0.184, -0.067, -0.246, -0.073, 0.32, 0.001,
 -0.008, 0.049, -0.041, 0.438, -0.016, -0.153, -0.208, 0.493, 0.381, -0.155} /* a */
};   /* in the order p, v, c, a */


FILE *frub, *flnf, *frst, *frst1, *fin;
char *ratef="rates";
enum {Fequal, F1x4, F3x4, Fcodon, F1x4MG, F3x4MG, FMutSel3x4, FMutSelCodon} CodonFreqs;
char *codonfreqs[]={"Fequal", "F1x4", "F3x4", "Fcodon", "F1x4MG", "F3x4MG", 
                    "FMutSel3x4", "FMutSelCodon"};
enum {NSbranchB=1, NSbranch2, NSbranch3} NSBranchModels;
char *NSbranchmodels[]={"One dN/dS ratio", 
     "free dN/dS Ratios for branches", "several dN/dS ratios for branches",
     "NSbranch3"};
enum {Poisson, EqualInput, Empirical, Empirical_F,
     FromCodon=6, REVaa_0=8, REVaa=9} AAModel;
char *aamodels[]={"Poisson", "EqualInput", "Empirical", "Empirical_F", "",
     "", "FromCodon", "", "REVaa_0", "REVaa"};
enum {NSneutral=1, NSselection, NSdiscrete, NSfreqs, NSgamma, NS2gamma, 
     NSbeta, NSbetaw, NSbetagamma, NSbeta1gamma, NSbeta1normal, NS02normal, 
     NS3normal, NSM8a, NSM8b} NSsitesModels;
char *NSsitesmodels[]={"one-ratio","neutral", "selection","discrete","freqs", 
     "gamma","2gamma","beta","beta&w","beta&gamma", "beta&gamma+1", 
     "beta&normal>1", "0&2normal>0", "3normal>0", "M8a beta&w=1", "M8b beta&w>=1"};
enum {FIT1=11, FIT2=12} SiteClassModels;
enum { AAClasses=7 } aaDistModels;
char *clockstr[]={"", "Global clock", "Local clock", "ClockCombined"};
enum {GlobalClock=1, LocalClock, ClockCombined} ClockModels;


#define CODEML 1
#include "treesub.c"
#include "treespace.c"

FILE *fout;
int ncatG0=10, insmodel=0, nnsmodels=1, nsmodels[14]={0};


int main (int argc, char *argv[])
{
   FILE *fseq=NULL;
   FILE *fpair[6]; 
   char pairfs[6][32]={"2NG.dS","2NG.dN","2NG.t", "2ML.dS","2ML.dN","2ML.t"};
   char ctlf[96]="codeml.ctl", *pmodel, timestr[64];
   char *seqtypestr[3]={"CODONML", "AAML", "CODON2AAML"};
   char *Mgenestr[]={"diff. rate", "separate data", "diff. rate & pi", 
                     "diff. rate & k&w", "diff. rate & pi & k&w"};
   int i, k, sconP0=0,s2=0, idata, nc, nUVR;

#ifdef NSSITESBandits
   atexit(finishup);
#endif
   starttime();
   com.ndata=1;
   noisy=0;           com.runmode=0;
   com.clock=0;       com.fix_rgene=0; /* 0: estimate rate factors for genes */
   com.cleandata=0;  /* 1: delete; 0:use missing data */
   com.seqtype=AAseq;
   com.model=Empirical_F;  strcpy(com.daafile, "jones.dat");
   com.icode=0;       com.nrate=0;
   com.fix_kappa=0;   com.kappa=1;    com.omega=2.1;
   com.fix_alpha=1;   com.alpha=0.;   com.ncatG=4;   /* alpha=0 := inf */
   com.fix_rho=1;     com.rho=0.;
   com.getSE=0;       com.print=0;    com.verbose=1;  com.fix_blength=0;
   com.method=0;      com.space=NULL;
   com.pf3x4MG=com.f3x4MG[0];

   frub=gfopen("rub","w");   frst=gfopen("rst","w");  frst1=gfopen("rst1","w");

/*
mergeSeqs(frst);  exit(0);
Ina();
*/
   SetSeed(4*(int)time(NULL)+1);
   /* SetSeed(123456789); */

#if (DSDN_MC || DSDN_MC_SITES)
   SimulateData2s61();
#endif


/*
{
int i;
double S, dSt, dNt, kappa=1, w=1, P[NCODE*NCODE];
com.seqtype=1; com.ls=1; FOR(i,64) com.pi[i]=1;
FOR(com.icode,12) {
   setmark_61_64 ();
   EigenQc(1,1,&S,&dSt,&dNt,NULL,NULL,NULL, &kappa,w,P);
   printf("%3d  %12.6f\n", com.icode, S);
}
exit(0);
}
*/


   if(argc>1) strcpy(ctlf,argv[1]);
   GetOptions(ctlf);
   fprintf(frst, "Supplemental results for CODEML (seqf: %s  treef: %s)\n", 
      com.seqf, com.treef);

   printf("%s in %s\n",seqtypestr[com.seqtype-1],VerStr);

   fout=gfopen(com.outf, "w");
   if((fseq=fopen(com.seqf,"r"))==NULL || com.seqf[0]=='\0') {
      printf ("\n\nSequence file %s not found!\n", com.seqf);
      exit (-1);
   }

   if(noisy && com.seqtype==CODONseq) 
      { printcu(F0,NULL,com.icode); puts("Nice code, uuh?"); }

   /* space for P&U&V&Root */
   nUVR=1; nc=20;
   if(com.seqtype==1) { nc=64; if(com.model>=1) nUVR=6; }
   else if (com.seqtype==CODONseq||com.model==FromCodon) nc=64;

   PMat=(double*)malloc((nc*nc+nUVR*nc*nc*2+nUVR*nc)*sizeof(double));
   if(PMat==NULL) error2("oom getting P&U&V&Root");
   U=_UU[0]=PMat+nc*nc;  V=_VV[0]=_UU[0]+nc*nc; Root=_Root[0]=_VV[0]+nc*nc;
   for(i=1; i<nUVR; i++) {
      _UU[i]=_UU[i-1]+nc*nc*2+nc; _VV[i]=_VV[i-1]+nc*nc*2+nc; 
      _Root[i]=_Root[i-1]+nc*nc*2+nc;
   }

   /* d4dSdN(fout); */
   if (com.aaDist==AAClasses) {
      SetAA1STEP();
      GetOmegaAA(OmegaAA);
   }
   else if (com.seqtype==AAseq && com.model==REVaa_0)
      SetAA1STEP();

   if(com.seqtype==1) {
      for(i=0; i<3; i++) 
         fpair[i]=(FILE*)gfopen(pairfs[i],"w");
      if(com.runmode==-2)
         for(; i<6;i++) fpair[i]=(FILE*)gfopen(pairfs[i],"w");
   }
   else if(com.runmode==-2)
      fpair[0]=(FILE*)gfopen("2AA.t","w");

   for (idata=0; idata<com.ndata; idata++) {

      if (com.ndata>1) {
         printf ("\n\nData set %d\n", idata+1);
         fprintf (fout, "\n\nData set %d\n", idata+1);
         fprintf(frst,"\t%d",idata+1);
      }
      if(idata)  GetOptions(ctlf); /* Is this necessary? */

      if(nnsmodels>1) {
         if(com.fix_omega) error2("fix omega during batch run?");
         if(com.model) error2("model should be 0 in the batch run?");
         if(com.runmode) error2("runmode?");
         com.NSsites=NSbetaw;  com.ncatG=ncatG0+1;
         printf("NSsites batch run (ncatG as in YNGP2000): ");
         for(i=0; i<nnsmodels; i++) printf(" %2d", nsmodels[i]); FPN(F0);
      }

      /* fprintf(frst1, "[%d]", idata+1); */
      if (com.readpattf) {
         if (com.Mgene) error2("Mgene & readpattf incompatible.");
         ReadPatternFreq (fout, com.seqf);
      }
      else
         ReadSeq((com.verbose?fout:NULL),fseq); /*may change seqtype*/
      if(com.ndata==1) fclose(fseq);

      i=(com.ns*2-1)*sizeof(struct TREEN);
      if((nodes=(struct TREEN*)malloc(i))==NULL) error2("oom nodes");

      if(com.ngene>1 && com.aaDist>=FIT1)  /* because of pcodon0[] */
         { error2("ngene for fitness models"); }

      pmodel=(com.seqtype==CODONseq?NSbranchmodels[com.model]:aamodels[com.model]);
      fprintf(fout,"%s (in %s)  ",seqtypestr[com.seqtype-1],VerStr);
      fprintf(fout,"  %s   Model: %s ",com.seqf,pmodel);
      if(com.clock) fprintf(fout," %s ",clockstr[com.clock]);
      if(com.seqtype==CODONseq||com.model==FromCodon) {
         if(com.fix_kappa) fprintf(fout, " kappa = %.3f fixed\n", com.kappa);
         if(com.fix_omega) fprintf(fout, " omega = %.3f fixed\n", com.omega);
      }
      if(com.seqtype==AAseq && (com.model==Empirical||com.model==Empirical_F))
         fprintf (fout, "(%s) ", com.daafile);
      if(com.seqtype==AAseq&&com.nrate) fprintf(fout,"(nrate:%d) ",com.nrate);
      if(com.alpha && com.rho) fprintf (fout, "Auto-");
      if(com.alpha) fprintf (fout,"dGamma (ncatG=%d) ", com.ncatG);
      if(com.ngene>1)
         fprintf (fout, " (%d genes: %s)  ", com.ngene, Mgenestr[com.Mgene]);

      if(com.alpha==0)  com.nalpha=0;
      else              com.nalpha=(com.nalpha?com.ngene:!com.fix_alpha);
      if(com.Mgene==1) com.nalpha=!com.fix_alpha;
      if(com.nalpha>1 && (!com.alpha || com.ngene==1 || com.fix_alpha))
         error2("Malpha");
      if(com.nalpha>1 && com.rho) error2("Malpha or rho");
      if(com.nalpha>1) fprintf (fout,"(%d gamma)", com.nalpha);
     
      if(com.Mgene && com.ngene==1) error2("Mgene for one gene.");
      if(com.seqtype==CODONseq) {
         fprintf (fout, "\nCodon frequencies: %s\n", codonfreqs[com.codonf]);
         if(com.alpha) 
            fputs("Warning: Gamma model for codons.  See documentation.",fout);
      }
      if((com.seqtype==CODONseq||com.model==FromCodon) 
         && (com.aaDist && com.aaDist<10 && com.aaDist!=AAClasses))
         fprintf(fout,"%s, %s\n",com.daafile,(com.aaDist>0?"geometric":"linear"));

      if(com.NSsites) {
         fprintf(fout,"Site-class models");
         if (nnsmodels==1) {
            fprintf(fout," %s",NSsitesmodels[com.NSsites]);
            if(com.NSsites>=NSdiscrete)fprintf(fout," (%d categories)",com.ncatG);
         }
         if(com.nparK) fprintf(fout," & HMM");
         FPN(fout);
         if(com.aaDist)
            fprintf(fout,"\nFitness models: aaDist: %d\n",com.aaDist);
      }
      com.sspace=max2(100000,3*com.ncode*com.ncode*(int)sizeof(double));
      if(com.NSsites) 
         com.sspace=max2(com.sspace,2*com.ncode*com.ncode+4*com.npatt)*(int)sizeof(double);
/*
      i=com.ns*(com.ns-1)/2;
      com.sspace=max2(com.sspace,
        (int)sizeof(double)*((com.ns*2-2)*(com.ns*2-2+4+i)+i));
*/
      if((com.space=(double*)realloc(com.space,com.sspace))==NULL) {
         printf("\nfailed to get %d bytes for space", com.sspace);
         error2("oom space");
      }
      SeqDistance=(double*)realloc(SeqDistance, i*sizeof(double));
      ancestor=(int*)realloc(ancestor, i*sizeof(int));
      if(SeqDistance==NULL||ancestor==NULL) error2("oom distance&ancestor");

      if(com.seqtype==AAseq) {
         Initialize (fout);
         if (com.model==FromCodon /* ||com.aaDist==AAClasses */)
            AA2Codonf(com.pi, com.fb61);  /* get codon freqs from aa freqs */ 
      }
      else {  /* codon sequences */
         if(com.sspace<max2(com.ngene+1,com.ns)*(64+12+4)*(int)sizeof(double)) {
            com.sspace=max2(com.ngene+1,com.ns)*(64+12+4)*sizeof(double);
            if((com.space=(double*)realloc(com.space,com.sspace))==NULL)
               error2("oom space for #c");
         }
         if (InitializeCodon(fout,com.space)) 
            error2("giving up on stop codons");
         if(com.Mgene==3) FOR (i,com.ngene) xtoy(com.pi,com.piG[i],com.ncode);
      }

      if(com.seqtype==CODONseq)
         DistanceMatNG86(fout,fpair[0],fpair[1],fpair[2],0);
      else
         DistanceMatAA(fout);

      fflush(fout);

      if(com.seqtype==AAseq && com.model==Poisson && !com.print) 
         PatternJC69like (fout);
      if(com.alpha || com.NSsites) {
         s2=com.npatt*com.ncatG*sizeof(double);
         if((com.fhK=(double*)realloc(com.fhK,s2))==NULL) error2("oom fhK");
      }

      if(com.runmode==-2 && com.Mgene!=1) {
         if(com.seqtype==CODONseq) 
            PairwiseCodon(fout,fpair[3],fpair[4],fpair[5],com.space);  
         else
            PairwiseAA(fout, fpair[0]);  
      }
      else {
         if(!com.cleandata) {
            sconP0=com.ns*com.ncode*com.npatt*sizeof(double);
            com.conP0=(double*)malloc(sconP0);
         }
         com.sconP1= 2 *com.ncode*com.npatt*sizeof(double);
         /* com.sconP1= (com.ns-1)*com.ncode*com.npatt*sizeof(double); */
         com.conP=(double*)realloc(com.conP,com.sconP1);

         printf("\n%9ld bytes for distance",com.ns*(com.ns-1)/2*sizeof(double));
         printf("\n%9d bytes for conP0\n%9d bytes for conP1\n",sconP0,com.sconP1);
         printf ("%9d bytes for fhK\n%9d bytes for space\n",s2,com.sspace);
         if((!com.cleandata&&com.conP0==NULL) || com.conP==NULL)
            error2("oom");

    /*
SlidingWindow(fout,com.space);
exit(0);
*/

         if (nnsmodels>1) {
            for(insmodel=0; insmodel<nnsmodels; insmodel++) {
               com.NSsites=nsmodels[insmodel];
               if(com.NSsites<=NSselection)  com.ncatG=com.NSsites+1;
               else if (com.NSsites==NSdiscrete) com.ncatG=3;
               else if (com.NSsites==NSfreqs) com.ncatG=5;
               else if (com.NSsites==NSbetaw||com.NSsites==NS02normal) 
                     com.ncatG=ncatG0+1;
               else  com.ncatG=ncatG0;

               com.nrate=com.nkappa=(com.hkyREV?5:!com.fix_kappa);
               if(com.NSsites==0 || com.NSsites==NSselection || com.NSsites==NSbetaw)
                  com.nrate+=!com.fix_omega;
               if(com.NSsites==NSdiscrete)   com.nrate+=com.ncatG;

               printf("\n\nModel %d: %s\n",com.NSsites, NSsitesmodels[com.NSsites]);
               fprintf(fout,"\n\nModel %d: %s",com.NSsites,NSsitesmodels[com.NSsites]);
               fprintf(frst,"\n\nModel %d: %s",com.NSsites,NSsitesmodels[com.NSsites]);
               fprintf(frub,"\n\nModel %d: %s",com.NSsites,NSsitesmodels[com.NSsites]);
               if(com.NSsites) fprintf(fout," (%d categories)",com.ncatG);
               FPN(fout);

#ifdef NSSITESBandits
            com.fix_blength = (com.NSsites>0 ? 2 : -1);
            if(com.NSsites>0)
               strcpy(com.treef,"M0tree");
#endif
               Forestry(fout);

               printf("\nTime used: %s\n", printtime(timestr));
               fprintf(fout,"\nTime used: %s\n", printtime(timestr));
            }
#ifdef NSSITESBandits
            fprintf(frst1, "\t%s\n", printtime(timestr));
#endif
            exit(0);
         }

         if (com.Mgene==1)        MultipleGenes(fout, fpair, com.space);
         else if (com.runmode==0) Forestry(fout);
         else if (com.runmode==3) StepwiseAddition(fout, com.space);
         else if (com.runmode>=4) Perturbation(fout,(com.runmode==4),com.space);
         else                     StarDecomposition(fout, com.space);
      }
      FPN(frst);  fflush(frst);  
      FPN(frst1); fflush(frst1);
      free(nodes);
   }  /* for (idata) */

   fclose(frst);   if(fin)fclose(fin);
   k=0;
   if(com.seqtype==1) k=(com.runmode==-2?6:3);
   else if (com.runmode==-2) k=1;
   FOR(i,k) fclose(fpair[i]);
   if(com.ndata>1 && fseq) fclose(fseq);  
   free(PMat);
   if(sitesNonsyn) free(sitesNonsyn);

   printf("\nTime used: %s\n", printtime(timestr));
   return (0);
}


/* x[]: t[ntime]; rgene[ngene-1]; kappa; p[](NSsites); omega[]; 
        { alpha(for NSsites) !! alpha, rho || rK[], fK[] || rK[], MK[] }
*/

int Forestry (FILE *fout)
{
   static int times=0;
   FILE *ftree, *frate=NULL;
   int  status=0, i,j=0,k, itree, ntree, np, iteration=1;
   int pauptree=0, haslength;
   double x[NP],xb[NP][2], xcom[NP-NBRANCH], lnL=0,lnL0=0, e=1e-6,tl=0;
   double *var=NULL, nchange=-1;
#ifdef NSSITESBandits
   FILE *fM0tree;
#endif

   if ((ftree=fopen(com.treef,"r"))==NULL) {
      printf("\ntree file %s not found.\n", com.treef);
      exit(-1);
   }
   GetTreeFileType(ftree, &ntree, &pauptree, 0);
   if (com.alpha)
      frate=(FILE*)gfopen(ratef,"w");
   if (ntree>10 && com.print) puts("\nlarge lnf file");
   flnf=gfopen("lnf","w+");
   fprintf(flnf,"%6d %6d %6d\n", ntree, com.ls, com.npatt);

   if(com.seqtype==1 && com.aaDist>=FIT1) {
      xtoy(com.pi,pcodon0,64);
      zero(paa0,20);
      FOR(i,com.ncode) paa0[GeneticCode[com.icode][FROM61[i]]]+=pcodon0[i];
      pcodonClass=(double*)malloc(com.ncatG*64*sizeof(double));
      if(pcodonClass==NULL) error2("oom pcodonClass");
   }
   if(!com.cleandata) InitConditionalPNode();

   for(itree=0; ntree==-1||itree<ntree; iteration=1,itree++) {
      if((pauptree && PaupTreeRubbish(ftree)) || 
         ReadaTreeN(ftree,&haslength, &i,0,1))
            { puts("end of tree file."); break; }
      printf("\nTREE # %2d\n", itree+1);
      fprintf(fout,"\nTREE # %2d:  ", itree+1);
      fprintf(flnf,"\n\n%2d\n", itree+1);
      if(com.print) fprintf (frst,"\n\nTREE # %2d\n", itree+1);
      fprintf(frub,"\n\nTREE #%2d\n", itree+1);

      if (com.fix_blength==2 && !haslength) error2("no branch lengths in tree");
      if (com.fix_blength>0 && !haslength) com.fix_blength=0;
      if (times++==0 && com.fix_blength>0 && haslength) {
         if(com.clock) puts("\nBranch lengths in tree are ignored");
         else {
            if(com.fix_blength==2) puts("\nBranch lengths in tree are fixed.");
            else if(com.fix_blength==1) puts("\nBranch lengths in tree used as initials.");
            if(com.fix_blength==1) {
               FOR(i,tree.nnode) 
                  if((x[nodes[i].ibranch]=nodes[i].branch)<0) 
                     x[nodes[i].ibranch]=1e-5;
            }
         }
      }
      LASTROUND=0;
      if(com.cleandata) nchange=MPScore(com.space);
      if(com.ns<40) { OutaTreeN(F0,0,0); printf("   MP score: %.0f",nchange); }
      OutaTreeN(fout,0,0); fprintf(fout,"   MP score: %.0f",nchange);

      if(!com.clock && nodes[tree.root].nson<=2 && com.ns>2) {
         puts("\nThis is a rooted tree, without clock.  Check.");
         fputs("\nThis is a rooted tree.  Please check!",fout);
      }
      fflush(flnf);

      GetInitials(x, &i);
      if(i==-1) iteration=0;
      if((np=com.np)>NP || np-com.ntime>NP-NBRANCH) error2("raise NP");
      if((i=spaceming2(np))>com.sspace) {
         com.sspace=i;
         printf ("\nspace adjusted to %9d bytes\n",com.sspace);
         if((com.space=(double*)realloc(com.space,com.sspace))==NULL) {
            printf("\ntrying to get %d bytes for ming2", com.sspace);
            error2("oom space");
         }
      }
      printf("\nntime & nrate & np:%6d%6d%6d\n",com.ntime,com.nrate,com.np);
      if(itree && !fin) for(i=0;i<np-com.ntime;i++) x[com.ntime+i]=xcom[i];

      if(noisy>=2&&com.npi0)
         printf("\n%d zero frequencies.\n", com.npi0);

      if(iteration && np) {
         SetxBound(np, xb);
         SetxInitials (np, x, xb); /* start within the feasible region */
      }
      if(com.NSsites==NSM8a || com.NSsites==NSM8b) /* bad for the next tree */
         com.NSsites=NSbetaw;
      PointconPnodes ();

      lnL = com.plfun (x,np);

      if(noisy) {
         printf("\nnp =%6d", np);
         if(noisy>2 && np<100) matout(F0,x,1,np);
         printf("\nlnL0 = %12.6f\n",-lnL);

         /* exit(0); */
/*
         gradient (np, x, lnL, com.space, com.plfun, com.space+np, 1);
         FOR(i,np) printf("%12.6f", com.space[i]);  FPN(F0);
*/
      }

      if(iteration && np) {
         if(com.method == 1)
            j = minB (noisy>2?frub:NULL, &lnL,x,xb, e, com.space);
         else if (com.method==3)
            j = minB2(noisy>2?frub:NULL, &lnL,x,xb, e, com.space);
         else
            j = ming2 (noisy>2?frub:NULL,&lnL,com.plfun,NULL,x,xb, com.space,e,np);

         if (j==-1 || lnL<=0 || lnL>1e7) status=-1;  else status=0;
         if(status) fprintf(fout,"\ncheck convergence..");

      }
      printf("Out..\nlnL  = %12.6f\n",-lnL);

      printf("%d lfun, %d eigenQc, %d P(t)\n",NFunCall, NEigenQ, NPMatUVRoot);
      if (itree==0)
         { lnL0=lnL;  FOR(i,np-com.ntime) xcom[i]=x[com.ntime+i]; }
      else if (!j)
         for (i=0; i<np-com.ntime; i++) xcom[i]=xcom[i]*.2+x[com.ntime+i]*0.8;

      if(!LASTROUND && (com.NSsites==NSselection||com.NSsites==NSdiscrete
        ||com.NSsites==NSfreqs||com.NSsites==NS3normal)) {
         /* transform back to p0, p1,... */
         k=com.ntime+com.nrgene+com.nkappa;

         if(com.nparK) {   /* HMM model for w */
            k+=com.ncatG;
            for(i=0; i<com.ncatG; i++,k+=com.ncatG-1) 
               f_and_x(x+k,x+k,com.ncatG,0,0);
         }
         else {
            j = (com.NSsites==NS3normal ? 3 : com.ncatG);
            if(com.model && com.model<=NSbranch2) j=3;
            f_and_x(x+k,x+k,j,0,0);
         }
      }
      LASTROUND=1;
      if(com.NSsites==NSdiscrete && com.aaDist==0 && com.model==0)
         sortwM3(x);

      if(com.clock) { /* move times into x[] */
         for(i=0,j=1; i<tree.nnode; i++) 
            if(i!=tree.root && nodes[i].nson && !nodes[i].fix_age) 
               x[j++]=nodes[i].age;
      }

      fprintf (fout,"\nlnL(ntime:%3d  np:%3d):%14.6f%+14.6f\n",
         com.ntime, np, -lnL, -lnL+lnL0);
      if(com.fix_blength<2) { OutaTreeB(fout);  FPN(fout); }
      FOR (i,np) fprintf(fout," %8.5f",x[i]); FPN(fout);
      if (com.getSE) {
         if(np>20) puts("Calculating SE's");
         if(com.sspace<np*(np+1)*(int)sizeof(double)) {
            com.sspace=np*(np+1)*sizeof(double);
            if((com.space=(double*)realloc(com.space,com.sspace))==NULL)
               error2("oom space for SE");
         }
         var=com.space+np;

         Hessian (np, x, lnL, com.space, var, com.plfun, var+np*np);
         matinv (var, np, np, var+np*np);
         fprintf(fout,"SEs for parameters:\n");
         FOR(i,np) fprintf(fout," %8.5f",(var[i*np+i]>0.?sqrt(var[i*np+i]):-1));
         FPN(fout);
         if (com.getSE==2) matout2(fout, var, np, np, 15, 10);
      }
      if(com.seqtype==1 && com.ntime && com.clock==0)
         fprintf(fout,"Note: Branch length is defined as number of nucleotide substitutions per codon (not per neucleotide site).\n");
      if(com.NSsites==NSselection && x[com.ntime+com.nrgene+com.nkappa+2]<1)
         fputs("\nNote: This model may have multiple optima. See doc.\n",fout);

      /* if (com.clock) SetBranch (x); */
      if(com.clock && com.nbtype>1)
         fputs("\nWarning: branch rates are not yet applied in tree length and branch lengths",fout);
      if(AbsoluteRate)
         fputs("\nNote: mutation rate is not applied to tree length.  Tree has times, for TreeView",fout);
      for(i=0,tl=0;i<tree.nnode;i++) if(i!=tree.root) tl+=nodes[i].branch;
      fprintf(fout,"\ntree length = %9.5f%s\n",tl,com.ngene>1?" (1st gene)":"");


#ifdef NSSITESBandits
if(com.NSsites==0) {
   fprintf(frst1,"%d\t%d\t%d",com.ns,com.ls,com.npatt);
   for(i=com.ntime; i<com.np; i++) fprintf(frst1,"\t%.3f",x[i]);
   fprintf(frst1,"\t%.2f\t%.3f",tl,-lnL);

   fM0tree=(FILE*)gfopen("M0tree", (insmodel==0?"w":"a"));
   fprintf(fM0tree, "%d  %d\n", com.ns, 1);
   OutaTreeN(fM0tree,1,1);  FPN(fM0tree);
   fclose(fM0tree);
}
else {
   for(i=com.ntime; i<com.np; i++) fprintf(frst1,"\t%.3f",x[i]);
   fprintf(frst1,"\t%.3f",-lnL);
}
#endif


for(i=0; i<com.np; i++) fprintf(frst1,"\t%.3f",x[i]);
fprintf(frst1,"\t%.3f",-lnL);
fflush(frst1);

      if(com.readpattf)  
         { FPN(fout); OutaTreeN(fout,0,1); FPN(fout); }
      else {
         FPN(fout); OutaTreeN(fout,0,1);  FPN(fout);
         FPN(fout); OutaTreeN(fout,1,1);  FPN(fout);
      }
      if(com.np-com.ntime||com.clock) DetailOutput(fout,x,var);
      if (com.seqtype==AAseq && com.model>=2 /* AAClasses ??? */)
         EigenQaa(fout, Root, U, V, x+com.ntime+com.nrgene); /* S & PAM */

      if (com.NSsites && !com.readpattf)  
         lfunNSsites_rate(frst,x,np);
      if (com.print) {
         if(com.rho==0 && com.nparK==0
            && com.clock<=1)
            AncestralSeqs(frst,x);
         if(!com.NSsites && com.plfun!=lfun)
            lfunRates(frate,x,np);
      }
      com.print-=9;  lnL=com.plfun(x,np);  com.print+=9;
      fflush(fout);
   }     /* for (itree) */

   fclose(ftree); 
   if(frate) fclose(frate);
   if (com.aaDist && com.aaDist<10 && com.aaDist!=AAClasses
      && (com.seqtype==CODONseq||com.model==FromCodon))
      printf("\n%s, %s.\n",com.daafile,(com.aaDist>0?"geometric":"linear"));
   if(com.seqtype==1 && com.aaDist>=FIT1) free(pcodonClass);

   if(ntree==-1)  ntree=itree;
   if(ntree>1) { rewind(flnf);  rell(flnf,fout,ntree); }
   fclose(flnf);

   return (0);
}


double *PointKappa (double xcom[], int igene)
{
/* This points to the kappa parameters in xcom[], by looking at com.model, 
   igene, et&c.
*/
   int k=com.nrgene;
   int nka=(com.hkyREV?5:1), nw=(com.aaDist==AAClasses?com.nOmegaType:1);

   if(com.Mgene>1 && com.Mgene>=3)
      k += igene*(nka + nw);

   if(com.fix_kappa) return(&com.kappa);

   return(xcom+k);
}

double *PointOmega (double xcom[], int igene, int inode, int isiteclass)
{
/* This points to the omega parameters in xcom[], by looking at com.model, 
   com.NSsites and igene.  This sometimes points to com.omega or com.rK[].
   This is called by SetParameters(), DetailOutput(), etc.
   
   Difficulties in using this with lfunt() etc.

   Trying to remove global variables com.pomega and com.KAPPA through 
   PointOmega and PointKappa, but was unsuccessful when too many changes were 
   made at the same time.  Perhaps look at this later.  Note that some 
   variables are passed over the head from lfunt() etc. to EigenQc().

   Ziheng Notes: 8 August 2003.

*/
   int k = com.nrgene+com.nkappa, backfore;
   int nka=(com.hkyREV?5:1), nw=(com.aaDist==AAClasses?com.nOmegaType:1);

   if (com.seqtype!=CODONseq && com.model!=FromCodon) 
      error2("should not be here.");

   if(com.NSsites==0 && com.model==0) { /* simple case: one ratio */
      if(com.ngene<=1) {
         if(com.fix_omega) return (&com.omega_fix);  /* fix_omega */
         else              ;
      }
      else if(com.Mgene>=3) 
         k += igene*(nka + nw) + nka;
   }
   else if(com.NSsites==0 && com.model) {  /* branch model */
      if (com.aaDist==0) {
         if(com.fix_omega && nodes[inode].label==com.nbtype-1) 
            return (&com.omega_fix);
         else k += (int)nodes[inode].label;
      }
      else if(com.aaDist==AAClasses)
         k += (int)nodes[inode].label*com.nOmegaType;
   }
   else if (com.NSsites && com.model==0) { /* site model */
      if(com.aaDist<10)
         k += com.ncatG-1+2*isiteclass;
      else if(com.aaDist==FIT1)
         k += com.ncatG-1+4*isiteclass;
      else if(com.aaDist==FIT2)
         k += com.ncatG-1+5*isiteclass;
      else 
         return (&com.rK[isiteclass]);
   }
   else if (com.NSsites && com.model<=NSbranch2) { /* branch&site models A&B */
      k += 2;   /* skip the frequencies. */
      backfore = (int)nodes[inode].label;
      if(isiteclass<2)
         return(&com.rK[isiteclass]);
      else if(isiteclass==2) {
         if(com.fix_omega && backfore) 
            return(&com.omega_fix);
         else
            k += 2 + (com.NSsites==NSselection?0:2) + backfore;
      }
   }
   else { /* NSbranch3: Bielawski's branch&site models C and D */
      k+=com.ncatG-1;   /* skip the frequencies. */
      backfore = (int)nodes[inode].label;
      if(isiteclass<com.ncatG-1)
         return(&com.rK[isiteclass]);
      else if(isiteclass==com.ncatG-1) {
         if(com.fix_omega && backfore) 
            return(&com.omega_fix);
         else
            k += 2 + (com.NSsites==NSselection?0:2) + backfore;
      }
   }
   return (xcom+k);
}


int sortwM3(double x[])
{
/* sort the w values for NSsites=NSdiscrete
   This assumes that com.freqK[] and com.rK[] have been initialized.
*/
   int i, k=com.ntime+com.nrgene+com.nkappa, rank[NCATG];
   double space[NCATG];

   if(com.NSsites!=NSdiscrete) error2("sortwM3");
   if(fabs(1-sum(com.freqK,com.ncatG))>1e-6) error2("sortwM3: freqK");

if(com.nparK) { puts("\asortwM3 for HMM not implemented yet.."); return(-1); }
    
   sort1(com.rK, com.ncatG, rank, 0, (int*)space);
   xtoy(com.rK,space,com.ncatG);
   FOR(i,com.ncatG) com.rK[i]=space[rank[i]];
   xtoy(com.freqK,space,com.ncatG);
   FOR(i,com.ncatG) com.freqK[i]=space[rank[i]];
   FOR(i,com.ncatG-1) x[k+i]=com.freqK[i];
   FOR(i,com.ncatG)   x[k+com.ncatG-1+i]=com.rK[i];
   return(0);
}


static int ijAAref=19*20+9; 
/* reference aa pair: VI (for REVaa, REVaa_0, AAClasses to estimate Grantham)
   The rate for this pair is set to 1, and other rates are relative to it.
*/
#define E1N(m,s) (s/sqrt(PI*2)*exp(-square((1-m)/s)/2)+m*(1-CDFNormal((1-m)/s)))


void DetailOutput (FILE *fout, double x[], double var[])
{
/* var[] is used for codon models if com.getSE=1 to calculate the variances 
   of dS and dN.
*/
   int i,j,k=com.ntime, np=com.np,npclass, ibtype;
   double om=-1,N=-1,S=0,dN=0,dS=0,dSt,dNt, vtw[4],vSN[4], omclass[NCATG];
   double phi1=0,phi2=0, t, *dSdNb=NULL;
   double mu[3]={0,1,2},sig[3]={-1}; /* 3normal: mu0=0 fixed. mu2 estimated */
   double *pkappa=com.KAPPA;

   fprintf(fout,"\nDetailed output identifying parameters\n");

   if(com.clock) OutputTimesRates(fout, x, var);
   k=com.ntime;
   if (com.nrgene) {
      fprintf (fout, "\nrates for %d genes:%6.0f", com.ngene, 1.);
      FOR (i,com.nrgene) fprintf (fout, " %8.5f", x[k++]);
      FPN(fout);
   }

   if (com.seqtype==CODONseq || com.model==FromCodon) {
      if(com.NSsites && com.model==0) {  /* dN/dS by averaged over classes */
         for(j=0,Qfactor_NS=0,dS=dN=0; j<com.ncatG; j++) {
            freqK_NS=com.freqK[j];
            if(com.aaDist) {
               k=com.ntime+com.nrgene+com.nkappa;
               
               if(com.aaDist<10)  com.pomega=x+k+com.ncatG-1+2*j;
               else if(com.aaDist>=FIT1) {
                  com.pomega=x+k+com.ncatG-1+j*(4+(com.aaDist==FIT2));
                  xtoy(pcodonClass+j*64, com.pi, com.ncode);
               }
            }
            EigenQc(1,1,&S,&dSt,&dNt,NULL,NULL,NULL, pkappa,com.rK[j],PMat);
            /* t=1 used here, and dS & dN used later for each branch */
            dS+=freqK_NS*dSt;  dN+=freqK_NS*dNt;
            omclass[j]=dNt/dSt;
         }
         Qfactor_NS=1/Qfactor_NS;
         om=dN/dS;  dS*=Qfactor_NS;  dN*=Qfactor_NS;  N=com.ls*3-S;
      }

      k=com.ntime+com.nrgene;
      if (com.hkyREV) {
         fprintf(fout,"a (TC) & b (TA) & c (TG) & d (CA) & e (CG): ");
         FOR(i,5) fprintf(fout,"%8.5f ", x[k++]);  FPN(fout);
      }
      else if (!com.fix_kappa)
         fprintf(fout,"\nkappa (ts/tv) = %8.5f\n", x[k++]);

      /* dN/dS rate ratios for classes */
      if (com.NSsites>=NSgamma) {
         fprintf(fout,"\nParameters in %s:\n ",NSsitesmodels[com.NSsites]);
         if(com.NSsites==NSgamma) 
            fprintf(fout,"  a=%9.5f  b=%9.5f\n",x[k],x[k+1]);
         else if(com.NSsites==NS2gamma)
            fprintf(fout," p0=%9.5f  a0=%9.5f  b0=%9.5f\n(p1=%9.5f) a1=%9.5f (b1=%9.5f)\n",
            x[k],x[k+1],x[k+2], 1-x[k], x[k+3],x[k+3]);
         else if(com.NSsites==NSbeta)
            fprintf(fout,"p=%9.5f  q=%9.5f\n",x[k],x[k+1]);
         else if(com.NSsites==NSbetaw)
            fprintf(fout," p0=%9.5f  p=%9.5f q=%9.5f\n (p1=%9.5f) w=%9.5f\n",
            x[k],x[k+1],x[k+2], 1-x[k], (com.fix_omega?com.omega:x[k+3]));
         else if(com.NSsites==NSbetagamma)
            fprintf(fout," p0=%9.5f  p=%9.5f  q=%9.5f\n(p1=%9.5f) a=%9.5f  b=%9.5f\n",
            x[k],x[k+1],x[k+2], 1-x[k], x[k+3],x[k+4]);
         else if(com.NSsites==NSbeta1gamma)
            fprintf(fout," p0=%9.5f  p=%9.5f  q=%9.5f\n(p1=%9.5f) a=%9.5f  b=%9.5f\n",
            x[k],x[k+1],x[k+2], 1-x[k], x[k+3],x[k+4]);
         else if(com.NSsites==NSbeta1normal)
            fprintf(fout," p0=%9.5f  p=%9.5f  q=%9.5f\n(p1=%9.5f) u=%9.5f  s=%9.5f\n",
            x[k],x[k+1],x[k+2], 1-x[k], x[k+3],x[k+4]);
         else if(com.NSsites==NS02normal)
            fprintf(fout,"p0=%9.5f  p1=%9.5f  u2=%9.5f  s1=%9.5f  s2=%9.5f\n",
            x[k],x[k+1],x[k+2],x[k+3],x[k+4]);
         else if(com.NSsites==NS3normal)
            fprintf(fout,"p0=%9.5f  p1=%9.5f (p2=%9.5f)\n u2=%9.5f  s1=%9.5f  s2=%9.5f  s3=%9.5f\n",
            x[k],x[k+1], 1-x[k]-x[k+1], x[k+2],x[k+3],x[k+4],x[k+5]);
      }

      if (com.NSsites==NSdiscrete && com.aaDist) { /* structural site classes */
         npclass=(com.aaDist<10 ? 2 : (com.aaDist==FIT1?4:5));
         fprintf(fout,"\nParameters in each class (%d)",npclass);
         fprintf(fout,"%s:\n\n",
            (com.aaDist<10 ? "(b, a)" : "(a_p, p*, a_v, v*, b)"));
         for(j=0,k+=com.ncatG-1; j<com.ncatG; j++,FPN(fout)) {
            fprintf(fout,"%d: f=%8.5f, ",j+1,com.freqK[j]);
            FOR(i,npclass) fprintf(fout,"%9.5f",x[k++]);
            fprintf(fout," dN/dS =%8.5f",omclass[j]);
         }
      }
      else if (com.NSsites && com.model) {
         fprintf(fout,"\n\ndN/dS for site classes (K=%d)\np: ",com.ncatG);
         FOR(i,com.ncatG) fprintf(fout," %8.5f",com.freqK[i]);
         fputs("\nw: ",fout);
         if(com.model<=NSbranch2) {
            FOR(i,com.ncatG-2) fprintf(fout," %8.5f",com.rK[i]);
            fprintf(fout, "  *        *\nw2 = %.5f\n", (com.fix_omega?com.omega_fix:x[com.np-1]));
         }
         else if (com.model==NSbranch3) {
            FOR(i,com.ncatG) fprintf(fout," %8.5f", com.rK[i]);
            k+=com.ncatG-1 + (com.NSsites==2?1:com.ncatG);
            fprintf(fout," %8.5f\n",x[k++]);
         }

      }
      else if (com.NSsites && com.aaDist==0) {
         fprintf(fout,"\n\ndN/dS for site classes (K=%d)\np: ",com.ncatG);
         FOR(i,com.ncatG) fprintf(fout," %8.5f",com.freqK[i]);
         fputs("\nw: ",fout);
         FOR(i,com.ncatG) fprintf(fout," %8.5f",com.rK[i]);  FPN(fout);
         if (com.nparK) {
            fprintf(fout,"\nTransition matrix M in HMM: M_ij=Prob(i->j):\n");
            matout(fout, com.MK, com.ncatG, com.ncatG);
         }
      }
      else if(com.aaDist && com.aaDist<=6) { /* one class (YNH98, Genetics) */
         k=com.ntime+com.nrgene+com.nkappa;
         fprintf (fout,"\nb = %9.5f", x[k++]);
         if (com.seqtype==CODONseq)  fprintf (fout,"\na = %9.5f\n", x[k++]);
      }
      else if(com.aaDist && com.aaDist>=11) { /* fitness, one class */
         fprintf (fout,"\nfitness model (a_p, p*, a_v, v*, (and w0 for FIT2):\n");
         k=com.ntime+com.nrgene+(com.hkyREV?5:!com.fix_kappa);
         FOR(i,4+(com.aaDist==FIT2)) fprintf(fout," %9.5f",x[k++]);  FPN(fout);
      }
      else if(com.model==0 && com.NSsites==0 && !com.fix_omega) 
         k++;
   }
   FOR(j,com.nalpha) {
      if (!com.fix_alpha)  
         fprintf(fout,"\nalpha (gamma) = %8.5f",(com.alpha=x[k++]));
      if(com.nalpha>1) 
         DiscreteGamma(com.freqK,com.rK,com.alpha,com.alpha,com.ncatG,0);
      fprintf (fout, "\nr (%2d):",com.ncatG);
      FOR (i,com.ncatG) fprintf (fout, " %8.5f", com.rK[i]);
      fprintf (fout, "\nf:     ");
      FOR (i,com.ncatG) fprintf (fout, " %8.5f", com.freqK[i]);
      FPN (fout);
   }

   if (com.rho) {
      if (!com.fix_rho) fprintf (fout, "rho (correlation) = %8.5f\n", x[k]);
      fprintf (fout, "transition probabilities between rate categories:\n");
      for(i=0;i<com.ncatG;i++,FPN(fout))  FOR(j,com.ncatG) 
         fprintf(fout," %8.5f",com.MK[i*com.ncatG+j]);
   }

#if 0
   /* phi1 & phi2 for NSsites models */
   /* NSgamma (2):       alpha, beta
      NS2gamma (4):      p0, alpha1, beta1, alpha2 (=beta2)
      NSbeta (2):        p_beta, q_beta
      NSbetaw (4):       p0, p_beta, q_beta, w (if !com.fix_omega, not used here)
      NSbetagamma (5):   p0, p_beta, q_beta, alpha, beta
      NSbeta1gamma (5):  p0, p_beta, q_beta, alpha, beta (1+gamma)
      NSbeta1normal (5): p0, p_beta, q_beta, mu, s (normal>1)
      NS02normal (5):    p0, p1, mu2, s1, s2 (s are sigma's)
      NS3normal (6):     p0, p1, mu2, s0, s1, s2 (s are sigma's)
   */
   if(com.seqtype==CODONseq && com.NSsites>=NSgamma) {
      puts("\nCalculting phi1 & phi2 from NSsites model.\n");
      k=com.ntime+com.nrgene+(com.hkyREV?5:!com.fix_kappa);
      if(com.NSsites!=NSbetaw && com.NSsites!=NSbeta&&com.NSsites!=NS02normal)
         phi1=1-CDFdN_dS(1,x+k);

      switch(com.NSsites) {
      case(NSgamma):  
         phi2=x[k]/x[k+1]*(1-CDFGamma(1,x[k]+1,x[k+1]));   break;
      case(NS2gamma): 
         phi2=x[k]*x[k+1]/x[k+2]*(1-CDFGamma(1,x[k+1]+1,x[k+2]))
             + (1-x[k])*(1-CDFGamma(1,x[k+3]+1,x[k+3]));    
         break;
      case(NSbetaw): 
         if(x[k+3]>1) { phi1=1-x[k]; phi2=(1-x[k])*x[k+3]; }
         break;
      case(NSbetagamma):
         phi2=(1-x[k])*x[k+3]/x[k+4]*(1-CDFGamma(1,x[k+3]+1,x[k+4]));
         break;
      case(NSbeta1gamma):  phi2=(1-x[k])*(1+x[k+3]/x[k+4]);  break;
      case(NSbeta1normal):
         phi2=(1-x[k])*E1N(x[k+3],x[k+4])/(1-CDFNormal((1-x[k+3])/x[k+4]));
         break;
      case(NS02normal):
         phi1=(1-x[k])*(1-CDFdN_dS(1,x+k));
         mu[2]=x[k+2]; sig[1]=x[k+3]; sig[2]=x[k+4];
         phi2=(1-x[k])*x[k+1]* (sig[1]/sqrt(PI*2)+.5)/CDFNormal(1/sig[1])
             +(1-x[k])*(1-x[k+1])*E1N(mu[2],sig[2])/CDFNormal(mu[2]/sig[2]);
         break;
      case(NS3normal):
         mu[2]=x[k+2]; sig[0]=x[k+3]; sig[1]=x[k+4]; sig[2]=x[k+5];

         phi2=x[k]*2*(sig[0]/sqrt(PI*2)*exp(-.5/square(sig[0]))
             +x[k+1]* (sig[1]/sqrt(PI*2)+.5)/CDFNormal(1/sig[1])
             +(1-x[k]-x[k+1])* E1N(mu[2],sig[2]))/CDFNormal(mu[2]/sig[2]);
         break;
      }

      printf("phi1=P(w>1) = %9.4f\tphi2=E(w|w>1) = %7.4f\n",phi1,phi2);
      fprintf(fout,"\nphi1=P(w>1) = %9.4f\nphi2=E(w|w>1)= %7.4f\n",phi1,phi2);
      fprintf(frst1," phi1&2:%9.4f%9.4f",phi1,phi2);
   }
#endif

   if (com.aaDist==AAClasses) {
      fprintf (fout, "\nw (dN/dS) classes for amino acid pairs:\n");
      FOR (k, com.nOmegaType) {
         for(ibtype=0; ibtype<com.nbtype; ibtype++)
            fprintf (fout, " %9.5f", 
               x[com.ntime+com.nrgene+com.nkappa+k+ibtype*com.nOmegaType]);
         fprintf (fout, ": ");
         FOR (i,20) FOR(j,i)
            if (OmegaAA[i*(i-1)/2+j]==k) fprintf(fout," %c%c", AAs[i],AAs[j]);
         if (k==0)  fprintf(fout, " (background ratio)");
         FPN(fout); 
      }
   }

   /* dN and dS for each branch in the tree */
   if(com.seqtype==CODONseq && com.ngene==1 && (com.model==0 || com.NSsites==0)
      /*||com.model==FromCodon||com.aaDist==AAClasses */){
      if(com.model) {
         dSdNb=(double*)malloc(tree.nnode*2*sizeof(double));
         if(dSdNb==NULL) error2("oom DetailOutput");
      }
      fputs("\ndN & dS for each branch\n\n",fout);
      fprintf(fout,"%7s%12s%9s%9s%9s%9s%9s %6s %6s\n\n",
              "branch","t","S","N","dN/dS","dN","dS","S*dS","N*dN");
      FOR (i,tree.nbranch) {
         fprintf(fout,"%4d..%-3d ",tree.branches[i][0]+1,tree.branches[i][1]+1);
         k=com.ntime+com.nrgene+com.nkappa;
         t=nodes[tree.branches[i][1]].branch;

         if(!com.NSsites) {
            if (com.aaDist) om=-1; /* not used in EigenQc() */
            else if (com.model==0 || com.model==FromCodon)
               om=(com.fix_omega?com.omega:x[k]);
            else if (com.model==NSbranchB) om=x[k+i];
            else if (com.model==NSbranch2) om=nodes[tree.branches[i][1]].omega;

            if(com.model && com.aaDist)
               com.pomega=x+com.ntime+com.nrgene+!com.fix_kappa
                         +(int)nodes[tree.branches[i][1]].label*com.nOmegaType;

            EigenQc(1,t,&S,&dS,&dN, NULL,NULL,NULL, pkappa,om,PMat); /* */
            if (com.aaDist) om=dN/dS;
/*
            if(dS<.01/com.ls) om=-1;
            else if(om==-1)   om=dN/dS;
            if(com.model==0)  om=com.omega;
*/
            N=com.ls*3-S;
            if(com.model) {  dSdNb[i]=dS; dSdNb[tree.nnode+i]=dN; }

            /* om not used in AAClasses model */
            if(com.getSE>1&&com.fix_blength<2&&!com.clock&&com.aaDist!=AAClasses){
               vtw[0]=var[i*np+i];  vtw[3]=var[k*np+k]; 
               vtw[1]=vtw[2]=var[i*np+k]; 
               VariancedSdN(t, com.omega, vtw, vSN);
               fprintf(fout,"dN = %7.5f +- %.5f   dS = %7.5f +- %.5f",
                  dN,(vSN[3]>0?sqrt(vSN[3]):-0),dS,(vSN[0]>0?sqrt(vSN[0]):-0));
               fprintf(fout," (by method 2)\n");
            }
            fprintf(fout,"%9.3f%9.1f%9.1f%9.4f%9.4f%9.4f %6.1f %6.1f\n",
                          t,S,N,om,dN,dS,S*dS,N*dN);
            /* fprintf(frst,"%8.1f%8.1f %9.5f%9.4f%9.4f",N,S,om,dN,dS); */
         }
         else if(com.model==0) {  /* NSsites & other site-class models */
            fprintf(fout,"%9.3f%9.1f%9.1f%9.4f%9.4f%9.4f %6.1f %6.1f\n",
                          t,S,N,om,dN*t,dS*t, S*dS*t,N*dN*t);
            /* fprintf(frst,"%8.1f%8.1f %9.5f%9.4f%9.4f",N,S,om,dN*t,dS*t); */
         }
         else {  /* NSbranchsites models */
            ;
            /*
            Qfactor_NS=Qfactor_NS_branch[(int)nodes[i].label];
            EigenQc(1,t,&S,&dS,&dN,NULL,NULL,NULL, pkappa,om,PMat);

            fprintf(fout,"%9.3f%9.1f%9.1f%9.4f%9.4f%9.4f %6.1f %6.1f\n",
                          t,S,N,om,dN*t,dS*t, S*dS*t,N*dN*t);
            */
         }
      }  /* for (i) */
      if(com.model && com.model!=7 && !com.NSsites) {
         FOR(i,tree.nbranch) nodes[tree.branches[i][1]].branch=dSdNb[i]; 
         fprintf(fout,"\ndS tree:\n");  OutaTreeN(fout,1,1);
         FOR(i,tree.nbranch) nodes[tree.branches[i][1]].branch=dSdNb[tree.nnode+i];
         fprintf(fout,"\ndN tree:\n");  OutaTreeN(fout,1,1);  FPN(fout);
         free(dSdNb);
      }
   }  /* if codonseqs */

   FPN(fout);
}



void GetNSsitesModels(char *line)
{
/* This reads the line  NSsites = 0 1 2 3 7 8  in codeml.ctl.
*/
   char *pline;
   int pop_digit, nNSsitesModels=16;

   if ((pline=strstr(line, "="))==NULL) error2(".ctl file error NSsites");
   pline++;
   for (nnsmodels=0; nnsmodels<nNSsitesModels; nnsmodels++) {
      if(sscanf(pline, "%d", &nsmodels[nnsmodels]) != 1) break;
      for(pop_digit=0; ; ) {
         if(isdigit(*pline)) { pline++; pop_digit=1; }
         else if(isspace(*pline)) {
            pline++;
            if(pop_digit) break;
         }
         else  error2(".ctl file NSsites line strange.");
      }
      if(nsmodels[nnsmodels]<0 || nsmodels[nnsmodels]>=nNSsitesModels)
         error2("NSsites model");
   }
   com.NSsites=nsmodels[0];
}


int GetOptions (char *ctlf)
{
   int i,j, nopt=35, lline=255;
   char line[255], *pline, opt[99], *comment="*#";
   char *optstr[] = {"seqfile", "outfile", "treefile", "seqtype", "noisy", 
        "cleandata", "runmode", "method", 
        "clock", "getSE", "RateAncestor", "CodonFreq", "verbose",
        "model", "hkyREV", "aaDist","aaRatefile",
        "NSsites", "NShmm", "icode", "Mgene", "fix_kappa", "kappa",
        "fix_omega", "omega", "fix_alpha", "alpha","Malpha", "ncatG", 
        "fix_rho", "rho", "ndata", "Small_Diff", "fix_blength", "readpattf"};
   double t;
   FILE  *fctl;
   char *daafiles[]={"", "grantham.dat", "miyata.dat", 
                     "g1974c.dat","g1974p.dat","g1974v.dat","g1974a.dat"};

   fctl=gfopen(ctlf,"r");
   if (noisy) printf ("\n\nReading options from %s..\n", ctlf);
   for (;;) {
      if (fgets (line, lline, fctl) == NULL) break;
      for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
         if (isalnum(line[i]))  { t=1; break; }
         else if (strchr(comment,line[i])) break;
      if (t==0) continue;
      sscanf (line, "%s%*s%lf", opt,&t);
      if ((pline=strstr(line, "="))==NULL) 
         error2("err: option file. add space around the equal sign?");
      for (i=0; i<nopt; i++) {
         if (strncmp(opt, optstr[i], 8)==0)  {
            if (noisy>=9)
               printf ("\n%3d %15s | %-20s %6.2f", i+1,optstr[i],opt,t);
            switch (i) {
               case ( 0): sscanf(pline+1, "%s", com.seqf);    break;
               case ( 1): sscanf(pline+1, "%s", com.outf);    break;
               case ( 2): sscanf(pline+1, "%s", com.treef);   break;
               case ( 3): com.seqtype=(int)t;     break;
               case ( 4): noisy=(int)t;           break;
               case ( 5): com.cleandata=(char)t;  break;
               case ( 6): com.runmode=(int)t;     break;
               case ( 7): com.method=(int)t;      break;
               case ( 8): com.clock=(int)t;       break;
               case ( 9): com.getSE=(int)t;       break;
               case (10): com.print=(int)t;       break;
               case (11): com.codonf=(int)t;      break;
               case (12): com.verbose=(int)t;     break;
               case (13): com.model=(int)t;       break;
               case (14): com.hkyREV=(int)t;      break;
               case (15): com.aaDist=(int)t;      break;
               case (16): sscanf(pline+2,"%s",com.daafile); break;
               case (17): GetNSsitesModels(line); break;
               case (18): com.nparK=(int)t;       break;
               case (19): com.icode=(int)t;       break;
               case (20): com.Mgene=(int)t;       break;
               case (21): com.fix_kappa=(int)t;   break;
               case (22): com.kappa=t;            break;
               case (23): com.fix_omega=(int)t;   break;
               case (24): com.omega=t;            break;
               case (25): com.fix_alpha=(int)t;   break;
               case (26): com.alpha=t;            break;
               case (27): com.nalpha=(int)t;      break;
               case (28): com.ncatG=(int)t;       break;
               case (29): com.fix_rho=(int)t;     break;
               case (30): com.rho=t;              break;
               case (31): com.ndata=(int)t;       break;
               case (32): Small_Diff=t;           break;
               case (33): com.fix_blength=(int)t; break;
               case (34): com.readpattf=(int)t;   break;
           }
           break;
         }
      }
      if (i==nopt)
        { printf ("\noption %s in %s not recognised\n", opt,ctlf); exit(-1); }
   }
   fclose (fctl);

   if (noisy) FPN(F0);
   setmark_61_64 ();
   if (com.ngene==1) com.Mgene=0;
   if (com.seqtype==AAseq || com.seqtype==CODON2AAseq) {
      if(com.NSsites) error2("use NSsites=0 for amino acids?");
      if(com.hkyREV)  error2("use hkyREV=0 for amino acids?");
      com.ncode=20;
      if(com.aaDist==AAClasses) 
         com.nrate=com.nkappa=(com.hkyREV ? 5 : !com.fix_kappa); 

      switch (com.model) {
      case (Poisson):  case (EqualInput): case (Empirical): case (Empirical_F):
         com.fix_kappa=1; com.kappa=0; com.nrate=0;   break;
      case (FromCodon): 
         com.nrate=com.nkappa=(com.hkyREV ? 5 : !com.fix_kappa);
         if(com.aaDist) com.nrate++;
         if(com.fix_omega) error2("fix_omega = 1");
         break;
      case (REVaa_0): com.fix_kappa=0; com.kappa=0; break; 
      case (REVaa):   com.fix_kappa=0; com.kappa=0; com.nrate=189; break;
      default: error2("model unavailable");
      }
      if (com.Mgene>2 || (com.Mgene==2 && (com.model==Fequal||com.model==2))) 
         error2 ("Mgene && model");
   }
   else if(com.seqtype==CODONseq) {
      if(com.nparK)
         if (com.model||com.aaDist||com.NSsites!=NSdiscrete||com.alpha||com.rho)
            error2("HMM model option");
      if(com.Mgene>1 && com.model) error2("Mgene & model?");
      if(com.fix_kappa && com.hkyREV)
         error2("can't fix kappa for the codon model you selected.");
      if(com.fix_kappa && com.codonf>=FMutSel3x4 && com.Mgene>=2)
         error2("option fix_kappa for Mgene models not implemented");
      if(com.hkyREV && (com.aaDist || com.Mgene>1))
         error2("hkyREV with aaDist or Mgene: check options?\a");
      if(com.aaDist && com.NSsites) 
         error2("aaDist & NSsites don't work together");
      if((com.model && com.aaDist)
         && (com.model>NSbranch2 || com.aaDist!=AAClasses))
            error2("model & aaDist");
      if(com.aaDist && com.fix_omega) 
         error2("can't fix_omega for aaDist models");
      
      com.nrate=com.nkappa=(com.hkyREV ? 5 : !com.fix_kappa);
      if(com.codonf>=FMutSel3x4)  com.nrate=com.nkappa += 3;
      if(com.NSsites==NSM8b)   com.fix_omega=0;
      if(com.NSsites==NSM8a) { com.fix_omega=1;  com.omega=1; }
      if (com.aaDist!=AAClasses) {
         if(com.fix_kappa>1) error2("fix_kappa>1, not tested.");  /** ???? */
         if (com.model>0 && (com.alpha || !com.fix_alpha)) 
            error2("dN/dS ratios among branches not implemented for gamma");
         if (com.model>0 && com.clock) 
            error2("model and clock don't work together");
         if (com.fix_omega) {
            com.omega_fix=com.omega;
            if((com.model==0 && com.NSsites==NSdiscrete)
               || (com.model && com.NSsites && com.NSsites!=NSselection
                   &&com.NSsites!=NSdiscrete&&com.NSsites!=NSbetaw))
               error2("\afix_omega?");
         }
         if (com.model>NSbranch3) error2("seqtype or model.");
/*
         if (com.model==NSbranch2 && com.clock==2) 
            error2("NSbranch & local clock.");
*/
         if (com.model==NSbranch2 && com.nbtype>1) 
            error2("com.nbtype & local clock models don't work together?");
         if (com.model==NSbranch3 && com.NSsites==2 && com.ncatG!=3) 
            { com.ncatG=3; puts("ncatG=3 reset."); }
         if(com.kappa<0)  error2("kappa..");
         if(com.runmode==-2 && (com.NSsites||com.alpha||com.aaDist))
            error2("err: incorrect model for pairwise comparison.\ncheck NSsites, alpha, aaDist.");
         if(com.runmode>0 && com.model==2) error2("tree search & model");
         if(com.aaDist && com.NSsites!=0 && com.NSsites!=NSdiscrete)
            error2("NSsites && aaDist.");

         if(com.aaDist==0) {
            if(!com.fix_omega || (com.Mgene && com.Mgene>=3)) com.nrate++;
         }
         else {
            if(com.aaDist<=6)          com.nrate+=2;   /* a & b, PSB2000 */
            else if(com.aaDist==FIT1)  com.nrate+=4; /* fitness models: */
            else if(com.aaDist==FIT2)  com.nrate+=5; /* ap, p*, av, v*, b */
            if(com.aaDist>=FIT1) FOR(i,2) FOR(j,20) AAchem[i][j]/=AAchem[i][20];
         }
         if (com.Mgene>=3 && com.nrate==0)  error2("Mgene");

         if(com.NSsites) {
            com.nrate=(com.hkyREV ? 5 : !com.fix_kappa);
            if(com.model && com.ngene>1) error2("NSbranchsites with ngene.");
            if(com.NSsites==NSfreqs && com.ncatG!=5)
               { puts("\nncatG changed to 5."); com.ncatG=5; }
            if(com.model && com.NSsites)
               if((com.model!=2  && com.model!=3) 
                  || (com.NSsites!=NSselection&&com.NSsites!=NSdiscrete))
               error2("only NSsites=2,3 & model=2,3 are compatible.");
            if(com.alpha || com.fix_alpha==0) error2("NSsites & gamma");
            if(com.NSsites<=NSselection) com.ncatG=com.NSsites+1;
            if(com.NSsites==NSbetaw || com.NSsites==NS02normal
               || com.NSsites==NSM8a || com.NSsites==NSM8b) com.ncatG++;

            if(com.model==2) { /* branchsite models A & B */
               com.ncatG=4; 
               com.nrate += (com.NSsites==2?1:3);
            }
            else if(com.model==3) { /* Joe's models C & D */
               if(com.NSsites==NSselection) {
                  com.ncatG=3;  com.nrate+=2;
               }
               if(com.NSsites==NSdiscrete) {
                  if(com.ncatG!=2  && com.ncatG!=3) 
                     error2("use 2 or 3 for ncatG for model=3?");
                  com.nrate += com.ncatG+1;
               }
            }
            else if(com.NSsites==NSselection || com.NSsites==NSbetaw || com.NSsites==NSM8b)
               { if(!com.fix_omega) com.nrate++; else com.omega_fix=com.omega; }
            else if(com.NSsites==NSdiscrete && com.aaDist) {
               if (com.aaDist<=6) com.nrate+=com.ncatG;  /* a&b PSB2000 */
               else {  /* fitness models */
                  com.nrate=!com.fix_kappa+4*com.ncatG;
                  if(com.aaDist==FIT2) com.nrate+=com.ncatG;
               }
            }
            else if(com.NSsites==NSdiscrete) 
               com.nrate+=com.ncatG;    /* omega's */
         }
      }
   }
   else
      error2 ("seqtype..");

   if(com.runmode==-2 && com.cleandata==0) {
      com.cleandata=1; puts("gaps are removed for pairwise comparison.");
   }
   if(com.method &&(com.clock||com.rho)) 
      { com.method=0; puts("\aiteration method reset: method = 0"); }

   if (com.runmode==3 && (com.clock)) error2("runmode+clock");
   if (com.aaDist<=6 && (com.seqtype==CODONseq || com.model==FromCodon))
      strcpy(com.daafile, daafiles[abs(com.aaDist)]);
   else if (com.seqtype==AAseq && com.model>=REVaa_0)
      strcpy(com.daafile,"mtmam.dat");

   if (com.fix_alpha && com.alpha==0) {
      if (com.rho) puts("rho set to 0.");  com.fix_rho=1; com.rho=0; 
   }
   if(!com.fix_alpha && com.alpha<=0) { com.alpha=0.5; puts("alpha reset"); }
   if(!com.fix_rho && com.rho==0) { com.rho=0.001;  puts("init rho reset"); }
   if(com.alpha||com.NSsites) 
      { if(com.ncatG<2 || com.ncatG>NCATG) error2("ncatG"); }
   else if (com.ncatG>1) com.ncatG=1;

#ifdef NSSITES_K1_K2_CLASSES
   if((com.NSsites==NSgamma||com.NSsites==NS2gamma||com.NSsites>=NSbetagamma)
      && com.ncatG<10) error2("need more categories for NSsites");
#endif

   fin=fopen("in.codeml","r");

   return(0);
}


int testx (double x[], int np)
{
/* This is used for LS branch length estimation by nls2, called only if(clock==0)
*/
   int i;
   double tb[]={.4e-6, 99};

   FOR (i,com.ntime)  
      if (x[i]<tb[0] || x[i]>tb[1]) 
         return (-1);
   return (0);
}



int SetxBound (int np, double xb[][2])
{
   int i=-1,j,k, K=com.ncatG;
   double tb[]={4e-6,50}, rgeneb[]={0.01,99}, rateb[]={1e-4,99};
   double alphab[]={0.005,99}, betab[]={0.005,99}, omegab[]={0.0001,999};
   double rhob[]={0.01,0.99}, pb[]={.000001,.999999};

   if(!com.clock)
      FOR (i,com.ntime)  FOR (j,2) xb[i][j]=tb[j];
   else {
      /* root relative age */
      xb[0][0]=(AbsoluteRate?AgeLow[tree.root]:tb[0]); xb[0][1]=tb[1];
      for(k=1; k<tree.nnode-com.ns-NFossils; k++)  /* prop times */
         { xb[k][0]=pb[0]; xb[k][1]=pb[1]; }
      for(; k<com.ntime; k++)                      /* rate and branch rates */
         FOR(j,2) xb[k][j]=rateb[j];
   }

   /* for(i=com.ntime;i<np;i++) { xb[i][0]=rateb[0]; xb[i][1]=rateb[1]; } */
   FOR (i,com.nrgene) FOR (j,2) xb[com.ntime+i][j]=rgeneb[j]; 
   FOR (i,com.nrate)  FOR (j,2) xb[com.ntime+com.nrgene+i][j]=rateb[j];
   k=com.ntime+com.nrgene+!com.fix_kappa; 
   if (com.NSsites) { /* p's before w's in xb[] */
      omegab[0]*=0.01;
      if(com.NSsites==NSM8b) omegab[0]=1;

      /* if(com.NSsites==NSselection) omegab[0]=1; */

      switch(com.NSsites) {
      case(NSneutral):  xb[k][0]=pb[0]; xb[k++][1]=pb[1];  break;  /* p0 */
      case(NSselection): /* for p0, p1, w2 */
         FOR(j,2) { xb[k][0]=-99; xb[k++][1]=99; }
         if(!com.fix_omega) /* for model=0 or 2 or 3 */
            { xb[k][0]=omegab[0]; xb[k++][1]=omegab[1]; } 
         if(com.model==3)  /* Joe's model C */
            { xb[k][0]=omegab[0];  xb[k++][1]=omegab[1]; }
         break;  
      case(NSdiscrete):  /* pK[] & rK[] */
         if(com.model==3) { /* Joe's model D */
            if(com.nparK) error2("model & NSsites & nparK");
            FOR(j,K-1) { xb[k][0]=-99; xb[k++][1]=99; }
            FOR(j,K+1) { xb[k][0]=omegab[0];  xb[k++][1]=omegab[1]; }
         }
         else if(com.model==2) {  /* branch-site model B */
            K=3;
            if(com.nparK==0) 
               FOR(j,K-1) { xb[k][0]=-99; xb[k++][1]=99; }
            FOR(j,K) { xb[k][0]=omegab[0];  xb[k++][1]=omegab[1]; }
            if(com.nparK) 
               FOR(j,K*(K-1)) { xb[k][0]=-99; xb[k++][1]=99; }
         }
         else  {                 /* NSsites models M3 */
            FOR(j,K-1) { xb[k][0]=-99; xb[k++][1]=99; }
            FOR(j,K) { xb[k][0]=omegab[0];  xb[k++][1]=omegab[1]; }
         }

         if(com.seqtype==CODONseq && com.aaDist)
            FOR(j,K) { xb[k][0]=omegab[0];  xb[k++][1]=omegab[1]; }
         break; 
      case(NSfreqs):     /* p0...pK */
         FOR(j,K-1) { xb[k][0]=-99; xb[k++][1]=99; }
         break; 
      case(NSgamma):
         FOR(j,2) { xb[k][0]=alphab[0]; xb[k++][1]=alphab[1]; } break;
      case(NS2gamma):    /* p0, alpha1,beta1,alpha2=beta2 */
         xb[k][0]=pb[0]; xb[k++][1]=pb[1];
         FOR(j,3) { xb[k][0]=alphab[0]; xb[k++][1]=alphab[1]; }
         break;
      case(NSbeta):       /* p_beta,q_beta */
         FOR(j,2) { xb[k][0]=betab[0]; xb[k++][1]=betab[1]; } 
         break;
      case(NSbetaw):
      case(NSM8a):
      case(NSM8b):
         /* p0, p_beta, q_beta, w */
         xb[k][0]=pb[0]; xb[k++][1]=pb[1]; /* p0 */
         FOR(j,2) { xb[k][0]=betab[0]; xb[k++][1]=betab[1]; }  /* p & q */
         if(!com.fix_omega) { xb[k][0]=omegab[0];  xb[k++][1]=omegab[1]; }
         break;
      case(NSbetagamma):  /* p0, p_beta, q_beta, alpha, beta */
         xb[k][0]=pb[0]; xb[k++][1]=pb[1]; /* p0 */
         FOR(j,4) { xb[k][0]=betab[0]; xb[k++][1]=betab[1]; }  /* p&q, a&b */
         break;
      case(NSbeta1gamma):  /* p0, p_beta, q_beta, alpha, beta */
         xb[k][0]=pb[0]; xb[k++][1]=pb[1]; /* p0 */
         FOR(j,4) { xb[k][0]=betab[0]; xb[k++][1]=betab[1]; }  /* p&q, a&b */
         break;
      case(NSbeta1normal):  /* p0, p_beta, q_beta, mu, s */
         xb[k][0]=pb[0]; xb[k++][1]=pb[1]; /* p0 */
         FOR(j,4) { xb[k][0]=betab[0]; xb[k++][1]=betab[1]; }  /* p&q, mu&s */
         xb[k-2][0]=1;  xb[k-2][1]=9;  /* mu */
         break;
      case(NS02normal):   /* p0, p1, mu2, s1, s2 */
         FOR(j,2) { xb[k][0]=pb[0];  xb[k++][1]=pb[1]; }  /* p0 & p1, */
         FOR(j,3) { xb[k][0]=.0001; xb[k++][1]=29; }  /* mu2,s1,s2 */
         break;
      case(NS3normal):    /* p0, p1, mu2, s0, s1, s2 */
         FOR(j,2) { xb[k][0]=-99;  xb[k++][1]=99; }  /* p0 & p1, tranformed */
         FOR(j,4) { xb[k][0]=.0001; xb[k++][1]=29; }  /* mu2,s0,s1,s2 */
         break;
      }
   }
   else if((com.seqtype==CODONseq||com.model==FromCodon)&&com.aaDist!=AAClasses)
     { if(!com.fix_omega) { xb[k][0]=omegab[0]; xb[k][1]=omegab[1]; } }
   if(com.seqtype==CODONseq&&com.model)
      FOR(j,com.nOmega-com.fix_omega) 
         { xb[k+j][0]=omegab[0]; xb[k+j][1]=omegab[1]; }

   if (com.aaDist<0 && (com.seqtype==1||com.model==FromCodon)) {
      /* linear relationship between d_ij and w_ij */
      if(com.nrate != !com.fix_kappa+1+(com.seqtype==1)) error2("in Setxbound");
      xb[com.ntime+com.nrgene+!com.fix_kappa][1]=1; /* 0<b<1 */
   }

   k=com.ntime+com.nrgene+com.nrate;
   for (i=0;i<com.nalpha;i++,k++)  FOR (j,2) xb[k][j]=alphab[j];
   if (!com.fix_rho)   FOR (j,2) xb[np-1][j]=rhob[j];

   if(noisy>=3 && np<100) {
      printf("\nBounds (np=%d):\n",np);
      FOR(i,np) printf(" %10.6f", xb[i][0]);  FPN(F0);
      FOR(i,np) printf(" %10.6f", xb[i][1]);  FPN(F0);
   }

   return(0);
}

void getpcodonClass(double x[], double pcodonClass[])
{
/* This uses pcodon0[], paa0[], and x[] to calculate pcodonclass[] and
   com.pi[] for the fitness models.
   pcodon0[] has the codon frequencies observed (codonFreq=3) or expected 
   (codonFreq=2 or 1 or 0) rootally.  Under the fitness models, the expected 
   codon frequencies pcodonClass[] differs among site classes and from the 
   rootal pi[] (pcodon0[]).
   This is called by SetParameters().
*/
   int i,iclass,iaa, k, nclass=(com.NSsites==0?1:com.ncatG);
   double paaClass[20], *w,fit;

   if(com.seqtype!=1 || com.aaDist<FIT1) error2("getpcodonClass");
   k=com.ntime+com.nrgene+!com.fix_kappa+nclass-1;
   FOR(iclass, nclass) {
      w=x+k+iclass*(4+(com.aaDist==FIT2));
      FOR(iaa,20) {
         fit=-w[0]*square(AAchem[0][iaa]-w[1])
              -w[2]*square(AAchem[1][iaa]-w[3]);
         paaClass[iaa]=exp(2*fit);
      }
      abyx(1/sum(paaClass,20), paaClass, 20);
      FOR(i,com.ncode) {
         iaa=GeneticCode[com.icode][FROM61[i]];
         pcodonClass[iclass*64+i]=pcodon0[i]/paa0[iaa]*paaClass[iaa];
      }

if(fabs(1-sum(pcodonClass+iclass*64,com.ncode))>1e-5) error2("pcodon!=1");
/*
fprintf(frst,"\nSite class %d: ",iclass+1);
matout (frst,paaClass,2, 10);
matout (frst,pcodonClass+iclass*64,16,4);
*/
   }
   if(nclass==1) FOR(i,com.ncode) com.pi[i]=pcodonClass[i];
}


int GetInitials (double x[], int* fromfile)
{
/* This caculates the number of parameters (com.np) and get initial values.
   This routine is too messy.  Perhaps try to restruct the code and make 
   two sections for amino acids and codons?
   com.nrate is initialised in getoptions().
*/
   static int times=0;
   int i, j,k=0, naa=20;
   int sconP1_new=(tree.nnode-com.ns)*com.ncode*com.npatt*sizeof(double);
   double t;
   
   NFunCall=NPMatUVRoot=NEigenQ=0;
   if(com.clock==ClockCombined && com.ngene<=1) 
      error2("Combined clock model requires mutliple genes.");
   GetInitialsTimes (x);

   com.plfun = (com.alpha==0 ? lfun : (com.rho==0?lfundG:lfunAdG));
   if(com.NSsites) com.plfun=lfundG;
   if(com.nparK) com.plfun=lfunAdG;

   if(com.method && com.fix_blength!=2 && com.plfun==lfundG) {
      com.conPSiteClass=1;
      sconP1_new*=com.ncatG;
   }
   if(com.sconP1<sconP1_new) {
      com.sconP1=sconP1_new;
      printf("\n%9d bytes for conP1, adjusted\n",com.sconP1);
      if((com.conP=(double*)realloc(com.conP,com.sconP1))==NULL) error2("oom conP1");
   }

   InitializeNodeScale();
   if(com.seqtype==1 && com.NSsites) /* mainly for NSsites = 1 & 2 (model = 0 or 2) */
      marksitesNonsyn();

   if(times++==0) {
      if((com.aaDist && com.aaDist<10 && com.aaDist!=AAClasses &&
          (com.seqtype==CODONseq||com.model==FromCodon)) ||
          (com.seqtype==AAseq &&
          (com.model==Empirical||com.model==Empirical_F||com.model>=REVaa_0))){
         GetDaa(NULL,com.daa);
      }
   }
   com.nrgene=(!com.fix_rgene)*(com.ngene-1);
   FOR(j,com.nrgene) x[com.ntime+j]=1;

/*
   if(com.seqtype==CODONseq && com.NSsites==0 && com.aaDist==0)
      com.nrate=!com.fix_kappa+!com.fix_omega;
*/
   if(com.seqtype==CODONseq && com.model && com.model<=NSbranch3) {
      if (com.model==NSbranchB) {
         com.nbtype=tree.nbranch;  
         FOR(i,tree.nbranch) nodes[(int)tree.branches[i][1]].label=i;
      }
      if(com.NSsites==0) {
         com.nOmega=com.nbtype;
         if(com.aaDist==0) 
            com.nrate=com.nkappa+!com.fix_omega+com.nbtype-1;
         else if (com.aaDist==AAClasses) 
            com.nrate=com.nkappa+com.nOmegaType*com.nbtype;
      }
   }

   if(com.Mgene>=3) com.nrate*=com.ngene;
   if(com.seqtype==1 && com.Mgene>=3 && com.fix_omega) com.nrate--; 
   com.np = com.ntime+com.nrgene+com.nrate;

   /* NSbranchsite models 
      com.nOmega=2 different w's at a site (three w's in the model: w0,w1,w2) */
   if(com.seqtype==CODONseq && com.model && com.NSsites) {
      if(com.nbtype!=2) error2("only two branch labels are allowed");
      if(com.model<=NSbranch2) { /* branch-site models A & B */
         com.ncatG=4;
         if(com.NSsites==NSdiscrete) 
            com.nrate = com.nkappa + 2 +!com.fix_omega+com.nbtype-1-1; /* add w0 and w1 */
         else 
            com.nrate= com.nkappa +!com.fix_omega+com.nbtype-1-1;
      }

      /* add p0 and p1.  check that this works for NSbranch2 */
      if(com.model==2) com.np = com.ntime+com.nrgene+com.nrate + com.ncatG-2;
      else             com.np = com.ntime+com.nrgene+com.nrate + com.ncatG-1;
      k=com.ntime+com.nrgene;
      if(com.hkyREV) {
         x[k++]=1+0.1*rndu(); 
         FOR(j,4) x[k++]=0.1+rndu();
      }
      else if(!com.fix_kappa)  x[k++]=com.kappa;

      if(com.model<=NSbranch2) {
         x[k++]=2.2;  if(com.ncatG==3) x[k++]=1.1;   /* p0 and p1 */
         if(com.NSsites==NSdiscrete) {
            x[k++]=0.1;  if(com.ncatG==3) x[k++]=0.8;   /* w0 and w1 */
         }
         if(!com.fix_omega)  x[k++]=3.5;  /* w2 */
      }
      else { /* NSbranch3: Joe's models C and D */
         x[k++]=2.2;  if(com.ncatG==3) x[k++]=1.1;   /* p0 and p1 */
         if(com.NSsites==NSdiscrete) {
            x[k++]=0.1;  if(com.ncatG==3) x[k++]=0.8;   /* w0 and w1 */
         }
         x[k++]=com.omega; x[k++]=com.omega;         /* additional w's */
      }
   }
   else if (com.NSsites) {        /* w's are counted in com.nrate */
      k=com.ntime+com.nrgene;  
      if(!com.fix_kappa)  x[k++]=com.kappa;
      switch(com.NSsites) {
      case(NSneutral):   com.np++; x[k++]=.7;  break;  /* p0 for w0=0 */
      case(NSselection): /* for p0, p1.  w is counted in nrate.  */
         com.np+=2; x[k++]=1.3; x[k++]=.3; 
         if(!com.fix_omega) { 
            x[k++]=com.omega;
            if(nnsmodels>1) x[k-1]=max2(2,com.omega);
         }
         break;  
      case(NSdiscrete):
         if(com.aaDist) {
            com.np+=com.ncatG-1; 
            FOR(i,com.ncatG-1) x[k++]=0.;
            if(com.aaDist<=6) FOR(i,com.ncatG) { x[k++]=1.1; x[k++]=1.2; }
            FOR(i,com.ncatG) /* ap,p*,av,v*, and b for each site class */
               FOR(j,4+(com.aaDist==FIT2)) x[k++]=rndu();
         }
         else if(com.nparK) { /* K*(K-1) paras in HMM of dN/dS over sites */
            zero(x+k,com.ncatG*(com.ncatG-1)); com.np+=com.ncatG*(com.ncatG-1);
         }
         else  {   /* p0...pK.  Note that w's are counted in nrate  */
            com.np+=com.ncatG-1; 
            FOR(i,com.ncatG-1) x[k++]=rndu();
            FOR(i,com.ncatG) 
               x[k++]=com.omega * (.5+i*2./com.ncatG*(0.8+0.4*rndu()));
         }
         break;
      case(NSfreqs):    /* p0...pK.  w's are fixed  */
         com.np+=com.ncatG-1; FOR(j,com.ncatG-1) x[k++]=(com.ncatG-j)/2.;
         break;
         break;
      case(NSgamma):  com.np+=2; x[k++]=1.1; x[k++]=1.1; break;
      case(NS2gamma):    /* p0, alpha1,beta1,alpha2=beta2 */
         com.np+=4; x[k++]=0.5; FOR(j,3) x[k++]=2*rndu()+j*0.1; break;
      case(NSbeta):       /* p_beta,q_beta */
         com.np+=2; x[k++]=.4; x[k++]=1.2+.1*rndu(); break; 
      case(NSbetaw):
      case(NSM8a):
      case(NSM8b):
         /* p0, p_beta, q_beta.  w is counted in nrate. */
         com.np+=3; x[k++]=.9; x[k++]=.2+rndu(); x[k++]=1.2+rndu();
         if(!com.fix_omega) {
            x[k++]=com.omega;
            if(nnsmodels>1) x[k-1]=rndu()*5;
         }
         break;
      case(NSbetagamma):  /* p0, p_beta, q_beta, alpha, beta */
         com.np+=5; x[k++]=.9; x[k++]=.4; x[k++]=1.2; x[k++]=1.1; x[k++]=1.1;
         break;
      case(NSbeta1gamma):  /* p0, p_beta, q_beta, alpha, beta */
         com.np+=5; x[k++]=.9; x[k++]=.4; x[k++]=1.2; x[k++]=.1; x[k++]=1.1;
         break;
      case(NSbeta1normal):  /* p0, p_beta, q_beta, alpha, beta */
         com.np+=5; x[k++]=.95; x[k++]=.4; x[k++]=1.2; x[k++]=1.1; x[k++]=1.1;
         break;
      case(NS02normal):    /* p0, p1, mu2, s1, s2 */
         com.np+=5; 
         x[k++]=.8; x[k++]=0.3;   /* p0 & p1, not transformed */
         x[k++]=.2; /* mu2 */ 
         x[k++]=5; x[k++]=1.1;  /* s1,s2 */
         break;
      case(NS3normal):    /* p0, p1, mu2, s0, s1, s2 */
         com.np+=6; 
         x[k++]=.77; x[k++]=0.22;   /* p0 & p1, transformed */
         x[k++]=.2; /* mu2 */ 
         x[k++]=0.5; x[k++]=5; x[k++]=1.1;  /* s0,s1,s2 */
         break;
      }
   }     /* if(com.NSsites) */

   k=com.ntime+com.nrgene;
   if (com.aaDist==AAClasses) { 
      if (!com.fix_kappa) x[k++]=com.kappa;
      FOR (i,com.nrate-!com.fix_kappa) x[k++]=com.omega;
      if (com.nOmegaType>65) 
         puts("\a\nget better initial values for AAclasses?");
   }
   else {
      if (com.seqtype==AAseq) {                     /* AAseq */
         if (com.nrate==0)  EigenQaa(NULL, Root, U, V, &t); /* once for all */
         if (com.model==REVaa_0) {
            FOR(i,naa) FOR(j,i) 
               if (AA1STEP[i*(i-1)/2+j] && i*naa+j!=ijAAref)
                  x[k++]=com.daa[i*naa+j]/com.daa[ijAAref];
         }
         else if (com.model==REVaa) { 
            for (i=1; i<naa; i++)  FOR(j,i)
               if(i*naa+j!=ijAAref) x[k++]=com.daa[i*naa+j]/com.daa[ijAAref];
         }
         else if (com.model==FromCodon) {
            FOR(j,com.nkappa)  x[k++]=com.kappa;
            FOR(j,com.nrate-com.nkappa)  x[k++]=com.omega; 
         }
      }
      else if (!com.NSsites) {                      /* CODONseq */
         if (com.nrate==0 && com.NSsites==0) 
            EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, 
            (com.nkappa>1?x+com.ntime+com.nrgene:&com.kappa), com.omega,PMat);
         else if (com.model==NSbranchB || com.model==NSbranch2) {
            if (!com.fix_kappa)  x[com.ntime+com.nrgene]=com.kappa; 
            FOR(i, (com.model==NSbranchB?tree.nbranch:com.nOmega-1+!com.fix_omega))
               x[com.ntime+com.nrgene+!com.fix_kappa+i]=com.omega;
         }
         else if(com.nrate) { /* either kappa, omega, or both for each gene */
            if(com.Mgene<=2) {
               if(com.hkyREV) { x[k++]=.5+rndu(); FOR(i,4) x[k++]=.1+rndu(); }
               else if (!com.fix_kappa) x[k++]=com.kappa;
               if(com.codonf>=FMutSel3x4) FOR(i,3) x[k++]=.5+rndu();
               if (!com.aaDist)
                  { if(!com.fix_omega)    x[k++]=com.omega; }
               else
                  { x[k++]=0.11; x[k++]=0.22; }
            }
            else { /* com.Mgene==3,4 */
               FOR(i,com.ngene) {
                  if(com.hkyREV) error2("hkyREV for ngene>1.  Fix me.");
                  if(!com.fix_kappa && !com.fix_omega)
                     { x[k+i*2]=com.kappa;  x[k+1+i*2]=com.omega; }
                  else if (com.fix_kappa) x[k+i]=com.omega;
                  else if (com.fix_omega) {
                     x[k+i*2]=com.kappa;  
                     if(i!=com.ngene-1) x[k+1+i*2]=com.omega;
                  }
               }
            }
         }
      }
   }
   for (i=0; i<com.nalpha; i++) x[com.np++]=com.alpha;

   if (!com.fix_rho) x[com.np++]=com.rho;
   if (com.rho)
      AutodGamma (com.MK, com.freqK, com.rK, &t, com.alpha, com.rho,com.ncatG);
   else if (com.alpha && com.fix_alpha && !com.NSsites)
      DiscreteGamma(com.freqK,com.rK,com.alpha,com.alpha,com.ncatG,0);

   if(com.fix_blength==-1)
      for(i=0; i<com.np; i++)  x[i] = (i<com.ntime ? .1+rndu() : .5+rndu());

   *fromfile=0;
   if(fin) {
      readx(x,fromfile);
      if(com.runmode>0 && fromfile && com.NSsites)  LASTROUND=1;
   }

   return (0);
}



int SetPGene (int igene, int _pi, int _UVRoot, int _alpha, double x[])
{
/* xcom[] does not contain time parameters
   Note that com.piG[][] have been homogeneized if (com.Mgene==3)
   Note calculation of nr1 for (com.Mgene>=3 && com.fix_omega), as only the 
   w for the last partition is fixed.
*/
   int nr1=(com.nrate+1)/com.ngene, k=com.nrgene+(com.Mgene>=3)*igene*nr1;
   double *xcom=x+com.ntime;

   if (_pi) {
      xtoy (com.piG[igene],com.pi,com.ncode);
      getpi_sqrt (com.pi, com.pi_sqrt, com.ncode, &com.npi0);
#if(defined(CODEML))
      if(com.codonf==F1x4MG || com.codonf==F3x4MG || com.codonf==FMutSel3x4)
         com.pf3x4MG=com.f3x4MG[igene];
#endif
   }
   if (_UVRoot) {
      if (com.seqtype==CODONseq) {
         if(!com.fix_kappa) com.kappa=xcom[k++];
         if(!com.fix_omega) com.omega=xcom[k++];
         else
            com.omega=(com.Mgene>2&&igene<com.ngene-1?xcom[k++]:com.omega_fix);
         if (!com.NSsites)
            EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, 
               (com.hkyREV||com.codonf>=FMutSel3x4?&xcom[com.nrgene]:&com.kappa),com.omega,PMat);
      }
      else
         EigenQaa(NULL, Root, U, V, xcom+k);
   }
   if (_alpha) {
      com.alpha=xcom[com.nrgene+com.nrate+igene];
      DiscreteGamma (com.freqK, com.rK, com.alpha, com.alpha, com.ncatG, 0);
   }
   return (0);
}


int SetParametersNSsites (double x[])
{
/* for NSsites and NSsitebranch models including HMM, NSbranchsite3 models
   p's are before w's in x[]
*/
   int k0=com.ntime+com.nrgene+(com.hkyREV?5:!com.fix_kappa), k=k0;
   int K=com.ncatG, i,j;
   double t, S,dS,dN, w0,w1,w2, spaceP2PI[NCATG*(NCATG+1)];
   double *pkappa=com.KAPPA;

   if(com.NSsites==0) error2("SetParametersNSsites : strange.");

   switch(com.NSsites) {
   case(NSneutral):  
      com.rK[0]=0; com.rK[1]=1;
      com.freqK[0]=x[k++]; com.freqK[1]=1-com.freqK[0]; break;
   case(NSselection): 
   case(NSdiscrete):
      if(com.model && com.model<=NSbranch2)   /* NSbranchsite A&B (Y&N2002) */
         K=3;
      if(com.nparK) {      /* HMM models, setting up p[] & w[] */
         FOR(j,K) com.rK[j]=x[k++];  /* w's for site classes */
         for (i=0; i<K; i++, k+=K-1) {
            if (!LASTROUND) f_and_x(x+k,com.MK+i*K,K,0,0);   /* x->f */
            else            xtoy  (x+k,com.MK+i*K,K-1);
            com.MK[i*K+K-1]=1-sum(com.MK+i*K,K-1);
         }
         PtoPi(com.MK, com.freqK, K, spaceP2PI);
         break;
      }

   case(NSfreqs):
      /* setting up p[] for NSselection, NSdiscrete, NSfreqs */
      if (!LASTROUND) {
         f_and_x(x+k,com.freqK,K,0,1);   /* x->f */
         k+=K-1;
      }
      else {
         for(j=0,com.freqK[K-1]=1;j<K-1;j++) 
            com.freqK[K-1]-=(com.freqK[j]=x[k++]);
         if(com.freqK[K-1]<0 || com.freqK[K-1]>1) 
            error2("freqK[]");
      }
      /* setting up w[] */
      if(com.NSsites==NSselection) {
         com.rK[0]=0; com.rK[1]=1; 
         com.rK[2]=(com.fix_omega?com.omega_fix:x[k++]); 
      }
      else if(com.NSsites==NSfreqs) {
         if(com.ncatG!=5) error2("NSfreqs, ncatG?");
         com.rK[0]=0; com.rK[1]=1./3; com.rK[2]=2./3; com.rK[3]=1; com.rK[4]=3;
      }
      else if(com.NSsites==NSdiscrete && com.aaDist==0)
         /* NSdiscrete: note this sets w0 & w1 for NSbranch2 & 3*/
         FOR(j,K) com.rK[j]=x[k++];
      break;

   case(NSgamma):
   case(NS2gamma):
   case(NSbeta):
   case(NSbetaw): 
   case(NSbetagamma):
   case(NSbeta1gamma):
   case(NSbeta1normal):
   case(NS02normal):
   case(NS3normal):
      DiscreteNSsites(x+k);  break;
   }
   /* calculates Qfactor_NS, to be used in eigenQc() for NSsites models */
   k=k0;
   if(com.model==0) {  /* NSsites models */
      for(j=0,Qfactor_NS=0; j<com.ncatG; j++) {
         freqK_NS=com.freqK[j];
         if(com.aaDist) {
            if(com.aaDist<10)         com.pomega=x+k+com.ncatG-1+2*j;
            else if(com.aaDist>=FIT1) {
               com.pomega=x+k+com.ncatG-1+j*(4+(com.aaDist==FIT2));
               xtoy(pcodonClass+j*64, com.pi, com.ncode);
            }
         }
         EigenQc(1,-1,&S,&dS,&dN,NULL,NULL,NULL, pkappa, com.rK[j], PMat);
      }
      Qfactor_NS=1/Qfactor_NS;
   }
   else if (com.model<=NSbranch2) { /* branch&site models */
      t=com.freqK[0]+com.freqK[1];
      com.freqK[2]=(1-t)*com.freqK[0]/(com.freqK[0]+com.freqK[1]);
      com.freqK[3]=(1-t)*com.freqK[1]/(com.freqK[0]+com.freqK[1]);
      com.rK[2]=com.rK[0]; com.rK[3]=com.rK[1];  /* w0 and w1 */
      /* calculates scale factors: background branches has two site classes
         while foreground branches has 3 site classes */
      FOR(i,2) {  /* i=0: background; i=1: foreground */
         for(j=0,Qfactor_NS=0; j<(i==0?2:3); j++) {
            w0=com.rK[j];
            if(i==0)       freqK_NS=com.freqK[j]/t;
            else if(j==2) {
               freqK_NS=1-t; 
               if(com.fix_omega) w0=com.omega_fix;
               else              w0=(com.NSsites==NSselection?x[k+2]:x[k+2+2]);
            }
            EigenQc(1,-1,&S,&dS,&dN,NULL,NULL,NULL, pkappa,w0,PMat);
         }
         Qfactor_NS_branch[i]=1/Qfactor_NS;
      }
      /* This calculates 5 sets of U&V&Root vectors (w0b,w0f,w1b,w1f,w2f), 
         which are used in GetPMatBranch():
            iclass=0: w0b,w0f
            iclass=1: w1b,w1f
            iclass=2: w0b,w2f
            iclass=3: w1b,w2f
         No eigenQc() calls are needed in ConditionalPNode() or minbranches().
      */
      k=k0+2;
      if(com.NSsites==NSselection) { w0=0; w1=1; }
      else                         { w0=x[k++]; w1=x[k++]; }
      w2=(!com.fix_omega?x[k++]:com.omega_fix);
      /* calculates UVRoot for 5 site classes. Note different Qfactor_NS */
      for(i=0; i<5; i++) {  /* (w0b,w0f,w1b,w1f,w2f) */
         Qfactor_NS = Qfactor_NS_branch[(i==1||i==3||i==4)];
         com.omega = (i<2 ? w0 : (i<4?w1:w2));
         EigenQc(0,-1,NULL,NULL,NULL,_Root[i],_UU[i],_VV[i], pkappa,com.omega,PMat);
      }
   }
   else { /* NSbranch3: Bielawski's branch&site models C and D */
      /* calculates scale factors: each branch has K=com.ncatG site classes */
      k=k0+K-1;   /* points to w[], by skipping the frequencies. */
      FOR(i,2) {  /* i=0: background branches; i=1: foreground */
         for(j=0,Qfactor_NS=0; j<K; j++) {
            freqK_NS=com.freqK[j];
            if(i==1 && j==K-1) w0 = (com.NSsites==2 ? x[k+1] : x[k+K-1+1]);
            else               w0 = com.rK[j];
            EigenQc(1,-1,&S,&dS,&dN,NULL,NULL,NULL,pkappa,w0,PMat);
         }
         Qfactor_NS_branch[i]=1/Qfactor_NS;
      }
      /* For K=3, this calculates 6 sets of U&V&Root vectors 
         (w0b,w0f,w1b,w1f,w2b,w2f), which are used in GetPMatBranch():
            iclass=0: w0b,w0f 
            iclass=1: w1b,w1f 
            iclass=2: w2b,w2f 
         no eigenQc() calls are needed in ConditionalPNode() or minbranches().
      */
      FOR(j,K) { /* (w0b,w0f,w1b,w1f,w2b,w2f) */
         FOR(i,2) {  /* i=0: background; i=1: foreground */
            Qfactor_NS = Qfactor_NS_branch[i];
            if(i==1 && j==K-1) w0 = (com.NSsites==2 ? x[k+1] : x[k+K-1+1]);
            else               w0 = com.rK[j];
            EigenQc(0,-1,NULL,NULL,NULL,_Root[j*2+i],_UU[j*2+i],_VV[j*2+i],pkappa,w0,PMat);
         }
      }
   }
   return(0);
}

int SetParameters (double x[])
{
/* Set com. variables and initialize U, V, Root etc. before each calculation 
   of the likelihood function.
   Is it a good idea to restruct this and/or Getinitials into two parts,
   one for aa's and another for codons?
   When (com.NSsites==NS02normal || NS3normal), p's are before w's in x[]; 
   see CDFdN_dS().
*/
   int i,j,k, ik=0, nUVR=6;
   double t,w0, *pkappa=com.KAPPA;

   if(com.fix_blength<2) SetBranch(x);
   if(com.np<=com.ntime) return(0);

   if(com.seqtype==1 || com.model==FromCodon || com.aaDist==AAClasses) {
      k=com.ntime+com.nrgene;
      if(com.hkyREV==0) {
         if(com.fix_kappa==1) { pkappa[0]=com.kappa; ik=1; }
         else                 com.kappa=x[k]; /* Is this necessary? */
      }

      for(i=0; i<(com.hkyREV?5:!com.fix_kappa); i++) pkappa[ik++]=x[k++];
      if(com.codonf>=FMutSel3x4) 
         for(i=0; i<3; i++) pkappa[ik++]=x[k++];

      com.pomega=x+com.ntime+com.nrgene+com.nkappa;
   }
   FOR(j,com.nrgene) com.rgene[j+1]=x[com.ntime+j];
   if(com.clock && AbsoluteRate) com.rgene[0]=x[0]; /* so that rgene are abs rates */

   if(com.seqtype==1 && com.aaDist>=FIT1) getpcodonClass(x, pcodonClass);

   k=com.ntime+com.nrgene+com.nkappa;
   if (com.nrate) {
      if(!com.model && !com.aaDist && !com.fix_omega && !com.NSsites) 
         com.omega=x[k];
      if(com.seqtype==AAseq)
         EigenQaa(NULL, Root, U, V, x+com.ntime+com.nrgene);
      else if(com.model==0 && com.NSsites==0 && com.Mgene<=1)
         /* CODONs, same dN/dS across branches & sites */ 
         EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, pkappa, com.omega,PMat);
      else if((com.model==NSbranchB || com.model==NSbranch2) && com.NSsites==0 
            && com.nbtype<=nUVR)
         FOR(i,com.nbtype) {
            if(com.aaDist==AAClasses)
               com.pomega=x+com.ntime+com.nrgene+com.nkappa+i*com.nOmegaType;
            else
               w0=(i==com.nOmega-1&&com.fix_omega?com.omega_fix:com.pomega[i]);
            EigenQc(0,-1,NULL,NULL,NULL,_Root[i],_UU[i],_VV[i], pkappa,w0,PMat);
         }
      k=com.ntime+com.nrgene+com.nrate;
   }
   if (com.seqtype==CODONseq && com.NSsites)
      SetParametersNSsites (x);
   if(com.model) com.omega=-1;  /* to force crash in case or error */
   /* branch models */
   if(com.seqtype==CODONseq && com.model && com.NSsites==0 && com.aaDist==0) {
      FOR(j,tree.nnode) {
         if(j==tree.root) continue;
         if (com.fix_omega && (int)nodes[j].label==com.nOmega-1)
            nodes[j].omega=com.omega_fix;
         else
            nodes[j].omega=com.pomega[(int)nodes[j].label];
      }
   }
   if (!com.fix_alpha && com.NSsites==0) {
      com.alpha=x[k++];
      if (com.fix_rho)
         DiscreteGamma(com.freqK,com.rK,com.alpha,com.alpha,com.ncatG,0);
   }
   if (!com.fix_rho) {
      com.rho=x[k++];
      AutodGamma(com.MK, com.freqK, com.rK, &t, com.alpha, com.rho, com.ncatG);
   }
   return (0);
}


int DiscreteNSsites(double par[])
{
/* This discretizes the continuous distribution for dN/dS ratios among sites
   and calculates freqK[] and rK[], using the median method.
   par[] contains all paras in the w distribution.  par[0] is the 
   proportion of beta if (com.NSsites==betaw), or the proportion of w=0 if 
   (com.NSsites=NS02normal).
   This routine uses com.NSsites, com.ncatG, com.freqK, com.rK.
   betaw has com.ncatG-1 site classes in the beta distribution, and 02normal 
   has com.ncatG-1 site classes in the mixed normal distribution.
   See the function CDFdN_dS() for definitions of parameters.
*/
   int status=0, j,off, K=com.ncatG-(com.NSsites==NSbetaw || com.NSsites==NS02normal);
   double xb[2]={1e-7,99};  /* bounds for omega.  */
   int K1=6, K2=4, UseK1K2=0;
   double p01=0, p,w0, lnbeta;

   if(com.NSsites==NSbeta || com.NSsites==NSbetaw) xb[1]=1;

#ifdef NSSITES_K1_K2_CLASSES
   if((com.NSsites==NSgamma||com.NSsites==NS2gamma||com.NSsites>=NSbetagamma)){
      K2=max2(K2,K/3);  K1=K-K2;  UseK1K2=1;
      p01=CDFdN_dS(1.,par);

      /* printf("\nK:%3d%3d\t p01=%9.5f\n",K1,K2,p01); */
      FOR(j,K) {
         if(j<K1) { p=p01*(j*2.+1)/(2.*K1); w0=p; }
         else     { p=p01+(1-p01)*((j-K1)*2.+1)/(2.*K2); w0=1.01+(j-K1)/K2; }
         com.rK[j]=InverseCDF(CDFdN_dS,p,w0,par,xb);
         com.freqK[j]=(j<K1 ? p01/K1 : (1-p01)/K2); thread
      }
   }
#endif

   if(!UseK1K2) { /* this is currently used */
      if(com.NSsites==NSbeta || com.NSsites==NSbetaw) {
         off=(com.NSsites==NSbetaw);  /* par[0] is proportion for beta for M8 */
         lnbeta=LnGamma(par[off])+LnGamma(par[off+1])-LnGamma(par[off]+par[off+1]);
         for(j=0; j<K; j++) {
            p=(j*2.+1)/(2.*K);
            com.rK[j]=InverseCDFBeta(p, par[off], par[off+1], lnbeta);
         }
      }
      else {
         FOR(j,K) {
            p=(j*2.+1)/(2.*K);
            w0=.01+j/K; if(com.rK[j]) w0=(w0+com.rK[j])/2;
            com.rK[j]=InverseCDF(CDFdN_dS,p,w0,par,xb);
         }
      }
      FOR(j,K) com.freqK[j]=1./K;
   }

   if(com.NSsites==NSbetaw) {
      if(!com.fix_omega) com.rK[com.ncatG-1]=par[3];
      else               com.rK[com.ncatG-1]=com.omega_fix;
      com.freqK[K]=1-par[0]; FOR(j,K) com.freqK[j]*=par[0];
   }
   if(com.NSsites==NS02normal) {
      for(j=K-1;j>=0;j--) /* shift to right by 1 to make room for spike at 0*/
         { com.rK[j+1]=com.rK[j]; com.freqK[j+1]=com.freqK[j];  }
      com.rK[0]=0;  com.freqK[0]=par[0];
      for(j=1;j<K+1;j++) com.freqK[j]*=(1-par[0]);
   }

   if(com.NSsites>=NSgamma){
      if(!status && com.NSsites==NSbeta) 
         for(j=1;j<com.ncatG;j++) if(com.rK[j]+1e-7<com.rK[j-1]) status=1;

      if(status) {
         printf("\nwarning: DiscreteNSsites\nparameters: ");
         FOR(j,(com.NSsites==7?2:4)) printf(" %12.6f", par[j]);  FPN(F0);
         FOR(j,com.ncatG)            printf("%13.5f", com.freqK[j]);  FPN(F0);
         FOR(j,com.ncatG)            printf("%13.5e", com.rK[j]);  FPN(F0);
      }
   }
   return(0);
}


double CDFdN_dS(double x,double p[])
{
/* This calculates the CDF of the continuous dN/dS distribution over sites, 
   to be used as argument to the routine InverseCDF().  When the distribution
   has spikes, the spikes are ignored in this routine, and the scaling
   is done outside this routine, for example, in DiscreteNSsites().
   All parameters (par) for the w distribution are passed to this routine, 
   although some (p0 for the spike at 0) are not used in this routine.  
   Parameters are arranged in the following order:

      NSgamma (2):       alpha, beta
      NS2gamma (4):      p0, alpha1, beta1, alpha2 (=beta2)
      NSbeta (2):        p_beta, q_beta
      NSbetaw (4):       p0, p_beta, q_beta, w (if !com.fix_omega, not used here)
      NSbetagamma (5):   p0, p_beta, q_beta, alpha, beta
      NSbeta1gamma (5):  p0, p_beta, q_beta, alpha, beta (1+gamma)
      NSbeta1normal (5): p0, p_beta, q_beta, mu, s (normal>1)
      NS02normal (5):    p0, p1, mu2, s1, s2 (s are sigma's)
      NS3normal (6):     p0, p1, mu2, s0, s1, s2 (s are sigma's)

   Parameters p0 & p1 are transformed if (!LASTROUND)

*/
   double cdf=-1;
   double z, f[3],mu[3]={0,1,2},sig[3]; /* 3normal: mu0=0 fixed. mu2 estimated */

   switch(com.NSsites) {
   case(NSgamma):  cdf=CDFGamma(x,p[0],p[1]);   break;
   case(NS2gamma): 
      cdf=p[0] *CDFGamma(x,p[1],p[2])+(1-p[0])*CDFGamma(x,p[3],p[3]);  break;
   case(NSbeta):   cdf=CDFBeta(x,p[0],p[1],0);  break;
   case(NSbetaw):  cdf=CDFBeta(x,p[1],p[2],0);  break;
   case(NSbetagamma):
      cdf=p[0]*CDFBeta(x,p[1],p[2],0)+(1-p[0])*CDFGamma(x,p[3],p[4]);  break;

   case(NSbeta1gamma):
      if(x<=1) cdf=p[0]*CDFBeta(x,p[1],p[2],0);
      else     cdf=p[0]+(1-p[0])*CDFGamma(x-1,p[3],p[4]);
      break;
   case(NSbeta1normal):
      if(x<=1) cdf=p[0]*CDFBeta(x,p[1],p[2],0);
      else {
         cdf=CDFNormal((p[3]-1)/p[4]);
         if(cdf<1e-9) {
            matout(F0,p,1,5);;
            printf("PHI(%.6f)=%.6f\n",(p[3]-1)/p[4],cdf);  getchar();
         }
         cdf=p[0]+(1-p[0])*(1- CDFNormal((p[3]-x)/p[4])/cdf);
      }
      break;
   case(NS02normal):
      mu[2]=p[2]; sig[1]=p[3]; sig[2]=p[4];
      f[1]=p[1];  f[2]=1-f[1];
      cdf = 1 - f[1]* CDFNormal(-(x-mu[1])/sig[1])/CDFNormal(mu[1]/sig[1])
              - f[2]* CDFNormal(-(x-mu[2])/sig[2])/CDFNormal(mu[2]/sig[2]);
      break;
   case(NS3normal):
      mu[2]=p[2]; sig[0]=p[3]; sig[1]=p[4]; sig[2]=p[5];

      if(LASTROUND) { f[0]=p[0]; f[1]=p[1]; }
      else          { z=(f[0]=exp(p[0]))+(f[1]=exp(p[1]))+1; f[0]/=z; f[1]/=z;}
      f[2]=1-f[0]-f[1];
      cdf = 1 - f[0]* 2*CDFNormal(-x/sig[0])
              - f[1]* CDFNormal(-(x-mu[1])/sig[1])/CDFNormal(mu[1]/sig[1])
              - f[2]* CDFNormal(-(x-mu[2])/sig[2])/CDFNormal(mu[2]/sig[2]);
      break;
   }
   return(cdf);
}


void GetSNnew(double pi[], double *Snew, double *Nnew, double *S4)
{
/* this calculates the synonymous and nonsynonymous sites according to the new 
   definition.
   S and N are sites per codon.
   It is not clear how to deal with stop codons.
*/
   int i,j,k, ic,b[3], aa0,aa1, *code=GeneticCode[com.icode];
   int by[3]={16,4,1}, nstop,s,n;
   double y;

   for(i=0,*Snew=*Nnew=*S4=0; i<com.ncode; i++) {
      ic=FROM61[i]; b[0]=ic/16; b[1]=(ic/4)%4; b[2]=ic%4;
      /* no need to check the first and second positions here */
      if(FourFold[b[0]][b[1]]) *S4 += pi[i];
      aa0=code[ic];

      for(j=0,s=n=nstop=0; j<3; j++) FOR(k,3) {
         aa1 = code[ic + ((b[j]+k+1)%4 - b[j])*by[j]];
         if(aa1==-1)        nstop++;
         else if(aa0==aa1)  s++;
         else               n++;
      }
      /* s + n ~= 9 */
      *Snew += pi[i]*s/9.*3.;
      *Nnew += pi[i]*n/9.*3.;
   }
   y = (*Snew + *Nnew)/3;
   *Snew /= y;  *Nnew /= y;
}


int EigenQc (int getstats, double blength, double *S, double *dS, double *dN,
    double Root[], double U[], double V[],
    double kappa[], double omega, double Q[])
{
/* This contructs the rate matrix Q for codon substitution and get the eigen
   values and vectors if getstats==0, or get statistics (dS & dN) if 
   getstats==1.  
   The routine is also called by Qcodon2aa for mechanistic amino acid 
   substitution models.
   Input parameters are kappa, omega and com.pi (or com.fb61).

   Statistics calculated include S, dS & dN.
   c0[0,1,2] and c[0,1,2] are rates for the 3 codon positions before and after 
   selection.  c4 is for 4-fold rates.  ts[3] and tv[3] are transition/
   transversion rates, not calculated.

   Under NSsites or other site-class models, this function does not scale 
   Q but calculates the Qfactor_NS.
   DetailOutput() uses this function to calculate dS & dN; for each omega 
   component, dS and dN are not scaled here but are scaled in DetailOutput()
   after Qfactor_NS is calculated.

   aaDist=FIT1 & FIT2:  ap,p*,av,v*, (and w0 for FIT2)
   The argument omega is used only if the model assumes one omega.  For 
   AAClasses, com.pomega is used instead.
*/
   int n=Nsensecodon, i,j,k, ic1,ic2,aa1,aa2, b1,b2;
   int ndiff,pos=0,from[3],to[3];
   double mr, rs0,ra0,rs,ra, y; /* rho's */
   double Snew, Nnew, S4;
   double d4=0, d0[3],d[3],ts[3],tv[3];  /* rates at positions and 4-fold sites */
   double *pi=(com.seqtype==AAseq?com.fb61:com.pi), w=-1, fit1,fit2, pijQij;
   double *pomega=com.pomega;

   NEigenQ++;
   if(blength>=0 && (S==NULL||dS==NULL||dN==NULL)) error2("EigenQc");
   FOR(i,3) d[i]=d0[i]=0;  FOR(i,3) ts[i]=tv[i]=0;
   FOR(i,n*n) Q[i]=0;
   for (i=0, rs0=ra0=rs=ra=0; i<n; i++) {
      ic1=FROM61[i]; from[0]=ic1/16; from[1]=(ic1/4)%4; from[2]=ic1%4;
      for(j=0; j<i; j++) {
         ic2=FROM61[j]; to[0]=ic2/16; to[1]=(ic2/4)%4; to[2]=ic2%4;
         for (k=0,ndiff=0; k<3; k++)  if (from[k]!=to[k]) { ndiff++; pos=k; }
         if (ndiff!=1)  continue;
         Q[i*n+j]=1;
         if(com.hkyREV) { /* REV model */
            b1=min2(from[pos],to[pos]); /* b1 and b2 are changed nucleotides */
            b2=max2(from[pos],to[pos]);
            if     (b1==0 && b2==1) Q[i*n+j]=kappa[0]; /* TC or CT, relative to AG */
            else if(b1==0 && b2==2) Q[i*n+j]=kappa[1]; /* TA */
            else if(b1==0 && b2==3) Q[i*n+j]=kappa[2]; /* TG */
            else if(b1==1 && b2==2) Q[i*n+j]=kappa[3]; /* CA */
            else if(b1==1 && b2==3) Q[i*n+j]=kappa[4]; /* CG */
         }
         else             /* HKY model */
            if(from[pos]+to[pos]==1 || from[pos]+to[pos]==5)
               Q[i*n+j]= kappa[0];

         /* MuseGaut-style models & mutation-selection model of Yang (2003) */
         if (com.codonf>=F1x4MG) {
            /* b1 and b2 are the 2 unchanged positions */
            if     (pos==0) { b1=1; b2=2; }
            else if(pos==1) { b1=2; b2=0; }
            else            { b1=0; b2=1; }
            if (com.codonf<=F3x4MG)
               Q[i*n+j] /= com.pf3x4MG[b1*4+to[b1]]*com.pf3x4MG[b2*4+to[b2]];
            else {  /* if(FMutSel3x4 || FMutSelCodon), 3 paras estimated */
               y=1;
               if(to[b1]!=3)   y  = kappa[(com.hkyREV?5:1)+to[b1]];
               if(to[b2]!=3)   y *= kappa[(com.hkyREV?5:1)+to[b2]];
               Q[i*n+j] /= y;
            }
         }

         Q[j*n+i]=Q[i*n+j];
         if(com.codonf)
            {  Q[i*n+j]*=pi[j];  Q[j*n+i]*=pi[i]; }

         pijQij=2*pi[i]*Q[i*n+j];
         if(getstats) {
            d0[pos] += pijQij;
            if(pos==2 && FourFold[to[0]][to[1]]) 
               d4 += pijQij;
         }

         aa1=GeneticCode[com.icode][ic1];  aa2=GeneticCode[com.icode][ic2];
         if(aa1==aa2)
            rs+=pijQij;
         else {
            ra0+=pijQij;  w=1;
            if (com.aaDist==AAClasses) {
               if (aa1<aa2)  { k=aa2; aa2=aa1; aa1=k; }
               k=aa1*(aa1-1)/2+aa2;
               if (pomega[OmegaAA[k]]<0) {
                  if (noisy)  printf("aa1 & aa2 & iw & w: %d %d %d %.5f\n", 
                                 aa1,aa2,OmegaAA[k],pomega[OmegaAA[k]]);
                  pomega[OmegaAA[k]]=0;
               }
               if (com.seqtype==AAseq && com.nrate>65 && aa1*20+aa2==ijAAref)
                   ;     /* if estimating grantham's matrix with aa sequences */
               else  w = pomega[OmegaAA[k]];
            }
            else if (com.aaDist==0)  w = omega; /* NSsites==0 or >0 */
            else if (com.aaDist<=6)  {          /* chemical properties: a & b */
               w = pomega[0]*com.daa[aa1*20+aa2];
               if(com.aaDist>0)           w = exp(-w);  /* geometric */
               else                       w = 1-w;      /* linear */
               if (com.seqtype==CODONseq) w *= pomega[1];
            }
            else if (com.aaDist>=FIT1) {   /* ap,p*,av,v* (and w0 for FIT2) */
               fit1=-pomega[0]*square(AAchem[0][aa1]-pomega[1])
                    -pomega[2]*square(AAchem[1][aa1]-pomega[3]);
               fit2=-pomega[0]*square(AAchem[0][aa2]-pomega[1])
                    -pomega[2]*square(AAchem[1][aa2]-pomega[3]);

               w = exp(-fit1-fit2);
               if(com.aaDist==FIT2) w *= pomega[4];
            }
            Q[i*n+j]*=w; Q[j*n+i]*=w;
            ra+=pijQij*w;
         }
         if(getstats) d[pos] += pijQij*(aa1==aa2?1:w);
      } /* for (j) */
   }    /* for (i) */

   mr=rs+ra;
   if(getstats) {
      if(com.NSsites) Qfactor_NS+=freqK_NS * (rs+ra);
      rs0=rs;
      w=(rs0+ra0);  rs0/=w;  ra0/=w;   *S=rs0*3*com.ls;
      if(!com.NSsites && blength>=0) {  /* calculates dS & dN */
         if(blength==0) *dS = *dN = 0;
         rs/=mr; ra/=mr;
         *dS=blength*rs/(3*rs0);  *dN=blength*ra/(3*ra0);
      }
      else if (com.NSsites)
         { *dS=rs/(rs0*3);  *dN=ra/(ra0*3); }
      if(blength>=0) {
         GetSNnew(com.pi, &Snew, &Nnew, &S4);

         FOR(i,3) { d[i]*=blength/mr;  d0[i]*=blength/mr; }
         d4*=blength/mr/S4;

/*
if(noisy>=3) {
   printf("\nd123[*] =%9.5f%9.5f%9.5f  average%9.5f\n", d[0],d[1],d[2], (d[0]+d[1]+d[2])/3);
   printf(  "    [B] =%9.5f%9.5f%9.5f  average%9.5f\n", d0[0],d0[1],d0[2], (d0[0]+d0[1]+d0[2])/3);
   printf("accept  =%9.5f%9.5f%9.5f\n\n", d[0]/d0[0],d[1]/d0[1],d[2]/d0[2]);
   printf("w =%9.5f dN =%9.5f dS =%9.5f d4 =%9.5f (%.1f four-fold sites)\n", *dN/ *dS, *dN,*dS, d4, S4*com.ls);
   printf("%12s dN*=%9.5f dS*=%9.5f S* =%7.2f N* =%7.2f\n", "",
      blength*ra/Nnew, blength*rs/Snew, Snew*com.ls, Nnew*com.ls);
}
*/

if(com.runmode==-2) {
fprintf(frst1, "\td123Ad123A\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\tw123\t%.4f\t%.4f\t%.4f\tdSdNd4w\t%.4f\t%.4f\t%.4f\t%.4f",
        d0[0],d0[1],d0[2], (d0[0]+d0[1]+d0[2])/3, 
        d[0],d[1],d[2], (d[0]+d[1]+d[2])/3,
        d[0]/d0[0],d[1]/d0[1],d[2]/d0[2],
        *dS, *dN, d4, omega);
fprintf(frst1, "\tkappa");
FOR(i,com.nkappa+(com.hkyREV==0&&com.fix_kappa)) 
fprintf(frst1, "\t%.4f", kappa[i]);
}
else
  fprintf(frst1, "\t%.4f\t%.4f", *dS, *dN);


      }
   }
   else {  /* get Root, U, & V */
      if (com.seqtype==AAseq) return (0);
      for (i=0; i<n; i++)
        Q[i*n+i]=-sum(Q+i*n,n);

      eigenQREV(Q, com.pi, com.pi_sqrt, n, com.npi0, Root, U, V);

      if(com.NSsites)  mr=1/Qfactor_NS;  /* Qfactor_NS calculated in SetParameter */
      FOR(i,n) Root[i]/=mr;
      UVRootChanged=1;
   }
   return (0);
}


int EigenQaa (FILE *fout, double Root[], double U[], double V[], double rate[])
{
/*  Codon-based AA model must use FromCodon, even if com.aaDist==AAClasses.
*/
   int naa=20, i,j,k;
   double Q[20*20], mr, t=0, space[NCODE*NCODE*2],*Qc=space+NCODE*NCODE;
   char aa3[4]="";

   FOR (i,naa) Q[i*naa+i]=0;
   switch (com.model) {
   case (Poisson)   : case (EqualInput) : 
      fillxc (Q, 1., naa*naa);  break;
   case (Empirical)   : case (Empirical_F):
      FOR(i,naa) FOR(j,i) Q[i*naa+j]=Q[j*naa+i]=com.daa[i*naa+j]/100;
      break;
   case (FromCodon):
      EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, 
         (com.hkyREV||com.codonf>=FMutSel3x4?rate:&com.kappa),com.omega,Qc);
      Qcodon2aa(Qc, com.fb61, Q, space);
      break;
   case (REVaa_0)  :
      for (i=1,k=0; i<naa; i++) for (j=0; j<i; j++)
         if (AA1STEP[i*(i-1)/2+j] && i*naa+j!=ijAAref)
            Q[i*naa+j]=Q[j*naa+i]=rate[k++];
      k=ijAAref;  Q[(k/naa)*naa+k%naa]=Q[(k%naa)*naa+k/naa]=1;
      break;
   case (REVaa)  : 
      for (i=0,k=0; i<naa; i++) FOR (j,i)
         if (i*naa+j!=ijAAref) Q[i*naa+j]=Q[j*naa+i]=rate[k++];
      Q[ijAAref]=Q[(ijAAref%naa)*naa+(ijAAref/naa)]=1; 
      break;
   }
   FOR (i,naa) FOR (j,naa) Q[i*naa+j]*=com.pi[j];
   for (i=0,mr=0; i<naa; i++) {
      Q[i*naa+i]=-sum(Q+i*naa,naa);  mr-=com.pi[i]*Q[i*naa+i]; 
   }
   if (fout && com.model>=REVaa_0) {
      fprintf (fout, "\n\nRate matrix (symmetrical part, Sij)\n");
      for(i=0,t=0;i<naa;i++) {
         if(com.pi[i]==0) error2("EigenQaa: do this now");
         FOR(j,i) t+=Q[i*naa+j]/com.pi[j]/(naa*(naa-1)/2.);
      }
      FOR (i,naa) {
         fprintf (fout, "\n%-5s", getAAstr(aa3,i));
         FOR (j,i) fprintf (fout, " %4.0f", Q[i*naa+j]/t/com.pi[j]*100);
      }
      fputs("\n     ",fout);  FOR(i,naa) fprintf(fout,"%5s",getAAstr(aa3,i));
      FPN(fout);  fflush(fout);
   }

   if (fout && frst1 && com.model>=REVaa_0) {
      fprintf(frst1, "\nRate matrix (symmetrical part, Sij) for bubble plot\n");
      FOR (i,naa)  FOR (j,i) 
         fprintf(frst1, "\t%d\t%d\t%.2f\n", i+1,j+1,Q[i*naa+j]/t/com.pi[j]*100);
   }

   eigenQREV(Q, com.pi, com.pi_sqrt, naa, com.npi0, Root, U, V);
   FOR(i,naa)  Root[i]=Root[i]/mr;
   UVRootChanged=1;

#if (DEBUG)
      if (fout) {
         fprintf (fout, "\nPI & Root\n");
         matout (fout, com.pi, 1, 20);
         matout (fout, Root, 1, 20);
         fprintf (fout, "\nPAM matrix, P(0.01)");
         PMatUVRoot (PMat, 0.01, naa, U, V, Root);

         FOR (i,naa) {
            fprintf (fout, "\n%5s", getAAstr(aa3,i));
            FOR (j,naa)  fprintf (fout, "%6.0f", PMat[i*naa+j]*100000);
         }
         FPN (fout);
      }
#endif
   return (0);
}


int Qcodon2aa (double Qc[], double pic[], double Qaa[], double piaa[])
{
/* Qc -> Qaa

   This routine constructs the rate matrix for amino acid replacement from
   the rate matrix for codon substitution, by congregating states in the
   Markov chain.  Both processes are time reversible, and only the
   symmetrical part of the rate matrix are constructed.  Codon frequencies 
   pic[] are used.  They should be filled with 1/61 if only amino acid
   sequences are available. 
   Qaa(aai,aaj) = SUMi SUMj (piC[i]*piC[j]]*Qc[i][j]) / (piAA[i]*piAA[j])
*/
   int i, j, aai, aaj, nc=Nsensecodon, naa=20;
   double ti, tij;

   zero(piaa,naa);  zero(Qaa,naa*naa);
   FOR (i,nc) piaa[GeneticCode[com.icode][FROM61[i]]] += pic[i];
   FOR (i,nc) {
      aai=GeneticCode[com.icode][FROM61[i]];
      ti=pic[i]/piaa[aai];
      FOR (j, i) {
         aaj=GeneticCode[com.icode][FROM61[j]];
         if (Qc[i*nc+j]==0 || aai==aaj) continue;
         tij=ti*pic[j]*Qc[i*nc+j]/piaa[aaj];
         Qaa[aai*naa+aaj]+=tij;
         Qaa[aaj*naa+aai]+=tij;
      }
   }

/*
for (i=0;i<naa;i++,FPN(frst)) FOR(j,i) {
   fprintf(frst,"%2d",(Qaa[i*naa+j]>0));
   if (Qaa[i*naa+j]<0) printf ("error2 Qij: %d %d \n", i,j);
}
fflush(frst);
*/
   return (0);
}



int ConditionalPNode (int inode, int igene, double x[])
{
   int n=com.ncode, i,j,k,h, ison, pos0=com.posG[igene],pos1=com.posG[igene+1];
   double t, smallw=1e-12;

   FOR(i,nodes[inode].nson)
      if(nodes[nodes[inode].sons[i]].nson>0 && !com.oldconP[nodes[inode].sons[i]])
         ConditionalPNode(nodes[inode].sons[i], igene, x);

   if(inode<com.ns)
      for(h=pos0*n; h<pos1*n; h++) nodes[inode].conP[h]=0; /* young ancestor */
   else 
      for(h=pos0*n; h<pos1*n; h++) nodes[inode].conP[h]=1;
   if (com.cleandata && inode<com.ns)
      for(h=pos0;h<pos1;h++) nodes[inode].conP[h*n+com.z[inode][h]]=1;

   if(com.NSsites && !IClass && com.rK[IClass]<smallw) {
      for(h=pos0; h<pos1; h++)
         if(sitesNonsyn[h])  FOR(j,n) nodes[inode].conP[h*n+j]=0;
   }
   else if(com.NSsites==2 && com.model==2 && IClass==2)
      for(h=pos0; h<pos1; h++)
         if(sitesNonsyn[com.npatt+h])  FOR(j,n) nodes[inode].conP[h*n+j]=0;

   FOR(i,nodes[inode].nson) {
      ison=nodes[inode].sons[i];
      t=nodes[ison].branch*_rateSite;

      if(com.clock)  t *= GetBranchRate(igene,(int)nodes[ison].label,x,NULL);
      else           t *= com.rgene[igene];
      GetPMatBranch(PMat, x, t, ison);

      if(com.cleandata && nodes[ison].nson<1)  /* tips */
         for(h=pos0; h<pos1; h++) {
            if((com.NSsites && !IClass && com.rK[IClass]<smallw && sitesNonsyn[h])
             ||(com.NSsites==2 && com.model==2 && IClass==2 && sitesNonsyn[com.npatt+h]))
               continue;
            FOR(j,n) nodes[inode].conP[h*n+j]*=PMat[j*n+com.z[ison][h]];
         }
      else
         for(h=pos0; h<pos1; h++) {
            if((com.NSsites && !IClass && com.rK[IClass]<smallw && sitesNonsyn[h])
             ||(com.NSsites==2 && com.model==2 && IClass==2 && sitesNonsyn[com.npatt+h]))
               continue;
            FOR(j,n) {
               for(k=0,t=0; k<n; k++)
                  t+=PMat[j*n+k]*nodes[ison].conP[h*n+k];
               nodes[inode].conP[h*n+j]*=t;
            }
         }
   }        /*  for (ison)  */
   if(com.nnodeScale && com.nodeScale[inode]) 
      NodeScale(inode, pos0, pos1);

   return (0);
}




int PMatJC69like (double P[], double t, int n)
{
   int i;
   double pii=1./n+(1.-1./n)*exp(-n/(n-1.)*t), pij=(1.-pii)/(n-1.);
   FOR(i,n*n) P[i]=pij;
   FOR(i,n) P[i*n+i]=pii;
   return (0);
}


int setmark_61_64 (void)
{
/* This sets two matrices FROM61[], and FROM64[], which translate between two 
   codings of codons.  In one coding, codons go from 0, 1, ..., 63 while in 
   the other codons range from 0, 1, ..., 61 with the three stop codons removed.
   FROM61[] translates from the 61-state coding to the 64-state coding, while 
   FROM64[] translates from the 64-state coding to the 61-state coding.

   This routine also sets up FourFold[4][4], which defines the 4-fold codon
   boxes.
*/
   int i,j,k, *code=GeneticCode[com.icode];
   int c[3],aa0,aa1, by[3]={16,4,1};
   double nSilent, nStop, nRepl;

   Nsensecodon=0;
   for (i=0; i<64; i++) {
      if (code[i]==-1)  FROM64[i]=-1; 
      else            { FROM61[Nsensecodon]=i; FROM64[i]=Nsensecodon++; }
   }
   com.ncode=Nsensecodon;

   FOR(i,4) FOR(j,4) {
      k=i*16+j*4;
      FourFold[i][j] = (code[k]==code[k+1] && code[k]==code[k+2] && code[k]==code[k+3]);
   }

   for (i=0,nSilent=nStop=nRepl=0; i<64; i++) {
      c[0]=i/16; c[1]=(i/4)%4; c[2]=i%4;
      if((aa0=code[i])==-1) continue;
      FOR(j,3) FOR(k,3) {
         aa1 = code[i + ((c[j]+k+1)%4 - c[j])*by[j]];
         if(aa1==-1)       nStop++;
         else if(aa0==aa1) nSilent++;
         else              nRepl++;
      }
   }
   printf("\ncode Stop Silent Replace\n");
   printf("%3d (%d)  %6.0f%6.0f%6.0f  %12.6f%12.6f\n", 
      com.icode, 64-com.ncode, nStop,nSilent,nRepl,nStop*3/(com.ncode*9),nSilent*3/(com.ncode*9));
   return (0);
}


int CanSkipThisSite (int h, char seqkept[])
{
/* Has site h have nonsynonymous differences among seqs kept?
   This is used for pick up sites which have prob=0 when w=0 
   under NSsites=1 and 2.
*/
   int i,j, nb[3], ib[3][4],ic, aa[2], allsyn=1, iAA=0;
   char codon[4]=" ";

   FOR(j,com.ns) {
      if(!seqkept[j]) continue;
      if(com.cleandata)
         aa[iAA]=GeneticCode[com.icode][FROM61[com.z[j][h]]];
      else {
         FOR(i,3) codon[i]=com.z[j][h*3+i];
         CodonListall(codon, nb, ib);
         if(nb[0]*nb[1]*nb[2]!=1) continue;
         ic=ib[0][0]*16+ib[1][0]*4+ib[2][0];
         aa[iAA]=GeneticCode[com.icode][ic];
      }
      if(iAA && aa[1]!=aa[0]) { allsyn=0; break; }
      iAA=1;
   }
   return (!allsyn);  /* 1 for yes, we can skip this site */
}


int InList (int key, int list[], int n)
{
/* Is key in list?
*/
   int i;
   for(i=0; i<n; i++) if(key==list[i]) return 1;
   return 0;
}

void marksitesNonsyn(void)
{
/* this marks sites (with 1) at which some differences are nonsynonymous.  
   Such sites have probability 0 when omega = 0 and may be skipped in conP 
   calculation under NSsites models.  For nforeground branches, the algorithm
   checks nforeground+1 sets of species: those that are descendent tips and 
   another set of non-descendents.
   com.z[] has coded data if (com.cleandata==1).
*/
   int  j,k, nclass, h, inode=-1;
   int  ifore, nforeground, nodesfore[NS*2]; /* foreground nodes */
   char seqkept[NS];

   nclass = (com.NSsites==2 && com.model==2 ? 2 : 1);
   if((sitesNonsyn=(char*)realloc(sitesNonsyn, com.npatt*nclass*sizeof(char)))==NULL)
      error2("oom SitesNonsyn");
   FOR(h,nclass*com.npatt) sitesNonsyn[h]=0;

   /* w0=0 */
   FOR(j,com.ns) seqkept[j]=1;
   FOR(h,com.npatt) sitesNonsyn[h]=(char)CanSkipThisSite(h,seqkept);

   /* w2=0 under branch-site model A */
   if(nclass==2) {
      for(j=0,nforeground=0;j<tree.nnode;j++) 
         if(j!=tree.root && nodes[j].label) nodesfore[nforeground++]=j;

      printf("\n%d foreground branches\n", nforeground);
      for(ifore=0; ifore<nforeground; ifore++) {
         inode=nodesfore[ifore];
         printf("\nforeground #%d: %d..%d\n", ifore+1,nodes[inode].father+1,inode+1);
         /* (A) descendent tips */
         FOR(j,com.ns) {
            seqkept[j]=0;
            for(k=j; k!=tree.root && k!=inode; k=nodes[k].father)
               if(InList(k,nodesfore,nforeground)) break;
            if(k==inode) seqkept[j]=1;
         }
         if(inode>=com.ns) FOR(h,com.npatt)  
            if(CanSkipThisSite(h,seqkept)) sitesNonsyn[com.npatt+h]=1;
      }
      /* (B) non-descendent tips, on the other side of the tree */
      FOR(j,com.ns) {
         seqkept[j]=0;
         for(k=j; k!=tree.root; k=nodes[k].father) 
            if(InList(k,nodesfore,nforeground)) break;
         if(k==tree.root) seqkept[j]=1;
      }
      FOR(h,com.npatt)  
         if(CanSkipThisSite(h,seqkept)) sitesNonsyn[com.npatt+h]=1;
   }
}


int printfcode (FILE *fout, double fb61[], double space[])
{
/* space[64*2]
*/
   int i, n=Nsensecodon;

   fprintf (fout, "\nCodon freq.,  x 10000\n");
   zero (space, 64);
   FOR(i,n) space[FROM61[i]] = fb61[i]*10000;
   printcu(fout, space, com.icode);
   return(0);
}


int printsmaCodon (FILE *fout,char * z[],int ns,int ls,int lline,int simple)
{
/* print, in blocks, multiple aligned and transformed codon sequences.
   indels removed.
   This is needed as codons are coded 0,1, 2, ..., 60, and 
   printsma won't work.
*/
   int ig, ngroup, lt, il,is, i,b, lspname=20;
   char equal='.',*pz, c0[4],c[4];

   ngroup = (ls-1)/lline + 1;
   for (ig=0,FPN(fout); ig<ngroup; ig++)  {
      fprintf (fout,"%-8d\n", ig*lline+1);
      FOR (is,ns) {
         fprintf(fout,"%-*s  ", lspname,com.spname[is]);
         lt=0; 
         for(il=ig*lline,pz=z[is]+il; lt<lline && il<ls; il++,lt++,pz++) {
            b=*pz;  b=FROM61[b]; 
            c[0]=(char)(b/16); c[1]=(char)((b%16)/4); c[2]=(char)(b%4); c[3]=0;
            FOR(i,3) c[i]=BASEs[c[i]];
            if (is && simple)  {
               b=z[0][il];  b=FROM61[b];
               c0[0]=(char)(b/16); c0[1]=(char)((b%16)/4); c0[2]=(char)(b%4);
               FOR(i,3) if (c[i]==BASEs[c0[i]]) c[i]=equal;
            }
           fprintf(fout,"%3s ", c);
         }
         FPN (fout);
      }
   }
   return (0);
}


double Fcodon_3x4 (double fcodon[], double fb3x4[]);
void OutFb3x4 (FILE*fout, double fb3x4[]);
void CountCodons (FILE *fout,double fcodonsg[],double fb3x4sg[],double fb4g[],int *miss);

double Fcodon_3x4(double fcodon[], double fb3x4[])
{
/* this converts the codon frequencies into a fb3x4 table. fcodon has 64 codons.
*/
   int b[3], k,j, nc=64, status=0;
   double t;

   zero(fb3x4,12);
   for(k=0; k<nc; k++) {
      b[0]=k/16; b[1]=(k%16)/4; b[2]=k%4;
      FOR(j,3)  fb3x4[j*4+b[j]]+=fcodon[k];
   }
   FOR(j,3)  {
      t=sum(fb3x4+j*4, 4);
      if(t<1e-15) status=-1;
      abyx(1/t, fb3x4+j*4, 4);
   }
   if(status)
      matout (F0, fcodon, 16, 4);
   return(status);
}

void OutFb3x4 (FILE*fout, double fb3x4[])
{
   int j,k;
   for(j=0; j<3; j++) {
      fprintf(fout, "\nposition %2d:", j+1);
      FOR(k,4) fprintf (fout,"%5c:%7.5f", BASEs[k],fb3x4[j*4+k]);
   }
}




int CodonListall(char codon[], int nb[3], int ib[3][4])
{
/* Changed from NucListall()
   Resolve a codon, possibly ambiguous into nb[k] nucleotides for 
   each codon position k.  
*/
   int k;

   FOR(k,3)  NucListall(codon[k], &nb[k], ib[k]);
   return(0);
}

void CountCodons (FILE *fout,double fcodonsg[],double fb3x4sg[],double fb4g[],int *miss)
{
/* Outputs codon counts and f3x4 tables, called from InitializeCodon(), where 
   more notes are found.
*/
   int h, j,k, nc=NCODE, ig, wname=15, nb[3],ib[3][4],ic;
   char codon[4]=" ";

   /* counts codons for output, species first, genes next */
   fputs("Codon usage in sequences\n",fout);
   zero(fcodonsg, com.ns*nc);
   FOR(j,com.ns) {
      for(h=0; h<com.npatt; h++) {
         if(com.cleandata) 
            getcodon(codon,FROM61[com.z[j][h]]);
         else 
            FOR(k,3) codon[k]=com.z[j][h*3+k];
         CodonListall(codon, nb, ib);

         k=nb[0]*nb[1]*nb[2];
         if(k>1)  { *miss=1; continue; }
         ic=ib[0][0]*16+ib[1][0]*4+ib[2][0];
         fcodonsg[j*nc+ic]+=com.fpatt[h];
      }
      Fcodon_3x4(fcodonsg+j*nc, fb3x4sg+j*12);
   }
   printcums(fout, com.ns, fcodonsg, com.icode);
   fputs("Codon position x base (3x4) table for each sequence.",fout);
   FOR (j,com.ns) {
      fprintf (fout,"\n\n#%d: %-*s", j+1,wname,com.spname[j]);
      OutFb3x4 (fout, fb3x4sg+j*12);
   }

   zero(fcodonsg,(com.ngene+1)*nc);  zero(fb4g,(com.ngene+1)*4);
   FOR(ig, com.ngene) {
      FOR(j, com.ns) {
         for(h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            if(com.cleandata) 
               getcodon(codon,FROM61[com.z[j][h]]);
            else 
               FOR(k,3) codon[k]=com.z[j][h*3+k];
            CodonListall(codon, nb, ib);
            k=nb[0]*nb[1]*nb[2];
            if(k>1) continue;
            ic=ib[0][0]*16+ib[1][0]*4+ib[2][0];
            fcodonsg[ig*nc+ic]+=com.fpatt[h];
         }
      }
      if(Fcodon_3x4(fcodonsg+ig*nc, fb3x4sg+ig*12)) {
         printf("Gene %d seems empty.", ig+1);
         exit(-1);
      }
   }
   if(com.ngene>1) {
      fputs("\n\nCodon usage in genes\n",fout);
      printcums(fout, com.ngene, fcodonsg, com.icode);
      fputs("Codon position x base (3x4) table for each gene.\n",fout);
      FOR (ig,com.ngene) {
         fprintf (fout,"\n\nGene #%d", ig+1);
         OutFb3x4 (fout, fb3x4sg+ig*12);
      }
   }
   
   FOR(ig,com.ngene)  FOR(k,nc) fcodonsg[com.ngene*nc+k]+=fcodonsg[ig*nc+k];
   Fcodon_3x4(fcodonsg+com.ngene*nc, fb3x4sg+com.ngene*12);
   FOR(ig,com.ngene+1)
      FOR(j,3) FOR(k,4) fb4g[ig*4+k]+=fb3x4sg[ig*12+j*4+k]/3;
   
   fputs("\n\nSums of codon usage counts",fout);
   printcu(fout, fcodonsg+com.ngene*nc, com.icode);
   if(*miss) fputs("\n(Ambiguity data are not used in the counts.)\n",fout);
   fputs("\n\nCodon position x base (3x4) table, overall\n",fout);
   OutFb3x4 (fout, fb3x4sg+com.ngene*12);
}


void AddCodonFreqSeqGene (int js, int ig, double fcodon0[], double fcodon[],
                    double fb3x40[], double fb3x4[], 
                    double fb40[], double fb4[]);

void AddCodonFreqSeqGene (int js, int ig, double fcodon0[], double fcodon[],
                    double fb3x40[], double fb3x4[], 
                    double fb40[], double fb4[])
{
/* This adds codon and nucleotide counts in sequence js in gene ig to fcodon,
   fb3x4, and fb4, using fcodon0, fb3x40, and fb40 to resolve ambiguities
   Similar to AddFreqSeqGene().
*/
   int h, k, i0,i1,i2, nc=NCODE;
   int nb[3],ib[3][4],ic=-1;
   double t,t1;
   char str[4]="   ", codon[4]=" ", ft[64];

   for(h=com.posG[ig]; h<com.posG[ig+1]; h++) {
      if(com.cleandata) 
         getcodon(codon,FROM61[com.z[js][h]]);
      else 
         FOR(k,3) codon[k]=com.z[js][h*3+k];
      CodonListall(codon, nb, ib);
      k=nb[0]*nb[1]*nb[2];
      FOR(k,3) {  /* f3x4 & f1x4, no regard for stop codons */
         for(i0=0,t=t1=0; i0<nb[k]; i0++) {
            t+=fb3x40[k*4+ib[k][i0]];
            t1+=fb40[ib[k][i0]];
         }
         FOR(i0, nb[k]) {
            fb3x4[k*4+ib[k][i0]] += 
               com.fpatt[h]* fb3x40[k*4+ib[k][i0]]/t;
               fb4[ib[k][i0]] += com.fpatt[h]* fb40[ib[k][i0]]/t1;
         }
      }
      FOR(i0,64) ft[i0]=0;
      for(i0=k=0,t=0; i0<nb[0]; i0++) FOR(i1,nb[1]) FOR(i2,nb[2]) {
         ic=ib[0][i0]*16+ib[1][i1]*4+ib[2][i2];         
         if(FROM64[ic]==-1) continue;
         ft[ic]=1;  k++;
         t+=fcodon0[ic];
      }
      if(k==0) printf("%s in seq. %d is stop (icode=%d)\n", 
         getcodon(str,ic),js+1,com.icode);
      if(t<1e-20)
         puts("difficulty in resolving codon..");
      FOR(ic,nc)  if(ft[ic]) 
         fcodon[ic] += (t>0 ? com.fpatt[h]*fcodon0[ic]/t : com.fpatt[h]/k);
   }
}


int InitializeCodon (FILE *fout, double space[])
{
/* Count codons for genes, calculate site patterns and fpatt.
   Sequences com.z[] are not coded and may contain ambiguity characters
   Space requirement for fcodonsg & fb3x4sg: max(ngene+1,ns)*(64+12+4).
   First we count codons for output, with ambiguity characters ignored.    
   Then we recount to resolve ambiguity characters, to be used for ML 
   calculation later on.
   set up com.pi[NCODE],com.piG[NGENE][64], according to com.codonf
   com.pi[] has freqs for all codon sites in the seqs if ngene>1.
   Space use is not economical as com.piG and fcodonsg do not overlap.
*/
   int j,k, nc=NCODE, ig, miss=0;
   int irf,nrf=20;
   double *fcodonsg=space, *fb3x4sg=space+max2((com.ngene+1),com.ns)*nc;
   double *fb4g=space+(com.ngene+1)*(64+12);
   double *ppi, fcodon0[64],fb3x40[12],fb40[4], d1,d2,d3;

   if(!com.readpattf) PatternWeight(fout);

   /* counts codons for output, species first, genes next */
   if(noisy) puts("Counting codons..");
   CountCodons(fout, fcodonsg, fb3x4sg, fb4g, &miss);

   /* Now to count fcodonsg, fb3x4sg, fb4g, to set up pi's for ML calculation.
      Three iterations are going on at the same time.
   */
   if (com.codonf!=Fequal && miss) { /* iteration to resolve ambiguities */
      FOR(ig, com.ngene) {    /* calculate com.piG[] */
         axtoy(1/sum(fcodonsg+ig*nc,nc), fcodonsg+ig*nc, fcodon0, nc);
         xtoy(fb3x4sg+ig*12,fb3x40, 12);
         xtoy(fb4g+ig*4,fb40, 4);

         for(irf=0; irf<nrf; irf++) {
            zero(fcodonsg+ig*nc,nc); zero(fb3x4sg+ig*12,12); zero(fb4g+ig*4,4);
            FOR(j,com.ns) {
               AddCodonFreqSeqGene (j, ig, fcodon0, fcodonsg+ig*nc, 
                  fb3x40, fb3x4sg+ig*12, fb40, fb4g+ig*4);
            }
            abyx(1/sum(fcodonsg+ig*nc,nc), fcodonsg+ig*nc, nc);
            FOR(k,3) abyx(1/sum(fb3x4sg+ig*12+k*4,4), fb3x4sg+ig*12+k*4, 4);
            abyx(1/sum(fb4g+ig*4,4), fb4g+ig*4, 4);
            d1=distance(fcodonsg+ig*nc, fcodon0, nc);
            d2=distance(fb3x4sg+ig*12, fb3x40, 12);
            d3=distance(fb4g+ig*4,fb40,4);
            if(d1<1e-8 && d2<1e-8 && d3<1e-8) 
               break;
            xtoy(fcodonsg+ig*nc, fcodon0, nc);
            xtoy(fb3x4sg+ig*12, fb3x40, 12);
            xtoy(fb4g+ig*4, fb40, 4);
         } /* for(irf) */
      }   /* for(ig) */

      axtoy(1/sum(fcodonsg+com.ngene*nc,nc),fcodonsg+com.ngene*nc,fcodon0,nc);
      xtoy(fb3x4sg+com.ngene*12,fb3x40, 12);
      xtoy(fb4g+com.ngene*4,fb40, 4);
      for(irf=0; irf<nrf; irf++) {  /* calculate com.pi[] */
         zero(fcodonsg+com.ngene*nc,nc); zero(fb3x4sg+com.ngene*12,12);
         zero(fb4g+com.ngene*4,4);
         FOR(ig, com.ngene)
            FOR(j,com.ns) {
               AddCodonFreqSeqGene(j, ig, fcodon0, fcodonsg+com.ngene*nc, 
                  fb3x40, fb3x4sg+com.ngene*12, fb40, fb4g+com.ngene*4);
            }
         abyx(1/sum(fcodonsg+com.ngene*nc,nc), fcodonsg+com.ngene*nc, nc);
         FOR(k,3) 
            abyx(1/sum(fb3x4sg+com.ngene*12+k*4,4), fb3x4sg+com.ngene*12+k*4, 4);
         abyx(1/sum(fb4g+com.ngene*4,4), fb4g+com.ngene*4, 4);
         d1=distance(fcodonsg+com.ngene*nc, fcodon0, nc);
         d2=distance(fb3x4sg+com.ngene*12, fb3x40, 12);
         d3=distance(fb4g+com.ngene*4,fb40,4);
         if(d1<1e-8 && d2<1e-8 && d3<1e-8)  break;
         xtoy(fcodonsg+com.ngene*nc, fcodon0, nc);
         xtoy(fb3x4sg+com.ngene*12, fb3x40, 12);
         xtoy(fb4g+com.ngene*4, fb40, 4);
      } /* for(irf) */
   }

   /* edit com.pi & com.piG according to com.codonf */
   FOR(ig,com.ngene+1) {
      ppi = (ig<com.ngene?com.piG[ig]:com.pi);
      zero(ppi,nc);
      if (com.codonf==Fequal)
         fillxc(ppi,1,com.ncode);
      else if (com.codonf==Fcodon || com.codonf==FMutSelCodon) {
         FOR(k,nc)  if(FROM64[k]>-1)  ppi[FROM64[k]]=fcodonsg[ig*nc+k]; 
      }
      else if (com.codonf==F3x4 || com.codonf==F3x4MG || com.codonf==FMutSel3x4) {
         FOR(k,nc)  if(FROM64[k]>-1)
            ppi[FROM64[k]]=fb3x4sg[ig*12+k/16]*fb3x4sg[ig*12+4+(k/4)%4]*fb3x4sg[ig*12+8+k%4];
         if(ig<com.ngene && com.codonf==F3x4MG)
            xtoy(fb3x4sg+ig*12, com.f3x4MG[ig], 12);
      }
      else if (com.codonf==F1x4 || com.codonf==F1x4MG) {
         FOR(k,nc)  if(FROM64[k]>-1)
            ppi[FROM64[k]]=fb4g[ig*4+k/16]*fb4g[ig*4+(k/4)%4]*fb4g[ig*4+k%4];
         if(ig<com.ngene && com.codonf==F1x4MG) {
            xtoy(fb3x4sg+ig*12, com.f3x4MG[ig], 12);
            for(k=0;k<4;k++)  {
               for(j=0,d1=0;j<3;j++) d1+=com.f3x4MG[ig][j*4+k];
               for(j=0;j<3;j++) com.f3x4MG[ig][j*4+k]=d1/3;;
            }
         }
      }
      abyx(1/sum(ppi,com.ncode),ppi,com.ncode);  /* ncode != nc */
      if(com.codonf==F1x4MG||com.codonf==F3x4MG) com.pf3x4MG=com.f3x4MG[0];
   }
   getpi_sqrt (com.pi, com.pi_sqrt, com.ncode, &com.npi0);

   if(com.verbose && com.ngene==1 && !miss) {
      fputs("\n\nCodon frequencies under model, for use in evolver:\n",fout); 
      FOR(k,64) {
        fprintf(fout,"%12.8f",GeneticCode[com.icode][k]==-1?0:com.pi[FROM64[k]]);
        if((k+1)%4==0) FPN(fout);
      }
   }

   return(0);
}



int AA2Codonf(double faa[20], double fcodon[])
{
/* get codon freqs from amino acid freqs, assuming equal freq. for each syn
   codon.  Used in codon-based amino acid substitution models.
*/
   int ic, iaa, i, NCsyn[20];

   FOR(i,20) NCsyn[i]=0;
   FOR(ic,64) if((iaa=GeneticCode[com.icode][ic])!=-1) NCsyn[iaa]++;
   zero(fcodon, 64);
   for(ic=0; ic<Nsensecodon; ic++) {
      iaa=GeneticCode[com.icode][FROM61[ic]];
      fcodon[ic]+=faa[iaa]/NCsyn[iaa];
   }
   if(fabs(1-sum(fcodon,64))>1e-6) printf("\n1 == %12.7f\n", sum(fcodon,64));

   return (0);
}


int DistanceMatAA (FILE *fout)
{
   int i,j, h;
   double p, lst;

   if(fout) fprintf(fout,"\nAA distances (raw proportions of different sites)\n");
   for(h=0,lst=0; h<com.npatt; h++)  lst+=com.fpatt[h];
   FOR(i, com.ns) {
      if(fout) fprintf(fout, "\n%-15s", com.spname[i]);
      FOR(j,i) {
         for(h=0,p=0; h<com.npatt; h++)  
            if (com.z[i][h]!=com.z[j][h]) p+=com.fpatt[h];
         p /= lst;
         SeqDistance[i*(i-1)/2+j]=p;
         if(fout) fprintf(fout, " %7.4f", p);
      }
   }
   if (fout) FPN(fout);
   return (0);
}

int DistanceMatNG86 (FILE *fout, FILE*fds, FILE*fdn, FILE*ft, double alpha)
{
/* Estimation of dS and dN by the method of Nei & Gojobori (1986)
   This works with both coded (com.cleandata==1) and uncoded data.
   In the latter case (com.cleandata==0), the method does pairwise delection.

   alpha for gamma rates is used for dN only.
*/
   char codon[2][4]={"   ", "   "};
   int is,js, k,i0,h, wname=20, status=0, ndiff,nsd[4];
   int nb[3],ib[3][4], missing;
   double ns,na, nst,nat, S,N, St,Nt, dS,dN,dN_dS,y, bigD=3, lst;
   double SEds, SEdn, p;

   fputs("\n\n\nNei & Gojobori 1986. dN/dS (dN, dS)",fout);
   if(com.cleandata==0) fputs("\n(Pairwise deletion)",fout);

   fputs("\n(Note: This matrix is not used in later m.l. analysis.\n",fout);
   fputs("Use runmode = -2 for ML pairwise comparison.)\n",fout);

   fputs("\nNumber of codon sites with 0,1,2,3 position differences\n",frst);

   fprintf(fds,"%6d\n",com.ns);  fprintf(fdn,"%6d\n",com.ns); 
   fprintf(ft,"%6d\n",com.ns);
   if(noisy>1 && com.ns>10)  puts("NG distances for seqs.:");
   FOR (is,com.ns) {
      fprintf(fout,"\n%-*s", wname,com.spname[is]);
      fprintf(fds,   "%-*s ",wname,com.spname[is]);
      fprintf(fdn,   "%-*s ",wname,com.spname[is]);
      fprintf(ft,    "%-*s ",wname,com.spname[is]);
      FOR (js,is) {
         FOR(k,4) nsd[k]=0;
         for (h=0,lst=0,nst=nat=S=N=0; h<com.npatt; h++)  {
            if(com.cleandata)
               FOR(i0,2) getcodon(codon[i0],FROM61[com.z[i0==0?is:js][h]]);
            else {
               FOR(i0,2) FOR(k,3) codon[i0][k]=com.z[i0==0?is:js][h*3+k];
               for(i0=0,missing=0;i0<2;i0++) {
                  CodonListall(codon[i0], nb, ib);
                  if(nb[0]*nb[1]*nb[2]!=1)  { missing=1; break; }
               }
               if(missing) continue;
            }
            lst+=com.fpatt[h];
            ndiff=difcodonNG(codon[0],codon[1],&St,&Nt,&ns,&na,0,com.icode);
            nsd[ndiff]+=(int)com.fpatt[h];
            S+=St*com.fpatt[h];
            N+=Nt*com.fpatt[h];
            nst+=ns*com.fpatt[h];
            nat+=na*com.fpatt[h];
         }  /* for(h) */
         y=lst*3./(S+N); S*=y; N*=y;       /* rescale for stop codons */

         if(noisy>=9)
           printf("\n%3d%3d:Sites %7.1f +%7.1f =%7.1f\tDiffs %7.1f +%7.1f =%7.1f",
             is+1,js+1,S,N,S+N,nst,nat, nst+nat);

         fprintf (frst, "%4d vs. %4d %6d %5d %5d %5d  ", 
            is+1,js+1,nsd[0],nsd[1],nsd[2],nsd[3]);

         dS=1-4./3*nst/S;  dN=1-4./3*nat/N;
         if(noisy>=9 && (dS<=0||dN<=0)) { puts("\nNG86 unusable."); status=-1;}
         if(dS==1) dS=0;
         else      dS=(dS<0?-1:3./4*(-log(dS)));
         if(dN==1) dN=0;
         else dN=(dN<0?-1:3./4*(alpha==0?-log(dN):alpha*(pow(dN,-1/alpha)-1)));

         dN_dS=(dS>0?dN/dS:-1);
         fprintf(fout,"%7.4f (%5.4f %5.4f)", dN_dS, dN, dS);
         fprintf(frst,"%7.4f (%6.4f %6.4f)\n", dN_dS, dN, dS);

         if(dN<0) dN=bigD; if(dS<0) dS=bigD;
         SeqDistance[is*(is-1)/2+js] = (S*dS+N*dN)*3/(S+N);
         fprintf(fds," %7.4f", dS);   fprintf(fdn," %7.4f", dN);
         fprintf(ft," %7.4f", (S*dS+N*dN)*3/(S+N));
         if(alpha==0 && dS<bigD) { p=nst/S; SEds=sqrt(9*p*(1-p)/(square(3-4*p)*S)); }
         if(alpha==0 && dN<bigD) { p=nat/N; SEdn=sqrt(9*p*(1-p)/(square(3-4*p)*N)); }
/*
 fprintf(frst1," NG %7.1f %7.1f%7.4f +- %6.4f %6.4f +- %6.4f %6.4f\n", 
    S,N, dN,SEdn, dS,SEds, dN_dS);
*/
      }
      FPN(fds); FPN(fdn); FPN(ft);
      if(noisy>1 && com.ns>10)  printf(" %3d",is+1);
   }    /* for(is) */
   FPN(F0); FPN(fout);
   if(status) fprintf (fout, "NOTE: -1 means that NG86 is inapplicable.\n");
   /* fputs("\n",frst); fflush (frst); */
   return (0);
}


int GetDaa (FILE* fout, double daa[])
{
/* Get the amino acid distance (or substitution rate) matrix 
   (grantham, dayhoff, jones, etc).
*/
   FILE * fdaa;
   char aa3[4]="";
   int i,j, naa=20;
   double dmax=0,dmin=1e40;

   fdaa=gfopen(com.daafile, "r");
   if(noisy>3) printf("\n\nReading matrix from %s..\n", com.daafile);
   if (com.model==REVaa_0||com.model==REVaa) puts("To get initial values.");

   for (i=0; i<naa; i++)  for (j=0,daa[i*naa+i]=0; j<i; j++)  {
      fscanf(fdaa, "%lf", &daa[i*naa+j]);
      daa[j*naa+i]=daa[i*naa+j];
      if (dmax<daa[i*naa+j]) dmax=daa[i*naa+j];
      if (dmin>daa[i*naa+j]) dmin=daa[i*naa+j];
   }
   if(com.aaDist && (com.seqtype==1||com.model==FromCodon)) { /* codon model */
      if(noisy) printf("\ndistance: %.2f --- %.2f\n", dmin, dmax);
      FOR (i,naa) FOR(j,naa) com.daa[i*naa+j]/=dmax;
   }
   else if (com.seqtype==AAseq && com.model==Empirical) {
      FOR(i,naa) if(fscanf(fdaa,"%lf",&com.pi[i])!=1) 
         error2("aaRatefile");
      if (fabs(1-sum(com.pi,20))>1e-6) {
         printf("\nSum of freq. = %.6f != 1 in aaRateFile\n",sum(com.pi,naa)); 
         exit(-1);
      }
      getpi_sqrt (com.pi, com.pi_sqrt, naa, &com.npi0);
   }
   fclose (fdaa);

   if(fout) {
      fprintf (fout, "\n%s\n", com.daafile);
      FOR (i,naa) {
         fprintf (fout, "\n%4s", getAAstr(aa3,i));
         FOR (j,i)  fprintf (fout, "%5.0f", daa[i*naa+j]); 
      }
      FPN (fout);
   }

/*
SetAA1STEP();
for(i=0,FPN(frst);i<naa;i++,FPN(frst))
   FOR(j,i) fprintf(frst,"%3d",AA1STEP[i*(i-1)/2+j]);

for(i=0,k=0;i<naa;i++) 
   FOR(j,i) if(AA1STEP[i*(i-1)/2+j]) {
      fprintf(frst,"%c%c\t%.2f\n",AAs[i],AAs[j],com.daa[i*naa+j]);
      k++;
   }
fprintf(frst,"\n%d one-step amino acid pairs\n", k);
exit (0);
*/

   return (0);
}


int SetAA1STEP (void)
{
/* Sets the global variable AA1STEP[19*20/2].
   Sets com.nrate for models like AAClasses and REVaa_0.
   AA1STEP[k] marks the k_th pair of amino acids that differ at one position, 
   Q[i*naa+j] is the k_th nonzero element if AA1STEP[k]=i*naa+j;
   Lower diagonal of Q is visited, with i>j.
*/
   int nc=Nsensecodon, naa=20, i,j,k, ic1,ic2, ndiff, from[3],to[3];
   int *Q=(int*)PMat;

   FOR (i, naa*naa) Q[i]=0;
   for (i=0; i<nc; i++) FOR (j,i) {
      ic1=FROM61[i]; from[0]=ic1/16; from[1]=(ic1/4)%4; from[2]=ic1%4;
      ic2=FROM61[j];   to[0]=ic2/16;   to[1]=(ic2/4)%4;   to[2]=ic2%4;
      for (k=0,ndiff=0; k<3; k++)  if (from[k]!=to[k]) ndiff++; 
      if (ndiff!=1)  continue; 
      ic1=GeneticCode[com.icode][ic1];
      ic2=GeneticCode[com.icode][ic2];
      Q[ic1*naa+ic2]++; 
      Q[ic2*naa+ic1]++;
   }
/*
#if DEBUG
      for (i=0,FPN(F0); i<naa; i++,FPN(F0)) FOR(j,i)printf("%3d",Q[i*naa+j]);
#endif
*/
   for (i=0,k=0; i<naa; i++) FOR(j,i) {
#if DEBUG
      if (Q[i*naa+j]!=Q[j*naa+i]) error2("strange in SetAA1STEP");
#endif
      if (Q[i*naa+j]>0) { AA1STEP[i*(i-1)/2+j]=1;  k++; }
      else                AA1STEP[i*(i-1)/2+j]=0;
   }
#if DEBUG
   for(i=0,FPN(F0);i<naa;i++,FPN(F0)) 
      FOR(j,i)printf("%3d",AA1STEP[i*(i-1)/2+j]);
#endif

   return(0);
}

int GetOmegaAA (int OmegaAA[])
{
/* This routine reads the file OmegaAA.dat to initialize the
   lower diagonal matrix OmegaAA, which specifies the aa substituion
   rate classes.  To be used with the codon substitution model
   AAClasses, which specifies several classes of the dN/dS ratio.

   OmegaAA[iaa*(iaa-1)/2+jaa]= -1 if no one-step change is possible; 
                             = 0 for the first, background, class
                             = i (1,..,nclass) if iaa and jaa are in class i 
*/
   char *OmegaAAf="OmegaAA.dat", line[1024];
   FILE *fin;
   int iomega, n1step=0, i,j,k, iaa,jaa, npair, naa=20, nline=1024;

   for(i=0,n1step=0; i<naa; i++) FOR(j,i)
      if (AA1STEP[i*(i-1)/2+j]) { OmegaAA[i*(i-1)/2+j]=0; n1step++; }
      else                        OmegaAA[i*(i-1)/2+j]=-1;
   if (noisy) {
       printf("\n\n%d one-step aa pairs.\n", n1step);
       printf("Reading omega class from %s.\n", OmegaAAf);
   }
   fin=gfopen(OmegaAAf,"r");
   fscanf(fin, "%d", &com.nOmegaType);  
   printf ("# of dN/dS classes requested: %d.\n", com.nOmegaType);
   if (com.nOmegaType<1 || com.nOmegaType>65-1) { 
      if (com.seqtype!=CODONseq) puts("\nTo be tested.\a");
      com.nOmegaType=0;
      if (com.seqtype==AAseq) {
         FOR(i,naa) FOR(j,i) if(i*naa+j!=ijAAref && AA1STEP[i*(i-1)/2+j])
             OmegaAA[i*(i-1)/2+j]=com.nOmegaType++;
      }
      else
         FOR(i,naa) FOR(j,i) 
           if(AA1STEP[i*(i-1)/2+j]) OmegaAA[i*(i-1)/2+j]=com.nOmegaType++;
      printf("%d dN/dS ratios estimated from data.\nCrtl-C if wrong. Enter to continue.",com.nOmegaType);
      getchar ();
   }
   else {
      FOR (iomega, com.nOmegaType-1) {
         fscanf(fin, "%d", &j);
         if (j!=iomega+1) { printf("err data file %s.", OmegaAAf); exit(-1); } 
         printf ("\nClass #%d: ", j);
         j=fgetc (fin);  if (j!=':') error2("err expecting :");
         fgets (line, nline, fin);
      
         printf ("%s\n", line);
         for (j=0,npair=0; j<nline-1&&line[j]&&line[j]!='\n'; j++) {
            iaa=line[j];
            if (!isalpha(iaa)) continue;
            jaa=line[++j];  if(!isalpha(jaa)) error2("err jaa");
            npair++;

            printf ("\npair %2d: |%c%c| ", npair, iaa,jaa);
            iaa=CodeChara((char)iaa,AAseq); jaa=CodeChara((char)jaa,AAseq);
            if(iaa<0||iaa>19||jaa<0||jaa>19) error2("aa not found");
            if (iaa<jaa)  { k=jaa, jaa=iaa; iaa=k; }
      
            printf ("|%c%c (%2d,%2d)| ", AAs[iaa], AAs[jaa],iaa,jaa);
            if (iaa==jaa) puts("This pair has no effect.");
            if (OmegaAA[iaa*(iaa-1)/2+jaa]==-1) {
               puts("\nThis pair cannot change in one step and is ignored!");
               continue;
            }
            else if (OmegaAA[iaa*(iaa-1)/2+jaa]) 
               error2("This pair has already been specified?");
            OmegaAA[iaa*(iaa-1)/2+jaa]=iomega+1;
            printf (" in class %d ",iomega+1);
         }
      }
   }
   fclose(fin);
   com.nrate=com.nkappa=(com.hkyREV?5:!com.fix_kappa);
   com.nrate += (com.nOmega=com.nOmegaType);
   printf ("\n%d AA substitution types (omega's)\n", com.nOmegaType);
/*
   for (i=0; i<naa; i++,FPN(F0)) 
       FOR(j,i) printf ("%3d",OmegaAA[i*(i-1)/2+j]);
*/
   return (0);
}



int GetCodonFreqs(double pi[])
{
/* Recalcualte the expected codon frequencies (com.pi[]) using the control
   variable com.codonf.
*/
   int n=com.ncode, i,j, ic,b[3];
   double fb3x4[12], fb4[4], GC[3]={0};

   if (com.codonf==Fequal) { fillxc(pi,1./n,n); return 0; }
   if (com.codonf!=Fcodon && com.codonf!=FMutSelCodon) {
      for (i=0,zero(fb3x4,12),zero(fb4,4); i<n; i++) {
         ic=FROM61[i];  b[0]=ic/16; b[1]=(ic/4)%4; b[2]=ic%4;
         FOR(j,3) { fb3x4[j*4+b[j]]+=pi[i];  fb4[b[j]]+=pi[i]/3.; }
      }
      for (i=0; i<n; i++) {
         ic=FROM61[i];  b[0]=ic/16; b[1]=(ic/4)%4; b[2]=ic%4;
         if (com.codonf==F3x4 || com.codonf==F3x4MG || com.codonf==FMutSel3x4)
            pi[i]=fb3x4[b[0]]*fb3x4[4+b[1]]*fb3x4[8+b[2]];
         else
            pi[i]=fb4[b[0]]*fb4[b[1]]*fb4[b[2]];
      }

      if(com.codonf==F1x4MG)
         FOR(j,3) xtoy(fb4, com.pf3x4MG+j*4, 4);
      else if(com.codonf==F3x4MG || com.codonf==FMutSel3x4)
         xtoy(fb3x4, com.pf3x4MG, 12);

      abyx (1./sum(pi,n), pi, n);

      GC[0]=(fb3x4[0+1]+fb3x4[0+3])*100;
      GC[1]=(fb3x4[4+1]+fb3x4[4+3])*100;
      GC[2]=(fb3x4[8+1]+fb3x4[8+3])*100;

fprintf(frst1, "\tGC123\t%.1f\t%.1f\t%.1f", GC[0],GC[1],GC[2]);

   }
   getpi_sqrt (com.pi, com.pi_sqrt, n, &com.npi0);

   return 0;
}


double lfun2dSdN (double x[], int np)
{
/* likelihood function for calculating dS and dN between 2 sequences,
   com.z[0] & com.z[1:
         prob(i,j) = PI_i * p(i,j,t)
   
   Data are clean and coded.
   Transition probability pijt is calculated for observed patterns only.
*/
   int n=com.ncode, h,i,k, ik;
   double lnL=0, fh,expt[NCODE];
   double *pkappa=com.KAPPA;

   NFunCall++;
   k=1, ik=0;
   if(com.hkyREV==0) {
      if(com.fix_kappa==1) { pkappa[0]=com.kappa; ik=1; }
      else                 com.kappa=x[k]; /* Is this necessary? */
   }
   for(i=0; i<(com.hkyREV?5:!com.fix_kappa); i++) pkappa[ik++]=x[k++];
   if(com.codonf>=FMutSel3x4) 
      for(i=0; i<3; i++) pkappa[ik++]=x[k++];

   if(!com.fix_omega) com.omega=x[1+com.nkappa];
   if(!com.fix_kappa || !com.fix_omega)
      EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, pkappa, com.omega,PMat);

   FOR(k,n) expt[k]=exp(x[0]*Root[k]);
   for (h=0; h<com.npatt; h++) {
      if(com.fpatt[h]<1e-20) continue;
      for(k=0,fh=0;k<n;k++) fh+=U[com.z[0][h]*n+k]*expt[k]*V[k*n+com.z[1][h]];
      fh*=com.pi[com.z[0][h]];
      if(fh<=0) {
         matout(F0,x,1,np); printf("lfun2dSdN: fh = %.9f\n",fh);
         fh=1e-70;
      }
      lnL-=log(fh)*com.fpatt[h];
   }
   return (lnL);
}


int VariancedSdN(double t, double omega, double vtw[2*2], double vdSdN[2*2])
{
/* This calculates the covariance matrix of dS & dN, using the 
   difference approximation, from the covariance matrix of t and 
   omega (vtw).  com.kappa and com.pi are used.  Sampling errors
   in parameters other than t and omega, such as kappa and pi[], 
   are ignored.
         JacobiSN = {{dS/dt, dS/dw}, {dN/dt,dN/dw}}
*/
   int np=2;
   double JacobiSN[2*2],T1[2*3],T2[2*3], S,dS,dN, dS1,dN1,dS2,dN2, eh;
   double *pkappa=com.KAPPA;

   if(vtw[0]<=0 || vtw[3]<=0)
      { puts("var(dS,dN) not calculable."); zero(vdSdN,4); return(-1); }

   /* printf("\nt & w: %.5f %.5f\n", t, omega);
      matout(F0,vtw, 2,2); */
   EigenQc(1,t,&S,&dS,&dN,NULL,NULL,NULL, pkappa,omega,PMat);

   eh=(t+1)*Small_Diff;
   EigenQc(1,t+eh,&S,&dS1,&dN1,NULL,NULL,NULL, pkappa,omega,PMat);
   EigenQc(1,t-eh,&S,&dS2,&dN2,NULL,NULL,NULL, pkappa,omega,PMat);
   JacobiSN[0*np+0]=(dS1-dS2)/(2*eh);
   JacobiSN[1*np+0]=(dN1-dN2)/(2*eh);
  
   eh=(omega+1)*Small_Diff;
   EigenQc(1,t,&S,&dS1,&dN1,NULL,NULL,NULL, pkappa,omega+eh,PMat);
   EigenQc(1,t,&S,&dS2,&dN2,NULL,NULL,NULL, pkappa,omega-eh,PMat);
   JacobiSN[0*np+1]=(dS1-dS2)/(2*eh);
   JacobiSN[1*np+1]=(dN1-dN2)/(2*eh);
  
   matby(JacobiSN,vtw,T1,2,2,2);
   mattransp2 (JacobiSN, T2, 2, 2);
   matby(T1,T2,vdSdN,2,2,2);

   /* matout(F0,vdSdN, 2,2); */

   return (0);
}

double distanceHKY85 (double x[], double *kappa, double alpha);
int distance3pos(double dHKY[], double kHKY[], int *sites4, char *z1, char *z2);

int distance3pos(double dHKY[], double kHKY[], int *sites4, char *z1, char *z2)
{
/* This calculates nucleotide-based distances between two protein-coding 
   DNA sequences z1 and z2, both of which are coded.  com.cleandata = 1 is 
   assumed.
*/
   int i,j, h, k, ic1, ic2, from[3], to[3];
   double fij[4][16]={{0}}, pi4[4]={0};
   /* [0,1,2] are for 3 positions, [3] is for 4-fold */

   for (h=0; h<com.npatt; h++) {
      ic1=FROM61[z1[h]]; from[0]=ic1/16; from[1]=(ic1/4)%4; from[2]=ic1%4;
      ic2=FROM61[z2[h]];   to[0]=ic2/16;   to[1]=(ic2/4)%4;   to[2]=ic2%4;
      for(k=0; k<3; k++) 
         fij[k][from[k]*4+to[k]] += com.fpatt[h]/com.ls;
      if(from[0]==to[0] && from[1]==to[1] && FourFold[to[0]][to[1]]) 
         fij[3][from[2]*4+to[2]] += com.fpatt[h];
   }
   *sites4 = (int) sum(fij[3], 16);
   FOR(k,16) fij[3][k] /= *sites4;
   FOR(i,4) FOR(j,4) pi4[i]+=fij[3][i*4+j]/2;
   FOR(i,4) FOR(j,4) pi4[j]+=fij[3][i*4+j]/2;

   for(k=0; k<4; k++)
      dHKY[k]=distanceHKY85 (fij[k], &kHKY[k], 0);

   return(0);
}




int PairwiseCodon (FILE *fout, FILE*fds, FILE*fdn, FILE*ft, double space[])
{
/* Calculates ds & dn for all pairwise codon sequence comparisons.
   It uses different npatt for different pairs.
   The data com.z[] should be encoded clean data, with ambiguity characters 
   removed.  Think of what to do with raw unclean data.
   JacobiSN has two columns, the 1st are deratives of dS (dS/dt, dS/dk, dS/dw)
   and the second of dN.
*/
   char *pz0[NS],codon[2][3];   /* pz0, npatt0, & fpatt0 hold the old information */
   int npatt0=com.npatt;
   double *fpatt0, ls0=com.ls;
   float fp[NCODE*NCODE];
   int n=com.ncode, is,js,j,k,h, i0,np, wname=15;
   int nb[3],ib[3][4],ic[2], missing=0, sites4;
   double x[10]={.9,1,.5,.5,.5,.5,.3}, xb[10][2], lnL, e=1e-5, *var=space+NP;
   double S,dS,dN;
   double JacobiSN[2*3],T1[2*3],T2[2*3],vSN[2*2], dS1,dN1,dS2,dN2,y[3],eh; 
          /* for calculating SEs of dS & dN */
   double dHKY[4], kHKY[4];
   double tb[2]={1e-5,39}, kappab[2]={.01,999}, omegab[2]={.001,99};
   double *pkappa=com.KAPPA;

   fpatt0=(double*)malloc(npatt0*3*sizeof(double));
   FOR(k,com.ns) pz0[k]=com.z[k];
   com.z[0]=(char*)(fpatt0+npatt0);  com.z[1]=com.z[0]+npatt0;
   FOR (k,npatt0) fpatt0[k]=(float)com.fpatt[k];

   if(!com.cleandata) puts("\nPairwiseCodon: pairwise deletion.");
   if (com.ngene>1 && com.Mgene==1) puts("ngene>1 to be tested.");
   if (noisy) printf("pairwise comparison (Goldman & Yang 1994).\n");
   fprintf(fout,"\npairwise comparison, codon frequencies: %s.\n",
      codonfreqs[com.codonf]);

   FOR(j,com.nkappa) { xb[1+j][0]=kappab[0]; xb[1+j][1]=kappab[1]; }
   if(!com.fix_omega)  { k=1+com.nkappa; xb[k][0]=omegab[0]; xb[k][1]=omegab[1]; }

   fprintf(fds,"%6d\n", com.ns);  fprintf(fdn,"%6d\n", com.ns);
   fprintf(ft,"%6d\n", com.ns);
   fprintf(frst, "\n\npairwise comparison (Goldman & Yang 1994)");
   fprintf(frst,
      "\nseq seq        N       S       dN       dS     dN/dS   Paras.\n");

   for(is=0;is<com.ns;is++) {
      fprintf(fds,"%-*s ", wname,com.spname[is]);
      fprintf(fdn,"%-*s ", wname,com.spname[is]);
      fprintf(ft,"%-*s ", wname,com.spname[is]);
      for(js=0; js<is; js++) {
         if(noisy>9) {
            puts("\ni & j (i>j)? "); scanf("%d%d",&is,&js);  
            is--; js--;
            if(is>com.ns || js<0 || is<js) error2("invalid pair");
         }
         printf ("\n\n%4d vs. %3d", is+1, js+1);
         fprintf(fout,"\n\n%d (%s) ... %d (%s)",
              is+1,com.spname[is], js+1,com.spname[js]);
         fprintf (frst, "%3d %3d ", is+1, js+1);
         if(noisy>2) fprintf(frub, "\n\n%d (%s) ... %d (%s)",
                  is+1,com.spname[is], js+1,com.spname[js]);
         FOR(k,n*n) fp[k]=0;
         if(com.cleandata) {
            for(h=0; h<npatt0; h++) {
               j=max2(pz0[is][h],pz0[js][h]); k=min2(pz0[is][h],pz0[js][h]);
               fp[j*n+k]+=(float)fpatt0[h];
            }
         }
         else {
            for(h=0,com.ls=0; h<npatt0; h++) {
               FOR(i0,2) FOR(k,3) codon[i0][k]=pz0[i0==0?is:js][h*3+k];
               for(i0=0,missing=0;i0<2;i0++) {
                  CodonListall(codon[i0], nb, ib);
                  if(nb[0]*nb[1]*nb[2]!=1)  { missing=1; break; }
                  else  ic[i0]=FROM64[ ib[0][0]*16+ib[1][0]*4+ib[2][0] ];
               }
               if(missing) continue;
               com.ls+=(int)fpatt0[h];

               j=max2(ic[0],ic[1]); k=min2(ic[0],ic[1]);
               fp[j*n+k]+=(float)fpatt0[h];
            }
         }

         for(j=0,com.npatt=0;j<n;j++) FOR(k,j+1)  if(fp[j*n+k]) {
            com.z[0][com.npatt]=(char)j;  com.z[1][com.npatt]=(char)k; 
            com.fpatt[com.npatt++]=fp[j*n+k];
         }
         if(noisy>2) printf("\n  npatt=%d ",com.npatt);
         for(j=0,zero(com.pi,n); j<com.npatt; j++) {
            com.pi[com.z[0][j]]+=com.fpatt[j]/(2.*com.ls);
            com.pi[com.z[1][j]]+=com.fpatt[j]/(2.*com.ls);

         }
         GetCodonFreqs(com.pi);

         distance3pos(dHKY, kHKY, &sites4, com.z[0], com.z[1]);

         np=com.np=(com.ntime=1)+com.nkappa+!com.fix_omega;  NFunCall=0;

         /* initial values and bounds */
         x[0]=SeqDistance[is*(is-1)/2+js]*(1+.2*rndu());  /* poor NG86 */
         if(x[0]<0.001)   x[0]=0.01*(0.5+rndu());
         else if(x[0]>3)  x[0]=1.5+rndu();
         xb[0][0]=max2(x[0]*0.1,tb[0]);
         xb[0][0]=min2(xb[0][0],0.1);
         xb[0][1]=min2(x[0]*5.0,tb[1]);
         if(com.nkappa==1) { /* HKY type model */
            if(is==0 && js==1)  x[1]=(com.icode==1?4:1.5)+rndu();
            else                x[1]=(x[1]*2+2+rndu())/3;
            if(x[1]>10) x[1]=5;
            xb[1][0]=0.4;
         }
         else         /* REV or FMutSel models, do something later */
            for(j=1,x[1]=.8+.4*rndu(); j<com.nkappa; j++) x[1+j]=.2+.4*rndu();

         if(!com.fix_omega) {
            k=1+com.nkappa;
            if(is==0 && js==0) x[k]=0.2+0.2*rndu();
            else               x[k]=(3*x[k]+0.6*rndu())/4;
            x[k]=max2(x[k],0.01);
            x[k]=min2(x[k],2);
         }

         if(noisy>=9) {
            FPN(F0);  FOR(k,np) printf(" %12.6f",x[k]); FPN(F0);
            FOR(k,np) printf(" %12.6f",xb[k][0]); FPN(F0);
            FOR(k,np) printf(" %12.6f",xb[k][1]); FPN(F0);
         }
         
         if(com.fix_kappa && com.fix_omega)  
            EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, pkappa,com.omega,PMat);
         if (np)
            ming2(noisy>3?frub:NULL,&lnL,lfun2dSdN,NULL,x,xb, space,e,np);
         else {  x[1]=x[2]=com.kappa=com.omega=0; lnL=0; }
         fprintf(fout,"\nlnL =%12.6f\n",-lnL);
         FOR(k,np) fprintf(fout," %8.5f",x[k]);  FPN(fout);

if(noisy>2) printf("\n\nt_NG = %.5f\tMLEs: %.5f  %.5f  %.5f\n",
               SeqDistance[is*(is-1)/2+js], x[0], x[1], x[np-1]);

         if (np&&com.getSE) {
            Hessian(np, x, lnL, space, var, lfun2dSdN, var+np*np);
            matinv(var, np, np, var+np*np);
            fprintf(fout,"SEs for parameters:\n");
            FOR(k,np) fprintf(fout," %8.5f",(var[k*np+k]>0.?sqrt(var[k*np+k]):-0));
            FPN(fout);
         }
         FPN(fout);
         EigenQc(1,x[0],&S,&dS,&dN, NULL,NULL,NULL, pkappa,com.omega,PMat);

         if(noisy>=3) {
            puts("\nNucleotide-based analysis (approximate MLEs under Q4x4):");
            printf("\ndHKY (123-4):");  FOR (k,4) printf(" %8.5f", dHKY[k]);
            printf("\nkHKY (123-4):");  FOR (k,4) printf(" %8.5f", kHKY[k]);
            printf(" (%d four-fold sites)\n", sites4);
         }
 
         fprintf(fds," %7.4f",dS);   fprintf(fdn," %7.4f",dN);
         fprintf(ft," %7.4f",x[0]);

         fprintf (fout,
             "t=%7.4f  S=%8.1f  N=%8.1f  dN/dS=%7.4f  dN=%7.4f  dS=%7.4f\n",
              x[0],S,com.ls*3-S,com.omega,dN,dS);

         fprintf(frst,"%8.1f %8.1f %8.4f %8.4f %8.4f",com.ls*3-S,S,dN,dS,com.omega);
         FOR(k,np) fprintf(frst," %8.4f",x[k]);  

fprintf(frst1,"\t%.3f", -lnL);

         k=np-1;
         if (com.getSE)
            fprintf(frst," +-%6.4f",(var[k*np+k]>0.?sqrt(var[k*np+k]):-1));
         fprintf(frst," %9.3f\n",-lnL);
         if(com.getSE && !com.fix_omega) {
            FOR(k,np) {
               FOR(j,np) y[j]=x[j];
               y[k] += (eh=(x[k]+1)*Small_Diff);
               if(!com.fix_kappa) com.kappa=y[1];
               com.omega=y[1+!com.fix_kappa];
               EigenQc(1,y[0],&S,&dS1,&dN1,NULL,NULL,NULL, pkappa,com.omega,PMat);
               y[k] -= 2*eh;
               if(!com.fix_kappa) com.kappa=y[1];
               com.omega=y[1+!com.fix_kappa];
               EigenQc(1,y[0],&S,&dS2,&dN2,NULL,NULL,NULL, pkappa, com.omega,PMat);

               JacobiSN[0*np+k]=(dS1-dS2)/(2*eh);
               JacobiSN[1*np+k]=(dN1-dN2)/(2*eh);
            }

            matby(JacobiSN,var,T1,2,np,np);
            mattransp2 (JacobiSN, T2, 2, np);
            matby(T1,T2,vSN,2,np,2);
/*
            fputs("\nvar(dS,dN):\n", fout);
            matout(fout,vSN,2,2);
*/
            fprintf(fout,"dN = %7.5f +- %.5f   dS = %7.5f +- %.5f",
                 dN,(vSN[3]>0?sqrt(vSN[3]):-0),dS,(vSN[0]>0?sqrt(vSN[0]):-0));
            fprintf(fout," (by method 1)\n");

            T1[0]=var[0]; T1[1]=T1[2]=var[0*np+np-1];
            T1[3]=var[(np-1)*np+(np-1)];
            if(com.getSE && !com.fix_omega)
               VariancedSdN(x[0], x[np-1], T1, vSN);

            fprintf(fout,"dN = %7.5f +- %.5f   dS = %7.5f +- %.5f",
               dN,(vSN[3]>0?sqrt(vSN[3]):-0),dS,(vSN[0]>0?sqrt(vSN[0]):-0));
            fprintf(fout," (by method 2)\n");

         }
         fflush(frst);  fflush(fout);
      }  /* for (js) */
      FPN(fds); FPN(fdn); FPN(ft);   fflush(fds); fflush(fdn); fflush(ft); 
   }     /* for (is) */

   com.ls=(int)ls0;  FOR(k,com.ns) com.z[k]=pz0[k];  
   com.npatt=npatt0;  FOR(h,npatt0) com.fpatt[h]=fpatt0[h]; free(fpatt0);
   return (0);
}


double lfun2AA (double t)
{
/* likelihood function for two amino acid sequences
         prob(i,j) = PI_i * p(i,j,t)
   
   The data are clean & coded (com.z[0] & com.z[1]).  
   Transition probability pijt is calculated for observed patterns only.
*/
   int n=20, h,k, aa0,aa1;
   double lnL=0, pijt,expt[20],al=com.alpha;

   if(al==0)  FOR(k,n) expt[k]=exp(t*Root[k]);
   else       FOR(k,n) expt[k]=pow(al/(al-t*Root[k]),al);
   FOR(h,com.npatt) {
      aa0=com.z[0][h]; aa1=com.z[1][h];
      for(k=0,pijt=0;k<n;k++) pijt+=U[aa0*n+k]*expt[k]*V[k*n+aa1];
      lnL-=log(com.pi[aa0]*pijt)*com.fpatt[h];
   }
   return (lnL);
}

int PairwiseAA (FILE *fout, FILE*f2AA)
{
/* Calculates pairwise distances using amino acid seqs.
   Data (com.z[]) are clean and coded.
   com.npatt for the whole data set is used which may be greater than 
   the number of patterns for each pair.
   SE is not calculated.
*/
   char *pz0[NS];
   int n=com.ncode, j, is,js;
   double x, xb[2]={0,19}, lnL, step;

   if (com.ngene>1 && com.Mgene==1) error2("ngene>1 to be tested.");
   if (noisy) printf("\npairwise ML distances of AA seqs.\n");

   if(com.model>Empirical_F)  error2("PairwiseAA: model wrong");
   if(com.model==0)  fillxc(com.pi,1./n, n);
   if(com.model>=Empirical)  GetDaa(NULL, com.daa);
   if(com.model==0 || com.model==Empirical)  EigenQaa(NULL, Root, U, V, NULL);

   FOR(j,com.ns) pz0[j]=com.z[j];
   fprintf(fout,"\nML distances of aa seqs.\n");
   if(com.alpha) 
      fprintf(fout,"\ncontinuous gamma with alpha = %.3f used.\n\n",com.alpha);

   fprintf(f2AA,"%6d\n", com.ns);
   for(is=0; is<com.ns; is++,FPN(F0),FPN(fout),FPN(f2AA)) {
      printf ("%4d vs", is+1);
      fprintf(f2AA,"%-14s ", com.spname[is]);
      fprintf(fout, "%-14s ", com.spname[is]);
      FOR(js,is) {
         com.z[0]=pz0[is]; com.z[1]=pz0[js]; 
         printf (" %2d", js+1);
         if(com.model==1||com.model==Empirical_F) {
            for (j=0,zero(com.pi,n); j<com.npatt; j++) {
               com.pi[com.z[0][j]]+=com.fpatt[j];
               com.pi[com.z[1][j]]+=com.fpatt[j];
            }
            abyx(1./sum(com.pi,n), com.pi, n);
            getpi_sqrt (com.pi, com.pi_sqrt, n, &com.npi0);
            EigenQaa(NULL,Root,U,V,NULL);
         }
         /* com.posG[1]=com.npatt; */

         xb[0]=SeqDistance[is*(is-1)/2+js];  x=xb[0]*1.5;  step=xb[0];
         LineSearch(lfun2AA, &lnL, &x, xb, step, 1e-7);
         fprintf(f2AA," %7.4f",x); fprintf(fout," %7.4f",x); 
         if (com.getSE) ;
      }  /* for (js) */
   }     /* for (is) */

   FOR(j,com.ns) com.z[j]=pz0[j];
   return (0);
}


int GetAASiteSpecies(int species, int sitepatt)
{
   int iaa;
   char str[4]="   ";

   if(com.cleandata) {
      iaa=com.z[species][sitepatt];
      iaa=GeneticCode[com.icode][FROM61[iaa]];
   }
   else {
      Codon2AA(com.z[species]+sitepatt*3, str,com.icode,&iaa);
      if(iaa==-1) iaa=20;
   }
   return (iaa);
}


void rainbowRGB (double temperature, int *R, int *G, int *B)
{
/* This returns the RGB values, each between 0 and 255, for given temperature 
   value in the range (0, 1) in the rainbow.  
   Curve fitting from the following data:

    T        R       G       B
    0        14      1       22
    0.1      56      25      57
    0.2      82      82      130
    0.3      93      120     60
    0.4      82      155     137
    0.5      68      185     156
    0.6      114     207     114
    0.7      223     228     70
    0.8      243     216     88
    0.9      251     47      37
    1        177     8       0

*/
   double T=temperature, maxT=1;

   *R = (int)fabs( -5157.3*T*T*T*T + 9681.4*T*T*T - 5491.9*T*T + 1137.7*T + 6.2168 );
   *G = (int)fabs( -1181.4*T*T*T + 964.8*T*T + 203.66*T + 1.2028 );
   *B = (int)fabs( 92.463*T*T*T - 595.92*T*T + 481.11*T + 21.769 );

   if(*R>255) *R=255;
   if(*G>255) *G=255;
   if(*B>255) *B=255;
}


int lfunNSsites_rate (FILE* frst, double x[], int np)
{
/* This calculates the dN/dS rates for sites under models with variabel dN/dS 
   ratios among sites (Nielsen and Yang 1998).  Modified from lfundG() 
   com.fhK[] holds the posterior probabilities.
*/
   int  h,hp, ir, i,it=0, refsp=0,iaa, k=com.ntime+com.nrgene+com.nkappa;
   double lnL=0, fh, cutoff=0.5;
   double w2=x[com.np-1],psel=0, *meanw, maxmw, minmw;
   char  *sig;

   FILE *fsites, *fras;
   int  continuous=0, R,G,B;
   int  ncolors=5;  /* continuous = 0 uses the specified colors */
   char sitelabel[96], *colors[5]={"darkblue", "lightblue", "purple", "pinkred", "red"};
   char *colorvalues[5]={"[2,2,120]", "[133,57,240]", "[186,60,200]", "[200,60,160]", "[250,5,5]"};

   if(com.nparK) error2("lfunNSsites_rate to be done for HMM.");

   if((meanw=(double*)malloc(com.npatt*sizeof(double)))==NULL) 
      error2("oom lfunNSsites_rate");  /* meanw useful for NSsites only */
   if(com.aaDist==0) {
      fprintf(frst,"\n\ndN/dS for site classes (K=%d)\np: ",com.ncatG);
      FOR(i,com.ncatG) fprintf(frst,"%9.5f",com.freqK[i]);
      fputs("\nw: ",frst);
      if(com.model==0) FOR(i,com.ncatG) fprintf(frst,"%9.5f",com.rK[i]);
      else if (com.model<=NSbranch2) {
         FOR(i,com.ncatG-2) fprintf(frst,"%9.5f",com.rK[i]);
         fprintf(frst, "  *        *\nw2 = %.5f", w2);
      }
      else if (com.model==NSbranch3) {
         FOR(i,com.ncatG+1) fprintf(frst,"%9.5f",x[k+com.ncatG-1+i]);
      }
      FPN(frst);
   }
   else  fputs("\nCheck main result file for parameter estimates\n",frst);
   FPN(frst);

   fx_r(x,np);
   if(com.nnodeScale)
      FOR(h,com.npatt) {
         for(ir=1,fh=com.fhK[h]; ir<com.ncatG; ir++)
            if(com.fhK[ir*com.npatt+h]>fh) fh=com.fhK[ir*com.npatt+h];
         for(ir=0; ir<com.ncatG; ir++)
            com.fhK[ir*com.npatt+h]=exp(com.fhK[ir*com.npatt+h]-fh);
         lnL-=fh*com.fpatt[h];
      }

   FOR(h,com.npatt) {
      for (ir=0,fh=meanw[h]=0; ir<com.ncatG; ir++) {
         fh += (com.fhK[ir*com.npatt+h]*=com.freqK[ir]);  /* prior=>posterior */
         meanw[h] += com.fhK[ir*com.npatt+h]*com.rK[ir];
      }
      for (ir=0,meanw[h]/=fh; ir<com.ncatG; ir++) com.fhK[ir*com.npatt+h]/=fh;
      lnL-=com.fpatt[h]*log(fh);
   }

   fprintf(frst,"\nPosterior P for %d classes (class)",com.ncatG);
   if(com.model==0) {
      fprintf(frst,"& postmean_w");
      if(com.rK[com.ncatG-1]>1)  fprintf(frst," & P(w>1)");
   }
   fprintf(frst,"\n%s used as reference\n\n",com.spname[refsp]);
   for(h=0; h<com.ls; h++,FPN(frst)) {
      hp=com.pose[h];
      iaa=GetAASiteSpecies(refsp, hp);
      fprintf(frst,"%4d %c  ",h+1,AAs[iaa]);
      for (ir=0,it=0,psel=0; ir<com.ncatG; ir++) {
         fprintf(frst," %9.4f", com.fhK[ir*com.npatt+hp]);
         if(com.fhK[ir*com.npatt+hp]>com.fhK[it*com.npatt+hp])  it=ir;
         if(com.model==0) {
            if(com.rK[ir]>1) psel+=com.fhK[ir*com.npatt+hp];
         }
      }
      fprintf(frst, " (%2d)", it+1);
      if(com.model==0) {
         fprintf(frst, "  %9.3f", meanw[hp]);
         if(psel) fprintf(frst, "  %9.3f", psel);
      }
   }

   /* list of positively selected sites */
   if(com.model==0) { /* NSsites models */
      for(ir=0,it=0; ir<com.ncatG; ir++) 
         if(com.freqK[ir]>1e-6 && com.rK[ir]>1) it=1;
      if(!com.aaDist && it) {
         fprintf(fout,"\nPositively selected sites\n\n\tProb(w>1)  mean w\n\n");
         fprintf(frst,"\nPositively selected sites\n\n\tProb(w>1)  mean w\n\n");
         FOR(h,com.ls) {
            hp=com.pose[h];
            for (ir=0,psel=0; ir<com.ncatG; ir++)
               if(com.rK[ir]>1) psel+=com.fhK[ir*com.npatt+hp];

            if(psel>cutoff) {
               sig="  ";  if(psel>.95) sig="* ";  if(psel>.99) sig="**";
               iaa=GetAASiteSpecies(refsp, hp);
               fprintf(fout,"%6d %c %.4f%s  %.2f\n",h+1,AAs[iaa],psel, sig, meanw[hp]);
               fprintf(frst,"%6d %c %.4f%s  %.2f\n",h+1,AAs[iaa],psel, sig, meanw[hp]);
            }
         }
         FPN(fout);
         if(com.rK[com.ncatG-2]>1)
            fputs("\nNote: more than one w>1.  Check rst for details\n",fout);
         if(com.rK[com.ncatG-1]<1.2)
            fputs("\nNote: w is small although >1.  Excise caution about the list.\n",fout);
      }
   }
   else if(com.model<=NSbranch2) {  /* branch&site models */
      if(com.rK[0]>1 || com.rK[1]>1) {  /* positive sites for all lineages */
         fputs("\n\nPositive sites for all lineages Prob(w>1):\n",fout);
         FOR(h,com.ls) {
            hp=com.pose[h]; iaa=GetAASiteSpecies(refsp, hp);
            psel=0;  if(com.rK[0]>1) psel=com.fhK[0*com.npatt+hp];
            if(com.rK[1]>1) psel+=com.fhK[1*com.npatt+hp];
            if(psel>cutoff) {
               sig="";  if(psel>.95) sig="*";  if(psel>.99) sig="**";
               fprintf(fout,"%6d %c %.4f%s\n",h+1,AAs[iaa],psel,sig);
            }
         }
      }
      if(w2>1) {  /* for foreground branches */
         fputs("\n\nPositive sites for foreground lineages Prob(w>1):\n",fout);
         FOR(h,com.ls) {
            hp=com.pose[h]; iaa=GetAASiteSpecies(refsp, hp);
            psel=com.fhK[2*com.npatt+hp]+com.fhK[3*com.npatt+hp];
            if(psel>cutoff) {
               sig="";  if(psel>.95) sig="*";  if(psel>.99) sig="**";
               fprintf(fout,"%6d %c %.4f%s\n",h+1,AAs[iaa],psel,sig);
            }
         }
      }
   }
   fprintf (frst,"\n\nlnL = %12.6f\n", -lnL);

   /* RasMol script for coloring structure */
   if(com.verbose && com.model==0) {
      fsites=(FILE*)fopen("SiteNumbering.txt", "r");
      if(fsites) {
         puts("\nCollecting RasMol commands for coloring structure into RasMol.txt");

         printf("Choose color scheme (0: %d colors, 1: white->red, 2: rainbow) ",ncolors);
         scanf("%d", &continuous);

         fras=(FILE*)gfopen("RasMol.txt", "w");
         for(h=0,maxmw=0,minmw=99; h<com.npatt; h++) {
            if(maxmw<meanw[h]) maxmw=meanw[h]; 
            if(minmw>meanw[h]) minmw=meanw[h]; 
         }
         if(continuous==0)
            for (it=0; it<ncolors; it++)
               printf("\t%-10s %-20s mean_w < %7.5f\n", 
                  colors[it], colorvalues[it], (it+1)*(maxmw-minmw)/ncolors);
         fprintf(fras, "cartoon\nset background white\n");
         FOR(h,com.ls) {
            fscanf(fsites, "%d%s", &it, sitelabel);
            if(it-1!=h)  { puts("Site number wrong.  Giving up."); break; }
            if(strchr(sitelabel, '?')) continue;
            hp=com.pose[h];

            if(continuous==0) {
               for (it=0; it<ncolors; it++)
                  if(meanw[hp]<minmw+(it+1.)*(maxmw-minmw)/ncolors+1e-9) break;
               fprintf(fras,"select %s\n\t\tcolor %s\n", sitelabel, colorvalues[it]);
            }
            else if (continuous==1) {
               it = 5+(int)(245*(meanw[hp]-minmw)/(maxmw-minmw+1e-9));
               fprintf(fras,"select %s\n\t\tcolor [250, %d, %d]\n", sitelabel, 255-it,255-it);
            }
            else {
               rainbowRGB((meanw[hp]-minmw)/(maxmw-minmw+1e-9), &R, &G, &B);
               fprintf(fras,"select %s\n\t\tcolor [%d, %d, %d]\n", sitelabel, R,G,B);
            }
         }
         fclose(fsites);  fclose(fras);
      }
   }
   free(meanw);
   return (0);
}


#ifdef NSSITESBandits
void finishup(void)
{
   FILE *fend=NULL;
   fend=(FILE*)gfopen("finished","w");
}
#endif


/*

(*) Codon models for variable dN/dS ratios among sites
    (com.nrate includes kappa & omega) (see also CDFdN_dS)

    NSsites          npara

    0  one-ratio     0:    one ratio for all sites
    1  neutral       1:    p0 (w0=0, w1=1)
    2  selection     3:    p0, p1, w2 (w0=0, w1=1)
    3  discrete      2K-1: p0,p1,..., and w0,w1,...
    4  freqs         K:    p's (w's are fixed)
    5  gamma         2:    alpha, beta
    6  2gamma        4:    p0, alpha1,beta1, alpha2=beta2
    7  beta          2:    p_beta, q_beta
    8  beta&w        4:    p0, p_beta, q_beta, w estimated
    9  beta&gamma    5:    p0, p_beta, q_beta, alpha, beta
   10  beta&1+gamma  5:    p0, p_beta, q_beta, alpha, beta (1+gamma used)
   11  beta&1>normal 5:    p0, p_beta, q_beta, mu, s    (normal truncated w>1)
   12  0&2normal     5:    p0, p1, mu2, s1, s2
   13  3normal       6:    p0, p1, mu2, s0, s1, s2
   14  M8a:beta&w=1  3:    p0, p_beta, q_beta, w=1 fixed
   15  M8a:beta&w>=1 4:    p0, p_beta, q_beta, w>=1 estimated

NSsites = 14 forces change to fix_omega so we can't have 2 models in one run.
NSsites = 15 would not set omegab[0] correctly for the next tree.


(*) Codon models for variable dN/dS ratios among both branches and sites
    (model=2, NSsites=3 or 2)
    (com.nrate includes kappa & omega)
    Parameters include branchlens, kappa, p0, p1, w0, w1, w2

    method = 0: SetPSiteClass copies w's to nodes[].omega and PMat is calculated
    in ConditionalPNode().  
    method = 1: PMat for branch of interest is calulated in lfuntdd_SiteClass().
    The two sets of branches have different Qfactor_NS: Qfactor_NS_branch[2].
    August 2000.


(*) Codon (perhaps aa as well) models for site-class models

    NSsites=3, ncatG=3 or 2 etc
  
    aaDist: 
       1-6 for G1974,Miyata,c,p,v,a
       FIT1 & FIT2 (11, 12): fitness model F_j = a_p*(p-p*)^2+a_v*(v-v*)^2
       FIT1:   w_ij = exp(F_j - F_i)
       FIT2:   w_ij = b*exp(F_j - F_i)

       FIT1 & FIT2 are also implemented for NSsites=0


(*) Amino acid models

    REVaa: The symmetrical part (S) of the rate matrix Q=S*PI are estimated, 
           making up 19*20/2-1=189 rate parameters for the matrix.  The aa 
           frequencies are estimated using the observed ones.  The Sij for 
           ijAAref=19*naa+9 (I-V) is set to one and others are relative rates;
    REVaa_0: AA1STEP[i*(i+1)+j] marks the aa pair i & j that are 
            interchangeable.  Sij for ijAAref=19*naa+9 (I-V) is set to one 
            and others are relative rates;


(*)
    Codon & amino acid models

    AAClasses: OmegaAA[i*(i-1)/2+j] marks the dN/dS ratio class for the pair 
            i & j.  Note kappa is before omega's in x[].
            OmegaAA[i*(i-1)/2+j]=-1, if AAs i & j are not interchangeable
                       =0,  for the background ratio
                       =1,...,nclass for AAs i & j specified in OmegaAA.dat.
            The total number of classes (com.nOmega) is one plus the number 
            specified in the file OmegaAAf.

   com.nOmega is the number of different dN/dS ratios in the NSbranchB, NSbranch2 models
      and in AAClasses.
   nodes[].label marks the dN/dS ratio for the node in the NSbranchB NSbranch2 models
   AA1STEP[i*(i-1)/2+j] =1 if AAs i & j differ at one codon position;
                        =0 otherwise.

(*) Codon and amino acid models

    aaDist = -5,-4,-3,-2,-1,1,2,3,4,5: 
    Geometric and linear relationships between amino acid distance and 
    substitution rate:
       wij = a*(1-b*dij/dmax)
       wij = a*exp(-b*dij/dmax)
    aaDist = 0:equal, +:geometric; -:linear, {1-5:Grantham,Miyata,c,p,v}

    aaDist = 11, 12: fitness models, see above.
*/


void Get4foldSites(void)
{
/* This collects the four-fold degerate sites into a file named 4fold.nuc.
   The data are not coded yet, and the routine is called from ReadSeq().
*/
   int ls4, j,k,h, ib[3][4], nb[3];
   char file4[12]="4fold.nuc", *mark4;
   FILE *f4;

   f4=gfopen(file4,"w");
   if ((mark4=(char*)malloc(com.ls*sizeof(char)))==NULL) error2("oom mark");
   FOR(h,com.ls) mark4[h]=0;

   for (h=0,ls4=0; h<com.ls; h++) {
      FOR(k,3)  NucListall(com.z[0][h*3+k], &nb[k], ib[k]);
      if(nb[0]==1 && nb[2]==1 && FourFold[ib[0][0]][ib[1][0]]) {
         for(j=1; j<com.ns; j++)
            FOR(k,2) if(com.z[j][h*3+k]!=com.z[0][h*3+k]) goto nextsite;
         mark4[h]=1;  ls4++;
      }
      nextsite: ;
   }     /* for(h) */

   fprintf (f4, "%6d  %6d\n", com.ns, ls4);
   for (j=0; j<com.ns; j++) {
      fprintf (f4, "\n%s\n", com.spname[j]);
      for (h=0; h<com.ls; h++)
         if(mark4[h]) fprintf (f4, "%c", com.z[j][h*3+2]);
      FPN (f4);
   }
   fprintf(f4, "\n\ncodons included\n");
   FOR(h,com.ls) if(mark4[h]) fprintf(f4, " %2d", h+1);  FPN(f4);

   fclose(f4);  free(mark4);
}


double distanceHKY85 (double x[], double *kappa, double alpha);

void d4dSdN(FILE* fout)
{
/* This looks at the 4-fold degerenate sites.
*/
   char str1[4]="   ", str2[4]="   ";
   int i,j,k, n=com.ncode, b[2][3], ic1,ic2,iaa;
   double pS4,d4,kappa4fold;
   double fij, fij4f[4*4], pi4f[4], pstop,t, S,dS,dN,dN_dS;
   double fb3x4[12]={.25, .25, .25, .25, 
                     .25, .25, .25, .25, 
                     .25, .25, .25, .25};

   int nii=18 /* 18 */,ii;
   double t0[]={0.001, 0.01,0.05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.2, 1.5, 2,2.5,3};


   com.ls=1; com.kappa=3; com.omega=1;

   fb3x4[0*4+0]=0.35;
   fb3x4[0*4+1]=0.15;
   fb3x4[0*4+2]=0.35;
   fb3x4[0*4+3]=0.15;
/*
   fb3x4[1*4+0]=0.35;  
   fb3x4[1*4+1]=0.15;
   fb3x4[1*4+2]=0.35;
   fb3x4[1*4+3]=0.15;
*/
   fb3x4[2*4+0]=0.35;  
   fb3x4[2*4+1]=0.15;
   fb3x4[2*4+2]=0.35;
   fb3x4[2*4+3]=0.15;


printf("\tt\tS\tdS\tdN\tdN/dS\tS4\td4\tk_4f\tpT_4f\n");

   zero(com.pi,64);
   FOR(k,64)  if(FROM64[k]>-1)
      com.pi[FROM64[k]]=fb3x4[k/16]*fb3x4[4+(k/4)%4]*fb3x4[8+k%4];
   pstop=1-sum(com.pi,n);
   abyx(1/(1-pstop),com.pi,n);

   getpi_sqrt (com.pi, com.pi_sqrt, com.ncode, &com.npi0);
   EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, &com.kappa,com.omega,PMat);

matout(frst,com.pi,16,4);

   FOR(ii,nii) {
      t=t0[ii];
      EigenQc(1,t,&S,&dS,&dN,NULL,NULL,NULL, &com.kappa,com.omega,PMat);      
      PMatUVRoot (PMat, t, n, U, V, Root);
      if(testTransP(PMat,n)) error2("testP");

matout(frst,PMat,n,n);

      for(i=0,zero(fij4f,16);i<n;i++) {
         ic1=FROM61[i]; b[0][0]=ic1/16; b[0][1]=(ic1/4)%4; b[0][2]=ic1%4;
         iaa=GeneticCode[com.icode][ic1];
         ic1-=b[0][2];
         FOR(k,4)  if(GeneticCode[com.icode][ic1+k]!=iaa)  break;
         if(k<4) continue;
         FOR(j,n) {
            fij=com.pi[i]*PMat[i*n+j];
            ic2=FROM61[j]; b[1][0]=ic2/16; b[1][1]=(ic2/4)%4; b[1][2]=ic2%4;

            if(b[0][0]!=b[1][0] || b[0][1]!=b[1][1]) continue;
            fij4f[b[0][2]*4+b[1][2]] += fij;

/* printf("%c %s %s  %.8f\n",AAs[iaa],getcodon(str1,ic1+b[0][2]),getcodon(str2,ic2),fij);
*/

         }
      }

      pS4=sum(fij4f,16)/3;
      abyx(1/sum(fij4f,16),fij4f,16);
      FOR(k,4) pi4f[k]=sum(fij4f+k*4,4);

/* matout(F0,fij4f,4,4); */

      d4 = distanceHKY85 (fij4f, &kappa4fold, 0);
      dN_dS = (dS>0 ? dN/dS : -1);
printf("\t%.4f\t%.5f\t%.5f\t%.5f\t%.5f\t%.3f\t%.5f\t%.3f\t%.4f\n", 
       t,S/3,dS,dN,dN_dS, pS4,d4,kappa4fold,pi4f[0]);

   }

printf("\nproportion of stop codons: %.4f\n", pstop);

   exit(0);
}


double distanceHKY85 (double x[], double *kappa, double alpha)
{
/* This is from SeqDivergence(), copied here to avoid linking to SeqDivergence.
*/
   int i,j;
   double p[4], Y,R, a1,a2,b, P1,P2,Q,fd,tc,ag;
   double largek=999, larged=9;

   if (testXMat(x)) 
      error2 ("X err..");
   for (i=0,fd=1,zero(p,4); i<4; i++) {
      fd -= x[i*4+i];
      FOR (j,4) { p[i]+=x[i*4+j]/2;  p[j]+=x[i*4+j]/2; }
   }
   P1=x[0*4+1]+x[1*4+0];
   P2=x[2*4+3]+x[3*4+2];
   Q = fd-P1-P2;
   Y=p[0]+p[1];    R=p[2]+p[3];  tc=p[0]*p[1]; ag=p[2]*p[3];

   a1=1-Y*P1/(2*tc)-Q/(2*Y);
   a2=1-R*P2/(2*ag)-Q/(2*R);
   b=1-Q/(2*Y*R);
   if (a1<=0 || a2<=0 || b<=0) return (larged);
   if (alpha<=0) { a1=-log(a1); a2=-log(a2); b=-log(b); }
   else   { a1=-gammap(a1,alpha); a2=-gammap(a2,alpha); b=-gammap(b,alpha);}
   a1 = -R/Y*b + a1/Y;
   a2 = -Y/R*b + a2/R;
   if (b>0) *kappa = min2((a1+a2)/(2*b), largek);
   return 2*(p[0]*p[1] + p[2]*p[3])*(a1+a2)/2 + 2*Y*R*b;
}
