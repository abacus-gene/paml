/* CODEML.c  (AAML.c & CODONML.c)
   version 3.56c. (c) Copyright 1993 by Ziheng Yang.
   Written by Ziheng Yang.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed.

   Maximum likelihood parameter estimation for codon sequences (seqtype=1) 
                    or amino-acid sequences (seqtype=2)

               cc -o codeml -fast codeml.c tools.o -lm
       cc -o codemlsites -fast -DNSSITESBATCH codeml.c tools.o -lm 
                         codeml <ControlFileName>
*/

/*
#define NSSITESBATCH
#define NSSITES_K1_K2_CLASSES
*/
/*
#define DSDN_MC  1
#define DSDN_MC_SITES  1
*/


#include "paml.h"

#ifdef macintosh
/* Added by Andrew Rambaut to accommodate Macs -
   (1) Brings up dialog box to allow command line parameters.
   (2) Reduces the size of statically allocated data.    */
#include <console.h>
#endif
#define NS            1000
#define NBRANCH       (NS*2-2)
#define NNODE         (NS*2-1)
#define NGENE         100
#define LSPNAME       30
#define NCODE         64
#define NCATG         40

#define NP            (NBRANCH*2+NGENE-1+2)
/*
#define NP            (NBRANCH+NGENE-1+189+2)
*/
extern char BASEs[],AAs[],Nsensecodon[];
extern int noisy, NFunCall, NEigenQ, NPMatUVRoot, *ancestor, GeneticCode[][64];
extern double *SeqDistance;

int Forestry (FILE *fout);
int sortwM3(double x[]);
void DetailOutput (FILE *fout, double x[], double var[]);
int GetOptions (char *ctlf);
int testx (double x[], int np);
int SetxBound (int np, double xb[][2]);
int SetxInitials (double x[]);
int GetInitials (double x[], int*fromfile);
int SetParameters (double x[]);
int SetPGene (int igene, int _pi, int _UVRoot, int _alpha, double xcom[]);
int SetPSiteClass(int iclass, double xcom[]);
int PMatJC69like (double P[], double t, int n);
int setmark_61_64 (void);
int printfcode (FILE *fout, double fb61[], double space[]);
int InitializeCodon (FILE *fout, double space[]);
int CodonListall(char codon[], int nb[3], int ib[3][4]);
int AA2Codonf (double faa[20], double fcodon[]);
int DistanceMatAA (FILE *fout);
int DistanceMatNG86 (FILE *fout, double alpha);
int GetDaa(FILE *fout, double daa[]);
void getpcodonClass(double x[], double pcodonClass[]);
int EigenQc(int getstats, double branchl, double *S, double *dS, double *dN,
    double Root[], double U[], double V[],
    double kappa[], double omega, double Q[]);
int EigenQaa(FILE *fout, double Root[], double U[], double V[],double rate[]);
int Qcodon2aa (double Qc[], double pic[], double Qaa[], double piaa[]);
int SetAA1STEP (void);
int GetOmegaAA(int OmegaAA[]);
int TestModelQc(FILE *fout, double x[]);
double lfun2dSdN (double x[], int np);
int VariancedSdN(double t, double omega, double vtw[2*2], double vdSdN[2*2]);
int GetCodonFreqs(double pi[]);
int PairwiseCodon (FILE *fout, double space[]);
int PairwiseAA (FILE *fout);
int lfunNSsites_rate (FILE* fout, double x[], int np);
int PartialLikelihood (int inode, int igene);
double CDFdN_dS(double x,double par[]);
int DiscreteNSsites(double par[]);
int mergeSeqs(FILE*fout);
int SlidingWindow(FILE*fout, double space[]);
void Get4foldSites(void);

int eigenQREV(double Q[], double pi[], double pi_sqrt[], int n, 
              double Root[], double U[], double V[]);

void SimulateData2s61(void);
void Ina(void);
void d4dSdN(FILE*fout);

struct common_info {
   char *z[NS], *spname[NS], seqf[96],outf[96],treef[96],daafile[96];
   int seqtype, ns, ls, ngene, posG[NGENE+1], lgene[NGENE], npatt,*pose;
   int runmode,clock,verbose,print, codonf,aaDist,model,NSsites, cleandata;
   int method, icode, ncode, Mgene, ndata;
   int fix_rgene,fix_kappa,fix_omega,fix_alpha,fix_rho,nparK,fix_branch,getSE;
   int np, ntime, nrgene, nkappa, nrate, nalpha, ncatG, sspace, slkl1;
   int npi0, hkyREV;
   double *fpatt, *space;
   double pi[NCODE],fb61[64],piG[NGENE][64],kappa,omega,alpha,rho,rgene[NGENE];
   double pi_sqrt[NCODE];
   double freqK[NCATG], rK[NCATG], MK[NCATG*NCATG],daa[20*20], *lkl,*lkl0,*fhK;
   double (*plfun)(double x[],int np);
}  com;
struct TREEB {
   int  nbranch, nnode, root, nlabel, branches[NBRANCH][2];
   double lnL;
}  tree;
struct TREEN {
   int father, nson, sons[NS], ibranch, label;
   double branch, divtime, omega, *lkl;
}  *nodes;


extern double Small_Diff;
int FROM61[64], FROM64[64];
int ChangedInIteration;  /* 1: t changed, update P(t); 2: paras changed, update UVRoot */
double *PMat,*U,*V,*Root, *_UU[5],*_VV[5],*_Root[5];
/* for 5 sets for branchsite models */

double xcom[NP-NBRANCH];
double *POMEGA, *PKAPPA;
  /* communication between SetParameters & PartialLikelihood & EigenQc.  
     Remove it? */
double pcodon0[64],paa0[20], *pcodonClass;  /* for aaDist=FIT1 or FIT2 */

int LASTROUND, UVRootChanged=0;
int IClass=-1;

double OMEGA_FIX=-1;  
/* fix the last dN/dS in the NSbranchB, NSbranch2 models with variable dN/dS ratios 
   for lineages.  Useful for testing whether w>1 for particular lineages. 
*/
int N_OMEGA=-1; /* number of omega for branches or in AAClass */
int OmegaAA[190], AA1STEP[190];

double _rateSite=1;
int *rateBranch=NULL, N_rateBranch=-1; /* rates for branches, local clocks */
double Qfactor_NS, Qfactor_NS_branch[2], freqK_NS;

/* is lkl memory allocated for each class? 
   =0 if (method==0) and =1 if (method==1) */
int lklSiteClass=0;
char _oldlkl[NNODE];       /* update lkl for nodes? to save computation */
int _nnodeScale;
char _nodeScale[NNODE];    /* nScale[ns-1] for interior nodes */
double *_nodeScaleF=NULL;  /* nScaleF[npatt] for scale factors */

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
enum {Fequal, F1x4, F3x4, Fcodon} CodonFreqs;
char *codonfreqs[]={"Fequal", "F1x4", "F3x4", "Fcodon"};
enum {NSbranchB=1, NSbranch2} NSBranchModels;   /* , AAClasses=7 */
char *NSbranchmodels[]={"One dN/dS ratio", 
     "free dN/dS Ratios for branches", "several dN/dS ratios for branches",
     "", "", "", "","AAClasses"};
enum {Poisson, EqualInput, Empirical, Empirical_F,
     FromCodon=6, AAClasses=7, REVaa_0=8, REVaa=9} AAModel;
char *aamodels[]={"Poisson", "EqualInput", "Empirical", "Empirical_F", "",
     "", "FromCodon", "AAClasses", "REVaa_0", "REVaa"};
enum {NSneutral=1, NSselection, NSdiscrete, NSfreqs, NSgamma, NS2gamma, 
     NSbeta, NSbetaw, NSbetagamma, NSbeta1gamma, NSbeta1normal, NS02normal, 
     NS3normal} NSsitesModels;
char *NSsitesmodels[]={"one-ratio","neutral", "selection","discrete","freqs", 
     "gamma","2gamma","beta","beta&w","beta&gamma", "beta&gamma+1", 
     "beta&normal>1", "0&2normal>0", "3normal>0"};
enum {FIT1=11, FIT2=12} SiteClassModels;


#ifdef SITELABELS
char *sitelabels[]={"-17", "-16", "-15"};
#endif


#define CODEML 1
#include "treesub.c"
#include "treespace.c"

FILE *fout;

int main(int argc, char *argv[])
{
   FILE *fseq=NULL;
   char ctlf[96]="codeml.ctl", *pmodel;
   char *seqtypestr[3]={"CODONML", "AAML", "CODON2AAML"};
   char *Mgenestr[]={"diff. rate", "separate data", "diff. rate & pi", 
                     "diff. rate & k&w", "diff. rate & pi & k&w"};
   char *clockstr[]={"", "global clock", "local clock", "clock with dated seqs"};
   int i, slkl0=0,s2=0, idata, nc, nUVR;
   double t;
   time_t tm;
   clock_t clock_start, clock_end;
   int ncatG0=10, im=0, nmodels=3, nsmodels[13]={0};


#ifdef __MWERKS__
/* Added by Andrew Rambaut to accommodate Macs -
   Brings up dialog box to allow command line parameters.
*/
   argc=ccommand(&argv);
#endif

   clock_start=clock();
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
   com.getSE=0;       com.print=0;    com.verbose=1;
   com.method=0;      com.space=NULL;

   frub=fopen("rub","w");   frst=fopen("rst","w");  frst1=fopen("rst1","w");
   if(frub==NULL || frst==NULL || frst1==NULL) error2("file open error");
/*
mergeSeqs(frst);  exit(0);
Ina();
*/
   SetSeed ((int)time(&tm));


#if (DSDN_MC || DSDN_MC_SITES)
SimulateData2s61();
#endif

   if(argc>1) strcpy(ctlf,argv[1]);
   GetOptions(ctlf);

#ifdef NSSITESBATCH
         if(com.fix_omega) error2("fix omega during batch run?");
         if(com.model) error2("model should be 0 in this batch run?");
         if(com.runmode) error2("runmode?");
         com.NSsites=NSbetaw;  com.ncatG=ncatG0+1;
#endif

   if ((fout=fopen (com.outf, "w"))==NULL) error2("outfile creation err.");

   if((fseq=fopen (com.seqf,"r"))==NULL) {
      printf ("\n\nSequence file %s not found!\n", com.seqf);
      exit (-1);
   }
   if(noisy && com.seqtype==CODONseq) 
      { printcu(F0,NULL,com.icode); puts("Nice code, uuh?"); }

   /* space for P&U&V&Root */
   if(com.seqtype==1) { nUVR=1; nc=64; if(com.model==2) nUVR=5; }
   else               { nUVR=1; nc=20; }
   PMat=(double*)malloc((nc*nc+nUVR*nc*nc*2+nUVR*nc)*sizeof(double));
   if(PMat==NULL) error2("oom getting P&U&V&Root");
   U=_UU[0]=PMat+nc*nc;  V=_VV[0]=_UU[0]+nc*nc; Root=_Root[0]=_VV[0]+nc*nc;
   for(i=1; i<nUVR; i++) {
      _UU[i]=_UU[i-1]+nc*nc*2+nc; _VV[i]=_VV[i-1]+nc*nc*2+nc; 
      _Root[i]=_Root[i-1]+nc*nc*2+nc;
   }

   /* d4dSdN(fout); */

   if (com.model==AAClasses) {
      SetAA1STEP();
      GetOmegaAA(OmegaAA);
   }
   else if (com.seqtype==AAseq && com.model==REVaa_0)
      SetAA1STEP();


   for (idata=0; idata<com.ndata; idata++) {

      if(com.ndata>1) fprintf(frst1,"%5d ",idata+1);

      if (com.ndata>1) {
         printf ("\n\nData set %d\n", idata+1);
         fprintf (fout, "\n\nData set %d\n", idata+1);
         fprintf(frst,"\t%d",idata+1);
      }
      if(idata)  GetOptions(ctlf); /* Is this necessary? */

      ReadSeq((com.verbose?fout:NULL),fseq); /*may change seqtype*/
      if(com.ndata==1) fclose(fseq);
      i=(com.ns*2-1)*sizeof(struct TREEN);
      if((nodes=(struct TREEN*)malloc(i))==NULL) error2("oom nodes");
      /* BootstrapSeq("boot");  exit(0); */

      if(com.ngene>1 && com.Mgene==1)  printSeqsMgenes ();
      if(com.ngene>1 && com.aaDist>=FIT1)  /* because of pcodon0[] */
         { error2("ngene for fitness models");}

      pmodel=(com.seqtype==CODONseq?NSbranchmodels[com.model]:aamodels[com.model]);
      fprintf(fout,"%s (in %s)  ",seqtypestr[com.seqtype-1],VerStr);
      fprintf(fout,"  %s   Model: %s ",com.seqf,pmodel);
      if(com.clock) fprintf(fout," %s ",clockstr[com.clock]);
      if(com.seqtype==CODONseq||com.model==FromCodon) {
         if (com.fix_kappa) fprintf(fout, " kappa = %.3f fixed\n", com.kappa);
         if (com.fix_omega) fprintf(fout, " omega = %.3f fixed\n", com.omega);
      }
      if (com.seqtype==AAseq && (com.model==Empirical||com.model==Empirical_F))
         fprintf (fout, "(%s) ", com.daafile);
      if(com.seqtype==AAseq&&com.nrate) fprintf(fout,"(nrate:%d) ",com.nrate);
      if (com.alpha && com.rho) fprintf (fout, "Auto-");
      if (com.alpha) fprintf (fout,"dGamma (ncatG=%d) ", com.ncatG);
      if (com.ngene>1)
         fprintf (fout, " (%d genes: %s)  ", com.ngene, Mgenestr[com.Mgene]);

      if (com.alpha==0) com.nalpha=0;
      else              com.nalpha=(com.nalpha?com.ngene:!com.fix_alpha);
      if (com.Mgene==1) com.nalpha=!com.fix_alpha;
      if(com.nalpha>1 && (!com.alpha || com.ngene==1 || com.fix_alpha))
         error2("Malpha");
      if (com.nalpha>1 && com.rho) error2("Malpha or rho");
      if (com.nalpha>1) fprintf (fout,"(%d gamma)", com.nalpha);
     
      if (com.Mgene && com.ngene==1) error2("Mgene for one gene.");
      if (com.seqtype==CODONseq) {
         fprintf (fout, "\nCodon frequencies: %s\n", codonfreqs[com.codonf]);
         if(com.alpha) 
            fputs("Warning: Gamma model for codons.  See documentation.",fout);
      }
      if ((com.seqtype==CODONseq||com.model==FromCodon) 
         && (com.aaDist &&com.aaDist<10))
         fprintf(fout,"%s, %s\n",com.daafile,(com.aaDist>0?"geometric":"linear"));

      if (com.NSsites) {
         fprintf(fout,"Site-class models");
#ifndef NSSITESBATCH
         fprintf(fout," %s",NSsitesmodels[com.NSsites]);
         if(com.NSsites>=NSdiscrete)fprintf(fout," (%d categories)",com.ncatG);
#endif
         if(com.nparK) fprintf(fout," & HMM");
         FPN(fout);
         if(com.aaDist)
            fprintf(fout,"\nFitness models: aaDist: %d\n",com.aaDist);
      }

      if(com.clock==2 && 
        (rateBranch=(int*)realloc(rateBranch,com.ns*2*sizeof(int)))==NULL) 
           error2("oom rateBranch");
      com.sspace=max2(800000,com.ls*(int)(com.ns*3*sizeof(char)+sizeof(int)));
      com.sspace=max2(com.sspace,3*com.ncode*com.ncode*(int)sizeof(double));
      i=com.ns*(com.ns-1)/2;
/*
      com.sspace=max2(com.sspace,
        (int)sizeof(double)*((com.ns*2-2)*(com.ns*2-2+4+i)+i));
*/
      if((com.space=(double*)realloc(com.space,com.sspace))==NULL)
         error2("oom space");

      SeqDistance=(double*)realloc(SeqDistance, i*sizeof(double));
      ancestor=(int*)realloc(ancestor, i*sizeof(int));
      if(SeqDistance==NULL||ancestor==NULL) error2("oom distance&ancestor");

      for(i=0;i<com.ncode;i++)  com.pi_sqrt[i]=-999;
      if(com.seqtype==AAseq) {
         Initialize (fout, com.space);
         if (com.model==FromCodon||com.model==AAClasses) 
            AA2Codonf(com.pi, com.fb61);  /* get codon freqs from aa freqs */ 
      }
      else {  /* codon sequences */
         if(InitializeCodon(fout,com.space)) error2("giving up on stop codons");
         if(com.Mgene==3) FOR (i,com.ngene) xtoy(com.pi,com.piG[i],com.ncode);
/*
matout(F0,com.pi,16,4); getchar();
FOR(i,com.ngene) { 
   matout(F0,com.piG[i],16,4); 
   getchar();
}
*/
      }

      if(com.seqtype==CODONseq && com.ns<200) DistanceMatNG86(fout,0);
      else                      DistanceMatAA(fout);

      fflush(fout);

      if(com.seqtype==AAseq && com.model==Poisson && !com.print) 
         PatternJC69like (NULL);
      if(com.alpha || com.NSsites) {
         s2=com.npatt*com.ncatG*sizeof(double);
         if((com.fhK=(double*)realloc(com.fhK,s2))==NULL) error2("oom fhK");
      }

      if(com.runmode==-2 && com.Mgene!=1) {
         if(com.seqtype==CODONseq) PairwiseCodon(fout,com.space);  
         else                      PairwiseAA(fout);  
      }
      else {
         if(!com.cleandata) {
            slkl0=com.ns*com.ncode*com.npatt*sizeof(double);
            com.lkl0=(double*)malloc(slkl0);
         }
         com.slkl1= 2 *com.ncode*com.npatt*sizeof(double);
         /* com.slkl1= (com.ns-1)*com.ncode*com.npatt*sizeof(double); */
         com.lkl=(double*)realloc(com.lkl,com.slkl1);

         printf("\n%9ld bytes for distance",com.ns*(com.ns-1)/2*sizeof(double));
         printf("\n%9d bytes for lkl0\n%9d bytes for lkl1\n",slkl0,com.slkl1);
         printf ("%9d bytes for fhK\n%9d bytes for space\n",s2,com.sspace);
         if((!com.cleandata&&com.lkl0==NULL) || com.lkl==NULL)
            error2("oom");

	 /*
SlidingWindow(fout,com.space);
exit(0);
*/
#ifdef NSSITESBATCH
         printf("\nNSsites batch run.\nUsing ncatG values as in YNGP2000.\n");
         FOR(i,NS3normal) printf("M% 2d: %s\n",i,NSsitesmodels[i]);
         printf("\nNumber of models and model codes (Suggestion: 6   0 1 2 3 7 8)? ");
         scanf("%d",&nmodels);  FOR(i,nmodels) scanf("%d",&nsmodels[i]);

         for(im=0; im<nmodels; im++) {
            com.NSsites=nsmodels[im];
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
            Forestry (fout);

            clock_end=clock(); i=(int)(t=(clock_end*1.0-clock_start)/CLOCKS_PER_SEC);
            fprintf(fout,"Time used: %02d:%02d:%02.1f.\n", i/3600,(i%3600)/60, t-(i/60)*60);
            clock_start=clock_end;

         }
         exit(0);
#endif
         if (com.Mgene==1)        MultipleGenes(fout, com.space);
         else if (com.runmode==0) Forestry(fout);
         else if (com.runmode==3) StepwiseAddition(fout, com.space);
         else if (com.runmode>=4) Perturbation(fout,(com.runmode==4),com.space);
         else                     StarDecomposition(fout, com.space);
      }
      FPN(frst);  fflush(frst);  fflush(frst1);
      free(nodes);
   }  /* for (idata) */
   if (noisy) putchar ('\a');
   fclose(frst);   if(fin)fclose(fin);
   if(com.ndata>1 && fseq) fclose(fseq);  
   free(PMat);
   clock_end=clock(); i=(int)(t=(clock_end*1.0-clock_start)/CLOCKS_PER_SEC);
   printf("\nTime used: %02d:%02d:%02.1f.\n", i/3600,(i%3600)/60, t-(i/60)*60);
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
   int pauptree=0, length_label;
   double x[NP],xb[NP][2], lnL=0,lnL0=0, e=1e-7,tl=0, *var=NULL, nchange=-1;

   if(com.clock==3) GetSeqTimes();

   if ((ftree=fopen(com.treef,"r"))==NULL) {
      printf("\ntree file %s not found.\n", com.treef);
      exit(-1);
   }
   GetTreeFileType(ftree, &ntree, &pauptree, 0);
   if (com.alpha)
      if ((frate=(FILE*)fopen(ratef,"w"))==NULL) error2("file err");
   if (ntree>10 && com.print) puts("\nlarge lnf file");
   FOR(i,NP-NBRANCH) xcom[i]=0.1;
   if((flnf=fopen("lnf","w+"))==NULL) error2("lnf file open error");
   fprintf(flnf,"%6d %6d %6d\n", ntree, com.ls, com.npatt);

   if(com.seqtype==1 && com.aaDist>=FIT1) {
      xtoy(com.pi,pcodon0,64);
      zero(paa0,20);
      FOR(i,com.ncode) paa0[GeneticCode[com.icode][FROM61[i]]]+=pcodon0[i];
      pcodonClass=(double*)malloc(com.ncatG*64*sizeof(double));
      if(pcodonClass==NULL) error2("oom pcodonClass");
   }
   if(!com.cleandata) InitPartialLikelihood();

   for(itree=0; ntree==-1||itree<ntree; iteration=1,itree++) {
      if((pauptree && PaupTreeRubbish(ftree)) || 
         ReadaTreeN(ftree,&length_label,1)) 
         { puts("err or end of tree file."); break; }

      printf("\nTREE # %2d\n", itree+1);
      fprintf(fout,"\nTREE # %2d:  ", itree+1);
      fprintf(flnf,"\n\n%2d\n", itree+1);
      if(com.print) fprintf (frst,"\n\nTREE # %2d\n", itree+1);
      fprintf(frub,"\n\nTREE #%2d\n", itree+1);

      LASTROUND=0;
      if(com.cleandata) nchange=MPScore(com.space);
      if(com.ns<40) { OutaTreeN(F0,0,0); printf("   MP score: %.0f",nchange); }
      OutaTreeN(fout,0,0); fprintf(fout,"   MP score: %.0f",nchange);

      if (times++==0 && (length_label==1 || length_label==3)) {
         if(com.clock) puts("\nBranch lengths in tree are ignored");
         else {
            puts("\ntree has branch lengths.\n0:ignore?\n1:initials?\n2:fixed?");
            com.fix_branch=1;
            scanf("%d",&com.fix_branch);

            if(com.fix_branch==1) {
               FOR(i,tree.nnode) 
                  if(i!=tree.root && (x[nodes[i].ibranch]=nodes[i].branch)<0)
                     error2("branch length in tree < 0");
            }
         }
      }
      if(!com.clock && nodes[tree.root].nson<=2 && com.ns>2) {
         puts("\nThis is a rooted tree, without clock.  Check.");
         fputs("\nThis is a rooted tree.  Please check!",fout);
      }
      fflush(fout),  fflush(flnf);

      GetInitials(x, &i);
      if(i==-1) iteration=0;
      if(iteration) SetxInitials (x); /* start within the feasible region */

      if((np=com.np)>NP || np-com.ntime>NP-NBRANCH) error2("raise NP");
      if((i=spaceming2(np))>com.sspace) {
         if((com.space=(double*)realloc(com.space,com.sspace=i))==NULL) 
            error2("oom space");
         printf ("\n%9d bytes for space, adjusted\n",com.sspace);
      }
      printf("\nntime & nrate & np:%6d%6d%6d\n",com.ntime,com.nrate,com.np);
      if(itree && !fin) for(i=0;i<np-com.ntime;i++) x[com.ntime+i]=xcom[i];
/*
      if(com.seqtype==CODONseq && com.NSsites==0 && com.model){
         printf("\n%d dN/dS ratios for branches assumed:\n",N_OMEGA);
         FOR(i,tree.nbranch) printf("%4d",OmegaBranch[i]); FPN(F0);
         fprintf (fout, "\n%d dN/dS ratios for branches assumed:\n",N_OMEGA);
         FOR(i,tree.nbranch) fprintf(fout,"%4d",OmegaBranch[i]); FPN(fout);
      }
*/
      if (com.clock==2) {
         printf("\n%d rates for branches assumed:\n",N_rateBranch);
         FOR (i,tree.nbranch) printf("%3d", rateBranch[i]); FPN(F0);
         FPN(fout); FOR(i,tree.nbranch)fprintf(fout,"%3d",rateBranch[i]); 
         FPN(fout);
      }

      PointLklnodes ();
      lnL = com.plfun (x,np);

      if(noisy) {
         printf("\nnp =%6d", np);
         if(noisy>2 && np<100) matout(F0,x,1,np);
         printf("\nlnL0 = %12.6f\n",-lnL);

/*
         gradient (np, x, lnL, com.space, com.plfun, space+np, 1);
         FOR(i,np) printf("%12.6f", com.space[i]);  FPN(F0);
*/
      }

      if(iteration) {
         SetxBound(np,xb);
         SetxInitials (x);
         if(com.method==1) {
            j=minB(noisy>2?frub:NULL, &lnL,x,xb, com.space);
            if(j==4) SetxBound(np,xb);
         }
         if(j==4 || com.method==0)
            j=ming2(noisy>2?frub:NULL,&lnL,com.plfun,NULL,x,xb, com.space,e,np); 

         if (j==-1 || lnL<=0 || lnL>1e7) status=-1;  else status=0;
         if(j) fprintf(fout,"\ncheck convergence..");

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

         if(!com.nparK) {
            j = (com.NSsites==NS3normal ? 3 : com.ncatG);
            if(com.model) j=3;
            f_to_x(x+k,x+k,j,0,0);
         }
         else {  /* HMM model for w */
            k+=com.ncatG;
            for(i=0; i<com.ncatG; i++,k+=com.ncatG-1) 
               f_to_x(x+k,x+k,com.ncatG,0,0);
         }
      }
      LASTROUND=1;
      if(com.NSsites==NSdiscrete && com.aaDist==0 && com.model==0)
         sortwM3(x);

      if(com.clock) {  /* this applies to all clock models (clock=1,2,3) */
         for(i=com.ns; i<tree.nnode; i++) 
            if(i!=tree.root) x[i-com.ns]=nodes[i].divtime;
      }

      fprintf (fout,"\nlnL(ntime:%3d  np:%3d):%14.6f%+14.6f\n",
         com.ntime, np, -lnL, -lnL+lnL0);
      if(com.fix_branch<2) { OutaTreeB(fout);  FPN(fout); }
      FOR (i,np) fprintf(fout," %8.5f",x[i]); FPN(fout); fflush(fout);

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
         FOR (i,np) var[i*np+i]=(var[i*np+i]>0.?sqrt(var[i*np+i]):-1);
         FOR (i,np) fprintf(fout," %8.5f", var[i*np+i]);  FPN (fout);
         /* if (com.getSE==2) matout2(fout,var,np,np,12,7); */
      }
      if(com.NSsites==NSselection && x[com.ntime+com.nrgene+com.nkappa+2]<1)
         fputs("\nNote: This model may have multiple optima. See doc.\n",fout);

      /* if (com.clock) SetBranch (x); */
      for(i=0,tl=0;i<tree.nnode;i++) if(i!=tree.root) tl+=nodes[i].branch;
      fprintf(fout,"\ntree length = %9.5f%s\n",tl,com.ngene>1?" (1st gene)":"");
/*
if(com.NSsites) {
   fprintf(frst1,"\n%-15s%4d %9.3f%9.2f ",
      NSsitesmodels[com.NSsites],com.ls,tl,-lnL);
   if(itree) fprintf(frst,"\t%.6f",x[0]);
   fprintf(frst,"\t%.3f",-lnL);
}
*/

fprintf(frst1,"\t%d",com.npatt);
for(i=com.ntime; i<com.np; i++) fprintf(frst1,"\t%.3f",x[i]);
fprintf(frst1,"\t%.3f\n",-lnL);
fflush(frst1);

      FPN(fout); OutaTreeN(fout,0,1);  FPN(fout);
      FPN(fout); OutaTreeN(fout,1,1);  FPN(fout);

      if(com.np-com.ntime||com.clock) DetailOutput(fout,x,var);
      if (com.seqtype==AAseq && com.model>=2 /* AAClasses */)
         EigenQaa(fout, Root, U, V, x+com.ntime+com.nrgene); /* S & PAM */


      if (com.print) {
         if(com.rho==0 && com.nparK==0)  AncestralSeqs(frst,x);
         if (com.NSsites)                lfunNSsites_rate(frst,x,np);
         else if(com.plfun!=lfun)        lfunRates(frate,x,np);
      }
      com.print-=9;  com.plfun(x,np);  com.print+=9;
      /*if(noisy)
           printf("%d eigenQc calls, %d lfun calls\n",N_eigenQc,NFunCall);*/
   }     /* for (itree) */

   fclose(ftree); 
   if(frate) fclose(frate);
   if (com.aaDist && com.aaDist<10 
      && (com.seqtype==CODONseq||com.model==FromCodon))
      printf("\n%s, %s.\n",com.daafile,(com.aaDist>0?"geometric":"linear"));
   if(com.seqtype==1 && com.aaDist>=FIT1) free(pcodonClass);

   if(ntree==-1) {
      ntree=itree;
      rewind(flnf);
      fprintf(flnf,"%6d", ntree);
   }
   if(noisy && ntree>1) {
      rewind(flnf);
      rell(flnf,fout,ntree);
   }
   fclose(flnf);

   return (0);
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
   int i,j,k=com.ntime, np=com.np,npclass;
   double om=-1,N=-1,S=0,dN=0,dS=0,dSt,dNt, vtw[4],vSN[4], omclass[NCATG];
   double phi1=0,phi2=0, tnode,t, *dSdNb=NULL;
   double mu[3]={0,1,2},sig[3]={-1}; /* 3normal: mu0=0 fixed. mu2 estimated */

   /* date estimation under molecular clocks */

   if((com.clock==1 || com.clock==2) && noisy>=9) { 
      /* SetBranch() called before this. */
      fputs("\nNode & distance to present\n",fout);
      for(i=0; i<tree.nnode-com.ns;i++) {
         printf("Node %3d distance %9.5f",com.ns+i+1,x[i]);
         if(com.getSE) printf(" +- %9.5f",var[i*com.np+i]);
         FPN(F0);
      }

      for (; ;) {
         printf("\nreference node & node time (-1 -1 to quit)? ");
         scanf("%d%lf", &j,&tnode);
         if(j<com.ns) break;
         FPN(F0);

         fprintf(fout,"\n\nNode %d Time %.3f\n\n",j,tnode);
         if(--j>=com.ns && tnode>0) {
            for(i=com.ns, tnode/=nodes[j].divtime; i<tree.nnode;i++) {
               printf("Node %3d Time %6.2f",i+1,nodes[i].divtime*tnode);
               fprintf(fout,"Node %3d Time %6.2f",i+1,nodes[i].divtime*tnode);
               if(com.getSE) { 
                  printf(" +- %6.2f",var[(i-com.ns)*com.np+(i-com.ns)]*tnode);
                  fprintf(fout," +- %6.2f",var[(i-com.ns)*com.np+(i-com.ns)]*tnode);
               }
               FPN(F0);  FPN(fout);
            }
         }
      }
   } 

   fprintf(fout,"\nDetailed output identifying parameters\n");
   if(com.clock==2) {
      fprintf (fout,"rates for branches:    1");
      for(k=tree.nnode-com.ns; k<com.ntime; k++) fprintf(fout," %8.5f",x[k]);
   }
   else if(com.clock==3) {
      fprintf(fout,"\nMutation rate = %.2e",x[com.ntime-1]/ScaleTimes_clock3);
      if(com.getSE) fprintf(fout," +- %.2e", var[(com.ntime-1)*com.np+(com.ntime-1)]/ScaleTimes_clock3);
      FPN(fout); FPN(fout);
      for(i=0; i<tree.nnode; i++) {
         fprintf(fout,"Node %3d Time %6.2f",
            i+1,YoungDate-nodes[i].divtime*ScaleTimes_clock3);
         if(com.getSE && i>=com.ns) 
            fprintf(fout," +- %6.2f\n", var[(i-com.ns)*com.np+(i-com.ns)]*ScaleTimes_clock3);
         else FPN(fout);
      }
   }

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
               
               if(com.aaDist<10)         POMEGA=x+k+com.ncatG-1+2*j;
               else if(com.aaDist>=FIT1) {
                  POMEGA=x+k+com.ncatG-1+j*(4+(com.aaDist==FIT2));
                  xtoy(pcodonClass+j*64, com.pi, com.ncode);
               }
            }
            EigenQc(1,1,&S,&dSt,&dNt,NULL,NULL,NULL,
               (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),com.rK[j],PMat);
            /* t=1 used here, and dS & dN used later for each branch */
            dS+=freqK_NS*dSt;  dN+=freqK_NS*dNt;
            omclass[j]=dNt/dSt;
         }
         Qfactor_NS=1/Qfactor_NS;
         om=dN/dS;  dS*=Qfactor_NS;  dN*=Qfactor_NS;  N=com.ls*3-S;
      }

      k=com.ntime+com.nrgene;
      if (com.hkyREV) {
         fprintf(fout,"a (TC) & b (TA) & c (TG) d (CA) & e (CG): ");
         FOR(i,5) fprintf(fout,"%8.5f ", x[k++]);  FPN(fout);
      }
      else if (!com.fix_kappa)
         fprintf(fout,"kappa (ts/tv) = %8.5f\n", x[k++]);

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
            x[k],x[k+1],x[k+2], 1-x[k], x[k+3]);
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
         FOR(i,com.ncatG) fprintf(fout,"%9.5f",com.freqK[i]);
         fputs("\nw: ",fout);
         FOR(i,com.ncatG-2) fprintf(fout,"%9.5f",com.rK[i]);
         fprintf(fout, "  *        *\nw2 = %.5f\n", (com.fix_omega?OMEGA_FIX:x[com.np-1]));
      }
      else if (com.NSsites && com.aaDist==0) {
         fprintf(fout,"\n\ndN/dS for site classes (K=%d)\np: ",com.ncatG);
         FOR(i,com.ncatG) fprintf(fout,"%9.5f",com.freqK[i]);
         fputs("\nw: ",fout);
         FOR(i,com.ncatG) fprintf(fout,"%9.5f",com.rK[i]);  FPN(fout);
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

   if (com.model==AAClasses) {
      fprintf (fout, "\nw (dN/dS) classes for amino acid pairs:\n");
      FOR (k,N_OMEGA) {
         fprintf (fout, "%9.5f:", x[com.ntime+com.nrgene+com.nkappa+k]);
         FOR (i,20) FOR(j,i)
            if (OmegaAA[i*(i-1)/2+j]==k) fprintf(fout," %c%c", AAs[i],AAs[j]);
         if (k==0)  fprintf(fout, " (background ratio)");
         FPN(fout); 
      }
      if (com.nrate>65) { /* estimating all acceptance rates */
         for(i=0,k=com.ntime+com.nrgene+com.nkappa; i<20; i++) FOR(j,i) {
            om=0;
            if (AA1STEP[i*(i-1)/2+j]) {
               om=x[k]*100;
               if (com.seqtype!=CODONseq && i*(i-1)/2+j==ijAAref)  om=100;
               else  k++;
            }
            fprintf(frst,"\t%d\t%d\t%.2f\n", i+1,j+1,om*AA1STEP[i*(i-1)/2+j]);
         }
      }
   }

   /* dN and dS for each branch in the tree */
   if(com.seqtype==CODONseq && com.ngene==1 && (com.model==0 || com.NSsites==0)
      /*||com.model==FromCodon||com.model==AAClasses */){
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
            if (com.model==AAClasses || com.aaDist) om=-1;
            else if (com.model==0 || com.model==FromCodon)
               om=(com.fix_omega?com.omega:x[k]);
            else if (com.model==NSbranchB) om=x[k+i];
            else if (com.model==NSbranch2) om=nodes[tree.branches[i][1]].omega;
/*
            else if (!com.fix_omega || OmegaBranch[i]<N_OMEGA-1)
               om=x[k+=OmegaBranch[i]];
            else  om=OMEGA_FIX;
*/
            if(com.aaDist) POMEGA=x+com.ntime+com.nrgene+!com.fix_kappa;

            EigenQc(1,t,&S,&dS,&dN, NULL,NULL,NULL,
               (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),om,PMat); /* */
            if(om==-1 && dS>1e-10) om=dN/dS;
            N=com.ls*3-S;
            if(com.model) {  dSdNb[i]=dS; dSdNb[tree.nnode+i]=dN; }

            /* om not used in AAClasses model */
            if(com.getSE>1&&com.fix_branch<2&&!com.clock&&com.model!=AAClasses){
               puts ("calculate the SEs for dN & dS for each branch");
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
            Qfactor_NS=Qfactor_NS_branch[nodes[i].label];
            EigenQc(1,t,&S,&dS,&dN,NULL,NULL,NULL,
            (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),om,PMat);

            fprintf(fout,"%9.3f%9.1f%9.1f%9.4f%9.4f%9.4f %6.1f %6.1f\n",
                          t,S,N,om,dN*t,dS*t, S*dS*t,N*dN*t);
            */
         }
      }  /* for (i) */
      if(com.model &&!com.NSsites) {
         FOR(i,tree.nbranch) nodes[tree.branches[i][1]].branch=dSdNb[i]; 
         fprintf(fout,"\ndS tree:\n");  OutaTreeN(fout,1,1);
         FOR(i,tree.nbranch) nodes[tree.branches[i][1]].branch=dSdNb[tree.nnode+i];
         fprintf(fout,"\ndN tree:\n");  OutaTreeN(fout,1,1);  FPN(fout);
         free(dSdNb);
      }
   }  /* if codonseqs */

   FPN(fout);
}



int GetOptions (char *ctlf)
{
   int i,j, nopt=33, lline=255;
   char line[255], *pline, opt[99], comment='*';
   char *optstr[] = {"seqfile", "outfile", "treefile", "seqtype", "noisy", 
        "cleandata", "runmode", "method", 
        "clock", "getSE", "RateAncestor", "CodonFreq", "verbose",
        "model", "hkyREV", "aaDist","aaRatefile",
        "NSsites", "NShmm", "icode", "Mgene", "fix_kappa", "kappa",
        "fix_omega", "omega", "fix_alpha", "alpha","Malpha", "ncatG", 
        "fix_rho", "rho", "ndata", "Small_Diff"};
   double t;
   FILE  *fctl;
   char *daafiles[]={"", "grantham.dat", "miyata.dat", 
                     "g1974c.dat","g1974p.dat","g1974v.dat","g1974a.dat"};

   if((fctl=fopen(ctlf,"r"))==NULL) 
      error2("can't find the control file .ctl\n");
   if (noisy) printf ("\n\nReading options from %s..\n", ctlf);
   for (;;) {
      if (fgets (line, lline, fctl) == NULL) break;
      for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
         if (isalnum(line[i]))  { t=1; break; }
         else if (line[i]==comment) break;
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
               case ( 5): com.cleandata=(int)t;   break;
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
               case (17): com.NSsites=(int)t;     break;
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
      switch (com.model) {
      case (Poisson):  case (EqualInput): case (Empirical): case (Empirical_F):
         com.fix_kappa=1; com.kappa=0; com.nrate=0;   break;
      case (FromCodon): 
         com.nrate=!com.fix_kappa+!com.fix_omega;
         if (com.kappa<0) error2("kappa.."); 
         break;
      case (AAClasses): break;
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
      if(com.Mgene && com.model) error2("Mgene & model?");
      if(com.hkyREV && com.fix_kappa) error2("hkyREV & fix_kappa?");
      if(com.hkyREV && (com.aaDist || com.Mgene)) 
         error2("hkyREV to do: check options?");

      com.ncode=Nsensecodon[com.icode];
      if (com.model!=AAClasses) {
         if(com.fix_kappa>1) error2("fix_kappa>1, not tested.");  /**** ???? */
         if (com.model>0 && (com.alpha || !com.fix_alpha)) 
            error2("dN/dS ratios among branches not implemented for gamma");
         if (com.fix_omega) {
            OMEGA_FIX=com.omega;
            if((com.model==0 && com.NSsites==NSdiscrete)
               || (com.model && com.NSsites && com.NSsites!=NSselection
                   &&com.NSsites!=NSdiscrete&&com.NSsites!=NSbetaw))
               error2("\afix_omega?");
               
            if(com.aaDist) error2("fix_omega & aaDist");

         }
         if (com.model>NSbranch2) error2("seqtype or model.");
/*
         if (com.model==NSbranch2 && com.clock==2) 
            error2("NSbranch & local clock.");
*/
         if(com.kappa<0)  error2("kappa..");
         if(com.runmode==-2 && (com.NSsites||com.alpha||com.aaDist))
            error2("err: incorrect model for pairwise comparison.");
         if(com.runmode>0 && com.model==2) error2("tree search & model");
         if(com.aaDist && com.NSsites!=0 && com.NSsites!=NSdiscrete)
            error2("NSsites && aaDist.");

         com.nrate=com.nkappa=(com.hkyREV ? 5 : !com.fix_kappa);
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
            if((com.model && com.model!=2)
               || (com.model&&com.NSsites!=NSselection&&com.NSsites!=NSdiscrete))
               error2("only NSsites=2,3 & model=2 are compatible.");
            if(com.alpha || com.fix_alpha==0) error2("NSsites & gamma");
            if(com.NSsites<=NSselection) com.ncatG=com.NSsites+1;
            if(com.NSsites==NSbetaw || com.NSsites==NS02normal) com.ncatG++;
            if(com.NSsites==NSdiscrete && com.model==2) com.ncatG=4;
            if(com.NSsites==NSselection && com.model==2) com.ncatG=4;
            if(com.NSsites==NSdiscrete && !com.model && com.ncatG>3) {
               puts("\nncatG may be too big for M3 (discrete)?");
               sleep(2000);
            }
            if(com.NSsites==NSfreqs && com.ncatG!=5) 
               { puts("\nncatG changed to 5."); com.ncatG=5; }
            if(com.NSsites==NSselection || com.NSsites==NSbetaw)
               { if(!com.fix_omega) com.nrate++; else OMEGA_FIX=com.omega; }
            else if(com.NSsites==NSdiscrete) {
               com.nrate+=com.ncatG;  /* omega's */
               if(com.aaDist) { 
                  if (com.aaDist<=6) com.nrate+=com.ncatG;  /* a&b PSB2000 */
                  else {  /* fitness models */
                     com.nrate=!com.fix_kappa+4*com.ncatG;
                     if(com.aaDist==FIT2) com.nrate+=com.ncatG;
                  }
               }
            }
         }
      }
   }
   else
      error2 ("seqtype..");

   if(com.runmode==-2 && com.cleandata==0) {
      com.cleandata=1; puts("gaps are removed for pairwise comparison.");
   }
   if(com.method==1 &&(com.clock||com.rho)) 
      { com.method=0; puts("\aiteration method reset"); }

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

   return (0);
}


int testx (double x[], int np)
{
   int i,k;
   double tb[]={.4e-5, 29}, rgeneb[]={0.1,20}, rateb[]={1e-5,999}, omegab=.005;
   double alphab[]={0.005,99}, rhob[]={0.01,0.99};

   FOR (i,(com.clock?1:com.ntime))   if (x[i]<tb[0] || x[i]>tb[1]) return (-1);
   if(com.clock) 
      for(i=1;i<com.ntime;i++) if (x[i]<1e-6 || x[i]>1-1e-6) return (-1);
   if (np==com.ntime) return (0);
   FOR (i,com.nrgene) 
      if (x[com.ntime+i]<rgeneb[0] || x[com.ntime+i]>rgeneb[1]) return (-1);

   FOR (i,com.nrate) if (x[com.ntime+com.nrgene+i]>rateb[1])    return (-1);
   if (com.seqtype==CODONseq && com.nrate==2) {
      if(x[com.ntime+com.nrgene]<rateb[0] || x[com.ntime+com.nrgene+1]<omegab) 
         return (-1);
   }
   else 
      FOR (i,com.nrate)  if (x[com.ntime+com.nrgene+i]<rateb[0]) return (-1);

   k=com.ntime+com.nrgene+com.nrate;

   for (i=0; i<com.nalpha; i++,k++)
      if (x[k]<alphab[0] || x[k]>alphab[1]) return (-1);
   if (!com.fix_rho) if (x[np-1]<rhob[0] || x[np-1]>rhob[1]) return (-1);
   return(0);
}

int SetxInitials (double x[])
{
/* This forces initial values into the boundary of the space
*/
   int i, k;
   double tb[]={.0001,25}, rgeneb[]={0.1,9}, rateb[]={-5,89};
   double alphab[]={0.05,9}, rhob[]={0.01,0.9};

   FOR (i,com.ntime)  
      { if (x[i]>tb[1]) x[i]=tb[1];  if (x[i]<tb[0]) x[i]=tb[0]; }
   FOR (i,com.nrgene) {
      if (x[com.ntime+i]>rgeneb[1]) x[com.ntime+i]=rgeneb[1];
      if (x[com.ntime+i]<rgeneb[0]) x[com.ntime+i]=rgeneb[0];
   }
   for (i=0,k=com.ntime+com.nrgene; i<com.nrate; i++,k++) {
      if (x[k]>rateb[1]) x[k]=rateb[1];
      if (x[k]<rateb[0]) x[k]=rateb[0];
   }
   for (i=0; i<com.nalpha; i++,k++)  {
      if (x[k]>alphab[1]) x[k]=alphab[1];
      if (x[k]<alphab[0]) x[k]=alphab[0];
   }
   if (!com.fix_rho)
      { if (x[k]>rhob[1]) x[k]=rhob[1]; if (x[k]<rhob[0]) x[k]=rhob[0]; }

/*
   FOR (i,com.np) if (x[i]<1e-4) printf("x[%d]=%.6f\n", i,x[i]);
*/
   return(0);
}


int SetxBound (int np, double xb[][2])
{
   int i,j,k, K=com.ncatG;
   double tb[]={4e-6,50}, rgeneb[]={0.01,99}, rateb[]={1e-4,99};
   double alphab[]={0.005,99}, betab[]={.05,99}, omegab[]={.0001,99};
   double rhob[]={0.01,0.99}, pb[]={.000001,.999999};

   if(com.NSsites) omegab[0]=0.00001;
   if(com.clock) {
      xb[0][0]=.00005;  xb[0][1]=tb[1];
      for(i=1;i<tree.nnode-com.ns;i++)  /* for node times, clock=1 or 2 */
         FOR(j,2) { xb[i][0]=.001; xb[i][1]=.9999; }
      for(; i<com.ntime; i++)           /* for rates for branches, clock=2 */
         FOR(j,2) { xb[i][0]=tb[0]; xb[i][1]=tb[1]; }
   }
   else 
      FOR (i,com.ntime)  FOR (j,2) xb[i][j]=tb[j];

   for(i=com.ntime;i<np;i++) { xb[i][0]=rateb[0]; xb[i][1]=rateb[1]; }
   FOR (i,com.nrgene) FOR (j,2) xb[com.ntime+i][j]=rgeneb[j]; 
   FOR (i,com.nrate)  FOR (j,2) xb[com.ntime+com.nrgene+i][j]=rateb[j];
   k=com.ntime+com.nrgene+!com.fix_kappa; 
   if (com.NSsites) { /* p's before w's in xb[] */
      switch(com.NSsites) {
      case(NSneutral):  xb[k][0]=pb[0]; xb[k++][1]=pb[1];  break;  /* p0 */
      case(NSselection): /* for p0, p1, w2 */
         FOR(j,2) { xb[k][0]=-99; xb[k++][1]=99; }
         if(!com.fix_omega) { xb[k][0]=omegab[0]; xb[k++][1]=omegab[1]; } 
         break;  
      case(NSdiscrete):  /* pK[] & rK[] */
         if(com.model) K=3;
         if(com.nparK==0) 
            FOR(j,K-1) { xb[k][0]=-99; xb[k++][1]=99; }
         FOR(j,K) { xb[k][0]=omegab[0];  xb[k++][1]=omegab[1]; }
         if(com.nparK) 
            FOR(j,K*(K-1)) { xb[k][0]=-99; xb[k++][1]=99; }
         
         if(com.seqtype==CODONseq&&com.aaDist)
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
      case(NSbetaw):      /* p0, p_beta, q_beta, w */
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
   else if((com.seqtype==CODONseq||com.model==FromCodon)&&com.model!=AAClasses)
     { if(!com.fix_omega) { xb[k][0]=omegab[0]; xb[k][1]=omegab[1]; } }
   if(com.seqtype==CODONseq&&com.model)
      FOR(j,N_OMEGA-com.fix_omega) 
         { xb[k+j][0]=omegab[0]; xb[k+j][1]=omegab[1]; }

   if (com.aaDist<0 && (com.seqtype==1||com.model==FromCodon)) {
      /* linear relationship between d_ij and w_ij */
      if(com.nrate != !com.fix_kappa+1+(com.seqtype==1)) error2("in Setxbound");
      xb[com.ntime+com.nrgene+!com.fix_kappa][1]=1; /* 0<b<1 */
   }

   k=com.ntime+com.nrgene+com.nrate;
   for (i=0;i<com.nalpha;i++,k++)  FOR (j,2) xb[k][j]=alphab[j];
   if (!com.fix_rho)   FOR (j,2) xb[np-1][j]=rhob[j];
/*
   if(noisy>=3 && np<20) {
      puts("\nBounds:\n");
      FOR(i,np) printf(" %10.6f", xb[i][0]);  FPN(F0);
      FOR(i,np) printf(" %10.6f", xb[i][1]);  FPN(F0);
   }
*/
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


int GetInitials(double x[], int* fromfile)
{
/* This caculates the number of parameters (com.np) and get initial values.
   Perhaps try to restruct the code and make two sections for amino acids 
   and codons?
   com.nrate is initialised in getoptions().
*/
   static int times=0;
   int i, j,k=0, naa=20;
   int slkl1_new=(tree.nnode-com.ns)*com.ncode*com.npatt*sizeof(double);
   double t;
   
   NFunCall=NPMatUVRoot=NEigenQ=0;
   com.plfun = (com.alpha==0 ? lfun : (com.rho==0?lfundG:lfunAdG));
   if(com.NSsites) com.plfun=lfundG;
   if(com.nparK) com.plfun=lfunAdG;

   if(com.method==1 && com.fix_branch!=2 && com.plfun==lfundG) {
      lklSiteClass=1;
      slkl1_new*=com.ncatG;
   }
   if(com.slkl1<slkl1_new) {
      com.slkl1=slkl1_new;
      printf("\n%9d bytes for lkl1, adjusted\n",com.slkl1);
      if((com.lkl=(double*)realloc(com.lkl,com.slkl1))==NULL) error2("oom lkl1");
   }
   
   InitializeNodeScale();

   if(times++==0) {
      if((com.aaDist && com.aaDist<10 &&
          (com.seqtype==CODONseq||com.model==FromCodon)) ||
          (com.seqtype==AAseq &&
          (com.model==Empirical||com.model==Empirical_F||com.model>=REVaa_0))){
         GetDaa(NULL,com.daa);
      }
   }
   if(com.fix_branch==2)  { com.ntime=0; com.method=0; }
   else if(com.clock==0)  com.ntime=tree.nbranch;
   else if(com.clock==1)  com.ntime=tree.nnode-com.ns+(tree.root<com.ns);
   else if(com.clock==3)  com.ntime=tree.nnode-com.ns+(tree.root<com.ns)+1;
   else {  /* if(com.clock==2) */
      for(i=0,k=0,N_rateBranch=0; i<tree.nbranch; i++) { 
         rateBranch[i]=j=nodes[tree.branches[i][1]].label;
         if(j+1>N_rateBranch)  N_rateBranch=j+1;
         if(j<0||j>tree.nbranch-1) error2("branch label in the tree.");
         else if (j)   k=1;
      }
      if (k)
         puts("\nUsing branch marks in the tree file. Stop if wrong"); 
      else {
         OutaTreeB(F0); FPN(F0);
         for (i=0,N_rateBranch=0; i<tree.nbranch; i++) {
            printf("Branch %2d: %2d..%-2d? ", i+1,
               tree.branches[i][0]+1, tree.branches[i][1]+1);
            scanf("%d", &rateBranch[i]);
            if(rateBranch[i]<0 || rateBranch[i]>=tree.nbranch) 
               error2("rates for branches.");
            if(rateBranch[i]+1>N_rateBranch) N_rateBranch=rateBranch[i]+1;
         }
      }
      printf("\n\n%d rates for branches. ",N_rateBranch);
      com.ntime=tree.nnode-com.ns + N_rateBranch-1;
      if(com.ntime>tree.nbranch) 
         printf("\a\nntime=%d, too many rates??\n",com.ntime);
   }

   if (com.clock) {
      for(j=1,x[0]=1; j<tree.nnode-com.ns; j++)
         x[j]=0.2+0.8*rndu();
      for(; j<com.ntime; j++) x[j]=1;   /* rates for branches */
      if(com.clock==3) SetDivTime_OldSon(tree.root);
   }
   else if(com.fix_branch==0) {
      FOR (j,com.ntime) x[j]=.2*rndu();
      if(com.ns<100) LSDistance(&t, x, testx);
   }
   com.nrgene=(!com.fix_rgene)*(com.ngene-1);
   FOR(j,com.nrgene) x[com.ntime+j]=1;

   /* if(com.seqtype==CODONseq && com.NSsites!=NSdiscrete) */
   /* problematic line: check carefully */
/*
   if(com.seqtype==CODONseq && com.NSsites==0 && com.aaDist==0)
      com.nrate=!com.fix_kappa+!com.fix_omega;
*/
   if(com.seqtype==CODONseq && com.model && com.model<=NSbranch2) {
      if (com.model==NSbranch2) {
         for(i=0,k=0,N_OMEGA=0; i<tree.nbranch; i++) { 
            j=nodes[tree.branches[i][1]].label;
            if(j+1>N_OMEGA)  N_OMEGA=j+1;
            if(j<0||j>tree.nbranch-1) 
               error2("branch label in the tree");
            else if (j)   k=1;
         }
         if (k)
            puts("\nUse branch marks in the tree file. Stop if wrong");
         else {
            printf("\nSpecify dN/dS ratios for all %d branches (0,1,2,etc.).\n",
                tree.nbranch);
            OutaTreeN(F0,0,0);  FPN(F0);  OutaTreeB(F0);  FPN(F0);
            for (i=0,N_OMEGA=0; i<tree.nbranch; i++) {
               printf("Branch %2d: %2d..%-2d? ", i+1,
                  tree.branches[i][0]+1,tree.branches[i][1]+1);
               scanf("%d",&j);

               if(com.NSsites && (j!=0 && j!=1))
                  error2("only labels 0 and 1 are allowed for the model.");
               nodes[tree.branches[i][1]].label=j;
               if (j+1>N_OMEGA)  N_OMEGA=j+1;
            }
         }
         if(N_OMEGA==1) error2("Only one w ratio, use model=0.");
         printf ("\n%d ratios are specified. Stop if wrong.", N_OMEGA);
      }
      else {
         N_OMEGA=tree.nbranch;  
         FOR(i,tree.nbranch) nodes[tree.branches[i][1]].label=i;
      }
      com.nrate=(com.hkyREV ? 5 : !com.fix_kappa)  +!com.fix_omega+N_OMEGA-1;
   }    /* if(dN/dS models for branches */

   if(com.Mgene>=3) com.nrate*=com.ngene;
   if(com.seqtype==1 && com.Mgene>=3 && com.fix_omega) com.nrate--; 
   com.np = com.ntime+com.nrgene+com.nrate;

   /* NSbranchsite models 
      N_OMEGA=2 different w's at a site (three w's in the model: w0,w1,w2) */
   if(com.seqtype==CODONseq && com.model && com.NSsites) {
      if(N_OMEGA!=2) error2("only two branch labels are allowed");
      com.ncatG=4;
      if(com.NSsites==NSdiscrete) 
         com.nrate = (com.hkyREV ? 5 : !com.fix_kappa) 
                   + 2 +!com.fix_omega+N_OMEGA-1-1; /* add w0 and w1 */
      else 
         com.nrate= (com.hkyREV ? 5 : !com.fix_kappa) +!com.fix_omega+N_OMEGA-1-1;
      com.np = com.ntime+com.nrgene+com.nrate + 2;  /* add p0 and p1 */

      k=com.ntime+com.nrgene;
      if(com.hkyREV) {
         x[k++]=1+0.1*rndu(); 
         x[k++]=rndu(); x[k++]=rndu(); x[k++]=rndu(); x[k++]=rndu();
      }
      else if(!com.fix_kappa)  x[k++]=com.kappa;
      x[k++]=2.2; x[k++]=1.1;      /* p0 and p1 */
      if(com.NSsites==NSdiscrete) { x[k++]=0.1; x[k++]=0.8; } /* w0 and w1 */
      if(!com.fix_omega)  x[k++]=3.5;  /* w2 */

   }
   else if (com.NSsites) {        /* w's are counted in com.nrate */
      k=com.ntime+com.nrgene;  
      if(!com.fix_kappa)  x[k++]=com.kappa;
      switch(com.NSsites) {
      case(NSneutral):   com.np++; x[k++]=.7;  break;  /* p0 for w0=0 */
      case(NSselection): /* for p0, p1.  w is counted in nrate.  */
         com.np+=2; x[k++]=1.3; x[k++]=.3; 
         if(!com.fix_omega) x[k++]=com.omega;
#ifdef NSSITESBATCH
         if(!com.fix_omega) x[k-1]=max2(2,com.omega);
#endif
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
            FOR(i,com.ncatG-1) x[k++]=0.;
            FOR(i,com.ncatG) x[k++]=com.omega;
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
      case(NSbetaw):        /* p0, p_beta, q_beta.  w is counted in nrate. */
         com.np+=3; x[k++]=.9; x[k++]=.4; x[k++]=1.2+rndu();
         if(!com.fix_omega) x[k++]=com.omega;
#ifdef NSSITESBATCH
         if(!com.fix_omega) x[k-1]=max2(2,com.omega);
#endif
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
   if (com.model==AAClasses) { 
      if (!com.fix_kappa) x[k++]=com.kappa;
      FOR (i,com.nrate-!com.fix_kappa) x[k++]=com.omega;
      if (com.nrate>65) puts("\a\nget better initial values for AAclasses?");
       /* if estimating all acceptance rates */
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
            if(!com.fix_kappa)  x[k++]=com.kappa;
            FOR(j,com.nrate-!com.fix_kappa)  x[k++]=com.omega; 
         }
      }
      else if (!com.NSsites) {                      /* CODONseq */
         if (com.nrate==0 && com.NSsites==0) 
            EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, 
               (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),com.omega,PMat);
         else if (com.model==NSbranchB || com.model==NSbranch2) {
            if (!com.fix_kappa)  x[com.ntime+com.nrgene]=com.kappa; 
            FOR(i, (com.model==NSbranchB?tree.nbranch:N_OMEGA-1+!com.fix_omega))
               x[com.ntime+com.nrgene+!com.fix_kappa+i]=com.omega;
         }
         else if(com.nrate) { /* either kappa, omega, or both for each gene */
            if(com.Mgene<=2) {
               if(com.hkyREV) { x[k++]=1; FOR(i,4) x[k++]=.5; }
               else if (!com.fix_kappa) x[k++]=com.kappa;
               if (!com.fix_omega)      x[k++]=com.omega;
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

   *fromfile=0;
   if(fin) {
      readx(x,fromfile);
      if(com.runmode>0 && fromfile && com.NSsites)  LASTROUND=1;
   }

   return (0);
}



int SetPGene (int igene, int _pi, int _UVRoot, int _alpha, double xcom[])
{
/* xcom[] does not contain time parameters
   Note that com.piG[][] have been homogeneized if (com.Mgene==3)
   Note calculation of nr1 for (com.Mgene>=3 && com.fix_omega), as only the 
   w for the last partition is fixed.
*/
   int j,nr1=(com.nrate+1)/com.ngene, k=com.nrgene+(com.Mgene>=3)*igene*nr1;
   double small=.1/(com.ns*com.ls);

   if (_pi) {
      xtoy (com.piG[igene],com.pi,com.ncode);

      for(j=0,com.npi0=0; j<com.ncode; j++)   
         if(com.pi[j]<small) com.npi0++;
      if(com.npi0==0)
         for(j=0; j<com.ncode; j++)  com.pi_sqrt[j]=sqrt(com.pi[j]);
   }
   if (_UVRoot) {
      if (com.seqtype==CODONseq) {
         if(!com.fix_kappa) com.kappa=xcom[k++];
         if(!com.fix_omega) com.omega=xcom[k++];
         else
            com.omega=(com.Mgene>2&&igene<com.ngene-1?xcom[k++]:OMEGA_FIX);
         if (!com.NSsites)
            EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, 
               (com.hkyREV?&xcom[com.nrgene]:&com.kappa),com.omega,PMat);
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



int SetParameters (double x[])
{
/* Set com. variables and initialize U, V, Root etc. before each calculation 
   of the likelihood function.
   Is it a good idea to restruct this and/or Getinitials into two parts,
   one for aa's and another for codons?
   When (com.NSsites==NS02normal || NS3normal), p's are before w's in x[]; 
   see CDFdN_dS().
*/
   int K=com.ncatG, i,j,k;
   double t, S,dS,dN, w0,w1,w2, spaceP2PI[NCATG*(NCATG+1)];

   if(com.fix_branch<2) SetBranch(x);
/*
OutaTreeN(F0,0,1);
*/

   if(com.np<=com.ntime) return(0);

   if(com.seqtype==1 || com.model==6)
      POMEGA=x+com.ntime+com.nrgene+!com.fix_kappa;
   FOR(j,com.nrgene) com.rgene[j+1]=x[com.ntime+j];
   if(com.seqtype==1 && com.aaDist>=FIT1) getpcodonClass(x, pcodonClass);

   k=com.ntime+com.nrgene;
   PKAPPA=x+k;
   if (com.nrate) {
      if(!com.fix_kappa && com.hkyREV==0) com.kappa=x[k];
      k+=com.nkappa;
      if(!com.model && !com.aaDist && !com.fix_omega) com.omega=x[k];
      if(com.seqtype==AAseq)
         EigenQaa(NULL, Root, U, V, x+com.ntime+com.nrgene);
      else if((com.model==0 && com.NSsites==0 && com.Mgene<=1) || com.model==AAClasses)
         /* CODONs, same dN/dS across branches & sites */ 
         EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, 
            (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),com.omega,PMat);
      else if (com.model==2 && com.NSsites==0 && N_OMEGA<=5)
         FOR(i,N_OMEGA) {
            w0=(i==N_OMEGA-1&&com.fix_omega?OMEGA_FIX:POMEGA[i]);
            EigenQc(0,-1,NULL,NULL,NULL,_Root[i],_UU[i],_VV[i],
               (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),w0,PMat);
         }

      k=com.ntime+com.nrgene+com.nrate;
   }

   if (com.NSsites) {   /* p's before w's in x[] */
      k=com.ntime+com.nrgene+!com.fix_kappa;
      switch(com.NSsites) {
      case(NSneutral):  
         com.rK[0]=0; com.rK[1]=1;
         com.freqK[0]=x[k++]; com.freqK[1]=1-com.freqK[0]; break;
      case(NSselection): 
      case(NSdiscrete):
         if(com.model) K=3;  /* note this. */
         if(com.nparK) {
            FOR(j,K) com.rK[j]=x[k++];  /* w's for site classes */
            for (i=0; i<K; i++, k+=K-1) {
               if (!LASTROUND) f_to_x(x+k,com.MK+i*K,K,0,0);
               else            xtoy  (x+k,com.MK+i*K,K-1);
               com.MK[i*K+K-1]=1-sum(com.MK+i*K,K-1);
            }
            PtoPi(com.MK, com.freqK, K, spaceP2PI);
         }
  
      case(NSfreqs):
         if (!LASTROUND && !com.nparK) {
            for(j=0,t=1;j<K-1;j++) t+=(com.freqK[j]=exp(x[k++]));
            FOR(j,K-1) com.freqK[j]/=t; com.freqK[K-1]=1/t;
         }
         else if(!com.nparK) {
            for(j=0,com.freqK[K-1]=1;j<K-1;j++) 
               com.freqK[K-1]-=(com.freqK[j]=x[k++]);
            if(com.freqK[K-1]<0 || com.freqK[K-1]>1) error2("freqK[]");
         }

         if(com.NSsites==NSselection) {
            com.rK[0]=0; com.rK[1]=1; 
            com.rK[2]=(com.fix_omega?OMEGA_FIX:x[k++]); 
         }
         else if(com.NSsites==NSfreqs) {
            if(com.ncatG!=5) error2("NSfreqs, ncatG?");
            com.rK[0]=0; com.rK[1]=1./3; com.rK[2]=2./3; com.rK[3]=1; com.rK[4]=3;
         }
         else if(com.aaDist==0 && com.nparK==0)
            /* NSdiscrete: note this sets w0 & w1 for NSbranchsite */
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
      /* calculates Qfactor_NS, to be used in eigenQc for NSsites models */
      if(com.model==0) {
         for(j=0,Qfactor_NS=0; j<com.ncatG; j++) {
            freqK_NS=com.freqK[j];
            if(com.aaDist) {
               k=com.ntime+com.nrgene+!com.fix_kappa;
               if(com.aaDist<10)         POMEGA=x+k+com.ncatG-1+2*j;
               else if(com.aaDist>=FIT1) {
                  POMEGA=x+k+com.ncatG-1+j*(4+(com.aaDist==FIT2));
                  xtoy(pcodonClass+j*64, com.pi, com.ncode);
               }
            }
            EigenQc(1,-1,&S,&dS,&dN,NULL,NULL,NULL,
               (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),com.rK[j],PMat);
         }
         Qfactor_NS=1/Qfactor_NS;
      }
   }  /* if(com.NSsites) */


   if (com.NSsites && com.model) { /* branch&site models */
      k=com.ntime+com.nrgene+!com.fix_kappa;
      t=com.freqK[0]+com.freqK[1];
      com.freqK[2]=(1-t)*com.freqK[0]/(com.freqK[0]+com.freqK[1]);
      com.freqK[3]=(1-t)*com.freqK[1]/(com.freqK[0]+com.freqK[1]);
      com.rK[2]=com.rK[0]; com.rK[3]=com.rK[1];  /* w0 and w1 */
      /* calculates scale factors: background branches has two site classes
         while foreground branches has 3 site classes */
      FOR(i,2) {
         for(j=0,Qfactor_NS=0; j<(i==0?2:3); j++) {
            w0=com.rK[j];
            if(i==0)       freqK_NS=com.freqK[j]/t;
            else if(j==2) {
               freqK_NS=1-t; 
               if(com.fix_omega) w0=OMEGA_FIX;
               else              w0=(com.NSsites==NSselection?x[k+2]:x[k+2+2]); 
            }
            EigenQc(1,-1,&S,&dS,&dN,NULL,NULL,NULL,
               (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),w0,PMat);
         }
         Qfactor_NS_branch[i]=1/Qfactor_NS;
      }
      /* This calculates 5 sets of U&V&Root vectors (w0b,w0f,w1b,w1f,w2f), 
         which are used in SetPSiteClass():
            iclass=0: w0b,w0f 
            iclass=1: w1b,w1f 
            iclass=2: w0b,w2f 
            iclass=3: w1b,w2f
         no eigenQc() calls are needed in PartialLikelihood() or minbranches().
      */
      k=com.ntime+com.nrgene+!com.fix_kappa+2;
      if(com.NSsites==NSselection) { w0=0; w1=1; }
      else                         { w0=x[k++]; w1=x[k++]; }
      w2=(!com.fix_omega?x[k++]:OMEGA_FIX);
      for(i=0; i<5; i++) {  /* (w0b,w0f,w1b,w1f,w2f) */
         Qfactor_NS = Qfactor_NS_branch[(i==1||i==3||i==4)];
         com.omega = (i<2 ? w0 : (i<4?w1:w2));
         EigenQc(0,-1,NULL,NULL,NULL,_Root[i],_UU[i],_VV[i],
            (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),com.omega,PMat);
      }
   }  /* branch&site models */

   if(com.model) com.omega=-1;  /* to force crash in case or error2 */
   if(com.seqtype==CODONseq && !com.NSsites && com.model) {
      FOR(j,tree.nnode) {
         if(j==tree.root) continue;
         if (com.fix_omega && nodes[j].label==N_OMEGA-1)
            nodes[j].omega=OMEGA_FIX;
         else
            nodes[j].omega=POMEGA[nodes[j].label];
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
      else               com.rK[com.ncatG-1]=OMEGA_FIX;
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
         FOR(j,(com.NSsites==7?2:4)) printf(" %12.5e", par[j]);  FPN(F0);
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


int EigenQc (int getstats, double branchl, double *S, double *dS, double *dN,
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
   c[0-2] are rates for the 3 codon positions, useful for calculating ratios
   like c[1]/c[0].
   k3[0,1,2] are the kappa's for the three codon positions expected under
   the codon model.

   Under NSsites or other site-class models, this function does not scale 
   Q but calculates the Qfactor_NS.
   DetailOutput() uses this function to calculate dS & dN; for each omega 
   component, dS and dN are not scaled but are scaled in DetailOutput() 
   after Qfactor_NS is calculated.

   aaDist=FIT1 & FIT2:  ap,p*,av,v*, (and w0 for FIT2)
*/
   int n=Nsensecodon[com.icode], i,j,k, ic1,ic2,aa1,aa2, b1,b2;
   int ndiff,pos=0,from[3],to[3];
   double rs0,ra0,rs,ra, c0[3],c[3],ts[3],tv[3], t,space[64*3], *ri=space;
   double *pi=(com.seqtype==AAseq?com.fb61:com.pi);
   double w, fit1,fit2;

   NEigenQ++;
   if(branchl>=0 && (S==NULL||dS==NULL||dN==NULL)) error2("EigenQc");
   FOR (i,3) c[i]=c0[i]=0;  FOR(i,3) ts[i]=tv[i]=0;
   FOR (i,n*n) Q[i]=0;
   for (i=0, rs0=ra0=rs=ra=0; i<n; i++) FOR (j,i) {
      ic1=FROM61[i]; from[0]=ic1/16; from[1]=(ic1/4)%4; from[2]=ic1%4;
      ic2=FROM61[j];   to[0]=ic2/16;   to[1]=(ic2/4)%4;   to[2]=ic2%4;
      for (k=0,ndiff=0; k<3; k++)  if (from[k]!=to[k]) { ndiff++; pos=k; }
      if (ndiff!=1)  continue;
      t=2*pi[i]*pi[j];
      Q[i*n+j]=1;

      if(com.hkyREV) {
         b1=min2(from[pos],to[pos]);  
         b2=max2(from[pos],to[pos]);
         if     (b1==0 && b2==1) Q[i*n+j]=kappa[0];
         else if(b1==0 && b2==2) Q[i*n+j]=kappa[1];
         else if(b1==0 && b2==3) Q[i*n+j]=kappa[2];
         else if(b1==1 && b2==2) Q[i*n+j]=kappa[3];
         else if(b1==1 && b2==3) Q[i*n+j]=kappa[4];
      }
      else
         if(from[pos]+to[pos]==1 || from[pos]+to[pos]==5)
            Q[i*n+j]= *kappa;

      c0[pos]+=t*Q[i*n+j];

      aa1=GeneticCode[com.icode][ic1]; 
      aa2=GeneticCode[com.icode][ic2];

      if(aa1!=aa2){
         ra0+=t*Q[i*n+j];
         if (com.model==AAClasses) {
            if (aa1<aa2)  { k=aa2; aa2=aa1; aa1=k; }
            k=aa1*(aa1-1)/2+aa2;
            if (POMEGA[OmegaAA[k]]<0) {
               if (noisy)  printf ("aa1 & aa2 & iw & w: %d %d %d %.5f\n", 
                              aa1,aa2,OmegaAA[k],POMEGA[OmegaAA[k]]);
               POMEGA[OmegaAA[k]]=0;
            }
            if (com.seqtype==AAseq && com.nrate>65 && aa1*20+aa2==ijAAref)
                ;     /* if estimating grantham's matrix with aa sequences */
            else  Q[i*n+j] *= POMEGA[OmegaAA[k]];
         }
         else if (com.aaDist==0)  Q[i*n+j] *= omega; /* NSsites==0 or >0 */
         else if (com.aaDist<=6)  {        /* chemical properties: a & b */
            /* w = POMEGA[0]*com.daa[aa1*20+aa2]*com.daa[aa1*20+aa2]; */
            w = POMEGA[0]*com.daa[aa1*20+aa2];
            if(com.aaDist>0)           Q[i*n+j] *= exp(-w);  /* geometric */
            else                       Q[i*n+j] *= 1-w;      /* linear */
            if (com.seqtype==CODONseq) Q[i*n+j] *= POMEGA[1];
         }
         else if (com.aaDist>=FIT1) {   /* ap,p*,av,v* (and w0 for FIT2) */
            fit1=-POMEGA[0]*square(AAchem[0][aa1]-POMEGA[1])
                 -POMEGA[2]*square(AAchem[1][aa1]-POMEGA[3]);
            fit2=-POMEGA[0]*square(AAchem[0][aa2]-POMEGA[1])
                 -POMEGA[2]*square(AAchem[1][aa2]-POMEGA[3]);

            Q[i*n+j] *= exp(-fit1-fit2);
            if(com.aaDist==FIT2) Q[i*n+j]*=POMEGA[4];
         }
         ra+=t*Q[i*n+j];
      }
      else 
         rs+=t*Q[i*n+j];

      if(from[pos]+to[pos]==1 || from[pos]+to[pos]==5) ts[pos]+=t*Q[i*n+j];
      else                                             tv[pos]+=t*Q[i*n+j];

      c[pos]+=t*Q[i*n+j];
      Q[j*n+i]=Q[i*n+j];

   }  /* for (i,j) */

   if(getstats) {
      if(com.NSsites) Qfactor_NS+=freqK_NS * (rs+ra);
      /* if(noisy>2) FOR(i,3) printf("\nts/tv[pos %d] = %.6f", i+1, ts[i]/tv[i]);*/
      rs0=rs;
      t=(rs0+ra0);  rs0/=t;  ra0/=t;   *S=rs0*3*com.ls;
      if(!com.NSsites && branchl>=0) {  /* calculates dS & dN */
         if(branchl==0) *dS = *dN = 0;
         t=(rs+ra); rs/=t; ra/=t;
         *dS=branchl*rs/(3*rs0);  *dN=branchl*ra/(3*ra0);
      }
      else if (com.NSsites)
         { *dS=rs/(3*rs0);  *dN=ra/(3*ra0); }

/*      
FPN(F0); FOR(i,3) printf(" c0[%d]=%.5f ",i+1,c0[i]); FPN(F0);
         FOR(i,3) printf("  c[%d]=%.5f ",i+1, c[i]); FPN(F0); FPN(F0);
*/

   }
   else {  /* get Root, U, & V */
      if (com.seqtype==AAseq) return (0);

      FOR(i,n) FOR(j,n) Q[i*n+j]*=pi[j];
      for (i=0,t=0; i<n; i++) 
        { Q[i*n+i]=-sum(Q+i*n,n); t-=pi[i]*Q[i*n+i]; }
      if(com.NSsites)  t=1/Qfactor_NS;  /* Qfactor_NS calculated in SetParameter */

      if(com.npi0==0)
         eigenQREV(Q,com.pi,com.pi_sqrt,n,Root,U,V);
      else {
         if (eigen (1,Q,n,Root,ri,U,V,space+n)) error2("eigenQc err.");
         xtoy (U, V, n*n);
         matinv (V, n, n, space);
      }
      FOR(i,n) Root[i]/=t;
      UVRootChanged=1;
   }
   return (0);
}


int EigenQaa (FILE *fout, double Root[], double U[], double V[], double rate[])
{
   int naa=20, i,j,k;
   double Q[20*20], mr, t=0, space[20*3];
   char aa3[4]="";

   FOR (i,naa) Q[i*naa+i]=0;
   switch (com.model) {
   case (Poisson)   : case (EqualInput) : 
      fillxc (Q, 1., naa*naa);  break;
   case (Empirical)   : case (Empirical_F):
      FOR(i,naa) FOR(j,i) Q[i*naa+j]=Q[j*naa+i]=com.daa[i*naa+j]/100;
      break;
   case (FromCodon): case (AAClasses): 
      EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, 
         (com.hkyREV?rate:&com.kappa),com.omega,PMat);
      Qcodon2aa(PMat, com.fb61, Q, space);
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

   if(com.npi0==0) 
      eigenQREV(Q, com.pi, com.pi_sqrt, naa, Root, U, V);
   else {
      if (eigen (1, Q, naa, Root, space, U, V, space+naa)) 
         error2("eigenQc err.");
      xtoy (U, V, naa*naa);
      matinv (V, naa, naa, space);
   }
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
   int i, j, aai, aaj, nc=Nsensecodon[com.icode], naa=20;
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



int PartialLikelihood (int inode, int igene)
{
   int n=com.ncode, i,j,k,h, ison, pos0=com.posG[igene],pos1=com.posG[igene+1];
   double t;

   FOR(i,nodes[inode].nson)
      if(nodes[nodes[inode].sons[i]].nson>0 && !_oldlkl[nodes[inode].sons[i]])
         PartialLikelihood(nodes[inode].sons[i], igene);

   fillxc (nodes[inode].lkl+pos0*n, (double)(inode>=com.ns), (pos1-pos0)*n);

   if (com.cleandata && inode<com.ns)
      for(h=pos0;h<pos1;h++) nodes[inode].lkl[h*n+com.z[inode][h]]=1;

   FOR(i,nodes[inode].nson) {
      ison=nodes[inode].sons[i];
      t=nodes[ison].branch*com.rgene[igene]*_rateSite;
      GetPMatBranch(PMat, PKAPPA, t, ison);

      if(com.cleandata && nodes[ison].nson<1)  /* tips */
         for(h=pos0; h<pos1; h++)
            FOR(j,n) nodes[inode].lkl[h*n+j]*=PMat[j*n+com.z[ison][h]];
      else {
         for(h=pos0; h<pos1; h++) 
            FOR(j,n) {
               for(k=0,t=0; k<n; k++)
                  t+=PMat[j*n+k]*nodes[ison].lkl[h*n+k];
               nodes[inode].lkl[h*n+j]*=t;
            }
      }
   }        /*  for (ison)  */
   if(_nnodeScale && _nodeScale[inode]) 
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
/* This sets FROM61[], which goes from 0, 1, ..., 61, and FROM64[], 
   which goes from 0, 1, ..., to 64.
*/
   int i,n;

   for (i=0,n=0; i<64; i++) {
      if (GeneticCode[com.icode][i]==-1)  FROM64[i]=-1; 
      else                              { FROM61[n]=i; FROM64[i]=n++; }
   }
   if (n!=Nsensecodon[com.icode])  error2 ("# sense codons?");
   return (0);
}


int printfcode (FILE *fout, double fb61[], double space[])
{
/* space[64*2]
*/
   int i, n=Nsensecodon[com.icode];

   fprintf (fout, "\nCodon freq.,  x 10000\n");
   zero (space, 64);
   FOR (i, n) space[FROM61[i]] = fb61[i]*10000;
   printcu (fout, space, com.icode);
   return (0);
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


double Fcodon_3x4(double fcodon[], double fb3x4[]);
void OutFb3x4 (FILE*fout, double fb3x4[]);
void OutCodonCounts (FILE *fout,double fcodonsg[],double fb3x4sg[],double fb4g[],int *miss);

double Fcodon_3x4(double fcodon[], double fb3x4[])
{
/* this converts the codon frequencies into a fb3x4 table. fcodon has 64 codons.
*/
   int b[3], k,j, nc=64;
   double t;

   zero(fb3x4,12);
   for(k=0; k<nc; k++) {
      b[0]=k/16; b[1]=(k%16)/4; b[2]=k%4;
      FOR(j,3)  fb3x4[j*4+b[j]]+=fcodon[k];
   }
   FOR(j,3)  {
      t=sum(fb3x4+j*4, 4);
      if(t<1e-15) error2("Fcodon_3x4: empty sequence?");
      abyx(1/t, fb3x4+j*4, 4);
   }
   return(0);
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
   Resolve an ambiguity codon (if !cleandata) into nb[k] nucleotides for 
   each codon position k.  
   If(cleandata) nb[0]=nb[1]=nb[2]=1, and the nucleotides are in ib[k][].
*/
   int k,ic;

   if(com.cleandata) {
      nb[0]=nb[1]=nb[2]=1;
      ic=FROM61[codon[0]];
      ib[0][0]=ic/16;  ib[1][0]=(ic/4)%4;   ib[2][0]=ic%4;
   }
   else
      FOR(k,3)
         NucListall(codon[k], &nb[k], ib[k]);
   return(0);
}

void OutCodonCounts (FILE *fout,double fcodonsg[],double fb3x4sg[],double fb4g[],int *miss)
{
/* Outputs codon counts and f3x4 tables, called from InitializeCodon(), where 
   more notes are found.
*/
   int h, j,k, nc=NCODE, ig, wname=15, nb[3],ib[3][4],ic;
   char n31=(char)(com.cleandata?1:3);

   /* counts codons for output, species first, genes next */
   fputs("Codon usage in sequences\n",fout);
   zero(fcodonsg, com.ns*nc);
   FOR(j,com.ns) {
      for(h=0; h<com.npatt; h++) {
         CodonListall(com.z[j]+h*n31, nb, ib);
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
            CodonListall(com.z[j]+h*n31, nb, ib);
            k=nb[0]*nb[1]*nb[2];
            if(k>1) continue;
            ic=ib[0][0]*16+ib[1][0]*4+ib[2][0];
            fcodonsg[ig*nc+ic]+=com.fpatt[h];
         }
      }
      Fcodon_3x4(fcodonsg+ig*nc, fb3x4sg+ig*12);
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
   char n31=(com.cleandata?1:3);
   double t,t1;
   char str[4]="   ", ft[64];

   for(h=com.posG[ig]; h<com.posG[ig+1]; h++) {
      CodonListall(com.z[js]+h*n31, nb, ib);
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
      if(k==0) printf("stop codon %s in seq. %d\n", getcodon(str,ic),js+1);
      if(t<1e-20) puts("difficulty in resolving codon..");
      FOR(k,nc) if(ft[k]) 
         fcodon[k]+=com.fpatt[h]*fcodon0[k]/t;
   }
}


int InitializeCodon (FILE *fout, double space[])
{
/* Count codons for genes, calculate site patterns and fpatts 
   sequences com.z[] are not coded and may contain ambiguity characters
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
   double small=.1/(double)com.ns/com.ls;

   PatternWeight (fout,space);

   /* counts codons for output, species first, genes next */
   if(noisy) puts("Counting codons..");
   OutCodonCounts (fout, fcodonsg, fb3x4sg, fb4g, &miss);

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
               AddCodonFreqSeqGene (j, ig, fcodon0, fcodonsg+com.ngene*nc, 
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
      else if (com.codonf==Fcodon) {
         FOR(k,nc)  if(FROM64[k]>-1)  ppi[FROM64[k]]=fcodonsg[ig*nc+k]; 
      }
      else if (com.codonf==F3x4) {
         FOR(k,nc)  if(FROM64[k]>-1)
            ppi[FROM64[k]]=fb3x4sg[ig*12+k/16]*fb3x4sg[ig*12+4+(k/4)%4]*fb3x4sg[ig*12+8+k%4];
      }
      else if (com.codonf==F1x4) {
         FOR(k,nc)  if(FROM64[k]>-1)
            ppi[FROM64[k]]=fb4g[ig*4+k/16]*fb4g[ig*4+(k/4)%4]*fb4g[ig*4+k%4];
      }
      abyx(1/sum(ppi,com.ncode),ppi,com.ncode);  /* ncode != nc */
   }
   for(k=0,com.npi0=0; k<com.ncode; k++)   if(com.pi[k]<small)  com.npi0++;
   if(com.npi0==0) 
      for(k=0; k<com.ncode; k++)  com.pi_sqrt[k]=sqrt(com.pi[k]);

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
   for(ic=0; ic<Nsensecodon[com.icode]; ic++) {
      iaa=GeneticCode[com.icode][FROM61[ic]];
      fcodon[ic]+=faa[iaa]/NCsyn[iaa];
   }
   if(fabs(1-sum(fcodon,64))>1e-6) printf("\n1 == %12.7f\n", sum(fcodon,64));

   return (0);
}


int DistanceMatAA (FILE *fout)
{
   int i,j, h;
   double p;

   if(fout) fprintf(fout,"\nAA distances (raw proportions of different sites)\n");
   FOR(i, com.ns) {
      if(fout) fprintf (fout, "\n%-15s", com.spname[i]);
      FOR(j,i) {
         for(h=0,p=0; h<com.npatt; h++)  
            if (com.z[i][h]!=com.z[j][h]) p+=com.fpatt[h];
         p /= com.ls;
         SeqDistance[i*(i-1)/2+j]=p;
         if(fout) fprintf(fout, " %7.4f", p);
      }
   }
   if (fout) FPN(fout);
   return (0);
}

int DistanceMatNG86 (FILE *fout, double alpha)
{
/* Estimation of dS and dN by the method of Nei & Gojobori (1986)
   This works with both coded (com.cleandata==1) and uncoded data.
   In the latter case (com.cleandata==0), the method does pairwise delection.

   alpha for gamma rates is used for dN only.
*/
   FILE *fds,*fdn,*ft;
   char dsf[32]="2NG.dS",dnf[32]="2NG.dN",tf[32]="2NG.t";
   char codon[2][3]={"   ", "   "};
   int is,js, k,i0,h, wname=20, status=0, ndiff,nsd[4];
   int nb[3],ib[3][4], missing;

   double ns,na, nst,nat, S,N, St,Nt, dS,dN,dN_dS,y, bigD=3, ls1;

   fputs("\n\n\nNei & Gojobori 1986. dN/dS (dN, dS)",fout);
   if(com.cleandata==0) fputs("\n(Pairwise deletion)",fout);

   fputs("\n(Note: This matrix is not used in later m.l. analysis.\n",fout);
   fputs("Use runmode = -2 for ML pairwise comparison.)\n",fout);

   fputs("\nNumber of codon sites with 0,1,2,3 position differences\n",frst);

   fds=(FILE*)fopen(dsf,"w");  fdn=(FILE*)fopen(dnf,"w");
   ft=(FILE*)fopen(tf,"w");
   if(fds==NULL || fdn==NULL || ft==NULL) error2("err DistanceMatNG86: file error2");
   fprintf(fds,"%6d\n",com.ns);  fprintf(fdn,"%6d\n",com.ns); 
   fprintf(ft,"%6d\n",com.ns);
   if(noisy)  puts("NG distances for seqs.:");
   FOR (is,com.ns) {
      fprintf(fout,"\n%-*s", wname,com.spname[is]);
      fprintf(fds,   "%-*s ",wname,com.spname[is]);
      fprintf(fdn,   "%-*s ",wname,com.spname[is]);
      fprintf(ft,    "%-*s ",wname,com.spname[is]);
      FOR (js,is) {
         FOR(k,4) nsd[k]=0;
         ls1=(com.cleandata?com.ls:0);
         for (h=0,nst=nat=S=N=0; h<com.npatt; h++)  {
            if(com.cleandata)
               FOR(i0,2) getcodon(codon[i0],FROM61[com.z[i0==0?is:js][h]]);
            else {
               FOR(i0,2) FOR(k,3) codon[i0][k]=com.z[i0==0?is:js][h*3+k];
               for(i0=0,missing=0;i0<2;i0++) {
                  CodonListall(codon[i0], nb, ib);
                  if(nb[0]*nb[1]*nb[2]!=1)  { missing=1; break; }
               }
               if(missing) continue;
               ls1+=com.fpatt[h];
            }
            ndiff=difcodonNG(codon[0],codon[1],&St,&Nt,&ns,&na,0,com.icode);
            nsd[ndiff]+=(int)com.fpatt[h];
            S+=St*com.fpatt[h];
            N+=Nt*com.fpatt[h];
            nst+=ns*com.fpatt[h];
            nat+=na*com.fpatt[h];
         }  /* for(h) */
         y=ls1*3./(S+N); S*=y; N*=y;       /* rescale for stop codons */

         if(noisy>=9)
           printf("\n%3d%3d:Sites%7.1f +%7.1f =%7.1f\tDiffs%7.1f +%7.1f =%7.1f",
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

/* fprintf(frst1,"NG %7.1f %7.1f%7.4f %6.4f %6.4f ", S,N,dN,dS,dN_dS);
*/
         if(dN<0) dN=bigD; if(dS<0) dS=bigD;
         SeqDistance[is*(is-1)/2+js] = (S*dS+N*dN)*3/(S+N);
         fprintf(fds," %7.4f", dS);   fprintf(fdn," %7.4f", dN);
         fprintf(ft," %7.4f", (S*dS+N*dN)*3/(S+N));
      }
      FPN(fds); FPN(fdn); FPN(ft);
      if(noisy>1 && noisy<9 && com.ns>10)  printf(" %3d",is+1);
   }    /* for(is) */
   FPN(F0); FPN(fout);
   if(status) fprintf (fout, "NOTE: -1 means that NG86 is inapplicable.\n");
   fputs("\n\n",frst); fflush (frst);
   fclose(fds);  fclose(fdn);  fclose(ft);
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

   if ((fdaa=fopen(com.daafile, "r"))==NULL) 
      { printf("\nAA dist file %s not found.", com.daafile); exit(-1); }
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
      FOR(i,naa) if(fscanf(fdaa,"%lf",&com.pi[i])!=1) error2("err aaRatefile");
      if (fabs(1-sum(com.pi,20))>1e-6) {
         printf("\nSum of freq. = %.6f != 1 in aaRateFile\n",sum(com.pi,naa)); 
         exit(-1);
      }

      for(j=0,com.npi0=0; j<naa; j++)   if(com.pi[j]==0)  com.npi0++;
      if(com.npi0==0)
         for(j=0; j<naa; j++)  com.pi_sqrt[j]=sqrt(com.pi[j]);
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
   int nc=Nsensecodon[com.icode], naa=20, i,j,k, ic1,ic2, ndiff, from[3],to[3];
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
   com.nrate=k-1;     /* one element (ijAAref) is fixed */
   return(0);
}

int GetOmegaAA (int OmegaAA[])
{
/* This routine reads the file OmegaAA.dat to initialize the
   lower diagonal matrix OmegaAA, which specifies the aa substituion
   rate classes.  To be used with the codon substitution model
   AAClasses, which specifies several classes of the dN/dS ratio.

   OmegaAA[iaa*(iaa-1)/2+jaa]=-1 if no change; 
                             =0 for the first, background, class
                             =i (1,..,nclass) if iaa and jaa are in class i 
*/
   char *OmegaAAf="OmegaAA.dat", line[1024];
   FILE *fin;
   int iomega, n1step=0, i,j,k, iaa,jaa, npair, naa=20, nline=1024;
   int fromfile=1;  /* if 0, estimate all elements in the distance matrix */

   for(i=0,n1step=0; i<naa; i++) FOR(j,i)
      if (AA1STEP[i*(i-1)/2+j]) { OmegaAA[i*(i-1)/2+j]=0; n1step++; }
      else                        OmegaAA[i*(i-1)/2+j]=-1;
   if (noisy) {
       printf("\n\n%d one-step aa pairs.\n", n1step);
       printf("Reading dN/dS class info from %s.\n", OmegaAAf);
   }
   if ((fin=fopen(OmegaAAf,"r"))==NULL) 
      { printf("file %s does not exist.", OmegaAAf);  fromfile=0; }
   else {
      fscanf(fin, "%d", &N_OMEGA);  
      N_OMEGA++;
      printf ("# of dN/dS classes requested: %d.\n", N_OMEGA);
      if (N_OMEGA<1 || N_OMEGA>65-1) { 
         fromfile=0;  fclose(fin); 
         printf ("This is reset to %d.\n", n1step);
      }
   }
   if (!fromfile) {
      if (com.seqtype!=CODONseq) puts("\nTo be tested.\a\a");
      N_OMEGA=0;
      if (com.seqtype==AAseq) {
         FOR(i,naa) FOR(j,i) if(i*naa+j!=ijAAref && AA1STEP[i*(i-1)/2+j])
             OmegaAA[i*(i-1)/2+j]=N_OMEGA++;
      }
      else
         FOR(i,naa) FOR(j,i) 
           if(AA1STEP[i*(i-1)/2+j]) OmegaAA[i*(i-1)/2+j]=N_OMEGA++;
      printf("%d dN/dS ratios estimated from data.\nStop if wrong.",N_OMEGA);
      getchar ();
   }
   else {
      FOR (iomega, N_OMEGA-1) {
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

            printf ("\npair %d: |%c%c| ", npair, iaa,jaa);
            iaa=CodeChara((char)iaa,AAseq); jaa=CodeChara((char)jaa,AAseq);
            if(iaa<0||iaa>19||jaa<0||jaa>19) error2("aa not found");
            if (iaa<jaa)  { k=jaa, jaa=iaa; iaa=k; }
      
            printf ("|%c%c (%2d,%2d)| ", AAs[iaa], AAs[jaa],iaa,jaa);
            if (iaa==jaa) printf ("This pair has no effect.");
            if (OmegaAA[iaa*(iaa-1)/2+jaa]) error2("This pair specified?");
            if (AA1STEP[iaa*(iaa-1)/2+jaa]==0) error2("This pair has rate 0!");
            OmegaAA[iaa*(iaa-1)/2+jaa]=iomega+1;
            printf (" in class %d ",iomega+1);
         }
      }
   }
   com.nrate = !com.fix_kappa + N_OMEGA;
   printf ("\nNw=%d\tnrate=%d\n", N_OMEGA, com.nrate);
   if (N_OMEGA>n1step-(com.seqtype==AAseq)) error2("too many classes.");
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
   double fb3x4[12], fb4[4];

   if (com.codonf==Fequal) { fillxc(pi,1./n,n); return 0; }
   if (com.codonf<Fcodon) {
      for (i=0,zero(fb3x4,12),zero(fb4,4); i<n; i++) {
         ic=FROM61[i];  b[0]=ic/16; b[1]=(ic/4)%4; b[2]=ic%4;
         FOR(j,3) { fb3x4[j*4+b[j]]+=pi[i];  fb4[b[j]]+=pi[i]/3.; }
      }
      for (i=0; i<n; i++) {
         ic=FROM61[i];  b[0]=ic/16; b[1]=(ic/4)%4; b[2]=ic%4;
         if (com.codonf==2)  pi[i]=fb3x4[b[0]]*fb3x4[4+b[1]]*fb3x4[8+b[2]];
         else                pi[i]=fb4[b[0]]*fb4[b[1]]*fb4[b[2]];
      }
      abyx (1./sum(pi,n), pi, n);
      for(j=0,com.npi0=0; j<com.ncode; j++)   if(com.pi[j]==0)  com.npi0++;
      if(com.npi0==0)
         for(j=0; j<com.ncode; j++)  com.pi_sqrt[j]=sqrt(com.pi[j]);
   }
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
   int n=com.ncode, h,k;
   double lnL=0, fh,expt[NCODE];

   NFunCall++;
   if(!com.fix_kappa && com.hkyREV==0) com.kappa=x[1];
   if(!com.fix_omega) com.omega=x[1+com.nkappa];
   if(!com.fix_kappa||!com.fix_omega)
      EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, 
         (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),com.omega,PMat);

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

   if(vtw[0]<=0 || vtw[3]<=0)
      { puts("var(dS,dN) not calculable."); zero(vdSdN,4); return(-1); }

   /* printf("\nt & w: %.5f %.5f\n", t, omega);
      matout(F0,vtw, 2,2); */
   EigenQc(1,t,&S,&dS,&dN,NULL,NULL,NULL, &com.kappa,omega,PMat);

   eh=(t+1)*Small_Diff;
   EigenQc(1,t+eh,&S,&dS1,&dN1,NULL,NULL,NULL, &com.kappa,omega,PMat);
   EigenQc(1,t-eh,&S,&dS2,&dN2,NULL,NULL,NULL, &com.kappa,omega,PMat);
   JacobiSN[0*np+0]=(dS1-dS2)/(2*eh);
   JacobiSN[1*np+0]=(dN1-dN2)/(2*eh);
  
   eh=(omega+1)*Small_Diff;
   EigenQc(1,t,&S,&dS1,&dN1,NULL,NULL,NULL, &com.kappa,omega+eh,PMat);
   EigenQc(1,t,&S,&dS2,&dN2,NULL,NULL,NULL, &com.kappa,omega-eh,PMat);
   JacobiSN[0*np+1]=(dS1-dS2)/(2*eh);
   JacobiSN[1*np+1]=(dN1-dN2)/(2*eh);
  
   matby(JacobiSN,vtw,T1,2,2,2);
   mattransp2 (JacobiSN, T2, 2, 2);
   matby(T1,T2,vdSdN,2,2,2);

   /* matout(F0,vdSdN, 2,2); */

   return (0);
}

int PairwiseCodon (FILE *fout, double space[])
{
/* Calculates ds & dn for all pairwise codon sequence comparisons.
   It uses different npatt for different pairs.
   The data com.z[] should be encoded clean data, with ambiguity characters 
   removed.  Think of what to do with raw unclean data.
   JacobiSN has two columns, the 1st are deratives of dS (dS/dt, dS/dk, dS/dw)
   and the second of dN.
*/
   FILE *fds,*fdn,*ft;
   char *pz0[NS];   /* pz0, npatt0, & fpatt0 hold the old information */
   int npatt0=com.npatt;
   double *fpatt0, ls0=com.ls;
   float fp[NCODE*NCODE];
   char dsf[32]="2ML.dS",dnf[32]="2ML.dN",tf[32]="2ML.t", codon[2][3];
   int n=com.ncode, is,js,j,k,h, i0,np, wname=15;
   int nb[3],ib[3][4],ic[2], missing=0;
   double x[7]={.9,1,.5,.5,.5,.5,.3}, xb[7][2], lnL, e=1e-7, *var=space+NP, S,dS,dN;
   double JacobiSN[2*3],T1[2*3],T2[2*3],vSN[2*2], dS1,dN1,dS2,dN2,y[3],eh; 
          /* for calculating SEs of dS & dN */
   double tb[2]={1e-5,39}, kappab[2]={.1,999}, omegab[2]={.001,99};

   fpatt0=(double*)malloc(npatt0*3*sizeof(double));
   FOR(k,com.ns) pz0[k]=com.z[k];
   com.z[0]=(char*)(fpatt0+npatt0);  com.z[1]=com.z[0]+npatt0;
   FOR (k,npatt0) fpatt0[k]=(float)com.fpatt[k];

   if(!com.cleandata) puts("\nPairwiseCodon: pairwise deletion.");
   if (com.ngene>1 && com.Mgene==1) puts("ngene>1 to be tested.");
   if (noisy) printf("\n\npairwise comparison (Goldman & Yang 1994)..\n");
   fprintf(fout,"\npairwise comparison, codon frequencies: %s.\n",
      codonfreqs[com.codonf]);
   fds=(FILE*)fopen(dsf,"w"); fdn=(FILE*)fopen(dnf,"w"); ft=(FILE*)fopen(tf,"w"); 
   if(fds==NULL || fdn==NULL || ft==NULL) error2("err PairwiseCodon: file error");

   xb[0][0]=tb[0]; xb[0][1]=tb[1];
   com.nkappa=(com.hkyREV ? 5 : !com.fix_kappa);
   FOR(j,com.nkappa) { xb[1+j][0]=kappab[0]; xb[1+j][1]=kappab[1]; }
   if(!com.fix_omega)  { k=1+com.nkappa; xb[k][0]=omegab[0]; xb[k][1]=omegab[1]; }

   fprintf(fds,"%6d\n", com.ns);  fprintf(fdn,"%6d\n", com.ns);
   fprintf(ft,"%6d\n", com.ns);
   fprintf(frst, "\n\npairwise comparison (Goldman & Yang 1994)");
   fprintf(frst,
      "\nseq seq        N       S       dN      dS   dN/dS    Paras.\n");

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
         if(com.codonf<Fcodon) GetCodonFreqs(com.pi);


printf("\n%s eigen routine\n", (com.npi0?"old":"new"));


         np=com.np=(com.ntime=1)+com.nkappa+!com.fix_omega;  NFunCall=0;
         x[0]=SeqDistance[is*(is-1)/2+js]*(.5+rndu());  /* NG86 as initial */
         if(x[0]<0.01||x[0]>3) x[0]=0.1*(.5+rndu()); 

         FOR(j,com.nkappa) x[1+j]=.5+.1*rndu();
         if(!com.fix_omega) x[k=1+com.nkappa]=.5*rndu();
         if(noisy>=9) {
            FPN(F0);  FOR(k,np) printf(" %12.6f",x[k]); FPN(F0);
            FOR(k,np) printf(" %12.6f",xb[k][0]); FPN(F0);
            FOR(k,np) printf(" %12.6f",xb[k][1]); FPN(F0);
         }
         if(com.fix_kappa && com.fix_omega)  
            EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, 
            (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),com.omega,PMat);
         if (x[0])
            ming2(noisy>2?frub:NULL,&lnL,lfun2dSdN,NULL,x,xb, space,e,np);
         else {  x[1]=x[2]=com.kappa=com.omega=0; lnL=0; }
         fprintf(fout,"\nlnL =%12.6f\n",-lnL);
         FOR(k,np) fprintf(fout," %8.5f",x[k]);  FPN(fout);

if(noisy>2) printf("\n\nt_NG = %.4f\tt_ML = %.4f w_ML = %.4f",
               SeqDistance[is*(is-1)/2+js]*1.1,x[0],x[np-1]);

         if (x[0]&&com.getSE) {
            Hessian(np, x, lnL, space, var, lfun2dSdN, var+np*np);
            matinv(var, np, np, var+np*np);
            fprintf(fout,"SEs for parameters:\n");
            FOR(k,np) fprintf(fout," %8.5f",(var[k*np+k]>0.?sqrt(var[k*np+k]):-0));
            FPN(fout);
         }
         FPN(fout);
         EigenQc(1,x[0],&S,&dS,&dN, NULL,NULL,NULL,
            (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),com.omega,PMat);
         fprintf(fds," %7.4f", dS);   fprintf(fdn," %7.4f", dN);
         fprintf(ft," %7.4f", x[0]);

         fprintf (fout,
             "t=%7.4f  S=%8.1f  N=%8.1f  dN/dS=%7.4f  dN=%7.4f  dS=%7.4f\n",
              x[0],S,com.ls*3-S,com.omega,dN,dS);

         fprintf(frst,"%8.1f%8.1f %8.4f%8.4f%8.4f",com.ls*3-S,S,dN,dS,com.omega);
         FOR(k,np) fprintf(frst," %8.4f",x[k]);  

/* 
fprintf(frst1,"ML: %6d%6d %8.4f%8.4f %8.4f\n",is+1,js+1,dS,dN,com.omega);
fprintf(frst1,"%8.4f%8.4f %8.4f\n", dS,dN,com.omega);
*/

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
               EigenQc(1,y[0],&S,&dS1,&dN1,NULL,NULL,NULL,
                  (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),com.omega,PMat);
               y[k] -= 2*eh;
               if(!com.fix_kappa) com.kappa=y[1];
               com.omega=y[1+!com.fix_kappa];
               EigenQc(1,y[0],&S,&dS2,&dN2,NULL,NULL,NULL,
                  (com.hkyREV?x+com.ntime+com.nrgene:&com.kappa),com.omega,PMat);

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
/*
         fprintf(frst, "%7.3f%7.3f%7.3f%10.2f", x[0],com.kappa,com.omega,-lnL);
         FPN (frst);
*/
         fflush(frst);  fflush(fout);
      }  /* for (js) */
      FPN(fds); FPN(fdn); FPN(ft);   fflush(fds); fflush(fdn); fflush(ft); 
   }     /* for (is) */
   fclose(fds); fclose(fdn); fclose(ft);

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

int PairwiseAA (FILE *fout)
{
/* Calculates pairwise distances using amino acid seqs.
   Data (com.z[]) are clean and coded.
   com.npatt for the whole data set is used which may be greater than 
   the number of patterns for each pair.
   SE is not calculated.
*/
   FILE *fdist;
   char distf[32]="2ML.aa", *pz0[NS];
   int n=com.ncode, i, is,js;
   double x, xb[2]={0,19}, lnL, step;

   if (com.ngene>1 && com.Mgene==1) error2("ngene>1 to be tested.");
   if (noisy) printf("\npairwise ML distances of AA seqs.\n");
   if((fdist=(FILE*)fopen(distf,"w"))==NULL) error2("PairwiseAA: file error");

   if(com.model>Empirical_F)  error2("PairwiseAA: model wrong");
   if(com.model==0)  fillxc(com.pi,1./n, n);
   if(com.model>=Empirical)  GetDaa(NULL, com.daa);
   if(com.model==0 || com.model==Empirical)  EigenQaa(NULL, Root, U, V, NULL);

   FOR(i,com.ns) pz0[i]=com.z[i];
   fprintf(fout,"\nML distances of aa seqs, also in file %s\n",distf);
   if(com.alpha) 
      fprintf(fout,"\ncontinuous gamma with alpha = %.3f used.\n\n",com.alpha);

   fprintf(fdist,"%6d\n", com.ns);
   for(is=0; is<com.ns; is++,FPN(F0),FPN(fout),FPN(fdist)) {
      printf ("%4d vs", is+1);
      fprintf(fdist,"%-14s ", com.spname[is]);
      fprintf(fout, "%-14s ", com.spname[is]);
      FOR(js,is) {
         com.z[0]=pz0[is]; com.z[1]=pz0[js]; 
         printf (" %2d", js+1);
         if(com.model==1||com.model==Empirical_F) {
            for (i=0,zero(com.pi,n); i<com.npatt; i++) {
               com.pi[com.z[0][i]]+=com.fpatt[i];
               com.pi[com.z[1][i]]+=com.fpatt[i];
            }
            abyx(1./sum(com.pi,n), com.pi, n);

            for(i=0,com.npi0=0; i<n; i++)  if(com.pi[i]==0)  com.npi0++;
            if(com.npi0==0)
               for(i=0; i<n; i++)  com.pi_sqrt[i]=sqrt(com.pi[i]);

            EigenQaa(NULL,Root,U,V,NULL);
         }
         /* com.posG[1]=com.npatt; */

         xb[0]=SeqDistance[is*(is-1)/2+js];  x=xb[0]*1.5;  step=xb[0];
         LineSearch(lfun2AA, &lnL, &x, xb, step);
         fprintf(fdist," %7.4f",x); fprintf(fout," %7.4f",x); 
         if (com.getSE) ;
      }  /* for (js) */
   }     /* for (is) */

   fclose(fdist);
   FOR(i,com.ns) com.z[i]=pz0[i];
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



int lfunNSsites_rate (FILE* frst, double x[], int np)
{
/* This calculates the dN/dS rates for sites under models with variabel dN/dS 
   ratios among sites (Nielsen and Yang 1998).  Modified from lfundG() 
   Only conditional mode is used?
   com.fhK[] holds the posterior probs.
*/
   int  h,hp, ir, i,it=0, refsp=0,iaa;
   double  lnL=0, fh, cutoff=0.5, mw=-1, w2=x[com.np-1],psel=0;
   char  *sig;

   int  ncolors=3;
   char *colors[3]={"gray", "blue", "red"};
   int nwanted=57, wanted=0;
   int wantedlist[]={5, 7, 9, 22, 24, 26, 57, 58, 59, 61, 62, 63, 64, 65, 66, 
                     67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 80, 81, 82, 
                     84, 95, 97, 99, 114, 116, 143, 145, 146, 147, 149, 150, 
                     151, 152, 154, 155, 156, 157, 158, 159, 161, 162, 163, 
                     165, 166, 167, 169, 171};

   if(com.nparK) error2("lfunNSsites_rate to be done for HMM.");

   if(!com.aaDist) {
      fprintf(frst,"\n\ndN/dS for site classes (K=%d)\np: ",com.ncatG);
      FOR(i,com.ncatG) fprintf(frst,"%9.5f",com.freqK[i]);
      fputs("\nw: ",frst);
      if(com.model==0) FOR(i,com.ncatG) fprintf(frst,"%9.5f",com.rK[i]);
      else {
         FOR(i,com.ncatG-2) fprintf(frst,"%9.5f",com.rK[i]);
         fprintf(frst, "  *        *\nw2 = %.5f", w2);
      }
      FPN(frst);
   }
   else  fputs("\nCheck main result file for parameter estimates\n",frst);
   FPN(frst);

   fx_r(x,np);
   if(_nnodeScale)
      FOR(h,com.npatt) {
         for(ir=1,fh=com.fhK[h]; ir<com.ncatG; ir++)
            if(com.fhK[ir*com.npatt+h]>fh) fh=com.fhK[ir*com.npatt+h];
         for(ir=0; ir<com.ncatG; ir++)
            com.fhK[ir*com.npatt+h]=exp(com.fhK[ir*com.npatt+h]-fh);
         lnL-=fh*com.fpatt[h];
      }

   FOR(h,com.npatt) {
      for (ir=0,fh=0; ir<com.ncatG; ir++)
         fh += (com.fhK[ir*com.npatt+h]*=com.freqK[ir]);
      FOR(ir, com.ncatG)  com.fhK[ir*com.npatt+h]/=fh;
      lnL-=com.fpatt[h]*log(fh);
   }

   fprintf(frst,"\nPosterior P for %d classes (class)",com.ncatG);
   if(com.model==0) fprintf(frst,"& mean w");
   fprintf(frst,"\n%s used as reference\n\n",com.spname[refsp]);
   for(h=0; h<com.ls; h++,FPN(frst)) {

/*
for(i=0,wanted=0; i<nwanted; i++)
   if(h+1==wantedlist[i]) { wanted=1; break; }
if(!wanted) continue;
*/

      hp=com.pose[h];
      iaa=GetAASiteSpecies(refsp, hp);
      fprintf(frst,"%4d %c  ",h+1,AAs[iaa]);
      for (ir=0,it=0,mw=psel=0; ir<com.ncatG; ir++) {
         fprintf(frst," %9.4f", com.fhK[ir*com.npatt+hp]);
         if(com.fhK[ir*com.npatt+hp]>com.fhK[it*com.npatt+hp])  it=ir;
         if(com.model==0) {
            mw+=com.fhK[ir*com.npatt+hp]*com.rK[ir];
            if(com.rK[ir]>1) psel+=com.fhK[ir*com.npatt+hp];
         }
      }

      fprintf(frst, " (%2d)", it+1);
      if(com.model==0) {
         fprintf(frst, "  %9.3f", mw);
         if(psel) fprintf(frst, "  %9.3f", psel);
      }
   }

   /* list of positively selected sites */
   if(com.model==0) { /* NSsites models */
      for(ir=0,it=0; ir<com.ncatG; ir++) if(com.rK[ir]>1) it=1;
      if(!com.aaDist && it) {
         fprintf(fout,"\n\nPositively selected sites Prob(w>1):\n\n");
         fprintf(frst,"\n\nPositively selected sites Prob(w>1):\n\n");
         FOR(h,com.ls) {
            hp=com.pose[h];
            for (ir=0,psel=0; ir<com.ncatG; ir++)
               if(com.rK[ir]>1) psel+=com.fhK[ir*com.npatt+hp];

            if(psel>cutoff) {
               sig="";  if(psel>.95) sig="*";  if(psel>.99) sig="**";
               iaa=GetAASiteSpecies(refsp, hp);
#ifdef SITELABELS
               fprintf(fout,"%6d  %4s %.4f%s\n",h+1,sitelabels[h],psel, sig);
#else
               fprintf(fout,"%6d %c %.4f%s\n",h+1,AAs[iaa],psel, sig);
               fprintf(frst,"%6d %c %.4f%s\n",h+1,AAs[iaa],psel, sig);
#endif
            }
	      }
         FPN(fout);
         if(com.rK[com.ncatG-2]>1)
            fputs("\nNote: more than one w>1.  Check rst for details\n",fout);
         if(com.rK[com.ncatG-1]<1.2)
            fputs("\nNote: w is small although >1.  Excise caution about the list.\n",fout);
      }
   }
   else {  /* branch&site models */
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

/* 
   RasMol script for coloring structure.  comment this out 
   */
   /*
   if(com.verbose && com.model==0 && com.NSsites==NSdiscrete && com.rK[2]>1) {
      FOR(h,com.ls) {
         hp=com.pose[h];
         for (ir=1,it=0; ir<com.ncatG; ir++)
            if(com.fhK[ir*com.npatt+hp]>com.fhK[it*com.npatt+hp]) it=ir;
         fprintf(frst1,"select %d\n color %s\n",h+1,colors[it]);
      }
   }
   */
   fprintf (frst,"\n\nlnL = %12.6f\n", -lnL);
   return (0);
}






int eigenQREV (double Q[], double pi[], double pi_sqrt[], int n, 
             double Root[], double U[], double V[])
{
/* 
   This finds the eigen solution of the rate matrix for a time-reversible 
   Markov process.
   Rate matrix Q = S * diag{pi} = U * diag{Root} * V, 
   where S is symmetrical, all elements of pi are positive, and U*V = I.

   This uses the algorithm for a real symmetric matrix and assumes 
   that all elements of pi[] are positive.

   Ziheng Yang, 23 May 2001
*/
   int status=0, i,j;

   /* construct symmetrical matrix S = sqrt(D) * Q * sqrt(-D) */
   if(pi_sqrt==NULL)
      for(j=0;j<n;j++) pi_sqrt[j]=sqrt(pi[j]);
   for(i=0;i<n;i++) {
      for(j=0;j<i;j++)
         U[i*n+j] = U[j*n+i] = (Q[i*n+j] * pi_sqrt[i]/pi_sqrt[j]);
      U[i*n+i] = Q[i*n+i];
   }
   status=eigenRealSym(U, n, Root, V);
   for(i=0;i<n;i++) for(j=0;j<n;j++)  V[i*n+j] = U[j*n+i] * pi_sqrt[j];
   for(i=0;i<n;i++) for(j=0;j<n;j++)  U[i*n+j] /= pi_sqrt[i];

   return(status);
}



/*

(*) Codon models for variable d_N/d_S ratios among lineages

    model=0: one d_N/d_S ratio for all lineages.
             fix_omega=1 will fix the single omega at the value provided
    model=1: each branch has its own omega, use fix_omega=0 only
    model=2: 2 to nbranch ratios.  Use consecutive branch marks 0,1,2,... 
             to specify omega.
             fix_omega=1 fixes the last omega at the value provided.
 
    Some computation can be saved under com.model==NSbranch2 if only 2 or a few more 
    ratios are used instead of all b ratios.  Check whether it is worth the 
    effort.

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


(*) Codon models for variable dN/dS ratios among both branches and sites
    (model=2, NSsites=3 or 2)
    (com.nrate includes kappa & omega)
    Parameters include branchlens, kappa, p0, p1, w0, w1, w2

    method = 0: SetPSiteClass copies w's to nodes[].omega and PMat is calculated
    in PartialLikelihood().  
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
            The total number of classes (N_OMEGA) is one plus the number 
            specified in the file OmegaAAf.

   N_OMEGA is the number of different dN/dS ratios in the NSbranchB, NSbranch2 models
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
