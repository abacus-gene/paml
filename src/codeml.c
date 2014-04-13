/* CODEML.c  (AAML.c & CODONML.c)
   Maximum likelihood parameter estimation for codon sequences (seqtype=1) 
                    or amino-acid sequences (seqtype=2)
               Copyright, Ziheng Yang (March 1993 onwards)

               cc -o codeml -fast codeml.c tools.o eigen.o -lm
       cc -o codemlsites -fast -DNSSITESBATCH codeml.c tools.o eigen.o -lm 
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


#include "tools.h"

#ifdef macintosh
/* Added by Andrew Rambaut to accommodate Macs -
   (1) Brings up dialog box to allow command line parameters.
   (2) Reduces the size of statically allocated data.    */
#include <console.h>
#define NS            200
#else
#define NS            500
#endif

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
extern int GenetCode[][64];
extern int noisy, NFunCall, *ancestor;
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
int printsmaCodon (FILE *fout,char * z[],int ns,int ls,int lline,int simple);
int InitializeCodon (FILE *fout, double space[]);
int AA2Codonf (double faa[20], double fcodon[]);
int DistanceMatAA (FILE *fout);
int DistanceMatNG86 (FILE *fout, double alpha);
int GetDaa(FILE *fout, double daa[]);
void getpcodonClass(double x[], double pcodonClass[]);
int EigenQc (int getstats, double branchl, double *S, double *dS, double *dN,
    double Root[], double U[], double V[],
    double kappa, double omega, double Q[]);
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

void SimulateData2s61(void);
void Ina(void);

struct common_info {
   char *z[NS], *spname[NS], seqf[96],outf[96],treef[96],daafile[96];
   int seqtype, ns, ls, ngene, posG[NGENE+1], lgene[NGENE], npatt,*pose;
   int runmode,clock,verbose,print, codonf,aaDist,model,NSsites, cleandata;
   int method, icode, ncode, Mgene, ndata;
   int fix_rgene,fix_kappa,fix_omega,fix_alpha,fix_rho,nparK,fix_branch,getSE;
   int np, ntime, nrgene, nrate, nalpha, ncatG, sspace, slkl1;
   double *fpatt, *space;
   double pi[NCODE],fb61[64],piG[NGENE][64],kappa,omega,alpha,rho,rgene[NGENE];
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
double *PMat,*U,*V,*Root, *_UU[5],*_VV[5],*_Root[5]; 
/* for 5 sets for branchsite models */

double xcom[NP-NBRANCH];
double *POMEGA;
  /* communication between SetParameters & PartialLikelihood & EigenQc.  
     Remove it? */
double pcodon0[64],paa0[20], *pcodonClass;  /* for aaDist=FIT1 or FIT2 */

int LASTROUND, N_eigenQc=0;
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

double AAchem[][20+1]={
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
char *NSbranchmodels[]={"One dN/dS ratio for branches", 
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
   SetSeed(13213);
/*
mergeSeqs(frst);  exit(0);
Ina();
*/

SetSeed(i=(int)clock()+1);
#if (DSDN_MC || DSDN_MC_SITES)
SimulateData2s61();
#endif

   if(argc>1) strcpy(ctlf,argv[1]);
   GetOptions(ctlf);

#ifdef NSSITESBATCH
         if(com.fix_omega) error2("fix omega during batch run?");
         if(com.runmode) error2("runmode?");
         com.NSsites=NSbetaw;  com.ncatG=ncatG0+1;
#endif

   if ((fout=fopen (com.outf, "w"))==NULL) error2("outfile creation err.");

   if((fseq=fopen (com.seqf,"r"))==NULL) {
      printf ("\n\nSequence file %s not found!\n", com.seqf);
      exit (-1);
   }
   if(noisy && com.seqtype==CODONseq) 
      { printcu(F0,NULL,com.icode); puts("The code is nice, uuh?"); }

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
      if((nodes=(struct TREEN*)malloc(i))==NULL) error2("oom");
      /* BootstrapSeq("boot");  exit(0); */

      if(com.ngene>1 && com.Mgene==1)  printSeqsMgenes ();
      if(com.Mgene==1 && com.print)
         { puts("Mgene==1 & RateAncestor not tested. Do not use?"); sleep(100);}
      if(com.ngene>1 && com.aaDist>=FIT1)  /* because of pcodon0[] */
         { error2("ngene for fitness models");}

      pmodel=(com.seqtype==CODONseq?NSbranchmodels[com.model]:aamodels[com.model]);
      fprintf(fout,"%s %s Model: %s ",seqtypestr[com.seqtype-1],com.seqf,pmodel);
      if(com.clock) fprintf(fout," %s ",clockstr[com.clock]);
      if(com.seqtype==CODONseq||com.model==FromCodon) {
         if (com.fix_kappa)  fprintf(fout, " kappa = %.3f\n", com.kappa);
         if (com.fix_omega) fprintf(fout, " omega? = %.3f fixed\n", com.omega);
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
         fprintf(fout,"Site-class models: %s",NSsitesmodels[com.NSsites]);
         if(com.NSsites>=NSdiscrete)fprintf(fout," (%d categories)",com.ncatG);
         if(com.nparK) fprintf(fout," & HMM");
         FPN(fout);
         if(com.aaDist)
            fprintf(fout,"\nFitness models: aaDist: %d\n",com.aaDist);
      }

      if(com.clock==2 && 
        (rateBranch=(int*)realloc(rateBranch,com.ns*2*sizeof(int)))==NULL) 
           error2("oom rateBranch");
      com.sspace=max2(100000,com.ls*(int)(com.ns*3*sizeof(char)+sizeof(int)));
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

      if(com.seqtype==AAseq) {
         Initialize (fout, com.space, com.seqtype);
         if (com.model==FromCodon||com.model==AAClasses) 
            AA2Codonf(com.pi, com.fb61);  /* get codon freqs from aa freqs */ 
      }
      else {
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

      if(com.seqtype==CODONseq) DistanceMatNG86(fout,0);
      else                      DistanceMatAA(fout);

      fflush(fout);

      if(com.seqtype==AAseq && com.model==Poisson && !com.print) 
         PatternJC69like (NULL);
      if(com.alpha || com.NSsites)
         if((com.fhK=(double*)realloc(com.fhK,
             s2=com.npatt*com.ncatG*sizeof(double)))==NULL) error2("oom fhK");

      if(com.runmode==-2 && com.Mgene!=1) {
         if(com.seqtype==CODONseq) PairwiseCodon(fout,com.space);  
         else                      PairwiseAA(fout);  
      }
      else {
         printf("\n%9ld bytes for distance",com.ns*(com.ns-1)/2*sizeof(double));

         if(!com.cleandata) {
            slkl0=com.ns*com.ncode*com.npatt*sizeof(double);
            if((com.lkl0=(double*)malloc(slkl0))==NULL) error2("oom lkl0");
         }
         com.slkl1=(com.ns-1)*com.ncode*com.npatt*sizeof(double);
         if((com.lkl=(double*)realloc(com.lkl,com.slkl1))==NULL) error2("oom lkl");

         printf("\n%9d bytes for lkl0\n%9d bytes for lkl1\n",slkl0,com.slkl1);
         printf ("%9d bytes for fhK\n%9d bytes for space\n",s2,com.sspace);

/*
SlidingWindow(fout,com.space);
exit(0);
*/
#ifdef NSSITESBATCH
         printf("\nNSsites batch run.\nUsing ncatG values as in YNGP2000.\n");
         FOR(i,NS3normal) printf("M% 2d: %s\n",i,NSsitesmodels[i]);
         printf("\nNumber of models and model codes (e.g., 6   0 1 2 3 7 8)? ");
         scanf("%d",&nmodels);  FOR(i,nmodels) scanf("%d",&nsmodels[i]);

         for(im=0; im<nmodels; im++) {
            com.NSsites=nsmodels[im];
            if(com.NSsites<=NSselection)  com.ncatG=com.NSsites+1;
            else if (com.NSsites==NSdiscrete) com.ncatG=3; 
            else if (com.NSsites==NSfreqs) com.ncatG=5;
            else if (com.NSsites==NSbetaw||com.NSsites==NS02normal) 
                  com.ncatG=ncatG0+1;
            else  com.ncatG=ncatG0;

            com.nrate=!com.fix_kappa;
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
         if (com.Mgene==1)         MultipleGenes(fout, com.space);
         else if (com.runmode==0)  Forestry(fout);
         else if (com.runmode==3)  StepwiseAddition(fout, com.space);
         else if (com.runmode>=4)  Perturbation(fout,(com.runmode==4),com.space);
         else                      StarDecomposition(fout, com.space);
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
   int  i,j=0,k, itree, ntree, np, iteration=1;
   int pauptree=0, length_label;
   double x[NP],xb[NP][2], lnL=0,lnL0=0, e=1e-7,tl=0, *var=NULL, nchange=-1;

   if(com.clock==3) GetSeqTimes();

   if ((ftree=fopen (com.treef,"r"))==NULL) error2 ("treefile not found");
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
      FOR(i,com.ncode) paa0[GenetCode[com.icode][FROM61[i]]]+=pcodon0[i];
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
            scanf("%d",&com.fix_branch);
            /* com.fix_branch=0; */

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
      NFunCall=0;
      lnL = com.plfun (x,np);
      if(noisy) {
         printf("\nnp =%6d", np);
         if(noisy>2 && np<100) matout(F0,x,1,np);
         printf("\nlnL0 = %12.6f\n",-lnL);


/*         
lnL = com.plfun (x,np);
printf("\nlnL0 = %12.6f\n",-lnL);
exit(0);         
*/         
         /* printf("%d eigenQc calls, %d lfun calls\n",N_eigenQc, NFunCall); */

/*
         gradient (np, x, lnL, com.space, com.plfun, space+np, 1);
         FOR(i,np) printf("%12.6f", com.space[i]);  FPN(F0);
*/
      }

      if(iteration) {
         SetxBound(np,xb);
         SetxInitials (x);

         if(com.method==0 || com.ntime==0 || com.clock)
            j=ming2(noisy>2?frub:NULL,&lnL,com.plfun,NULL,x,xb,com.space,e,np);
         else
            minB(noisy>2?frub:NULL, &lnL,x,xb, com.space);
         if(j) fprintf(fout,"\ncheck convergence..");

      }
      printf("Out..\nlnL  = %12.6f\n",-lnL);
      if (itree==0)
         { lnL0=lnL;  FOR(i,np-com.ntime) xcom[i]=x[com.ntime+i]; }
      else if (!j)
         for (i=0; i<np-com.ntime; i++) xcom[i]=xcom[i]*.2+x[com.ntime+i]*0.8;

      if(!LASTROUND && (com.NSsites==NSselection||com.NSsites==NSdiscrete
        ||com.NSsites==NSfreqs||com.NSsites==NS3normal)) {
         /* transform back to p0, p1,... */
         k=com.ntime+com.nrgene+!com.fix_kappa;
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
         FOR (i,np) var[i*np+i]=(var[i*np+i]>0.?sqrt(var[i*np+i]):-0);
         FOR (i,np) fprintf(fout," %8.5f", var[i*np+i]);  FPN (fout);
         /* if (com.getSE==2) matout2(fout,var,np,np,12,7); */
      }
      if(com.NSsites==NSselection && x[com.ntime+com.nrgene+!com.fix_kappa+2]<1)
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

      if (com.NSsites)
         lfunNSsites_rate(frst,x,np);
      else if (com.print) {
         if(com.rho==0 && com.nparK==0)        AncestralSeqs(frst,x);
         if(com.NSsites==0&& com.plfun!=lfun)  lfunRates(frate,x,np);
      }
      if(ntree>1 || ntree==-1) {
         com.print-=9;  com.plfun(x,np);  com.print+=9;
      }
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
   int i, k=com.ntime+com.nrgene+!com.fix_kappa, rank[NCATG];
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
      if(com.NSsites && com.model==0) {  /* calculate dN/dS by averaging over classes */
         for(j=0,Qfactor_NS=0,dS=dN=0; j<com.ncatG; j++) {
            freqK_NS=com.freqK[j];
            if(com.aaDist) {
               k=com.ntime+com.nrgene+!com.fix_kappa;
               if(com.aaDist<10)         POMEGA=x+k+com.ncatG-1+2*j;
               else if(com.aaDist>=FIT1) {
                  POMEGA=x+k+com.ncatG-1+j*(4+(com.aaDist==FIT2));
                  xtoy(pcodonClass+j*64, com.pi, com.ncode);
               }
            }
            EigenQc(1,1,&S,&dSt,&dNt,NULL,NULL,NULL,com.kappa,com.rK[j],PMat);
            /* t=1 used here, and dS & dN used later for each branch */
            dS+=freqK_NS*dSt;  dN+=freqK_NS*dNt;
            omclass[j]=dNt/dSt;
         }
         Qfactor_NS=1/Qfactor_NS;
         om=dN/dS;  dS*=Qfactor_NS;  dN*=Qfactor_NS;  N=com.ls*3-S;
      }

      k=com.ntime+com.nrgene;
      if (!com.fix_kappa) fprintf(fout,"kappa (ts/tv) = %8.5f\n", x[k++]);

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
         k=com.ntime+com.nrgene+!com.fix_kappa;
         fprintf (fout,"\nb = %9.5f", x[k++]);
         if (com.seqtype==CODONseq)  fprintf (fout,"\na = %9.5f\n", x[k++]);
      }
      else if(com.aaDist && com.aaDist>=11) { /* fitness, one class */
         fprintf (fout,"\nfitness model (a_p, p*, a_v, v*, (and w0 for FIT2):\n");
         k=com.ntime+com.nrgene+!com.fix_kappa;
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
      k=com.ntime+com.nrgene+!com.fix_kappa;
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
         fprintf (fout, "%9.5f:", x[com.ntime+com.nrgene+!com.fix_kappa+k]);
         FOR (i,20) FOR(j,i)
            if (OmegaAA[i*(i-1)/2+j]==k) fprintf(fout," %c%c", AAs[i],AAs[j]);
         if (k==0)  fprintf(fout, " (background ratio)");
         FPN(fout); 
      }
      if (com.nrate>65) { /* estimating all acceptance rates */
         for(i=0,k=com.ntime+com.nrgene+!com.fix_kappa; i<20; i++) FOR(j,i) {
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

   /* Calculate dN and dS for each branch in the tree */
   if(com.seqtype==CODONseq && (!com.model || !com.NSsites)
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
         k=com.ntime+com.nrgene+!com.fix_kappa;
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

            EigenQc(1,t,&S,&dS,&dN, NULL,NULL,NULL,com.kappa,om,PMat); /* */
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

            /*
            Qfactor_NS=Qfactor_NS_branch[nodes[i].label];
            EigenQc(1,t,&S,&dS,&dN,NULL,NULL,NULL,com.kappa,om,PMat);

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
   int i,j, nopt=32, lline=255;
   char line[255], *pline, opt[99], comment='*';
   char *optstr[] = {"seqfile", "outfile", "treefile", "seqtype", "noisy", 
        "cleandata", "runmode", "method", 
        "clock", "getSE", "RateAncestor", "CodonFreq", "verbose",
        "model","aaDist","aaRatefile",
        "NSsites", "NShmm", "icode", "Mgene", "fix_kappa", "kappa",
        "fix_omega", "omega", "fix_alpha", "alpha","Malpha", "ncatG", 
        "fix_rho", "rho", "ndata", "Small_Diff"};
   double t;
   FILE  *fctl;
   char *daafiles[]={"", "grantham.dat", "miyata.dat", 
                     "g1974c.dat","g1974p.dat","g1974v.dat","g1974a.dat"};

   if((fctl=fopen(ctlf,"r"))==NULL) error2(".ctl file open error.\n");
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
            if (noisy>2)
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
               case (14): com.aaDist=(int)t;      break;
               case (15): sscanf(pline+2,"%s",com.daafile); break;
               case (16): com.NSsites=(int)t;     break;
               case (17): com.nparK=(int)t;       break;
               case (18): com.icode=(int)t;       break;
               case (19): com.Mgene=(int)t;       break;
               case (20): com.fix_kappa=(int)t;   break;
               case (21): com.kappa=t;            break;
               case (22): com.fix_omega=(int)t;   break;
               case (23): com.omega=t;            break;
               case (24): com.fix_alpha=(int)t;   break;
               case (25): com.alpha=t;            break;
               case (26): com.nalpha=(int)t;      break;
               case (27): com.ncatG=(int)t;       break;
               case (28): com.fix_rho=(int)t;     break;
               case (29): com.rho=t;              break;
               case (30): com.ndata=(int)t;       break;
               case (31): Small_Diff=t;           break;
           }
           break;
         }
      }
      if (i==nopt)
         { printf ("\noption %s in %s\n", opt,ctlf);  exit(-1); }
   }
   fclose (fctl);

   if (noisy) FPN(F0);
   setmark_61_64 ();
   if (com.ngene==1) com.Mgene=0;
   if (com.seqtype==AAseq || com.seqtype==CODON2AAseq) {
      if(com.NSsites) error2("use NSsites=0 for amino acids?");
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
      com.ncode=Nsensecodon[com.icode];
      if (com.model!=AAClasses) {
         if(com.fix_kappa>1) error2("fix_kappa>1, not tested.");
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
         if (com.kappa<0)  error2("kappa..");
         if(com.runmode==-2 && (com.NSsites||com.alpha||com.aaDist))
            error2("err: incorrect model for pairwise comparison.");
         if(com.runmode>0 && com.model==2) error2("tree search & model");
         if(com.aaDist && com.NSsites!=0 && com.NSsites!=NSdiscrete)
            error2("NSsites && aaDist.");

         com.nrate=!com.fix_kappa;
         if(com.aaDist==0)  com.nrate+=!com.fix_omega;
         else {
            if(com.aaDist<=6)          com.nrate+=2;   /* a & b, PSB2000 */
            else if(com.aaDist==FIT1)  com.nrate+=4; /* fitness models: */
            else if(com.aaDist==FIT2)  com.nrate+=5; /* ap, p*, av, v*, b */
            if(com.aaDist>=FIT1) FOR(i,2) FOR(j,20) AAchem[i][j]/=AAchem[i][20];
         }
         if (com.Mgene>=3 && com.nrate==0)  error2("Mgene");

         if(com.NSsites) {
            com.nrate=!com.fix_kappa;
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
   if(com.alpha||com.NSsites) if(com.ncatG<2 || com.ncatG>NCATG) error2("ncatG");

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
   double tb[]={.0001,25}, rgeneb[]={0.1,9}, rateb[]={.0001,89};
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
   double tb[]={4e-6,50}, rgeneb[]={0.1,99}, rateb[]={1e-4,99};
   double alphab[]={0.005,99}, betab[]={.005,99}, omegab[]={.001,89};
   double rhob[]={0.01,0.99}, pb[]={.000001,.999999};

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
         iaa=GenetCode[com.icode][FROM61[i]];
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
   int slkl1_new=(tree.nnode-com.ns)*com.ncode*com.npatt*com.ncatG*sizeof(double);
   double t;

   
   com.plfun = (com.alpha==0 ? lfun : (com.rho==0?lfundG:lfunAdG));
   if(com.NSsites) com.plfun=lfundG;
   if(com.nparK) com.plfun=lfunAdG;

   if(com.method==1 && com.fix_branch!=2 && com.plfun==lfundG) {
      lklSiteClass=1;
      if(com.slkl1!=slkl1_new) {
         com.slkl1=slkl1_new;
         printf("\n%9d bytes for lkl1, adjusted\n",com.slkl1);
         if((com.lkl=(double*)realloc(com.lkl,com.slkl1))==NULL)
            error2("oom lkl1");
      }
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

   if(com.fix_branch==2)  com.ntime=0;
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
      for(j=1,x[0]=1; j<tree.nnode-com.ns; j++) x[j]=0.8;
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
      com.nrate=!com.fix_kappa+!com.fix_omega+N_OMEGA-1;
   }    /* if(dN/dS models for branches */

   if (!com.fix_kappa && com.Mgene>=3) com.nrate*=com.ngene;
   com.np = com.ntime+com.nrgene+com.nrate;

   /* NSbranchsite models 
      N_OMEGA=2 different w's at a site (three w's in the model: w0,w1,w2) */
   if(com.seqtype==CODONseq && com.model && com.NSsites) {
      if(N_OMEGA!=2) error2("only two branch labels are allowed");
      com.ncatG=4;
      if(com.NSsites==NSdiscrete) 
         com.nrate=!com.fix_kappa+ 2 +!com.fix_omega+N_OMEGA-1-1;  /* add w0 and w1 */
      else 
         com.nrate=!com.fix_kappa+!com.fix_omega+N_OMEGA-1-1;
      com.np = com.ntime+com.nrgene+com.nrate + 2;  /* add p0 and p1 */

      k=com.ntime+com.nrgene;
      if(!com.fix_kappa)  x[k++]=com.kappa;
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
            EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, com.kappa,com.omega,PMat);
         else if (com.model==NSbranchB || com.model==NSbranch2) {
            if (!com.fix_kappa)  x[com.ntime+com.nrgene]=com.kappa; 
            FOR(i, (com.model==NSbranchB?tree.nbranch:N_OMEGA-1+!com.fix_omega))
               x[com.ntime+com.nrgene+!com.fix_kappa+i]=com.omega;
         }
         else if(com.nrate) { /* either kappa, omega, or both for each gene */
            FOR(i, (com.Mgene>=3?com.ngene:1)) {
               if(!com.fix_kappa && !com.fix_omega)
                  { x[k+i*2]=com.kappa;  x[k+1+i*2]=com.omega; }
               else if (com.fix_kappa) x[k+i]=com.omega;
               else if (com.fix_omega) x[k+i]=com.kappa;
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
*/
   int nr1=com.nrate/com.ngene, k=com.nrgene+(com.Mgene>=3)*igene*nr1;

   if (_pi) xtoy (com.piG[igene],com.pi,com.ncode);
   if (_UVRoot) {
      if (com.seqtype==CODONseq) {
         if(!com.fix_kappa) com.kappa=xcom[k++];
         if(!com.fix_omega) com.omega=xcom[k++];
         if (!com.NSsites)
            EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, com.kappa,com.omega,PMat);
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
   if (com.nrate) {
      if(!com.fix_kappa) com.kappa=x[k++];
      if(!com.model && !com.aaDist && !com.fix_omega) com.omega=x[k];
      if(com.seqtype==AAseq)
         EigenQaa(NULL, Root, U, V, x+com.ntime+com.nrgene);
      else if((com.model==0 && com.NSsites==0 && com.Mgene<=1) || com.model==AAClasses)
         /* CODONs, same dN/dS across branches & sites */ 
         EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, com.kappa,com.omega,PMat);
      else if (com.model==2 && com.NSsites==0 && N_OMEGA<=5)
         FOR(i,N_OMEGA) {
            w0=(i==N_OMEGA-1&&com.fix_omega?OMEGA_FIX:POMEGA[i]);
            EigenQc(0,-1,NULL,NULL,NULL,_Root[i],_UU[i],_VV[i],com.kappa,w0,PMat);
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
            EigenQc(1,-1,&S,&dS,&dN,NULL,NULL,NULL,com.kappa,com.rK[j],PMat);
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
            EigenQc(1,-1,&S,&dS,&dN,NULL,NULL,NULL,com.kappa,w0,PMat);
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
         EigenQc(0,-1,NULL,NULL,NULL,_Root[i],_UU[i],_VV[i],com.kappa,com.omega,PMat);
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
   int status, j,off, K=com.ncatG-(com.NSsites==NSbetaw || com.NSsites==NS02normal);
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
      for(j=status=0;j<com.ncatG;j++) if(com.rK[j]<-1e-9) status=1;
      if(!status && com.NSsites==NSbeta) 
         for(j=1;j<com.ncatG;j++) if(com.rK[j]+1e-7<com.rK[j-1]) status=1;
      if(status) {
         printf("\nwarning: DiscreteNSsites\nw's");
         matout(F0,par,1,(com.NSsites==7?2:4));
         matout(F0,com.rK,1,com.ncatG);
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
    double kappa, double omega, double Q[])
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
   int n=Nsensecodon[com.icode], i,j,k, ic1,ic2,aa1,aa2;
   int ndiff,pos=0,from[3],to[3];
   double rs0,ra0,rs,ra, c0[3],c[3],ts[3],tv[3], t,space[64*3], *ri=space;
   double *pi=(com.seqtype==AAseq?com.fb61:com.pi);
   double w, fit1,fit2;

   N_eigenQc++;
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
      if ((from[pos]+to[pos]-1)*(from[pos]+to[pos]-5)==0)  Q[i*n+j]=kappa;
      c0[pos]+=t*Q[i*n+j];

      aa1=GenetCode[com.icode][ic1]; 
      aa2=GenetCode[com.icode][ic2];

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

      if ((from[pos]+to[pos]-1)*(from[pos]+to[pos]-5)==0) ts[pos]+=t*Q[i*n+j];
      else                                                tv[pos]+=t*Q[i*n+j];

      c[pos]+=t*Q[i*n+j];
      Q[j*n+i]=Q[i*n+j];

   }  /* for (i,j) */

   if(getstats) {
      if(com.NSsites) Qfactor_NS+=freqK_NS * (rs+ra);
      /* if(noisy>2) FOR(i,3) printf("\nts/tv[pos %d] = %.6f", i+1, ts[i]/tv[i]);*/
      rs0=rs;
      t=(rs0+ra0);  rs0/=t;  ra0/=t;   *S=rs0*3*com.ls; 
      if(!com.NSsites && branchl>0) {  /* calculates dS & dN */
         t=(rs+ra); rs/=t; ra/=t;
         *dS=branchl*rs/(3*rs0);  *dN=branchl*ra/(3*ra0);
      }
      else if (com.NSsites) { *dS=rs/(3*rs0);  *dN=ra/(3*ra0); }
   }
   else {  /* get Root, U, & V */
      if (com.seqtype==AAseq) return (0);

      FOR(i,n) FOR(j,n) Q[i*n+j]*=pi[j];
      
      for (i=0,t=0; i<n; i++) 
        { Q[i*n+i]=-sum(Q+i*n,n); t-=pi[i]*Q[i*n+i]; }
      if(com.NSsites)  t=1/Qfactor_NS;  /* Qfactor_NS calculated in SetParameter */
      FOR(i,n*n) Q[i]/=t;
      if (eigen (1,Q,n,Root,ri,U,V,space+n)) error2("eigenQc err.");
      xtoy (U, V, n*n);
      matinv (V, n, n, space);
   }
   return (0);
}


int EigenQaa (FILE *fout, double Root[], double U[], double V[], double rate[])
{
   int naa=20, i,j,k;
   double Q[20*20], mr, t=0, space[20*3];
   char aa3[4]="";

   FOR (i,naa*naa) Q[i]=0;
   switch (com.model) {
   case (Poisson)   : case (EqualInput) : 
      fillxc (Q, 1., naa*naa);  break;
   case (Empirical)   : case (Empirical_F):
      FOR(i,naa) FOR(j,i) Q[i*naa+j]=Q[j*naa+i]=com.daa[i*naa+j]/100;
      break;
   case (FromCodon): case (AAClasses): 
      EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, com.kappa,com.omega,PMat);
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
      Q[i*naa+i]=0; Q[i*naa+i]=-sum(Q+i*naa,naa);  mr-=com.pi[i]*Q[i*naa+i]; 
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

   if (eigen(1,Q,naa,Root,space,U,V,space+com.ncode))  error2("err: eigenQaa");
   xtoy(U,V,naa*naa);
   matinv(V,naa,naa,space);

   FOR(i,naa)  Root[i]=Root[i]/mr;

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
   FOR (i,nc) piaa[GenetCode[com.icode][FROM61[i]]] += pic[i];
   FOR (i,nc) {
      aai=GenetCode[com.icode][FROM61[i]];
      ti=pic[i]/piaa[aai];
      FOR (j, i) {
         aaj=GenetCode[com.icode][FROM61[j]];
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
      GetPMatBranch(PMat, t, ison);

      if(com.cleandata && nodes[ison].nson<1)       /* end node */
         for(h=pos0; h<pos1; h++)
            FOR(j,n) nodes[inode].lkl[h*n+j]*=PMat[j*n+com.z[ison][h]];
      else
         for(h=pos0; h<pos1; h++) 
            FOR(j,n) {
               for(k=0,t=0; k<n; k++)
                  t+=PMat[j*n+k]*nodes[ison].lkl[h*n+k];
               nodes[inode].lkl[h*n+j]*=t;
            }
   }        /*  for (ison)  */

   if(_nnodeScale && _nodeScale[inode]) {  /* scaling to avoid underflow */
      for(i=0,k=0; i<tree.nnode; i++)   /* k-th node for scaling */
         if(i==inode) break;  else if(_nodeScale[i]) k++;

      for(h=pos0; h<pos1; h++) {
         for(j=0,t=0;j<n;j++)
            if(nodes[inode].lkl[h*n+j]>t) t=nodes[inode].lkl[h*n+j];
         
         if(t<1e-300)  _nodeScaleF[k*com.npatt+h]=-500;
         else {  
            for(j=0;j<n;j++)  nodes[inode].lkl[h*n+j]/=t;
            _nodeScaleF[k*com.npatt+h]=log(t);
         }
      }
   }
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
      if (GenetCode[com.icode][i]==-1)  FROM64[i]=-1; 
      else                            { FROM61[n]=i; FROM64[i]=n++; }
   }
   if (n!=Nsensecodon[com.icode])  error2 ("# sense codons?");
   return (0);
}

int printfcode (FILE *fout, double fb61[], double space[])
{
/* space[64*2]
*/
   int i, n=Nsensecodon[com.icode];

   fprintf (fout, "\nCodon freq. x10000\n");
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
   int ig, ngroup, lt, il,is, i,b;
   char equal='.',*pz, c0[4],c[4];

   ngroup = (ls-1)/lline + 1;
   for (ig=0,FPN(fout); ig<ngroup; ig++)  {
      fprintf (fout,"%-8d\n", ig*lline+1);
      FOR (is,ns) {
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


int InitializeCodon (FILE *fout, double space[])
{
/* Count codons for genes, calculate site patterns and fpatts 
   sequences com.z[] are not coded and may contain ambiguity characters
   fbsg[ngene*ns*64] (space)
   set up com.pi[NCODE],com.piG[NGENE][64], according to com.codonf
   com.pi[] has freqs for all codon sites in the seqs if ngene>1.
   Ambiguity characters are not used in counting codon freqs.
*/
   int h, j,k, nc=NCODE, ig,lt, wname=20,miss=0;
   int nb[3],ib[3][4],ic, b[3], i0, status=0;
   double *fbsg=space, *fbs, fb3x4[3*4],fb3x4t[3*4],fb4[4],fb4t[4];
   double t=-1,sump1, *ppi;
   char str[4]="";

   PatternWeight (fout,space);
   for(j=0,zero(fbsg,com.ngene*com.ns*nc); j<com.ns; j++) {
      for(h=0,ig=0,lt=0; h<com.npatt; h++) {
         for(k=0;k<3;k++)  NucListall(com.z[j][h*3+k], &nb[k], ib[k]);
         k=nb[0]*nb[1]*nb[2];       /* k = no. compatible codons */
         if(k>1)  { miss=1; continue; }

         ic=ib[0][0]*16+ib[1][0]*4+ib[2][0];
         fbsg[j*com.ngene*nc+ig*nc+ic]+=com.fpatt[h];

/*       ambiguity codons are ignored
         FOR(i0,nb[0])  FOR(i1,nb[1])  FOR(i2,nb[2]) {
            ic=ib[0][i0]*16+ib[1][i1]*4+ib[2][i2];
            fbsg[j*com.ngene*nc+ig*nc+ic] += com.fpatt[h]/(double)k;
            if(GenetCode[com.icode][ic]==-1) {
               printf("stop codon in seq %3d: ",j+1);
               for(k=0;k<3;k++)  printf("%c",com.z[j][h*3+k]);  FPN(F0);
               status=-1; 
            }
         }
*/
         if ((lt+=(int)com.fpatt[h])==com.lgene[ig]) ig++;
      }
   }

   if((fbs=(double*)malloc((com.ns+1)*nc*sizeof(double)))==NULL) error2("fbs?");
   zero(fbs,(com.ns+1)*nc);
   FOR(j,com.ns) FOR(ig,com.ngene)
      FOR(k,nc) fbs[j*nc+k]+=fbsg[j*com.ngene*nc+ig*nc+k];
   FOR(j,com.ns) FOR(k,nc) fbs[com.ns*nc+k]+=fbs[j*nc+k];
   fputs("\nCodon usage in species\n",fout);
   printcums(fout, com.ns, fbs, com.icode);
   fputs("\nSums of codon usage counts across species\n",fout);
   printcu(fout, fbs+com.ns*64, com.icode);
   FOR(ig,com.ngene) for(j=0,zero(com.piG[ig],nc); j<com.ns; j++) 
      FOR(k,nc) com.piG[ig][k]+=fbsg[j*com.ngene*nc+ig*nc+k];
   if(com.ngene>1) {
      fputs("\nCodon usage in genes\n",fout);
      printcums (fout, com.ngene, &com.piG[0][0], com.icode);
   }

   fputs("\n\nCodon position x base (3x4) table for each sequence.",fout);
   for(j=0,zero(fb3x4t,12); j<com.ns; j++) {
      for(k=0,zero(fb3x4,12); k<nc; k++) {
         b[0]=k/16; b[1]=(k%16)/4; b[2]=k%4;
         FOR(i0,3)  fb3x4[i0*4+b[i0]]+=fbs[j*nc+k];
      }
      FOR(i0,3) abyx(1/sum(fb3x4+i0*4,4), fb3x4+i0*4, 4);
      FOR(k,12) fb3x4t[k]+=fb3x4[k]/com.ns;

      fprintf (fout,"\n\n#%d: %-*s", j+1,wname,com.spname[j]);
      FOR (i0,3) {
         fprintf (fout, "\nposition %2d:", i0+1);
         FOR(k,4) fprintf (fout,"%5c:%7.5f", BASEs[k], fb3x4[i0*4+k]);
      }
   }
   FOR(k,4) fb4t[k]=(fb3x4t[k]+fb3x4t[4+k]+fb3x4t[4*2+k])/3;
   for(j=0,fputs("\n\nAverage",fout); j<3; j++) {
      fprintf(fout, "\nposition %2d:", j+1);
      FOR(k,4) fprintf(fout,"%5c:%7.5f", BASEs[k], fb3x4t[j*4+k]);
   }
   fprintf(fout,"\n\n%-12s","Mean");
   FOR(k,4) fprintf(fout,"%5c:%7.5f", BASEs[k],fb4t[k]);  FPN(fout);

   if(com.ngene>1) {
      fputs("\n\nCodon position x base (3x4) table for each gene.",fout);
      for(ig=0; ig<com.ngene; ig++) {
         t=(ig?com.lgene[ig]-com.lgene[ig-1]:com.lgene[ig])*com.ns;
         for(k=0,zero(fb3x4,12); k<nc; k++) {
            b[0]=k/16; b[1]=(k%16)/4; b[2]=k%4;
            FOR(j,3) fb3x4[j*4+b[j]]+=com.piG[ig][k]/t;
         }
         for(j=0,fprintf(fout,"\n\nGene #%d:",ig+1); j<3; j++) {
            fprintf(fout, "\nposition %2d:", j+1);
            FOR(k,4) fprintf (fout,"%5c:%7.5f", BASEs[k],fb3x4[j*4+k]);
         }
      }
   }

/* edit com.pi & com.piG according to com.codonf */
   if (com.codonf==Fequal) {
      FOR (j,com.ngene) fillxc(com.piG[j],1./com.ncode,com.ncode);
      xtoy(com.piG[0],com.pi,com.ncode);
   }
   else if (com.codonf>=Fcodon) {
      FOR(ig,com.ngene) FOR(k,nc) if(FROM64[k]>-1)
         com.piG[ig][FROM64[k]]=com.piG[ig][k];
      FOR(k,nc) if(FROM64[k]>-1)
         com.pi[FROM64[k]]=fbs[com.ns*nc+k]/com.ls;
   }
   else if (com.codonf==F1x4 || com.codonf==F3x4) {
      for(ig=0; ig<com.ngene; ig++) {
         t=(ig?com.lgene[ig]-com.lgene[ig-1]:com.lgene[ig])*com.ns;
         for(k=0,zero(fb3x4,12); k<nc; k++) {
            b[0]=k/16; b[1]=(k%16)/4; b[2]=k%4;
            FOR(j,3) fb3x4[j*4+b[j]]+=com.piG[ig][k]/t;
         }
         FOR(k,4) fb4[k]=(fb3x4[k]+fb3x4[4+k]+fb3x4[4*2+k])/3;
         if (com.codonf==F1x4) {
            FOR(k,nc)  if(FROM64[k]>-1)
               com.piG[ig][FROM64[k]]=fb4[k/16]*fb4[(k/4)%4]*fb4[k%4];
       }
         else if (com.codonf==F3x4) {
           FOR(k,nc)  if(FROM64[k]>-1) 
              com.piG[ig][FROM64[k]]=fb3x4[k/16]*fb3x4[4+(k/4)%4]*fb3x4[8+k%4];
       }
   }  /* for(ig) */
   if (com.codonf==F1x4) {
      FOR(k,nc)  if(FROM64[k]>-1)
         com.pi[FROM64[k]]=fb4t[k/16]*fb4t[(k/4)%4]*fb4t[k%4];
      }
      else if (com.codonf==F3x4) {
        FOR(k,nc)  if(FROM64[k]>-1) 
           com.pi[FROM64[k]]=fb3x4t[k/16]*fb3x4t[4+(k/4)%4]*fb3x4t[8+k%4];
      }
   }
   for(k=com.ncode;k<nc;k++) { com.pi[k]=0;FOR(ig,com.ngene) com.piG[ig][k]=0;}
   FOR(ig,com.ngene) abyx(1/sum(com.piG[ig],nc),com.piG[ig],nc);
   abyx(1/sum(com.pi,nc),com.pi,nc);

   /* data manipulation if come codons are not observed in the data
      codonf=4: Nick's idea of adding 1 to each codon count; 
      codonf=5: Ziheng's idea of adding 0.5 to each missing codon count (0)
   */
   FOR(ig,com.ngene+1) {
      ppi=(ig==0?com.pi:com.piG[ig-1]);
      for(k=0,j=0; k<com.ncode; k++)
         if(ppi[k]==0) { j++; t=1./(com.ns*com.ls); }
      if(j && noisy) printf("\n%d codons are not observed in the data\n",j);
      if (j && com.codonf==4) /* Nick's idea */
         FOR(k,com.ncode) ppi[k]=(ppi[k]+t)/(1.+com.ncode*t);
      else if (j && com.codonf==5) { /* Ziheng's idea */
         for(k=0,sump1=1; k<com.ncode; k++)  
            if(ppi[k]==0) { ppi[k]=0.5*t; sump1-=0.5*t; }
         FOR(k,com.ncode) if(ppi[k]>0.5*t) ppi[k]*=sump1;
      }

if(fabs(1-sum(ppi,com.ncode))>1e-6) error2("ppi!=1.");

   }

   
   free(fbs);
   if(miss) fputs("\n(Ambiguity codons are not used to calculate freqs.)\n",fout);

   if(com.verbose && com.ngene==1) {
      fputs("\n*Codon frequencies under model, for evolver.\n",fout); 
      FOR(j,64) {
        fprintf(fout,"%12.6f",GenetCode[com.icode][j]==-1?0:com.pi[FROM64[j]]);
        if((j+1)%4==0) FPN(fout);
      }
   }

   if(com.cleandata) {  /* transform the sequences if cleandata */
      FOR(j,com.ns) {
         if(transform(com.z[j],com.npatt*3,1,0)) error2("strange??");
         FOR(h,com.npatt) {
            b[0]=com.z[j][h*3]; b[1]=com.z[j][h*3+1]; b[2]=com.z[j][h*3+2];
            if(FROM64[k=b[0]*16+b[1]*4+b[2]]==-1) {
               printf("\nstop codon %s in seq. %d\n", getcodon(str,k),j+1);
               exit(-1);
           }
           com.z[j][h]=(char)FROM64[k];
         }
         com.z[j]=(char*)realloc(com.z[j],com.npatt);
      }
   }
   return(status);
}


int AA2Codonf(double faa[20], double fcodon[])
{
/* get codon freqs from amino acid freqs, assuming equal freq. for each syn
   codon.  Used in codon-based amino acid substitution models.
*/
   int ic, iaa, i, NCsyn[20];

   FOR(i,20) NCsyn[i]=0;
   FOR(ic,64) if((iaa=GenetCode[com.icode][ic])!=-1) NCsyn[iaa]++;
   zero(fcodon, 64);
   for(ic=0; ic<Nsensecodon[com.icode]; ic++) {
      iaa=GenetCode[com.icode][FROM61[ic]];
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
   char codon[2][3];
   int is,js, k,i0,h, wname=20, status=0, ndiff,nsd[4];
   int nb[3],ib[3][4], missing;
   double ns,na, nst,nat, S,N, St,Nt, dS,dN,dN_dS,y, bigD=3, ls1;

   fputs("\n\nNei & Gojobori 1986. dN/dS (dN, dS)",fout);
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
      fprintf(ft,   "%-*s ",wname,com.spname[is]);
      FOR (js,is) {
         FOR(k,4) nsd[k]=0;
         ls1=(com.cleandata?com.ls:0);
         for (h=0,nst=nat=S=N=0; h<com.npatt; h++)  {
            if(com.cleandata)
               FOR(i0,2) getcodon(codon[i0],FROM61[com.z[i0==0?is:js][h]]);
            else {
               FOR(i0,2) FOR(k,3) codon[i0][k]=com.z[i0==0?is:js][h*3+k];
               for(i0=0,missing=0;i0<2;i0++) {
                  FOR(k,3) NucListall(codon[i0][k],&nb[k],ib[k]);
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
   if(status) fprintf (fout, "NOTE: -1 means the method is inapplicable.\n");
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
      if (fabs(1-sum(com.pi,20))>1e-4) {
         printf("\nSum of freq. = %.6f != 1 in aaRateFile\n",sum(com.pi,20)); 
         exit(-1);
      }
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
      ic1=GenetCode[com.icode][ic1];
      ic2=GenetCode[com.icode][ic2];
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
   if (com.codonf==Fcodon) return 0;

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
   if(!com.fix_kappa) com.kappa=x[1];
   if(!com.fix_omega) com.omega=x[1+!com.fix_kappa];
   if(!com.fix_kappa||!com.fix_omega)
      EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, com.kappa,com.omega,PMat);

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
   EigenQc(1,t,&S,&dS,&dN,NULL,NULL,NULL,com.kappa,omega,PMat);

   eh=(t+1)*Small_Diff;
   EigenQc(1,t+eh,&S,&dS1,&dN1,NULL,NULL,NULL,com.kappa,omega,PMat);
   EigenQc(1,t-eh,&S,&dS2,&dN2,NULL,NULL,NULL,com.kappa,omega,PMat);
   JacobiSN[0*np+0]=(dS1-dS2)/(2*eh);
   JacobiSN[1*np+0]=(dN1-dN2)/(2*eh);
  
   eh=(omega+1)*Small_Diff;
   EigenQc(1,t,&S,&dS1,&dN1,NULL,NULL,NULL,com.kappa,omega+eh,PMat);
   EigenQc(1,t,&S,&dS2,&dN2,NULL,NULL,NULL,com.kappa,omega-eh,PMat);
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
   int n=com.ncode, is,js,j,k,h, i0,np, wname=20;
   int nb[3],ib[3][4],ic[2], missing=0;
   double x[3]={1,1,1}, xb[3][2], lnL, e=1e-6, *var=space+NP, S,dS,dN;
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
   if(!com.fix_kappa)  { xb[1][0]=kappab[0]; xb[1][1]=kappab[1]; }
   if(!com.fix_omega)  { k=1+!com.fix_kappa; xb[k][0]=omegab[0]; xb[k][1]=omegab[1]; }

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
                  FOR(k,3) NucListall(codon[i0][k],&nb[k],ib[k]);
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

         np=com.np=(com.ntime=1)+!com.fix_kappa+!com.fix_omega;  NFunCall=0;
         x[0]=SeqDistance[is*(is-1)/2+js];  /* NG86 estimate as initial */
         if(x[0]<0.01||x[0]>3) x[0]=0.1; 

         if(!com.fix_kappa) x[1]=2;
         if(!com.fix_omega) 
            { k=1+!com.fix_kappa; if(x[k]<0.05||x[k]>5) x[k]=.5;}
         if(noisy>=9) {
            FPN(F0);  FOR(k,np) printf(" %12.6f",x[k]); FPN(F0);
            FOR(k,np) printf(" %12.6f",xb[k][0]); FPN(F0);
            FOR(k,np) printf(" %12.6f",xb[k][1]); FPN(F0);
         }
         if(com.fix_kappa && com.fix_omega)  
            EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, com.kappa,com.omega,PMat);
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
         EigenQc(1,x[0],&S,&dS,&dN, NULL,NULL,NULL,com.kappa,com.omega,PMat);
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
               EigenQc(1,y[0],&S,&dS1,&dN1,NULL,NULL,NULL,com.kappa,com.omega,PMat);
               y[k] -= 2*eh;
               if(!com.fix_kappa) com.kappa=y[1];
               com.omega=y[1+!com.fix_kappa];
               EigenQc(1,y[0],&S,&dS2,&dN2,NULL,NULL,NULL,com.kappa,com.omega,PMat);

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
      iaa=GenetCode[com.icode][FROM61[iaa]];
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
   int  h,hp, ir, i,it=0, k=-1, refsp=0,iaa;
   double  lnL=0, fh, cutoff=0.5, mw=-1, w2=x[com.np-1],psel=0;
   char  *sig;
   int  ncolors=3;
   char *colors[3]={"gray", "blue", "red"};

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
      fprintf(frst, " (%1d)", it+1);
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
#endif
            }
	    }
         FPN(fout);
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
   fprintf (frst,"\n\nlnL=%14.6f\n", -lnL);
   return (0);
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





#if 0  /* routines for testing codon-models */

int GetCategoryQc (char z[NS])
{
/* the category ID for a codon site with z[NS], transformed
   classified into 19 categories 
*/
   int i,j, icat, ncat=19, it, b[NS][3], nbase[3], markb[4];

   puts("\nDo not work with missing data, GetCategoryQc.");
   for (j=0; j<com.ns; j++) {
      it=FROM61[(int)z[j]];  b[j][0]=it/16; b[j][1]=(it/4)%4; b[j][2]=it%4;
   }
   FOR (i,3) {
      FOR (j,4) markb[j]=0;
      FOR (j,com.ns) markb[b[j][i]]=1;
      nbase[i]=markb[0]+markb[1]+markb[2]+markb[3]-1;
   }
   if(nbase[1]>=2) icat=ncat-1;
   else {
      if(nbase[0]>2) nbase[0]=2;  if(nbase[2]>2) nbase[2]=2;
      icat = nbase[1]*9+nbase[0]*3+nbase[2];
   }
   return (icat);
}

int TestModelQc (FILE * fout, double x[])
{
/* Test the Qc model, slower than simulations
*/
   char z[NS];
   int h, npatt, it, icat, j, nodeb[NS], imposs;
   int n=Nsensecodon[com.icode], isum, nsum, ncat=19;
   double  fh, y, nobs[19], pexp[19], Pt[8][NCODE*NCODE];

   puts("\nDo not work with missing data, GetCategoryQc.");
   puts("\ntest Qc..\n");
   for (h=0,zero(nobs,ncat); h<com.npatt; h++) {
      for (j=0; j<com.ns; j++) z[j]=com.z[j][h]-1;
      icat = GetCategoryQc(z);
      nobs[icat]+=com.fpatt[h];
   }
   FOR (j,ncat) 
      printf("cat #%4d: %4d%4d%4d%6.0f\n", j+1,j/9+1,(j/3)%3+1,j%3+1,nobs[j]);

   if (com.ns>5 || com.alpha || com.ngene>1)
      error2 ("TestModelQc: ns>5 || alpha>0.");
   if (SetParameters (x)) puts ("\npar err..");
   for (j=0,npatt=1; j<com.ns; j++)  npatt*=n;
   for (isum=0,nsum=1; isum<tree.nnode-com.ns; nsum*=n,isum++) ;
   printf("\nTest Qc: npatt = %d\n", npatt);
   FOR (j, tree.nbranch) 
      PMatUVRoot (Pt[j], nodes[tree.branches[j][1]].branch, n, U, V, Root);

   for (h=0,zero(pexp,ncat); h<npatt; h++) {
      for (j=0,it=h; j<com.ns; nodeb[com.ns-1-j]=it%n,it/=n,j++) ;
      for (j=0,imposs=0; j<com.ns; j++) 
         { z[j]=nodeb[j];  if (com.pi[(int)z[j]]==0) imposs=1; }
      if (imposs) continue;
      
      if ((icat=GetCategoryQc(z)) == ncat-1) continue;
      if ((h+1)%100==0) 
         printf("\rTest Qc:%9d%4d%9.2f%%", h+1, icat, 100.*(h+1.)/npatt);

      for (isum=0,fh=0; isum<nsum; isum++) {
         for (j=0,it=isum; j<tree.nbranch-com.ns+1; j++)
            { nodeb[com.ns+j]=it%n; it/=n; }
         for (j=0,y=com.pi[nodeb[tree.root]]; j<tree.nbranch; j++) 
            y*=Pt[j][nodeb[tree.branches[j][0]]*n+nodeb[tree.branches[j][1]]];
         fh += y;
      }
      if (fh<=0) {
         matout (F0, x, 1, com.np);
         printf ("\a\ntest Qc: h=%4d  fh=%9.4f \n", h, fh);
      }
      pexp[icat]+=fh;
   }    
   pexp[ncat-1]=1-sum(pexp,ncat-1);

   FOR (j,ncat) 
      fprintf(fout, "\ncat # %4d%4d%4d%4d%6.0f%10.5f%10.2f", 
         j+1, j/9+1, (j/3)%3+1, j%3+1, nobs[j], pexp[j], com.ls*pexp[j]);
   return (0);
}

#endif

#if (DSDN_MC || DSDN_MC_SITES)

void SimulateData2s61(void)
{
/* This simulates two codon sequences and analyze using ML (GY94).
   It generates site pattern freqs and then samples from them
   to generate the seq data.  Codons are coded as 0,1,...,60.  There
   is another routine of a similar name in the file dsdn.c where the
   codons are coded as 0,1,...,63.  The two routines should not be
   mixed.
   Note that com.pi[] is changed in the analysis but not reused to 
   calculate Efij[]
   Ifdef (DSDN_MC_SITES), the data will be simulated with the NSsites models
   but analysed assuming one omega for all sites, so the model is wrong.
*/
   char infile[32]="in.codon2s", seqfile[32]="codonseq.tmp",str[4]="";
   FILE *fseq, *fin;
   int ir,nr=100, ip, i,j,k,h, n=Nsensecodon[com.icode];
   int npatt0=n*(n+1)/2, nobs[NCODE*NCODE];
   int il,nil, ls[50]={0,200,500};
   double y, x[6]={1,1,1},xb[6][2], S,dN,dS,dNt,dSt,om,lnL;
   double t0,kappa0,omega0=.5,pi0[NCODE], mx[6],vx[6],mse[6]; /* t,k,w,dS,dN */
   double Efij[NCODE*(NCODE+1)/2], space[50000];

   com.icode=0; com.seqtype=1; com.ns=2;
   com.ncode=n; com.cleandata=1; setmark_61_64 ();
   FOR(j,com.ns) com.z[j]=(char*) malloc(npatt0*sizeof(char));
   if(com.z[com.ns-1]==NULL) error2 ("oom z");
   if((com.fpatt=(double*)malloc(npatt0*sizeof(double)))==NULL)
   error2("oom fpatt");
   FOR(j,3) { xb[j][0]=.0001; xb[j][1]=99; }

#if (DSDN_MC_SITES)
   strcpy(infile,"in.codon2sSites");
#endif
   printf("\nTwo codon seq. simulation for ML (GY94), input from %s\n",infile);
   if((fin=fopen(infile,"r"))==NULL) error2("data file read error");
   
   fscanf (fin,"%d%d%d%d", &k,&nr, &com.codonf,&nil);  SetSeed(k);
   printf("\n%d replicates, %s model for analysis\nLc:",
      nr,codonfreqs[com.codonf]);
   FOR(il,nil) fscanf(fin, "%d", &ls[il+1]);
   matIout(F0,ls+1,1,nil);
   for(i=0,k=0; i<NCODE; i++) {
      fscanf(fin,"%lf",&y);
      if(GenetCode[com.icode][i]>-1) pi0[k++]=y;
      else if(y!=0)                  error2("stop codon freq !=0");
   }
   printf("sum pi = 1 = %.6f\n", sum(pi0,n));

   for(ip=0; ip<99; ip++) {
      fscanf(fin, "%lf%lf", &t0,&kappa0);
      if(t0<0) exit(0);
      printf("\n\nParameter set %d\nt0 =%.2f  kappa0 =%.2f\n",ip+1,t0,kappa0);
      fprintf(frst,"\n\nParameter set %d\nt0 =%.2f  kappa0 =%.2f\n",ip+1,t0,kappa0);

      FOR(j,n) com.pi[j]=pi0[j];  com.ls=1;
#if (DSDN_MC_SITES)
      com.NSsites=3;
      fscanf(fin,"%d", &com.ncatG);
      FOR (i,com.ncatG) fscanf(fin,"%lf", &com.freqK[i]);
      FOR (i,com.ncatG) fscanf(fin,"%lf", &com.rK[i]);

      printf("\nSite classe model (K=%d)\np: ",com.ncatG);
      FOR(i,com.ncatG) printf("%7.4f",com.freqK[i]);
      printf("\nw: "); FOR(i,com.ncatG) printf("%7.4f",com.rK[i]); FPN(F0);
      fprintf(frst,"\nSite classe model (K=%d)\np: ",com.ncatG);
      FOR(i,com.ncatG) fprintf(frst,"%7.4f",com.freqK[i]);
      fputs("\nw: ",frst); FOR(i,com.ncatG) fprintf(frst,"%7.4f",com.rK[i]); FPN(frst);

      if(1-sum(com.freqK,com.ncatG)) error2("freqs do not sum to 1");
      for(j=0,Qfactor_NS=0,dS=dN=0; j<com.ncatG; j++) {
         freqK_NS=com.freqK[j];
         EigenQc(1,1,&S,&dSt,&dNt,NULL,NULL,NULL,kappa0,com.rK[j],PMat);
         dS+=freqK_NS*dSt;  dN+=freqK_NS*dNt;
      }
      Qfactor_NS=1/Qfactor_NS;
      om=(dS>0?dN/dS:-1);  dS*=t0*Qfactor_NS;  dN*=t0*Qfactor_NS;

#else
      fscanf(fin,"%lf", &omega0);
      EigenQc(1,t0,&S,&dS,&dN, NULL,NULL,NULL,kappa0,omega0,space);
      om=omega0;
#endif
      printf("\nCorrect values"); 
      printf("\nS%%=%7.4f  dS=%7.4f  dN=%7.4f  w=%7.4f\n",S/3,dS,dN,om);
      fprintf(frst,"\nCorrect values");
      fprintf(frst,"\nS%%=%7.4f  dS=%7.4f  dN=%7.4f  w=%7.4f\n",S/3,dS,dN,om);
      
      /* calculate Efij[], the site pattern probabilities */
      FOR(j,n) com.pi[j]=pi0[j];
#if (DSDN_MC_SITES)
      com.NSsites=3;
      for(k=0,zero(Efij,npatt0); k<com.ncatG; k++) {
         EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, kappa0,com.rK[k],PMat);
         PMatUVRoot(PMat, t0, n, U, V, Root);
         for(i=0,h=0;i<n;i++) for(j=0;j<=i;j++) {
            y=com.pi[i]*PMat[i*n+j];
            Efij[h++] += (i==j?y:2*y) * com.freqK[k];
         }
      }
      com.NSsites=0;
#else
      EigenQc(0,-1,NULL,NULL,NULL,Root, U, V, kappa0, omega0, PMat);
      PMatUVRoot (PMat, t0, n, U, V, Root);
      for(i=0,h=0;i<n;i++) for(j=0;j<=i;j++) { /* why for each il? */
         y=com.pi[i]*PMat[i*n+j];
         Efij[h++]=(i==j?y:2*y);
      }
#endif
      for(i=h=0,com.ls=1,com.npatt=npatt0;i<n;i++) for(j=0;j<=i;j++) {
         com.z[0][h]=(char)i; com.z[1][h]=(char)j;
         com.fpatt[h]=Efij[h];  h++;
      }
      if(fabs(1-sum(Efij,npatt0))>1e-6) error2("sum Efij != 1");

      for(il=0; il<nil+1; il++) {
         com.ls=ls[il];
         if(com.ls==0) {
            puts("\nML estimates from infinite data"); 
            com.ls=1;
            x[0]=t0*rndu(); x[1]=kappa0; x[2]=omega0*rndu();
            if(com.codonf<Fcodon) GetCodonFreqs(com.pi);
            ming2(NULL,&lnL,lfun2dSdN,NULL,x,xb, space,1e-10,3);
            printf("lnL = %.6f\n",-lnL);
            EigenQc(1,x[0],&S,&dS,&dN, NULL,NULL,NULL,x[1],x[2],space);
            printf("S%%=%7.4f  dS=%7.4f  dN=%7.4f  w=%7.4f\n",S/3,dS,dN,x[2]);
            fprintf(frst,"ML estimates from infinite data\nt=%7.4f  k=%7.4f",x[0],x[1]);
            fprintf(frst,"  S%%=%7.4f  dS=%7.4f  dN=%7.4f  w=%7.4f\n",S/3,dS,dN,x[2]);

            for(h=1;h<npatt0; h++) Efij[h]+=Efij[h-1];
            puts("\nt & k & w & dS & dN");
            fputs("\nLc & t & k & w & dS & dN\n",frst);  fflush(frst);
            continue;
         }

         printf("\nls = %d\n", com.ls);
         for(ir=0,zero(mx,6),zero(vx,6),zero(mse,6); ir<nr; ir++) {
            MultiNomial(com.ls, npatt0, Efij, nobs, NULL);
            for(i=0,com.npatt=0,zero(com.pi,n);i<n;i++) for(j=0;j<=i;j++)
               if(nobs[k=i*(i+1)/2+j]) {
                  com.z[0][com.npatt]=i; com.z[1][com.npatt]=j;
                  com.fpatt[com.npatt++]=nobs[k];
               }
            for(i=0,zero(com.pi,n); i<com.npatt; i++) {
               com.pi[com.z[0][i]]+=com.fpatt[i]/(2.*com.ls);
               com.pi[com.z[1][i]]+=com.fpatt[i]/(2.*com.ls);
            }
            if(com.codonf<Fcodon) GetCodonFreqs(com.pi);

            x[0]=t0;  x[1]=kappa0; x[2]=omega0;
            /* printf("\nlnL=%9.6f\n",-lfun2dSdN(x,3)); */
            ming2((noisy?F0:NULL),&lnL,lfun2dSdN,NULL,x,xb, space,1e-7,3);
            EigenQc(1,x[0],&S,&x[3],&x[4], NULL,NULL,NULL,x[1],x[2],space);
            FOR(j,5) {
               vx[j]+=(x[j]-mx[j])*(x[j]-mx[j]);
               mx[j]=(mx[j]*ir+x[j])/(ir+1.);
            }
            mse[0]+=square(x[2]-omega0);
            printf("\r%4d%8.4f%8.4f%8.4f  %8.4f%8.4f%8.4f%8.4f%8.4f",
                   ir+1,x[0],x[1],x[2],mx[0],mx[1],mx[2],mx[3],mx[4]);
#if 0
            if(ir==9) {
               if((fseq=fopen(seqfile,"w"))==NULL) error2("seq file err");
               fprintf(fseq,"%6d %6d\n", com.ns,com.ls*3);
               for(i=0;i<2;i++,FPN(fseq),fflush(fseq)) {
                  fprintf(fseq,"seq.%-5d  ", i+1);
                  FOR(h,com.npatt) FOR(k,(int)com.fpatt[h]) 
                     fprintf(fseq,"%s", getcodon(str,FROM61[com.z[i][h]]));
               }
               fclose(fseq); exit(0);
           }
#endif
         }
         if(nr>1) { FOR(j,5) vx[j]=sqrt(vx[j]/(nr-1.)/nr); mse[0]=sqrt(mse[0]/nr); }
         fprintf(frst,"%4d ", com.ls);
         FOR(i,5) fprintf(frst,"%7.4f +%7.4f", mx[i],vx[i]);  FPN(frst);
      }  /* for (ii) */
   }   /* for(ip) */
   exit(0);
}


void Ina(void)
{
/* This simulates two codon sequences and analyze them using Ina's method.  
   Ina's program is modified to output result in Ina1.tmp.  Consistency
   analysis is done by generating long sequences.
   Note that com.pi[] is not changed in the analysis, which is done outside
   the program.  
*/
   char seqfile[32]="codonseq.tmp",tmpfile[32]="Ina1.tmp",str[4]="";
   FILE *fseq, *ftmp;
   int ip,ir,nr=500, i,j,k,h, n=Nsensecodon[com.icode];
   int npatt0=n*(n+1)/2, nobs[NCODE*NCODE];
   int il,nil=1, ls[]={500,100,200,300,400,500,600,800,1000}, fcodon=1;
   double y, t=.5,f3x4[12], x[3]={1,1,1}, S,dS,dN;
   double t0=1,kappa0=1,omega0=1, mx[6],vx[6],mse[6]; /* t,k,w,dS,dN */
   double Efij[NCODE*NCODE], space[50000];
   double f3x4_data[][3*4]={
                            {0.25, 0.25, 0.25, 0.25,
                             0.25, 0.25, 0.25, 0.25,
                             0.25, 0.25, 0.25, 0.25},

                            {0.20517, 0.28293, 0.30784, 0.20406, /* mt pri */
                             0.40979, 0.27911, 0.18995, 0.12116,
                             0.15105, 0.43290, 0.37123, 0.04482},

                            {0.19020, 0.16201, 0.36655, 0.28124, /* hiv */
                             0.28889, 0.18805, 0.30179, 0.22127,
                             0.24875, 0.16894, 0.36822, 0.21410},

                            {0.14568, 0.24519, 0.33827, 0.27086,
                             0.35556, 0.18765, 0.24049, 0.21630, 
                             0.26444, 0.25728, 0.21012, 0.26815} /* lysinNew*/
                           };

   puts("\nThis simulates data and analyses them using Ina95.");

printf ("fcodon? ");
scanf ("%d", &fcodon);

   FOR(j,12) f3x4[j]=f3x4_data[fcodon][j];
   for(j=0,h=0,y=0; j<64; j++) {
      if (GenetCode[com.icode][j]==-1) continue;
      com.pi[h]=f3x4[j/16]*f3x4[4+(j%16)/4]*f3x4[8+j%4];
      y+=com.pi[h++];
   }
   FOR(j,n) com.pi[j]/=y;
   printf("fcodon: %d\n",fcodon);
   matout(frst,f3x4,3,4);
   com.icode=0; com.seqtype=1; com.ns=2; com.ls=1; npatt0=n*(n+1)/2;
   com.ncode=n; setmark_61_64 ();
   FOR(j,com.ns) com.z[j]=(char*) malloc(npatt0*sizeof(char));
   if(com.z[com.ns-1]==NULL) error2 ("oom z");
   if((com.fpatt=(double*)malloc(npatt0*sizeof(double)))==NULL)
      error2("oom fpatt");

   printf("\nInfinite sequences.\nsum pi=1=%.6f\n",sum(com.pi,NCODE));
   noisy=0;  FOR(i,6) x[i]=0;

   FOR(ip,99) {
      printf("\nt0 & kappa0 & omega0? ");
      scanf("%lf%lf%lf", &t0,&kappa0,&omega0);
      if(t0<0) exit(0);
      printf("t0 =%.2f & kappa0 =%.2f & omega0 =%.2f\n",t0,kappa0,omega0);
      fprintf(frst, "\nt & k & w: %8.2f%8.2f%8.2f\n\n", t0,kappa0,omega0);
      EigenQc(1,t0,&S,&dS,&dN, NULL,NULL,NULL,kappa0,omega0,space); 
      fprintf(frst,"\nS/(S+N)=%7.4f  dS=%7.4f  dN=%7.4f\n",S/3,dS,dN);
      fflush(frst);
      fputs("Lc & t & k & w & dS & dN\n",frst);
   
      EigenQc(0,-1,NULL,NULL,NULL,Root, U, V, kappa0, omega0, PMat);
      PMatUVRoot (PMat, t0, n, U, V, Root);
      for(i=0,h=0;i<n;i++) for(j=0;j<=i;j++) {
         y=com.pi[i]*PMat[i*n+j];
         Efij[h++]=(i==j?y:2*y);
      }
      for(i=h=0,com.ls=1,com.npatt=npatt0;i<n;i++) for(j=0;j<=i;j++) {
         com.z[0][h]=(char)i; com.z[1][h]=(char)j;
         com.fpatt[h]=Efij[h];  h++;
      }
      for(h=1;h<npatt0; h++) Efij[h]+=Efij[h-1];
      if(fabs(1-Efij[npatt0-1])>1e-6) puts("Sum p_ij != 1.");
      for(il=0; il<nil; il++) {
   
         com.ls=ls[il];
         printf("\nls = %d\n", com.ls);
         for(ir=0,zero(mx,6),zero(vx,6),zero(mse,6); ir<nr; ir++) {
            printf("\r%4d", ir+1);
            MultiNomial (com.ls, npatt0, Efij, nobs, NULL);
            for(i=0,com.npatt=0;i<n;i++) for(j=0;j<=i;j++)
               if(nobs[k=i*(i+1)/2+j]) {
                  com.z[0][com.npatt]=i; com.z[1][com.npatt]=j; 
                  com.fpatt[com.npatt++]=nobs[k];
               }
            if((fseq=fopen(seqfile,"w"))==NULL) error2("seq file err");
            fprintf(fseq,"> %6d %6d\n", com.ns,com.ls*3);
            for(i=0;i<2;i++,FPN(fseq),fflush(fseq)) {
               fprintf(fseq,"seq.%-5d  ", i+1);
               FOR(h,com.npatt) FOR(k,(int)com.fpatt[h]) 
                  fprintf(fseq,"%s", getcodon(str,FROM61[com.z[i][h]]));
            }
            fclose(fseq);
            if(com.ls>2000) system("Ina1Large codonseq.tmp >t");
            else            system("Ina1 codonseq.tmp >t");
            if((ftmp=fopen(tmpfile,"r"))==NULL) error2("tmp file err");
            if(fscanf(ftmp,"%lf%lf%lf",&x[0],&x[1],&x[2]) !=3) 
               error2("reading tmpf");
            fclose(ftmp);
            FOR(j,5) {
               vx[j]+=(x[j]-mx[j])*(x[j]-mx[j]);
               mx[j]=(mx[j]*ir+x[j])/(ir+1.);
            }
            mse[0]+=square(x[0]-omega0);

            printf("%7.4f%7.4f%7.4f  %7.4f%7.4f%7.4f%7.4f%7.4f",
                   x[0],x[1],x[2],mx[0],mx[1],mx[2],mx[3],mx[4]);

/*            fprintf(frst1,"%7.4f%7.4f%7.4f  %7.4f%7.4f%7.4f%7.4f%7.4f\n",
                   x[0],x[1],x[2],mx[0],mx[1],mx[2],mx[3],mx[4]);
*/
         }
         if(nr>1) { FOR(j,5) vx[j]=sqrt(vx[j]/(nr-1.)/nr); mse[0]=sqrt(mse[0]/nr); }
         fprintf(frst,"%4d ", com.ls);
         FOR(i,5) fprintf(frst,"%7.3f +%7.4f", mx[i],vx[i]);
         FPN(frst); fflush(frst);
         fprintf(frst1,"%6d%6d %7.2f%7.2f%7.2f: %8.3f +%7.3f\n",
             com.ls,nr, t0,kappa0,omega0, mx[0],mse[0]);  fflush(frst1);
      }    /* for (il) */
   }       /* for (ip) */
   exit(0);
}

#endif

#if 0

int mergeSeqs(FILE*fout)
{
/* This concatenate multiple genes (data sets) for the same set of species
   into one file of a long gene.  Used to process Anne Yoders' alignment.
*/

   char *filenames[12]={"NADH1.fin","NADH2.fin","COI.fin","COII.fin","ATPase8.fin",
        "ATPase6.fin","COIII.fin","NADH3.fin","NADH4L.fin","NADH4.fin",
        "NADH5.fin", "Cytb.fin"};

   int ns0=32, nfile=12, ifile, ls0, lswhole=20000, i,h, lgene0[32];
   char *z0[32], *spname0[32]={"Artibeus", "B.musculus", "B.physalus", "Bos", 
      "Canis", "Cavia", "Ceratother", "Dasypus", "Didelphis", "E.asinus", 
      "E.caballus","Erinaceus", "Felis", "Glis", "Gorilla", "Halichoeru", "Homo",
      "Hylobates", "Macropus", "Mus", "Ornithorhy", "Oryctolagu", "Ovis",
      "P.paniscus", "P.troglody", "Papio", "Phoca", "P.abelii",
      "P.pygmaeus", "Rattus", "Rhinoceros", "Sus"};
   FILE *fseq;

   noisy=0;
   FOR(i,ns0) if((z0[i]=(char*)malloc(lswhole*sizeof(char)))==NULL) 
      error2("oom z");
   for(ifile=0,ls0=0; ifile<nfile; ifile++) {
      printf("Reading data set %2d/%2d (%s)", ifile+1,nfile,filenames[ifile]);
      if((fseq=fopen (filenames[ifile],"r"))==NULL)  error2("file not found");
      ReadSeq(NULL,fseq,1);
      lgene0[ifile]=com.ls;  com.ls*=3;
      FOR(i,ns0) if(strcmp(spname0[i],com.spname[i])) error2("spname different"); 
      FOR(i,ns0)  FOR(h,com.ls) z0[i][ls0+h]=com.z[i][h];
      ls0+=com.ls;
      printf(" + %5d = %5d\n", com.ls, ls0);
   }
   fprintf(fout,"%6d %6d  G\nG %4d ", ns0,ls0,nfile);
   FOR(ifile,nfile) fprintf(fout, " %4d", lgene0[ifile]);  FPN(fout);

   for(i=0;i<ns0;i++,FPN(fout)) {
      fprintf(fout,"%-12s  ", spname0[i]);
      FOR(h,ls0) {
         fprintf(fout,"%c", z0[i][h]);
         if((h+1)%3==0) fprintf(fout," ");
      }
   }
   return(0);
}

int SlidingWindow(FILE*fout,double space[])
{
/* sliding window analysis, clean data only
*/
   int wlen=100, wstart, n=com.ncode, j, h;
   int npatt0=com.npatt; 
   char *z0[NS];
   double *fpatt0, pi0[NCODE*NCODE];

   if(!com.cleandata || com.ngene>1 || com.runmode)
      error2("clean data & one gene only for sliding window analysis");   
   FOR(j,com.ns) z0[j]=com.z[j];
   FOR(j,com.ns) if((com.z[j]=malloc(npatt0*sizeof(char)))==NULL) error2("oom z");
   if((fpatt0=(double*)malloc(com.npatt*sizeof(double)))==NULL) error2("oom fp");
   FOR(h,com.npatt) fpatt0[h]=com.fpatt[h];
   FOR(j,n*n) pi0[j]=com.pi[j];

   FOR(j,com.ns) FOR(h,com.npatt) com.z[j][h]=-1;  /* to crash if in error2 */

   puts("\nSliding window analysis.\nwindow size (# codons or amino acids)?" );
   scanf("%d",&wlen);
   if(wlen<10 || wlen>com.ls) error2("strange win size");
   for (wstart=0; wstart+wlen<com.ls; wstart++) {
      FOR(h,npatt0) com.fpatt[h]=0;
      for(h=wstart; h<wstart+wlen; h++)
         com.fpatt[com.pose[h]]++;

      for(h=0,com.npatt=0,zero(com.pi,n); h<npatt0;h++) if(com.fpatt[h]>0) {
         FOR(j,com.ns) com.z[j][com.npatt]=z0[j][h];
         com.fpatt[com.npatt]=com.fpatt[h];
         FOR(j,com.ns) com.pi[z0[j][h]]+=com.fpatt[h];
         com.npatt++;
      }
      com.posG[0]=0; com.posG[1]=com.npatt;
      abyx(1/sum(com.pi,n),com.pi,n);
      if(com.seqtype==CODONseq && com.codonf<3) 
         GetCodonFreqs(com.pi);

      printf("\nsites %3d -- %3d  npatt:%4d",wstart+1,wstart+wlen,com.npatt);
      fprintf(fout,"\nsites %3d -- %3d  %4d",wstart+1,wstart+wlen,com.npatt);
      fprintf(frst1,"sites %3d -- %3d  %4d",wstart+1,wstart+wlen,com.npatt);

      Forestry(fout);

   }
   FOR(h,com.npatt) com.fpatt[h]=fpatt0[h];
   xtoy(pi0,com.pi,n);
   free(fpatt0);  
   FOR(j,com.ns) { free(com.z[j]); com.z[j]=z0[j]; }

   return (0);
}

#endif
