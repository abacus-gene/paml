/* CODEML.c  (AAML.c & CODONML.c)
   Maximum likelihood parameter estimation for codon sequences (seqtype=1) 
                    or amino-acid sequences (seqtype=2)
                 Copyright, Ziheng Yang (March 1993 onwards)

                 cc -o codeml -fast codeml.c tools.o eigen.o -lm
                        codeml <ControlFileName>
*/
#include "tools.h"
#define NS            100
#define NBRANCH       (NS*2-2)
#define NNODE         (NS*2-1)
#define NGENE         100
#define LSPNAME       30
#define NCODE         64
#define NCATG         40
/*
#define NP            (NBRANCH*2+NGENE-1+2)
*/
#define NP            (NBRANCH+NGENE-1+189+2)

extern char NUCs[],AAs[],GenetCode[][64],Nsensecodon[];
extern int noisy, NFunCall, *ancestor;
extern double *SeqDistance;

int Forestry (FILE *fout, double space[]);
void DetailOutput (FILE *fout, double x[], double space[]);
int GetOptions (char *ctlf);
int testx (double x[], int np);
int SetxBound (int np, double xb[][2]);
int SetxInitials (double x[]);
int GetInitials (double x[]);
int SetParameters (double x[]);
int SetPGene (int igene, int _pi, int _UVRoot, int _alpha, double x[]);
double lfunAdG (double x[], int np);
double lfundG (double x[], int np);
double lfun (double x[], int np);
int PartialLikelihood (int inode, int igene);
int PMatJC69like (double P[], double t, int n);
int setmark_61_64 (void);
int printfcode (FILE *fout, double fb61[], double space[]);
int InitializeCodon (FILE *fout, double space[]);
int printsmaCodon (FILE *fout,char * z[],int ns,int ls,int simple);
int Getfcodon (double faa[20], double fcodon[]);
int DistanceMatAA (FILE *fout, double alpha);
int DistanceMatNG86 (FILE *fout, double alpha);
int GetDaa(FILE *fout, double daa[]);
int EigenQ20(FILE *fout, double Root[], double U[], double V[],double rate[]);
int EigenQ61(FILE *fout, double Root[], double U[], double V[], 
    double kappa, double omega, double Q[64*64]);
int Qcodon2aa (double Qc[], double pic[], double Qaa[], double piaa[]);
int SetAA1STEP (void);
int GetOmegaAA(int OmegaAA[]);
int TestModelQ61(FILE *fout, double x[]);
int Pairwise (FILE *fout, double space[]);
double lfunNSsites (double x[], int np);
int lfunNSsites_rate (FILE* fout, double x[], int np);

int NSQ61(double kappa, double pi[]);

struct common_info {
   char *z[NS], spname[NS][LSPNAME+1], seqf[64],outf[64],treef[64],daafile[64];
   int seqtype, ns, ls, ngene, posG[NGENE+1], lgene[NGENE], npatt,*pose;
   int runmode,clock,verbose,print, codonf,aaDist,model,NSsites;
   int icode, ncode, Mgene;
   int fix_rgene, fix_kappa, fix_omega, fix_alpha, fix_rho, getSE;
   int np, ntime, nrgene, nrate, nalpha, ncatG, maxnp;
   float *fpatt;
   double pi[NCODE], fb61[64], piG[NGENE][64];
   double kappa, omega, alpha, rho, rgene[NGENE];
   double freqK[NCATG], rK[NCATG], MK[NCATG*NCATG], daa[20*20], *chunk, *fhK;
   double (*plfun)(double x[],int np);
}  com;
struct TREEB {
   int  nbranch, nnode, origin, branches[NBRANCH][2];
   double lnL;
}  tree; 
struct TREEN {
   int father, nson, sons[NS], ibranch;
   double branch, divtime, *lkl;
}  nodes[2*NS-1];

int FROM61[64], FROM64[64];
double PMat[NCODE*NCODE],U[NCODE*NCODE],V[NCODE*NCODE],Root[NCODE],_BRANCHLEN;
double *XCOPY;
  /* communication between SetParameters & PartialLikelihood & EigenQ61.  
     Remove it? */
int LASTROUND=0;

double OMEGA_FIX=-1;  
/* fix the last dN/dS in the NSb, NS2 models with variable dN/dS ratios 
   for lineages.  Useful for testing whether w>1 for particular lineages. */

int N_OMEGA=-99, *OmegaBranch, OmegaAA[190], AA1STEP[190];

FILE *frub, *flfh, *frst;
enum {Fequal, F1x4, F3x4, Fcodon} CodonFreqs;
enum {NS1, NSb, NS2} CodonModels;   /* , AAClasses=7 */
enum {Poisson, EqualInput, Empirical, Empirical_F,
     FromCodon=6, AAClasses, REVaa_0=8, REVaa=9} AAModel;

char *codonfreqs[]={"Fequal", "F1x4", "F3x4", "Fcodon"};
char *codonmodels[]={"One dN/dS ratio", "free dN/dS Ratios", 
     "several dN/dS ratios", "", "", "", "","AAClasses"};
char *NSsitesmodels[]={"one-ratio","neutral model","selection model"};
char *aamodels[]={"Poisson", "EqualInput", "Empirical", "Empirical_F", "",
     "", "FromCodon", "AAClasses", "REVaa_0", "REVaa"};

#define CODEML
#include "treesub.c"
#include "treespace.c"

int main(int argc, char *argv[])
{
   FILE *fout, *fseq;
   char ctlf[32]="codeml.ctl", *pmodel;
   char *seqtypestr[3]={"CODONML", "AAML", "CODON2AAML"};
   char *Mgenestr[]={"diff. rate", "separate data", "diff. rate & pi", 
                     "diff. rate & kappa", "diff. rate & pi & kappa"};
   int i, s2=0, s3=0;
   double *space;

/*
printf ("\nrandom number seed? ");
scanf ("%d", &i);
SetSeed(i);
*/
   noisy=2;
   com.runmode=0;
   com.clock=0;       /* 1: clock, rooted tree;  0: no clock, unrooted tree  */
   com.fix_rgene=0; /* 0: estimate rate factors for genes */

   com.seqtype=AAseq;
   com.model=Empirical_F;  strcpy(com.daafile, "jones.dat");
   com.icode=0;
   com.nrate=0;
   com.fix_kappa=0;   com.kappa=1;    com.omega=2.1;
   com.fix_alpha=1;   com.alpha=0.;   com.ncatG=4;   /* alpha=0 := inf */
   com.fix_rho=0;    com.rho=0.4; 

   com.getSE=0;
   com.print=0;        com.verbose=1; 

   frub=fopen ("rub","w");   frst=fopen ("rst","w");
   flfh=fopen ("lfh", "w");

   if(argc>1) { strcpy(ctlf,argv[1]); printf("\nctlfile set to %s.\n",ctlf); }
   GetOptions (ctlf);
   if ((fout=fopen (com.outf, "w"))==NULL) error("outfile creation err.");
   if((fseq=fopen (com.seqf,"r"))==NULL)  {
       printf ("\n\nSequence file %s not found!\n", com.seqf);
       exit (-1);
   }
   ReadSeq ((com.verbose?fout:NULL), fseq, com.seqtype);
   if (com.ngene>1 && com.Mgene==1)  OutSeqsMgenes ();
   if (com.Mgene==1 && com.print) puts("Mgene & RateAncestor not well tested");

   if (com.seqtype==CODONseq && (com.model==NSb||com.model==NS2)) 
      if((OmegaBranch=(int*)malloc(NBRANCH*sizeof(int)))==NULL) error("oom");
   
   pmodel=(com.seqtype==CODONseq?codonmodels[com.model]:aamodels[com.model]);
   fprintf(fout,"%s %s Model: %s ",seqtypestr[com.seqtype-1],com.seqf,pmodel);
   if (com.seqtype==CODONseq||com.model==FromCodon) {
      if (com.fix_kappa)  fprintf (fout, " kappa = %.3f\n", com.kappa);
      if (com.fix_omega) fprintf (fout, " omega? = %.3f fixed\n", com.omega);
   }
   if ((com.seqtype==AAseq||com.seqtype==CODON2AAseq)
     &&(com.model==Empirical||com.model==Empirical_F))
      fprintf (fout, "(%s) ", com.daafile);
      
   /* if (com.nrate) fprintf(fout, "(nrate:%d)  ", com.nrate); */
   if (com.clock) fprintf (fout, " Clock  ");
   if (com.alpha && com.rho) fprintf (fout, "Auto-");
   if (com.alpha) fprintf (fout,"dGamma (ncatG=%d) ", com.ncatG);
   if (com.ngene>1) 
      fprintf (fout, " (%d genes: %s)  ", com.ngene, Mgenestr[com.Mgene]);
   if (com.alpha==0) com.nalpha=0;
   else             com.nalpha=(com.nalpha?com.ngene:!com.fix_alpha);
   if (com.nalpha>1 && com.rho) error ("Malpha or rho");
   if (com.nalpha>1) fprintf (fout,"(%d gamma)", com.nalpha);
   if (com.Mgene>=3) com.nrate*=com.ngene;
   if (com.Mgene && com.ngene==1) error("Mgene with ngene=1.");
   if (com.seqtype==CODONseq)
      fprintf (fout, "\nCodon frequencies: %s\n", codonfreqs[com.codonf]);
   if (com.aaDist && (com.seqtype==CODONseq||com.model==FromCodon))
      fprintf(fout,"%s, %s\n",
         com.daafile,(com.aaDist>0?"geometric":"linear"));

   if (com.NSsites) fprintf(fout, "%s, with %d dN/dS ratios for sites\n",
                   NSsitesmodels[com.NSsites], com.ncatG);

   com.maxnp = max((2*com.ns-2)*2+com.nrate+com.ngene-1+2, 50);
/*
printf ("\nmaxnp = %d\n", com.maxnp);
getchar ();
*/
   s3 = com.maxnp*(com.maxnp*2+2+12)*sizeof(double);
   s3 = max (s3, com.ls*(com.ns*sizeof(char)+sizeof(int)));
   i = com.ns*(com.ns-1)/2;
   s3 = max (s3, sizeof(double)*((com.ns*2-2)*(com.ns*2-2 + 4 + i) + i));
   if ((space=(double*)malloc(s3))==NULL)  error ("oom space");


      SeqDistance=(double*)malloc(com.ns*(com.ns-1)/2*sizeof(double));
      ancestor=(int*)malloc(com.ns*(com.ns-1)/2*sizeof(int));
      if (SeqDistance==NULL||ancestor==NULL) error("oom");

   if(com.seqtype==AAseq) { 
      Initialize (fout, space, com.seqtype);
      Getfcodon (com.pi, com.fb61);   /* get codon freqs from aa freqs */ 
   }
   else
      InitializeCodon (fout, space); /* seqtype change!*/

   if (com.seqtype==CODONseq) DistanceMatNG86 (fout, 0);
   else                       DistanceMatAA (fout, com.alpha);
   fflush(fout);

   if (com.Mgene==3)  FOR (i,com.ngene) xtoy(com.pi, com.piG[i], com.ncode);

   if (com.seqtype==AAseq && com.model==Poisson && !com.print) 
      PatternJC69like (NULL);

   if (com.alpha || com.NSsites)
     if((com.fhK=(double*)malloc(s2=com.npatt*com.ncatG*sizeof(double)))==NULL)
         error ("oom");

   if (com.runmode==-2)
      if (com.seqtype==CODONseq) { Pairwise(fout, space);  exit(0); }
      else error ("runmode");

   printf ("\n\n%10d bytes for distance ",com.ns*(com.ns-1)/2*sizeof(double));
   printf ("\n%10d bytes for partial likelihood\n",
           (com.ns-1)*com.ncode*com.npatt*sizeof(double));
   com.chunk = 
      (double*) malloc ((com.ns-1)*com.ncode*com.npatt*sizeof(double));
   if (com.chunk==NULL) error ("oom chunk");
   printf ("%10d bytes for fhK\n%10d bytes for working space\n", s2, s3);

   if (com.Mgene==1)         MultipleGenes (fout, space);
   else if (com.runmode==0)  Forestry (fout, space);
   else if (com.runmode==3)  StepwiseAddition (fout, space);
   else if (com.runmode>=4)  Perturbation(fout, (com.runmode==4), space);
   else                      StarDecomposition (fout, space);
   FPN(F0);
   if (noisy) putchar ('\a');
   return (0);
}


int Forestry (FILE *fout, double space[])
{
   int  i,j,k, itree, ntree, np;
   double x[NP],xcom[NP-NBRANCH],xb[NP][2], lnL,lnL0=0,lnLm=0, e=1e-7,tl=0;
   double *var=space+NP, nchange;
   FILE *ftree, *fin=(FILE*)fopen("in.codeml","r"), *frates=NULL;

   if (com.alpha)
      if ((frates=(FILE*)fopen("Rates.out","w"))==NULL) error("file err");
   if ((ftree=fopen (com.treef,"r"))==NULL)  error ("no tree file.");
   fscanf (ftree, "%d%d", &i, &ntree);
   if (i!=com.ns) error ("# seqs do not match."); 
   if (ntree>10 && com.print) puts("\nlarge lfh file");
   fprintf (flfh, "%6d%6d%6d\n", ntree, com.ls, com.npatt);
   FOR (itree, ntree) {
      printf ("\nTREE # %2d\n", itree+1);
      fprintf (fout,"\nTREE # %2d:  ", itree+1);
      fprintf (flfh,"\n\n%2d\n", itree+1);
      fprintf (frub,"\n\nTREE #%2d", itree+1);

      LASTROUND=0;
      if (ReadaTreeN(ftree, &i, 1)) error ("err tree..");

      nchange=MPScore (space);
      OutaTreeN (F0, 0, 0);     printf ("   MP score: %.0f", nchange);
      OutaTreeN (fout, 0, 0);   fprintf (fout, "   MP score: %.0f", nchange);
      fflush (fout),  fflush (flfh);  fflush(frst);

      GetInitials (x);
      if ((np=com.np)>NP || np-com.ntime>NP-NBRANCH) error("raise NP");
      printf ("\nntime & nrate & np: %d  %d  %d ", com.ntime,com.nrate,com.np);
      if (itree)  for (i=0; i<np-com.ntime; i++) x[com.ntime+i]=xcom[i];
      SetxInitials (x);

      if (com.seqtype==CODONseq && (com.model==NSb || com.model==NS2)) {
         printf ("\n%d dN/dS ratios for branches assumed:\n", N_OMEGA);
         FOR (i,tree.nbranch) printf("%4d", OmegaBranch[i]); FPN(F0);
         fprintf (fout, "\n%d dN/dS ratios for branches assumed:\n", N_OMEGA);
         FOR (i,tree.nbranch) fprintf(fout, "%4d", OmegaBranch[i]); FPN(fout);
      }
      if (fin) {
         puts("\nInitials from in.codeml. Break if not correct");
         getchar ();
         FOR (i,np) if (fscanf(fin,"%lf",&x[i])!=1) break;
         if (i<np)  {
            printf("err at #%d in in.codeml. Edit or remove it.",i+1);
            exit (-1);
         }
      }
/*
      printf("\nInitials (np = %d)?\n", np);
      FOR (i,np) scanf("%lf", &x[i]);
*/

      PointLklnodes ();
      NFunCall=0;

      if (noisy) {
         printf ("\nnp%6d", np);       
         if (noisy>2) matout (F0, x, 1, np); 
         lnL = com.plfun (x, np);
         printf ("\nlnL0 = %12.6f\n", -lnL);
      }

      if (com.clock==0 /* && com.model>REVaa_0  && com.NSsites<2 */ ) {
         SetxBound (np, xb);
         j=ming2(noisy>2?frub:NULL,&lnL,com.plfun,NULL,x,xb, space,e,np);
      }
      else 
         j=ming1(noisy>2?frub:NULL,&lnL,com.plfun,NULL,testx,x,space,e,np);

      printf ("\nOut..\nlnL  = %12.6f\n", -lnL);
      if (j) fprintf (fout,"\ncheck convergence..");
      if (itree==0) lnL0=lnLm=lnL;
      else if (lnL<lnLm) lnLm=lnL;
      if (itree==0) 
         for (i=0; i<np-com.ntime; i++) xcom[i]=x[com.ntime+i];
      else if (!j)  
         for (i=0; i<np-com.ntime; i++) xcom[i]=xcom[i]*.2+x[com.ntime+i]*0.8;

      if (com.seqtype==1 && com.NSsites==2) {
          if (com.getSE>0) LASTROUND=1;
          k=com.ntime+com.nrgene+com.nrate;
          tl=1+exp(x[k])+exp(x[k+1]);
          x[k]=exp(x[k])/tl;  x[k+1]=exp(x[k+1])/tl;
      }
      fprintf (fout,"\nlnL(ntime:%3d  np:%3d):%14.6f%+14.6f\n",
         com.ntime, np, -lnL, -lnL+lnL0);
      OutaTreeB (fout);  FPN (fout);
      FOR (i,np) fprintf (fout," %8.5f", x[i]);  FPN (fout);  fflush(fout);
/*
      fprintf(frst," %10.3f%7.3f", -lnL,x[com.ntime]);

for(i=com.ntime;i<com.np;i++) fprintf(frst," %8.4f", x[i]);  
fprintf(frst,"%12.3f\n", -lnL); 
*/

      if (com.getSE>0) {
         Hessian (np, x, lnL, space, var, com.plfun, var+np*np);
         matinv (var, np, np, var+np*np);
         fprintf (fout,"SEs for parameters:\n");
         FOR(i,np) fprintf(fout," %8.5f",var[i*np+i]>0.?sqrt(var[i*np+i]):0.);
         FPN (fout);
         if (com.getSE==2) matout2(fout, var, np, np, 12,7);
      }
      if (com.clock) {
         SetBranch (x, 0);
         for (i=0,tl=0;i<tree.nnode;i++) if(i-tree.origin) tl+=nodes[i].branch;
      }
      else   tl=sum(x, com.ntime);       
      fprintf (fout,"\ntree length: %9.5f\n", tl);
      OutaTreeN (fout, 1, 1);  FPN(fout);
      if (com.np>com.ntime) DetailOutput (fout, x, space);

      if (com.NSsites)            lfunNSsites_rate (frst, x, np);
      else if (com.print) {
         AncestralSeqs (frst, x, space);
         if (com.plfun!=lfun)     lfunAdG_rate (frates, x, np);
      }
      if (com.seqtype==AAseq && com.model>=2 /* AAClasses */)
         EigenQ20(fout, Root, U, V, x+com.ntime+com.nrgene); /* S & PAM */

      com.print-=9;
      com.plfun (x, np);
      com.print+=9;
      if (com.clock) {   /* collaps zero-length internal branch? */
         for (i=com.ns,j=0; i<tree.nnode; i++)
            if (i!=tree.origin && nodes[i].branch<5e-5)
               { CollapsNode (i, x);  j=1; break; }
         if (j) {
            fputs("\nAn internal branch may be collapsed to give\n", fout); 
            OutaTreeN (fout, 0, 0);  FPN(fout);
            OutaTreeN (fout, 1, 0);  FPN(fout);
            fputs("You may repeat the analysis using this tree.\n", fout); 
	 }
      }
   }         /* for (itree) */
   if (ntree>2) fprintf (fout, "\nBest likelihood:%12.6f\n", -lnLm);
   fclose (ftree);   if (fin) fclose(fin);

   if (com.aaDist && (com.seqtype==CODONseq||com.model==FromCodon))
      printf("\n%s, %s.\n",com.daafile,(com.aaDist>0?"geometric":"linear"));

   return (0);
}


static int ijAAref=19*20+9; 
/* reference aa pair: VI (for REVaa, REVaa_0, AAClasses to estimate Grantham)
   The rate for this pair is set to 1, and other rates are relative to it.
*/

void DetailOutput (FILE *fout, double x[], double space[])
{
   int i,j,k=com.ntime;
   double om;

   fprintf (fout, "\nDetailed output identifying parameters\n");
   if (com.nrgene) {
      fprintf (fout, "\nrates for %d genes:%6.0f", com.ngene, 1.);
      FOR (i,com.nrgene) fprintf (fout, " %8.5f", x[k++]);
      FPN(fout);
   }
   if (com.seqtype==CODONseq || com.model==FromCodon) {
      if (!com.fix_kappa) fprintf(fout,"kappa (ts/tv) = %8.5f\n", x[k++]);
      if (com.NSsites) {
         fprintf (fout,"\nFrequencies for categories (K=%d)", com.ncatG);
         matout (fout, com.freqK, 1, com.ncatG);
         fprintf (fout,"\ndN/dS rate ratios for categories (K=%d)", com.ncatG);
         matout (fout, com.rK, 1, com.ncatG);
      }
   }
   else  k+=com.nrate;
   if (!com.fix_alpha)  fprintf (fout, "alpha (gamma) = %8.5f\n", x[k++]);
   if (!com.fix_rho)   fprintf (fout, "rho (correlation) = %8.5f\n", x[k]);

   if (com.model==AAClasses) {
      fprintf (fout, "\nw (dN/dS) classes for amino acid pairs:\n");
      FOR (k,N_OMEGA) { 
         fprintf (fout, "%9.5f:", x[com.ntime+com.nrgene+!com.fix_kappa+k]);
         FOR (i,20) FOR(j,i)
            if (OmegaAA[i*(i-1)/2+j]==k) fprintf(fout," %c%c", AAs[i],AAs[j]);
         if (k==0)  fprintf(fout, " (background ratio)");
         FPN(fout); 
      }
      if (com.nrate>65) { /* estimating grantham matrix */
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
   if((com.seqtype==CODONseq||com.model==FromCodon||com.model==AAClasses) 
      && com.NSsites==0) {
      fprintf (fout, "\ndN & dS for each branch\n");
      XCOPY=x;
      FOR (i,tree.nbranch) {
         _BRANCHLEN=x[i];
         fprintf(fout," %3d..%-3d  ",
              tree.branches[i][0]+1,tree.branches[i][1]+1);
         /* fprintf (frst, "\n%9.4f", _BRANCHLEN=x[j]); */
         j=com.ntime+com.nrgene+!com.fix_kappa;
         if (com.model==AAClasses) om=-1;
         else if (com.model==NS1 ||com.model==FromCodon)
                                   om=(com.fix_omega?com.omega:x[j]);
         else if (com.model==NSb)  om=x[j+i];
         else if (!com.fix_omega || OmegaBranch[i]<N_OMEGA-1)
                              om=x[j+=OmegaBranch[i]];
         else                 om=OMEGA_FIX;
         EigenQ61(fout,Root,U,V,com.kappa,om,space); 
                         /* om not used in AAClasses model */
      }
   }
   if (com.alpha) {
      fprintf (fout, "rates for categories:");
      FOR (i,com.ncatG) fprintf (fout, " %8.5f", com.rK[i]);
      fprintf (fout, "\nfreqs for categories:");
      FOR (i,com.ncatG) fprintf (fout, " %8.5f", com.freqK[i]);
      FPN (fout);
   }
   if (com.rho) {
      fprintf (fout, "transition probabilities between rate categories:\n");
      for(i=0;i<com.ncatG;i++,FPN(fout))  FOR(j,com.ncatG) 
         fprintf (fout, " %8.5f", com.MK[i*com.ncatG+j]);
   }
   FPN (fout);
}


extern double SmallDiff;

int GetOptions (char *ctlf)
{
   int i, nopt=28, lline=255;
   char line[255], *pline, opt[99], comment='*';
   char *optstr[] = {"seqfile", "outfile", "treefile", "seqtype", "noisy", 
        "runmode", "clock", "getSE", "RateAncestor", "CodonFreq", "verbose",
        "model","aaDist","aaRatefile",
        "NSsites", "icode", "Mgene", "fix_kappa", "kappa",
        "fix_omega", "omega", "fix_alpha", "alpha","Malpha", "ncatG", 
        "fix_rho", "rho", "Small_Diff"};
   double t;
   FILE  *fctl=fopen (ctlf, "r");
   char *daafiles[]={"", "grantham.dat", "miyata.dat", 
                     "g1974c.dat","g1974p.dat","g1974v.dat"};

   if (fctl) {
      if (noisy) printf ("\n\nReading options from %s..\n", ctlf);
      for (;;) {
         if (fgets (line, lline, fctl) == NULL) break;
         for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
            if (isalnum(line[i]))  { t=1; break; }
            else if (line[i]==comment) break;
         if (t==0) continue;
         sscanf (line, "%s%*s%lf", opt,&t);
         if ((pline=strstr(line, "= "))==NULL) 
            error("err: option file. space around = ?");

         for (i=0; i<nopt; i++) {
            if (strncmp(opt, optstr[i], 8)==0)  {
               if (noisy>2)
                  printf ("\n%3d %15s | %-20s %6.2f", i+1,optstr[i],opt,t);
               switch (i) {
                  case ( 0): sscanf(pline+2, "%s", com.seqf);    break;
                  case ( 1): sscanf(pline+2, "%s", com.outf);    break;
                  case ( 2): sscanf(pline+2, "%s", com.treef);   break;
                  case ( 3): com.seqtype=(int)t;     break;
                  case ( 4): noisy=(int)t;           break;
                  case ( 5): com.runmode=(int)t;     break;
                  case ( 6): com.clock=(int)t;       break;
                  case ( 7): com.getSE=(int)t;       break;
                  case ( 8): com.print=(int)t;       break;
                  case ( 9): com.codonf=(int)t;      break;
                  case (10): com.verbose=(int)t;     break;
                  case (11): com.model=(int)t;       break;
                  case (12): com.aaDist=(int)t;      break;
                  case (13): sscanf(pline+2,"%s",com.daafile); break;
                  case (14): com.NSsites=(int)t;     break;
                  case (15): com.icode=(int)t;       break;
                  case (16): com.Mgene=(int)t;       break;
                  case (17): com.fix_kappa=(int)t;   break;
                  case (18): com.kappa=t;            break;
                  case (19): com.fix_omega=(int)t;   break;
                  case (20): com.omega=t;            break;
                  case (21): com.fix_alpha=(int)t;   break;
                  case (22): com.alpha=t;            break;
                  case (23): com.nalpha=(int)t;      break;
                  case (24): com.ncatG=(int)t;       break;
                  case (25): com.fix_rho=(int)t;     break;
                  case (26): com.rho=t;              break;
                  case (27): SmallDiff=t;            break;
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

   setmark_61_64 ();
   if (com.ngene==1) com.Mgene=0;
   if (com.model==AAClasses) {
         SetAA1STEP();  
         GetOmegaAA(OmegaAA);
   }  
   if (com.seqtype==AAseq || com.seqtype==CODON2AAseq) {
      com.ncode=20;
      switch (com.model) {
      case (Poisson):  case (EqualInput): case (Empirical): case (Empirical_F):
         com.fix_kappa=1; com.kappa=0; com.nrate=0;   break;
      case (FromCodon): 
         com.nrate=!com.fix_kappa+!com.fix_omega;
         if (com.kappa<0) error("kappa.."); 
         break;
      case (AAClasses): break;
      case (REVaa_0): com.fix_kappa=0; com.kappa=0; SetAA1STEP(); break; 
      case (REVaa):  com.fix_kappa=0; com.kappa=0; com.nrate=189; break;
      default: error ("model unavailable");
      }
      if ((com.Mgene>=2 && com.model==Fequal) 
        || (com.Mgene>=3 && com.model!=FromCodon)) error ("Mgene");
   }
   else if (com.seqtype==CODONseq) {
      com.ncode=Nsensecodon[com.icode];
      if (com.model!=AAClasses) {
         if (com.model>NS1 && com.NSsites) error ("Diff NS + NSsites.");
         if (com.model>NS1 || com.NSsites) {
            if(com.fix_kappa>1) error("fix_kappa, not implemented.");
            if(com.Mgene>1) error("model: Diff NS+Mgene");
         }
         if (com.fix_omega) {
            OMEGA_FIX=com.omega;
            if (com.model==NSb || com.NSsites) error("fix_omega?");
         }
         if (com.model>NS2) error ("seqtype or model.");
         if (com.kappa<0)  error("kappa..");
         if (com.Mgene>=3 && com.nrate==0)  error("Mgene");

         com.nrate=!com.fix_kappa+!com.fix_omega;
         if (com.aaDist)  com.nrate++;
         if (com.NSsites) {
            com.nrate=!com.fix_kappa;
            if(com.runmode==-2) error ("NSsites models..");
            if(com.NSsites==3)
                 { com.nrate++; if(com.alpha<=0) error("alpha<=0");}
            else               com.ncatG=com.NSsites+1;
         }
      }
   }
   else  
      error ("seqtype..");
   if (com.aaDist && (com.seqtype==CODONseq || com.model==FromCodon))
      strcpy(com.daafile, daafiles[abs(com.aaDist)]);

   else if (com.seqtype==AAseq && com.model>=REVaa_0) 
      strcpy(com.daafile,"mtmam.dat");

   if (com.fix_alpha && com.alpha==0) {
      if (com.rho) puts("rho set to 0.");  com.fix_rho=1; com.rho=0; 
   }
   if(!com.fix_alpha && com.alpha<=0) { com.alpha=0.5; puts("alpha reset"); }
   if(!com.fix_rho && com.rho==0) { com.rho=0.001;  puts("init rho reset"); }
   if(com.alpha)  if (com.ncatG<2 || com.ncatG>NCATG) error ("ncatG");

   if (com.runmode==3 && (com.clock))  error("runmode+clock");

   return (0);
}


int testx (double x[], int np)
{
   int i,k;
   double tb[]={.4e-5, 9}, rgeneb[]={0.1,20}, rateb[]={1e-5,999}, omegab=.005;
   double alphab[]={0.005,99}, rhob[]={0.01,0.99};

   if (com.clock && SetBranch (x, 0))  return (-1);
   FOR (i,com.ntime)  if (x[i]<tb[0] || x[i]>tb[1])  return (-1);
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
   FOR (i,com.NSsites) if (x[k+i]<0 || x[k+i]>1) return (-1);
   if (com.NSsites==2) { 
      if (1-x[k]-x[k+1]<0) return (-1);
      if (x[k+2]<rateb[0] || x[k+2]>rateb[1]) return (-1);
   }
   for (i=0; i<com.nalpha; i++,k++)
      if (x[k]<alphab[0] || x[k]>alphab[1]) return (-1);
   if (!com.fix_rho) if (x[np-1]<rhob[0] || x[np-1]>rhob[1]) return (-1);
   return(0);
}

int SetxBound (int np, double xb[][2])
{
   int i, j, k;
   double tb[]={.0001,9}, rgeneb[]={0.1,99}, rateb[]={1e-4,999};
   double alphab[]={0.005,99}, rhob[]={0.01,0.99}, omegab[]={.005,99}; 

   FOR (i,com.np)  FOR (j,2) xb[i][j]=tb[j];
   FOR (i,com.nrgene) FOR (j,2) xb[com.ntime+i][j]=rgeneb[j]; 
   FOR (i,com.nrate)  FOR (j,2) xb[com.ntime+com.nrgene+i][j]=rateb[j];
   k=com.ntime+com.nrgene+!com.fix_kappa;
   if ((com.seqtype==CODONseq||com.model==FromCodon) && com.model!=AAClasses) 
      { xb[k][0]=omegab[0]; xb[k][1]=omegab[1]; }

   if (com.aaDist<0 && (com.seqtype==1||com.model==FromCodon)) {
      /* linear relationship between d_ij and w_ij */
      if(com.nrate != !com.fix_kappa+1+(com.seqtype==1)) error("in Setxbound");
      xb[com.ntime+com.nrgene+!com.fix_kappa][1]=1; /* 0<b<1 */
   }

   k=com.ntime+com.nrgene+com.nrate;
   if (com.NSsites==1) { xb[k][0]=1e-4; xb[k][1]=.9999; } /*prob for class 1*/
   else if (com.NSsites==2) { xb[k][0]=xb[k+1][0]=-99; xb[k][1]=xb[k+1][1]=99;}

   k=com.ntime+com.nrgene+com.nrate;
   for (i=0;i<com.nalpha;i++,k++)  FOR (j,2) xb[k][j]=alphab[j];
   if (!com.fix_rho)   FOR (j,2) xb[np-1][j]=rhob[j];

   return(0);
}


int SetxInitials (double x[])
{
/* This forces initial values into the boundary of the space
*/
   int i, j, k;
   double tb[]={.0002,3}, rgeneb[]={0.1,9}, rateb[]={.01,29};
   double alphab[]={0.05,5}, rhob[]={0.01,0.9};

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

int GetInitials (double x[])
{
/* 
   perhaps try to restruct the code and make two sections for amino acids 
   and codons?
*/
   static int times=0;
   int i, j,k=0, nc=com.ncode, naa=20;
   double t;

   com.plfun = (com.alpha==0 ? lfun : (com.rho==0?lfundG:lfunAdG));
   if (com.NSsites)  com.plfun=lfunNSsites; 

   if (times++==0) {
      if ((com.aaDist && (com.seqtype==CODONseq||com.model==FromCodon)) ||
          (com.seqtype==AAseq &&
          (com.model==Empirical||com.model==Empirical_F||com.model>=REVaa_0))){
         GetDaa(NULL, com.daa);
      }
   }
   if (com.seqtype==CODONseq && (com.model==NSb||com.model==NS2)) {
      if (com.model==NS2) {
         puts ("\nSpecify the dN/dS ratio for each branch (use 0,1,2,etc.).");
         OutaTreeN(F0,0,0);  FPN(F0);  OutaTreeB(F0);  FPN(F0);
         for (i=0,N_OMEGA=0; i<tree.nbranch; i++) {
            printf("Branch %2d: %2d..%-2d? ", i+1,
               tree.branches[i][0]+1, tree.branches[i][1]+1);
            scanf ("%d", &OmegaBranch[i]);
            if (OmegaBranch[i]+1>N_OMEGA) N_OMEGA=OmegaBranch[i]+1;
         }
         if (N_OMEGA==1) error("Only one ratio, use model=0.");
         printf ("\n\n%d ratios are specified. Stop if incorrect.", N_OMEGA);
      }
      else {
         N_OMEGA=tree.nbranch;  FOR(i,tree.nbranch) OmegaBranch[i]=i;
      }
      com.nrate=!com.fix_kappa+!com.fix_omega+N_OMEGA-1;
   }
   com.nrgene=(!com.fix_rgene)*(com.ngene-1);
   if (!com.fix_kappa && com.Mgene>=3) com.nrate*=com.ngene;

   com.np = com.ntime+com.nrgene+com.nrate+(!com.fix_alpha)+(!com.fix_rho);
   if (com.NSsites==1) com.np++;
   else if (com.NSsites==2) com.np+=3;
   FOR (j, com.ntime)  x[j]=0.05+0.01*(com.ntime-j);

   LSDistance (&t, x, testx);
/*
   FOR (j, com.ntime) x[j]=0.1;
*/
   FOR (j, com.nrgene)  x[com.ntime+j]=1;

   k=com.ntime+com.nrgene;
   if (com.model==AAClasses) { 
      if (!com.fix_kappa) x[k++]=com.kappa;
      FOR (i,com.nrate-!com.fix_kappa) x[k++]=com.omega;
      if (com.nrate>65) puts("\a\nget better initial values?");
       /* if estimating all acceptance rates */
   }
   else {
      if (com.seqtype==AAseq) {                     /* AAseq */
         if (com.nrate==0)  EigenQ20(NULL, Root, U, V, &t); /* once for all */
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
            if (!com.fix_kappa)  x[k++]=com.kappa;
            if (j,com.nrate-!com.fix_kappa)  x[k++]=com.omega; 
         }
      }
      else {                                        /* CODONseq */
         if (com.nrate==0 && com.NSsites==0) 
            EigenQ61(NULL, Root, U, V, com.kappa, com.omega, PMat);
         else if (com.model==NSb || com.model==NS2) {
            if (!com.fix_kappa)  x[com.ntime+com.nrgene]=com.kappa; 
            FOR(i, (com.model==NSb?tree.nbranch:N_OMEGA-1+!com.fix_omega))
               x[com.ntime+com.nrgene+!com.fix_kappa+i]=com.omega;
         }
         else {
            FOR (i, com.nrate*(com.Mgene>=3?com.ngene:1)) {
               if (com.fix_kappa)      x[k+i*com.nrate]=com.omega; 
               else if (com.fix_omega) x[k+i*com.nrate]=com.kappa; 
               else {
                  x[k+i*com.nrate]=com.kappa; 
                  FOR(j,com.nrate-com.fix_kappa) x[++k+i*com.nrate]=com.omega;
               }
            }
            /* to be corrected later */
            k=com.ntime+com.nrgene+!com.fix_kappa;
            if (com.NSsites==1) x[k]=.5;
            else if (com.NSsites==2) { x[k]=.1; x[k+1]=.7; x[k+2]=com.omega; }
         }
      }
   }
   if (!com.fix_alpha) x[com.ntime+com.nrgene+com.nrate]=com.alpha;
   if (!com.fix_rho) x[com.np-1]=com.rho;
   if (com.rho)
      AutodGamma (com.MK, com.freqK, com.rK, &t, com.alpha, com.rho,com.ncatG);
   else if (com.alpha) {
      DiscreteGamma (com.freqK, com.rK, com.alpha, com.alpha, com.ncatG, 0);
      if (com.NSsites==3)  FOR (j,com.ncatG) com.rK[j]*=com.omega;
   }

   return (0);
}



int SetPGene (int igene, int _pi, int _UVRoot, int _alpha, double x[])
{
   int k, nr1=com.nrate/com.ngene;

   if (_pi) xtoy (com.piG[igene], com.pi, com.ncode);
   if (_UVRoot) {
      k=com.ntime+com.nrgene+(com.Mgene>=3)*igene*nr1;
      if (!com.fix_kappa)  com.kappa=x[k++];
      if (!com.fix_omega) com.omega=x[k++];
      if (com.seqtype==CODONseq)
         if (!com.NSsites) EigenQ61(NULL,Root,U,V,com.kappa,com.omega,PMat);
      else 
         EigenQ20(NULL, Root, U, V, x+k);
   }
   if (_alpha) {
      com.alpha=x[com.ntime+com.nrgene+com.nrate+igene];
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
*/
   int i, k, naa=20;
   double t;

   XCOPY=x;

   k=SetBranch (x, 1);
#if (DEBUG)
   if(k) puts ("\nbranch len err..");
#endif
   FOR (i, com.nrgene) com.rgene[i+1]=x[com.ntime+i];

   k=com.ntime+com.nrgene;
   if (com.nrate) {
      if (!com.fix_kappa) com.kappa=x[k++]; 
      if (com.model!=AAClasses && !com.fix_omega) com.omega=x[k++]; 
      if (com.seqtype==AAseq)
         EigenQ20(NULL, Root, U, V, x+com.ntime+com.nrgene);
      else if ((com.model==0 && com.NSsites==0) || com.model==AAClasses)
         /* CODONs, same S/N for branch site */ 
         EigenQ61(NULL, Root, U, V, com.kappa, com.omega, PMat);
   }
   k=com.ntime+com.nrgene+com.nrate;
   if (com.NSsites==1) {        /* only one parameter, p of neutral sites */
      com.freqK[0]=x[k];    com.freqK[1]=1-x[k];
      com.rK[0]=1;          com.rK[1]=0;
   }
   else if (com.NSsites==2) {
      if (LASTROUND)
         for (i=0,com.freqK[2]=1; i<2;i++) com.freqK[2]-=(com.freqK[i]=x[k+i]);
      else {
         t=1+exp(x[k])+exp(x[k+1]);
         com.freqK[0]=exp(x[k])/t;com.freqK[1]=exp(x[k+1])/t; com.freqK[2]=1/t;
      }
      com.rK[0]=1;  com.rK[1]=0;  com.rK[2]=x[k+2];
   }
   if (!com.fix_alpha) {
      com.alpha=x[k++];
      if (com.fix_rho) {
         DiscreteGamma (com.freqK, com.rK, com.alpha, com.alpha, com.ncatG, 0);
       if (com.NSsites==3)  FOR (i,com.ncatG) com.rK[i]*=com.omega;
      }
   }
   if (!com.fix_rho) {
      com.rho=x[k++];
      AutodGamma (com.MK, com.freqK, com.rK, &t, com.alpha, com.rho, com.ncatG);
   }

   return (0);
}

int EigenQ20 (FILE *fout, double Root[], double U[], double V[], double rate[])
{
   int naa=20, i,j,k;
   double Q[20*20], mr, t=0, space[20*3];

   FOR (i,naa*naa) Q[i]=0;
   switch (com.model) {
   case (Poisson)   : case (EqualInput) : 
      fillxc (Q, 1., naa*naa);  break;
   case (Empirical)   : case (Empirical_F):
      FOR(i,naa) FOR(j,i) Q[i*naa+j]=Q[j*naa+i]=com.daa[i*naa+j]/100;
      break;
   case (FromCodon): case (AAClasses): 
      EigenQ61(NULL, Root, U, V, com.kappa, com.omega, PMat);
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
   if (fout) {
      fprintf (fout, "\n\nRate matrix (symmetrical part, Sij)\n");
      for(i=0,t=0;i<naa;i++) FOR(j,i) t+=Q[i*naa+j]/com.pi[j]/(naa*(naa-1)/2.);
      FOR (i,naa) {
         fprintf (fout, "\n%-5s", getaa(0, i+1));
         FOR (j,i) fprintf (fout, " %4.0f", Q[i*naa+j]/t/com.pi[j]*100);
      }
      fprintf(fout, "\n     ");  FOR(i,naa) fprintf(fout,"%5s",getaa(0,i+1));
      fprintf(fout, "\n     ");  FOR(i,naa) fprintf(fout,"%5s",getaa(1,i+1));
      FPN (fout);
      fflush(fout);
   }
/*
   if (fout && frst) {
      fprintf(frst, "\nRate matrix (symmetrical part, Sij) for bubble plot\n");
      FOR (i,naa)  FOR (j,i) 
         fprintf(frst, "\t%d\t%d\t%.2f\n", i+1,j+1,Q[i*naa+j]/t/com.pi[j]*100);
   }
*/
   if (eigen (1, Q, naa, Root, space, U, V, space+com.ncode))
      error ("eigenQ20 err");
   xtoy (U, V, naa*naa);
   matinv (V, naa, naa, space);

   FOR (i,naa)  Root[i]=Root[i]/mr;

#if (DEBUG)
      if (fout) {
         fprintf (fout, "\nPI & Root\n");
         matout (fout, com.pi, 1, 20);
         matout (fout, Root, 1, 20);
         fprintf (fout, "\nPAM matrix, P(0.01)");
         PMatUVRoot (PMat, 0.01, naa, U, V, Root);

         FOR (i,naa) {
            fprintf (fout, "\n%5s", getaa(0, i+1));
            FOR (j,naa)  fprintf (fout, "%6.0f", PMat[i*naa+j]*100000);
         }
         FPN (fout);
      }
#endif
   return (0);
}


int EigenQ61 (FILE* fout, double Root[], double U[], double V[], 
    double kappa, double omega, double Q[64*64])
{
   int n=Nsensecodon[com.icode], i,j,k, ic1,ic2,ndiff,pos=0,from[3],to[3];
   double rs0,ra0,rs,ra, S,N,dS,dN,dN_dS,c0[3],c[3],t,space[64*3], *ri=space;
   double *pi=(com.seqtype==AAseq?com.fb61:com.pi);
   double *pomega=XCOPY+com.ntime+com.nrgene+!com.fix_kappa, w;

   FOR (i,3) c[i]=c0[i]=0;  
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
      if((ic1=GenetCode[com.icode][ic1]-1)!=(ic2=GenetCode[com.icode][ic2]-1)){
         ra0+=t*Q[i*n+j];
         if (com.model==AAClasses) {
            if (ic1<ic2)  { k=ic2; ic2=ic1; ic1=k; }
            k=ic1*(ic1-1)/2+ic2;
            if (pomega[OmegaAA[k]]<0) {
               if (noisy)  printf ("ic1 & ic2 & iw & w: %d %d %d %.5f\n", 
                              ic1,ic2,OmegaAA[k], pomega[OmegaAA[k]]);
               pomega[OmegaAA[k]]=0;
            }

            if (com.seqtype==AAseq && com.nrate>65 && ic1*20+ic2==ijAAref)
                ;     /* if estimating grantham's matrix with aa sequences */
            else  Q[i*n+j] *= pomega[OmegaAA[k]];
         }
         else if (com.aaDist==0)  Q[i*n+j] *= omega;
         else {
            w = pomega[0]*com.daa[ic1*20+ic2];
            if(com.aaDist>0)           Q[i*n+j] *= exp(-w);  /* geometric */
            else                       Q[i*n+j] *= 1-w;      /* linear */
            if (com.seqtype==CODONseq) Q[i*n+j]*=pomega[1];
/*
            Q[i*n+j]*=pomega[1]*(1-fabs(com.daa[ic1*20+ic2]-pomega[0]));
*/
         }
         ra+=t*Q[i*n+j];
      }
      else 
         rs+=t*Q[i*n+j];
      c[pos]+=t*Q[i*n+j];
      Q[j*n+i]=Q[i*n+j];
   }
   if (com.seqtype==AAseq && !fout) return (0);
   if (com.NSsites==0)  FOR(i,n) FOR(j,n) Q[i*n+j]*=pi[j]/(t=rs+ra);
   else                 FOR(i,n) FOR(j,n) Q[i*n+j]*=pi[j];

   FOR (i,n) Q[i*n+i]=-sum(Q+i*n,n);


   if (fout) {
      rs0=rs;
      t=rs0+ra0; rs0/=t;  ra0/=t;   S=rs0*3*com.ls; N=ra0*3*com.ls;
      t=rs+ra;   rs/=t;   ra/=t;
      dN_dS= (ra/rs)/(ra0/rs0);
      if (_BRANCHLEN)  {
         dS=_BRANCHLEN*rs/3/rs0;  dN=_BRANCHLEN*ra/3/ra0;
         fprintf (fout,
             "t=%7.3f  S=%7.3f  N=%7.3f  dN/dS=%7.3f  dN=%7.4f  dS=%7.4f\n",
              _BRANCHLEN,S,N,dN_dS,dN,dS);
/*
         fprintf (frst, "%8.1f%8.1f %9.4f%9.4f%9.3f", N, S, dN, dS, dN_dS);
*/
      }
      else 
         fprintf (fout, " N=%7.3f  S=%7.3f  dN/dS=%7.3f\n", N,S, dN_dS);
/*
      fprintf(fout,"\nc1:2:3 [0] = 1 :%6.3f :%6.3f  [*] = 1 :%6.3f :%6.3f\n",
           c0[1]/c0[0], c0[2]/c0[0], c[1]/c[0], c[2]/c[0]);
*/
   }
   if (com.seqtype==AAseq) return (0);

   if (eigen (1,Q,n,Root,ri,U,V,space+n)) error("eigenQ61 err.");
   xtoy (U, V, n*n);
   matinv (V, n, n, space);

   return (0);
}


int Qcodon2aa (double Qc[], double pic[], double Qaa[], double piaa[])
{
/* Q61 -> Q20

   This routine constructs the rate matrix for amino acid replacement from
   the rate matrix for codon substitution, by congregating states in the
   Markov chain.  Both processes are time reversible, and only the
   symmetrical part of the rate matrix are constructed.  Codon frequencies 
   pic[] are used.  They should be filled with 1/61 if only amino acid
   sequences are available. 
   Q20(aai,aaj) = SUMi SUMj (piC[i]*piC[j]]*Q61[i][j]) / (piAA[i]*piAA[j])
*/
   int i, j, aai, aaj, nc=Nsensecodon[com.icode], naa=20;
   double ti, tij;

   zero(piaa,naa);  zero(Qaa,naa*naa);
   FOR (i,nc) piaa[GenetCode[com.icode][FROM61[i]]-1] += pic[i];
   FOR (i,nc) {
      aai=GenetCode[com.icode][FROM61[i]]-1;
      ti=pic[i]/piaa[aai];
      FOR (j, i) {
         aaj=GenetCode[com.icode][FROM61[j]]-1;
         if (Qc[i*nc+j]==0 || aai==aaj) continue;
         tij=ti*pic[j]*Qc[i*nc+j]/piaa[aaj];
         Qaa[aai*naa+aaj]+=tij;    
         Qaa[aaj*naa+aai]+=tij; 
      }
   }

/*
for (i=0;i<naa;i++,FPN(frst)) FOR(j,i) {
   fprintf(frst,"%2d",(Qaa[i*naa+j]>0));
   if (Qaa[i*naa+j]<0) printf ("error Qij: %d %d \n", i,j);
}
fflush(frst);
*/
   return (0);
}

int PartialLikelihood (int inode, int igene)
{
   int n=com.ncode, i,j,k,h, ison, pos0=com.posG[igene],pos1=com.posG[igene+1];
   double t,om;

   FOR (i,nodes[inode].nson)
      if (nodes[nodes[inode].sons[i]].nson>0)
         PartialLikelihood (nodes[inode].sons[i], igene);
   fillxc (nodes[inode].lkl+pos0*n, (double)(inode>=com.ns), (pos1-pos0)*n);
   if (inode<com.ns) 
      for (h=pos0; h<pos1; h++) nodes[inode].lkl[h*n+com.z[inode][h]-1]=1;
   FOR (i, nodes[inode].nson) {
      ison=nodes[inode].sons[i];
      if (com.seqtype==CODONseq && (com.model==NSb||com.model==NS2)) {
         j=com.ntime+com.nrgene+!com.fix_kappa;
         om=XCOPY[j+OmegaBranch[k=nodes[ison].ibranch]];
         if (com.fix_omega && OmegaBranch[k]==N_OMEGA-1)  om=OMEGA_FIX;
         EigenQ61(NULL, Root, U, V, com.kappa, om, PMat);
      }
      t=nodes[ison].branch*com.rgene[igene];
      if (com.seqtype==AAseq && com.model==Poisson)
         PMatJC69like (PMat, t, n);
      else 
         PMatUVRoot (PMat, t, n, U, V, Root);
      if (nodes[ison].nson<1)               /* end node */
         for (h=pos0; h<pos1; h++) 
            FOR (j,n) nodes[inode].lkl[h*n+j]*=PMat[j*n+com.z[ison][h]-1];
      else
         for (h=pos0; h<pos1; h++) 
            FOR (j,n) {
               for (k=0,t=0; k<n; k++) t+=PMat[j*n+k]*nodes[ison].lkl[h*n+k];
               nodes[inode].lkl[h*n+j]*=t;
            }
   }        /*  for (ison)  */
   return (0);
}

int PMatJC69like (double P[], double t, int n)
{
   int i;
   double pii=1./n+(1.-1./n)*exp(-n/(n-1.)*t), pij=(1.-pii)/(n-1.);
   FOR (i, n*n) P[i]=pij;
   FOR (i, n) P[i*n+i]=pii;
   return (0);
}

int setmark_61_64 (void)
{
/* This sets FROM61[], which goes from 0, 1, ..., 61, and FROM64[], 
   which goes from 0, 1, ..., to 64.
*/
   int i, n;

   for (i=0,n=0; i<64; i++) {
      if (GenetCode[com.icode][i]==0) 
         FROM64[i]=-1; 
      else 
         { FROM61[n]=i; FROM64[i]=n++; }
   }
   if (n!=Nsensecodon[com.icode])  error ("# sense codons?");
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
   printcu (fout, space, com.icode, space+64);
   return (0);
}


int printsmaCodon (FILE *fout,char * z[],int ns,int ls,int simple)
{
/*  print, in blocks, multiple aligned and transformed codon sequences.
    indels removed.
    This is needed as codons are coded 1, 2, ..., 61, and printsma won't work.
 */
   int ig, ngroup, lt, il,is, i,b, lline=15;
   char equal='.',*pz, c0[4],c[4];

   ngroup = (ls-1)/lline + 1;
   for (ig=0,FPN(fout); ig<ngroup; ig++)  {
      fprintf (fout,"%-8d\n", ig*lline+1);
      FOR (is,ns)     {
	 lt=0; 
	 for(il=ig*lline,pz=z[is]+il; lt<lline && il<ls; il++,lt++,pz++) {
   	    b=(int) *pz;
            b=FROM61[b-1]; c[0]=b/16; c[1]=(b%16)/4; c[2]=b%4; c[3]=0;
            FOR(i,3) c[i]=NUCs[c[i]];
            if (is && simple)  {
     	         b=(int)z[0][il];
               b=FROM61[b-1]; c0[0]=b/16; c0[1]=(b%16)/4; c0[2]=b%4; c0[3]=0;
               FOR(i,3) if (c[i]==NUCs[c0[i]]) c[i]=equal;
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
/* for genes, site patterns and fpatts 
*/
   int h, i,j, it, ic[3], nc=com.ncode, nconstp, ig, lt, wname=20;
   double t, fb3x4[3*4], *fb64s, fb4t[4], fb3x4t[12], fb[64], *fbg;
   double *fpi[NGENE+1], lmax;

   if (com.ls%3) error ("Codon seq?");
   FOR (j, com.ns) 
      if (transform (com.z[j], com.ls, 1, 0)) printf("in seq #%d\n", j+1);

   fprintf (fout, "\n\nns =%4d\tls =%4d", com.ns, com.ls);
   FOR (j,com.ns) fprintf (fout,"\n%s", com.spname[j]);
/*
   FOR (j,com.ns) prints (fout, com.z[j], com.ls, 60, 3, 1);
   printsma (fout, com.z, com.ns, com.ls, 60, 3, 1, 1, CODONseq);
*/
   fbg=(fb64s=(double*)malloc((max(com.ns,com.ngene)+1)*64*sizeof(double)));
   if(fb64s==NULL) error("oom InitCodon");
   zero (fb3x4t, 3*4);   zero (fb64s, (com.ns+1)*64);   zero (com.fb61, 64);
   fprintf (fout,"\n\nCodon position x base (3x4) table for each sequence.");
   for (j=0; j<com.ns; j++) {
      for (h=0,zero(fb3x4, 3*4); h<com.ls; h+=3) {
         FOR (i,3) {
            ic[i]=(char)(com.z[j][h+i]-1);
            fb3x4[i*4+ic[i]] += 3./(double)com.ls;
            fb3x4t[i*4+ic[i]] += 3./((double)com.ns*com.ls);
         }
         it=ic[0]*16+ic[1]*4+ic[2];
         if (GenetCode[com.icode][it]==0) {
            printf ("\n%s: %s at codon #%2d", com.spname[j],getcode(it),h/3+1);
            error ("stop codon inside.");
         }
         fb64s[j*64+it] ++;
         com.z[j][h/3]=(char)(FROM64[it]+1); /* recode codons into 1,2,...61 */
      }
      fprintf (fout,"\n\n%-*s", wname, com.spname[j]);
      FOR (i,3) {
         fprintf (fout, "\nposition %2d:", i+1);
         FOR (h,4) fprintf (fout,"%5c:%7.5f", NUCs[h], fb3x4[i*4+h]);
      }
   }
   fprintf (fout, "\n\nAverage");
   FOR (i,3) {
      fprintf (fout, "\nposition %2d:", i+1);
      FOR (h,4) fprintf (fout,"%5c:%7.5f", NUCs[h], fb3x4t[i*4+h]);
   }
   FOR (i,4) fb4t[i]=(fb3x4t[i]+fb3x4t[4+i]+fb3x4t[4*2+i])/3;
   com.ls/=3;
   FOR (j,com.ns) FOR (i, 64) fb64s[com.ns*64+i]+=fb64s[j*64+i];
   FOR (i, 64) com.fb61[FROM64[i]]=fb64s[com.ns*64+i]/(com.ls*com.ns);
   fprintf (fout, "\n\nCodon usage for each species and their sums\n");
   printcums (fout, com.ns+1, fb64s, com.icode);
/*
printcu (fout, fb64s+64*com.ns, com.icode, space);
*/
/*
   fprintf (fout, "\nAverage sense codon frequencies, see the table below.");
   matout (fout, com.fb61, 16, 4);
   printfcode (fout, com.fb61, space);
*/
   if (com.seqtype==CODON2AAseq) {
      FOR (j,com.ns)  FOR (h,com.ls) 
         com.z[j][h]=GenetCode[com.icode][FROM61[com.z[j][h]-1]];
      com.seqtype=AAseq;
      fprintf (fout, "\nTranslated amino acid sequences..\n");
      fprintf (fout,"%10d %10d\n", com.ns, com.ls);
      FOR (j,com.ns)  {
         fprintf (fout,"%-*s ", LSPNAME, com.spname[j]);
         FOR (h,com.ls) fprintf(fout,"%1c", AAs[com.z[j][h]-1]);
         FPN (fout);
      }
/*
      printsma (fout, com.z, com.ns, com.ls, 60, 10, 1, 0, AAseq);
*/
   }
   PatternWeight (fout, space);

   if (fout && com.verbose) {
      fprintf(fout,"\nCodon site pattern freqs (patterns shown below)\n");
      FOR (h,com.npatt) {
         fprintf (fout," %3.0f", com.fpatt[h]);
         if ((h+1)%15==0) FPN (fout);
     }
     printsmaCodon (fout, com.z, com.ns, com.npatt, 1);
   }
   t=1./((double)com.lgene[0]*com.ns);  
   zero (com.pi, nc);  zero (fb, nc);  zero (fbg, (com.ngene+1)*64);
   for (h=0,ig=0,lt=0,nconstp=0; h<com.npatt; ) {
      for (j=1; j<com.ns; j++)  if(com.z[j][h]!=com.z[0][h]) break;
      if (j==com.ns) nconstp += (int)com.fpatt[h];
      FOR (j,com.ns) fb[com.z[j][h]-1]+=(double)com.fpatt[h]*t;
      if (com.seqtype==CODONseq) 
         FOR (j,com.ns) fbg[ig*64+FROM61[com.z[j][h]-1]]+=com.fpatt[h];
      else
         FOR (j,com.ns) fbg[ig*nc+com.z[j][h]-1]+=com.fpatt[h];

      lt+=(int)com.fpatt[h++];
      if (lt==com.lgene[ig]) {
         FOR (j,nc) com.pi[j]+=fb[j]/(t*com.ls*com.ns);
         xtoy (fb, com.piG[ig], nc);
         ig++;
         zero (fb,nc);
         if (ig<com.ngene) t=1./((com.lgene[ig]-com.lgene[ig-1])*com.ns);
      }
   }
   if (com.seqtype==CODONseq && com.ngene>1) {
      fprintf (fout, "\n\nSummed codon usage for each gene\n");
      printcums (fout, com.ngene, fbg, com.icode);
   }

   if (com.seqtype==AAseq) {      /* CODON2AAseq */
      printsma (fout, com.z, com.ns, com.npatt, 60, 10, 1, 0, AAseq);
      fprintf (fout, "\n\n%-12s", "AA freq for genes\n");
      FOR (ig, com.ngene) { 
         FOR (j,20) fprintf (fout,"%3c:%7.0f", AAs[j], fbg[ig*nc+j]);
         FPN (fout);
      }
      fprintf (fout, "\n\n%-12s", "Average aa freq.\n");
      FOR (j,20) fprintf (fout,"%3c:%7.5f", AAs[j], com.pi[j]);
   }
/* edit com.pi & com.piG for the model??? */
   if ((com.seqtype==AAseq && com.model==Fequal) 
       || (com.seqtype==CODONseq && com.codonf==Fequal)) {
      fillxc (com.pi, 1./nc, nc);
      FOR (j,com.ngene) fillxc (com.piG[j], 1./nc, nc);
   }
   if (com.seqtype==CODONseq && (com.codonf==F1x4 || com.codonf==F3x4)) {
      if (com.ngene>1) 
         fprintf (fout,"\nCodon position x base (3x4) table for each gene.\n");

      fpi[0]=com.pi;
      FOR (it, com.ngene) fpi[it+1]=com.piG[it];
      FOR (it, com.ngene+1) {
         zero (fb4t, 4); zero(fb3x4t, 12);         
         FOR (i, nc) {
            j=FROM61[i]; ic[0]=j/16; ic[1]=(j/4)%4; ic[2]=j%4; 
            FOR (j,3) fb4t[ic[j]]+=fpi[it][i];
            FOR (j,3) fb3x4t[j*4+ic[j]]+=fpi[it][i];
         }
         if (com.codonf==F1x4) 
            FOR (i,nc) {
               j=FROM61[i];
               fpi[it][i]=fb4t[j/16]*fb4t[(j/4)%4]*fb4t[j%4];
            }
         else if (com.codonf==F3x4) 
            FOR (i,nc) {
               j=FROM61[i];
               fpi[it][i]=fb3x4t[j/16]*fb3x4t[4+(j/4)%4]*fb3x4t[2*4+j%4];
            }
         abyx (1/sum(fpi[it],nc), fpi[it], nc);

         if (com.ngene>1) {
            if (it==0)  fprintf (fout,"\n\nAverage");
            else        fprintf (fout,"\n\nGene # %d", it);
            FOR (i,3) {
              fprintf (fout, "\nposition %2d:", i+1);
              FOR (j,4) fprintf (fout,"%5c:%7.5f", NUCs[j], fb3x4t[i*4+j]);
            }
         }
      }
   }
   fprintf (fout,"\n\n# constant sites: %6d (%6.2f%%)",
      nconstp, (double)nconstp*100./com.ls);
   for (h=0,lmax=-(double)com.ls*log((double)com.ls); h<com.npatt; h++)
      if (com.fpatt[h]>1)
         lmax+=com.fpatt[h]*log((double)com.fpatt[h]);
   fprintf (fout,"\nmax{ln(L)}:%16.6f\n", lmax);
   free(fb64s);
   return(0);
}

int Getfcodon (double faa[20], double fcodon[])
{
/* get codon freqs from amino acid freqs, assuming equal freq. for each syn
   codon.  Used in codon-based amino acid substitution models.
*/
   int ic, iaa, i, NCsyn[20];

   FOR(i,20) NCsyn[i]=0;
   FOR(ic,64) if((iaa=GenetCode[com.icode][ic])>0) NCsyn[iaa-1]++;
   zero (fcodon, 64);
   for (ic=0; ic<Nsensecodon[com.icode]; ic++) {
      iaa=GenetCode[com.icode][FROM61[ic]]-1;
      fcodon[ic]+=faa[iaa]/NCsyn[iaa];
   }
   if(fabs(1-sum(fcodon,64))>1e-6) printf("\n1 == %12.7f\n", sum(fcodon,64));

   return (0);
}


int OutSeqsCodon (FILE* fseq, int *pose)
{
/* Print codon sequence data into fseq.  Use pose=NULL if called before site 
   patterns are collapsed.  Sequences have been coded into 1 to 61.
*/ 
   char *pc;
   int h, j,ic;

   fprintf (fseq, "%8d%8d\n", com.ns, com.ls*3);
   for (j=0; j<com.ns; FPN(fseq),j++) {
      if (com.spname[j]) fprintf (fseq, "\n%s\n", com.spname[j]);
      else               fprintf (fseq, "\nseq.%d\n", j+1);
      for (h=0; h<com.ls; h++) {
         ic=com.z[j][pose?pose[h]:h];
         if(ic<1 || ic>Nsensecodon[com.icode])  error("err: OutSeqsCodon");
         ic=FROM61[ic-1]; pc=getcode(ic);
         fprintf (fseq, "%c%c%c ", pc[0],pc[1],pc[2]);
      }
   }
   fflush (fseq);
   return (0);
}



int DistanceMatAA (FILE *fout, double alpha)
{
   int i,j, h;
   double p;

   if (fout)  fprintf(fout,"\ndistances (alpha set at %.2f)\n", alpha);
   FOR (i, com.ns) {
      if (fout) fprintf (fout, "\n%-15s", com.spname[i]);
      FOR (j, i) {
         for (h=0,p=0; h<com.npatt; h++)  
            if (com.z[i][h]!=com.z[j][h]) p+=com.fpatt[h];
         p = p/com.ls;
         SeqDistance[i*(i-1)/2+j]=p= (alpha==0 ? p : alpha*(pow(1-p,-1/alpha)-1));
         if (fout) fprintf (fout, "%7.3f", p);
      }
   }
   if (fout) FPN(fout);
   return (0);
}

int DistanceMatNG86 (FILE *fout, double alpha)
{
/* Nei & Gojobori (1986)
   alpha used for nonsynonymous rates only
*/
   int i,j, h, wname=15, status=0, ndiff,nsd[4];
   char codon1[4], codon2[4];
   double ns, na, nst, nat, S, N, St, Nt, bigD=-1;

   if (fout)  {
      fprintf(fout,"\nNei & Gojobori 1986. dN/dS (dN, dS)");
      fprintf(fout,"\n(This matrix is not used in later m.l. analysis.)\n");
   }
   fputs("\nNumber of codon sites with 0, 1, 2, 3 position differences\n",frst);
   FOR (i,com.ns) {
      if (fout)  fprintf (fout, "\n%-*s", wname, com.spname[i]);
      FOR (j,i) {
         FOR (h,4) nsd[h]=0;
         for (h=0,nst=nat=St=Nt=0; h<com.npatt; h++)  {
            strcpy (codon1, getcode(FROM61[com.z[i][h]-1]));
            strcpy (codon2, getcode(FROM61[com.z[j][h]-1]));
            ndiff=difcodon (codon1, codon2, &S, &N, &ns, &na, 0, com.icode);
            nsd[ndiff]+=(int)com.fpatt[h];

            St+=S*com.fpatt[h];
            Nt+=N*com.fpatt[h];
            nst+=ns*com.fpatt[h];
            nat+=na*com.fpatt[h];
         }

         fprintf (frst, "\n%4d%4d %6d%5d%5d%5d", 
            i+1, j+1,nsd[0],nsd[1],nsd[2],nsd[3]);

         if (noisy>=9)
           printf("\n%3d%3d:Site%8.2f +%8.2f =%8.2f\tDiff%8.2f +%8.2f =%8.2f",
             i+1,j+1,St,Nt,St+Nt,nst,nat, nst+nat);

         S = 1-4./3*nst/St;
         N = 1-4./3*nat/Nt;
         if (S<=0 || N<=0)  { puts ("\nlarge distance.."); status=-1; }
         if (S==1) nst=0;
         else      nst=(S<0?bigD:3./4*(-log(S)));
         if (N==1) nat=0;
         else nat=(N<0?bigD:3./4*(alpha==0?-log(N):alpha*(pow(N,-1/alpha)-1)));

         SeqDistance [i*(i-1)/2+j] = (0.21*nst+0.79*nat)*3;
         if (nst>0) S=nat/nst;  if(S<0) S=-1;
         if (fout) fprintf (fout, "%7.3f(%6.4f %6.4f)", S, nat, nst);
      }
   }
   if (fout) FPN (fout);
   if (status) fprintf (fout, "NOTE: -1 means the method is inapplicable.\n");
   FPN(frst);  FPN(frst); fflush (frst);
   return (0);
}


int GetDaa (FILE* fout, double daa[])
{
/* Get the amino acid distance (or substitution rate) matrix 
   (grantham, dayhoff, jones, etc).
*/
   FILE * fdaa;
   int i,j, naa=20;
   double dmax=0,dmin=1e40;

int ncat=6, k, ncount;
double d=0;

   if ((fdaa=fopen(com.daafile, "r"))==NULL) 
      { printf("\nAA dist file %s not found.", com.daafile); exit(-1); }
   printf ("\nReading matrix from %s..\n", com.daafile);
   if (com.model==REVaa_0||com.model==REVaa) puts("To get initial values.");

   for (i=0; i<naa; i++)  for (j=0,daa[i*naa+i]=0; j<i; j++)  {
      fscanf(fdaa, "%lf", &daa[i*naa+j]);
      daa[j*naa+i]=daa[i*naa+j];
      if (dmax<daa[i*naa+j]) dmax=daa[i*naa+j];
      if (dmin>daa[i*naa+j]) dmin=daa[i*naa+j];
   }
   if(com.aaDist && (com.seqtype==1||com.model==FromCodon)) { /* codon model */
      if(noisy) printf("distance: %.2f --- %.2f\n", dmin, dmax);
      FOR (i,naa) FOR(j,naa) com.daa[i*naa+j]/=dmax;
   }
   else if (com.seqtype==AAseq && com.model==Empirical) {
      FOR (i,naa)  fscanf(fdaa, "%lf", &com.pi[i]);
      if (fabs(1-sum(com.pi,20))>1e-4) 
         error("sum of freqs != 1 in the aa rate file");
   }
   fclose (fdaa);

   if (fout) {
      fprintf (fout, "\n%s\n", com.daafile);
      FOR (i,naa) {
         fprintf (fout, "\n%4s", getaa (0, i+1));
         FOR (j,i)  fprintf (fout, "%5.0f", daa[i*naa+j]); 
      }
      FPN (fout);
   }

#if 0
printf("distance: %.2f --- %.2f\n", dmin, dmax);
printf("\ncat? ");
scanf("%d", &ncat);
dmax*=1.000000001;
d=(dmax-dmin)/ncat;

SetAA1STEP();  
for(i=0,FPN(F0);i<naa;i++,FPN(F0)) 
   FOR(j,i)printf("%3d",AA1STEP[i*(i-1)/2+j]);
for (k=0,ncount=0;k<ncat; k++) {
   printf ("\nd (%7.1f, %7.1f): ", dmin+k*d, dmin+(k+1)*d);
   for(i=0;i<naa;i++) FOR(j,i) {
      if (AA1STEP[i*(i-1)/2+j]==0) continue;
      if (com.daa[i*naa+j]>=dmin+k*d && com.daa[i*naa+j]<dmin+d*(k+1.)) {
         printf (" %c%c", AAs[i],AAs[j]);
         ncount++;
      }
   }
}

printf("\n\n%d pairs in all classes.\n", ncount);
exit (0);
#endif

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
      ic1=GenetCode[com.icode][ic1]-1;
      ic2=GenetCode[com.icode][ic2]-1;
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
      if (Q[i*naa+j]!=Q[j*naa+i]) error("strange in SetAA1STEP");
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
/* This routine reads from the file OmegaAA.dat to initialize the
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
       printf("\n\n%d one-step aa pairs\n", n1step);
       printf("\nReading AA classes from %s\n", OmegaAAf);
   }
   if ((fin=fopen(OmegaAAf,"r"))==NULL) 
      { printf("file %s does not exist.", OmegaAAf);  fromfile=0; }
   else { 
      fscanf(fin, "%d", &N_OMEGA);  
      N_OMEGA++;
      printf ("\n%d classes of dN/dS, including background class.\n", N_OMEGA);
      if (N_OMEGA>65-1) { fromfile=0;  fclose(fin); }
   }
   if (!fromfile) {
      puts("Estimating the aa fixation rate matrix.  Stop if not.");
      getchar ();
      if (com.seqtype!=CODONseq) puts("\nTo be tested.\a\a");
      N_OMEGA=0;
      if (com.seqtype==AAseq) {
         FOR(i,naa) FOR(j,i) if(i*naa+j!=ijAAref && AA1STEP[i*(i-1)/2+j])
             OmegaAA[i*(i-1)/2+j]=N_OMEGA++;
      }
      else
         FOR(i,naa) FOR(j,i) 
           if(AA1STEP[i*(i-1)/2+j]) OmegaAA[i*(i-1)/2+j]=N_OMEGA++;
   }
   else { 
      FOR (iomega, N_OMEGA-1) {
         fscanf(fin, "%d", &j);
         if (j!=iomega+1) { printf("err data file %s.", OmegaAAf); exit(-1); } 
         printf ("\nClass #%d: ", j);
         j=fgetc (fin);  if (j!=':') error("err expecting :");
         fgets (line, nline, fin);
      
         /* printf ("%s\n", line); */
         for (j=0,npair=0; j<nline-1&&line[j]&&line[j]!='\n'; j++) {
            iaa=line[j];
            if (!isalpha(iaa)) continue;
            jaa=line[++j];  if(!isalpha(jaa)) error("err jaa");
            npair++;

            /* printf ("\npair %d: |%c%c| ", npair, iaa,jaa); */
            iaa=CodeChara((char)iaa,AAseq)-1; jaa=CodeChara((char)jaa,AAseq)-1;
            if (iaa<jaa)  { k=jaa, jaa=iaa; iaa=k; }
      
            /* printf ("|%c%c (%2d,%2d)| ", AAs[iaa], AAs[jaa],iaa,jaa); */
            if (iaa==jaa) printf ("This pair has no effect.");
            if (OmegaAA[iaa*(iaa-1)/2+jaa]) error("This pair specified?");
            if (AA1STEP[iaa*(iaa-1)/2+jaa]==0) error("This pair has rate 0!");
            OmegaAA[iaa*(iaa-1)/2+jaa]=iomega+1;
            /*  printf (" in class %d ", iomega+1); */
         }
      }
   }
   com.nrate = !com.fix_kappa + N_OMEGA;
   printf ("\nNw=%d\tnrate=%d\n", N_OMEGA, com.nrate);
   if (N_OMEGA>n1step-(com.seqtype==AAseq)) error("too many classes.");
/*
   for (i=0; i<naa; i++,FPN(F0)) 
       FOR(j,i) printf ("%3d",OmegaAA[i*(i-1)/2+j]);
*/
   return (0);
}


#ifdef UNDEFINED

int GetCategoryQ61 (char z[NS])
{
/* the category ID for a codon site with z[NS], transformed
   classified into 19 categories 
*/
   int i,j, icat, ncat=19, it, b[NS][3], nbase[3], markb[4];

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

int TestModelQ61 (FILE * fout, double x[])
{
/* Test the Q61 model, slower than simulations
*/
   char z[NS];
   int h, npatt, it, icat, j, nodeb[NS], imposs;
   int n=Nsensecodon[com.icode], isum, nsum, ncat=19;
   double  fh, y, nobs[19], pexp[19], Pt[8][NCODE*NCODE];

   printf ("\ntest Q61..\n");
   for (h=0,zero(nobs,ncat); h<com.npatt; h++) {
      for (j=0; j<com.ns; j++) z[j]=com.z[j][h]-1;
      icat = GetCategoryQ61(z);
      nobs[icat]+=com.fpatt[h];
   }
   FOR (j,ncat) 
      printf("cat #%4d: %4d%4d%4d%6.0f\n", j+1,j/9+1,(j/3)%3+1,j%3+1,nobs[j]);

   if (com.ns>5 || com.alpha || com.ngene>1)
      error ("TestModelQ61: ns>5 || alpha>0.");
   if (SetParameters (x)) puts ("\npar err..");
   for (j=0,npatt=1; j<com.ns; j++)  npatt*=n;
   for (isum=0,nsum=1; isum<tree.nnode-com.ns; nsum*=n,isum++) ;
   printf("\nTest Q61: npatt = %d\n", npatt);
   FOR (j, tree.nbranch) 
      PMatUVRoot (Pt[j], nodes[tree.branches[j][1]].branch, n, U, V, Root);

   for (h=0,zero(pexp,ncat); h<npatt; h++) {
      for (j=0,it=h; j<com.ns; nodeb[com.ns-1-j]=it%n,it/=n,j++) ;
      for (j=0,imposs=0; j<com.ns; j++) 
         { z[j]=nodeb[j];  if (com.pi[(int)z[j]]==0) imposs=1; }
      if (imposs) continue;
      
      if ((icat=GetCategoryQ61(z)) == ncat-1) continue;
      if ((h+1)%100==0) 
         printf("\rTest Q61:%9d%4d%9.2f%%", h+1, icat, 100.*(h+1.)/npatt);

      for (isum=0,fh=0; isum<nsum; isum++) {
         for (j=0,it=isum; j<tree.nbranch-com.ns+1; j++)
            { nodeb[com.ns+j]=it%n; it/=n; }
         for (j=0,y=com.pi[nodeb[tree.origin]]; j<tree.nbranch; j++) 
            y*=Pt[j][nodeb[tree.branches[j][0]]*n+nodeb[tree.branches[j][1]]];
         fh += y;
      }
      if (fh<=0) {
         matout (F0, x, 1, com.np);
         printf ("\a\ntest Q61: h=%4d  fh=%9.4f \n", h, fh);
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


void GetFreqPair(int codonf, char *z1,char *z2,int npatt,float fpatt[],
     double pi[])
{
/* Get the expected codon frequencies for a pair of transformed codon seqs.
*/
   int n=com.ncode, i,j, ic,c[3];
   double fb3x4[12], fb4[4];

   for(i=0,zero(pi,n); i<npatt; i++) {
      pi[z1[i]-1]+=fpatt[i]/(2.*com.ls);  
      pi[z2[i]-1]+=fpatt[i]/(2.*com.ls);  
   }
   if (codonf==3) return;
   else if (codonf==0) { fillxc(pi,1./n,n); return; } 

   for (i=0,zero(fb3x4,12),zero(fb4,4); i<n; i++) {
      ic=FROM61[i];  c[0]=ic/16; c[1]=(ic/4)%4; c[2]=ic%4;
      FOR (j,3) { fb3x4[j*4+c[j]]+=pi[i];  fb4[c[j]]+=pi[i]/3.; }
   }
   for (i=0; i<n; i++) {
      ic=FROM61[i];  c[0]=ic/16; c[1]=(ic/4)%4; c[2]=ic%4;
      if (codonf==2)  pi[i]=fb3x4[c[0]]*fb3x4[4+c[1]]*fb3x4[8+c[2]];
      else            pi[i]=fb4[c[0]]*fb4[c[1]]*fb4[c[2]];
   }
   abyx (1./sum(pi,n), pi, n);
   return;
}

int Pairwise (FILE *fout, double space[])
{
/* calculate ds & dn for all pairwise codon sequence comparisons
   use different npatt for different pairs
*/
   char *pz0[NS];
   int  n=com.ncode, ii, jj,k,h, np;
   int  ns0=com.ns, npatt0=com.npatt;
   float *fpatt0, fp[NCODE*NCODE];
   double x[NP], xb[NP][2], lnL, e=1e-8, *var=space+NP;

   if (com.ngene>1) error("ngene>1 not supported");
   if (noisy) printf("\n\npairwise comparison (Goldman & Yang 1994)..\n");
   fprintf(fout,"\npairwise comparison, codon frequencies: %s.\n",
      codonfreqs[com.codonf]);

   fprintf(frst, "\n\npairwise comparison (Goldman & Yang 1994)");

   fprintf(frst,
      "\nseq seq        N       S        dN       dS    dN/dS    Paras.\n");
   tree.branches[0][0]=0;  tree.branches[0][1]=1;
   tree.nnode = (com.ntime=tree.nbranch=1) + 1;
   tree.origin=0;
   com.ns=2;  com.clock=0;
   BranchToNode ();

   fpatt0=(float*) malloc(npatt0*3*sizeof(float));
   FOR (k, ns0) pz0[k]=com.z[k];
   com.z[0]=(char*)(fpatt0+npatt0);   com.z[1]=com.z[0]+npatt0;

   com.chunk =(double*) malloc ((com.ns-1)*n*com.npatt*sizeof(double));
   if (com.chunk==NULL||fpatt0==NULL) error ("oom");
   FOR (h, npatt0) fpatt0[h]=com.fpatt[h];
   PointLklnodes ();

   FOR (ii, ns0) {
      FOR (jj, ii) {
         printf ("\n\n%4d vs. %3d", ii+1, jj+1);
         fprintf (fout, "\n\n%d (%s) ... %d (%s )",
              ii+1,com.spname[ii], jj+1,com.spname[jj]);
         fprintf (frst, "%3d %3d ", ii+1, jj+1);

         FOR (k,n*n) fp[k]=0;
         GetFreqPair (com.codonf, pz0[ii], pz0[jj], npatt0, fpatt0, com.pi);
         for(h=0; h<npatt0; h++) {
            fp[(pz0[ii][h]-1)*n+pz0[jj][h]-1]+=fpatt0[h];
            fp[(pz0[jj][h]-1)*n+pz0[ii][h]-1]+=fpatt0[h];
         }
         for (k=0,com.npatt=0; k<n; k++) for (h=0; h<=k; h++) {
            if (fp[k*n+h]) { 
               com.z[0][com.npatt]=(char)(k+1); 
               com.z[1][com.npatt]=(char)(h+1);
               if (h==k) com.fpatt[com.npatt++]=fp[k*n+h]/2;
               else      com.fpatt[com.npatt++]=fp[k*n+h];
            }
         }
         com.posG[1]=com.npatt;
         GetInitials (x);   np=com.np;  NFunCall=0;
         SetxBound (np, xb);
         ming2(noisy>2?frub:NULL,&lnL,com.plfun,NULL,x,xb, space,e,np);

         fprintf (fout,"\nlnL%14.6f\n", -lnL);
         FOR (k, np)  fprintf (fout," %8.5f", x[k]);  FPN (fout);

         if (com.getSE) {
            Hessian (np, x, lnL, space, var, com.plfun, var+np*np);
            matinv (var, np, np, var+np*np);
            fprintf (fout,"SEs for parameters:\n");
            FOR (k, np) var[k*np+k] = (var[k*np+k]>0.?sqrt(var[k*np+k]):-0);
            FOR (k,np) fprintf (fout," %8.5f", var[k*np+k]);  FPN (fout);
         }
         FPN (fout);
         _BRANCHLEN=x[0];
         EigenQ61(fout, Root, U, V, com.kappa, com.omega, PMat);

         FOR (k, np)  fprintf (frst," %8.3f", x[k]);  
         if (com.getSE) fprintf (frst, " +-%7.3f", var[np*np-1]);
         FPN (frst);

/*
         fprintf(frst, "%7.3f%7.3f%7.3f%10.2f", x[0],com.kappa,com.omega,-lnL);
         FPN (frst);
*/
         fflush(frst);  fflush(fout);
/*
         printf ("%14.6f%6d%6d", -lnL, com.npatt, NFunCall);
         FOR (k,com.np) printf ("%10.5f", x[k]);  FPN(F0);
*/
      }  /* for (jj) */
   }     /* for (ii) */
   free(fpatt0);  free(com.chunk);
   FPN(F0);   return (0);
}

double TreeScore(double x[], double space[])
{
   double xb[NP][2], e=1e-6, lnL;

   PointLklnodes ();
   com.ntime = com.clock ? tree.nnode-com.ns : tree.nbranch;
   GetInitials (x);

   SetxBound (com.np, xb);
   ming2(NULL,&lnL,com.plfun,NULL,x,xb, space,e,com.np);
   return (lnL);
}


/* modify this to calculate S and N for a pairwise method

int NSQ61(double kappa, double pi[])
{
   int n=Nsensecodon[com.icode], i,j,k, ic1,ic2, ndiff, pos=0, from[3],to[3];
   double rs0, ra0, t, Qij;

   for (i=0,rs0=ra0=0; i<n; i++) FOR (j,n) {
      ic1=FROM61[i]; from[0]=ic1/16; from[1]=(ic1/4)%4; from[2]=ic1%4;
      ic2=FROM61[j];   to[0]=ic2/16;   to[1]=(ic2/4)%4;   to[2]=ic2%4;
      for (k=0,ndiff=0; k<3; k++)  if (from[k]!=to[k]) { ndiff++; pos=k; }
      if (ndiff!=1)  continue; 
      t=2*pi[i];
      if ((from[pos]+to[pos]-1)*(from[pos]+to[pos]-5)==0)  Qij=kappa;
      else   Qij=1;
      if((ic1=GenetCode[com.icode][ic1]-1)!=(ic2=GenetCode[com.icode][ic2]-1))
         ra0+=t*Qij;
      else 
         rs0+=t*Qij;
   }
   t=rs0+ra0; rs0/=t;  ra0/=t;
   matout(F0, pi, 1, n);
   printf ("\n%9.3f%9.3f%9.1f%9.1f\n", rs0, ra0, 3*com.ls*rs0, 3*com.ls*ra0);
   return (0);
}
*/

double lfunNSsites (double x[], int np)
{
/* variable d_N/d_S ratios among sites
*/
   int  h, ir, i, ig;
   double lnL, fh=0;

   NFunCall++;
   if (com.seqtype!=1) error("seqtype..");
   if (SetParameters (x)) puts ("\npar err..");
   zero (com.fhK, com.npatt);
   for (ig=0; ig<com.ngene; ig++) {
      SetPGene(ig, com.Mgene>1, com.Mgene>1, com.nalpha>1, x);
      for (ir=0; ir<com.ncatG; ir++) {
         EigenQ61(NULL, Root, U, V, com.kappa, com.rK[ir], PMat);
         PartialLikelihood (tree.origin, ig);
         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            if (com.fpatt[h]==0) continue;
            for (i=0,fh=0; i<com.ncode; i++)
               fh += com.pi[i]*nodes[tree.origin].lkl[h*com.ncode+i];
            com.fhK[h] += com.freqK[ir]*fh;
         }
      }
   }
   for (h=0,lnL=0; h<com.npatt; h++) {
      if (com.fhK[h]<=0) {
         matout (F0, x, 1, np);
         printf ("\nlfunNSsites: h=%4d  fhK=%9.4f", h, com.fhK[h]);
         getchar ();
      }
      lnL-=log(com.fhK[h])*com.fpatt[h];
      if (com.print<0)
         fprintf (flfh, "\n%6d%6.0f%16.10f%9.2f  ", 
            h+1,com.fpatt[h],log(com.fhK[h]),com.ls*com.fhK[h]);
   }
   return (lnL);
}

int lfunNSsites_rate (FILE* fout, double x[], int np)
{
/* This calculates the dN/dS rates for sites under models with variabel dN/dS 
   ratios among sites (Nielsen and Yang 1998).  Modified from lfundG() 
   Only conditional mode is used?
*/
   int  ig, h, ir, it=0, i;
   double lnL, fh, t=0;

   if (SetParameters (x)) puts ("par err. lfunNSsites_rate");
   fprintf (fout,"\nFrequencies for categories (K=%d)", com.ncatG);
   matout (fout, com.freqK, 1, com.ncatG);
   fprintf (fout,"\ndN/dS rate ratios for categories (K=%d)", com.ncatG);
   matout (fout, com.rK, 1, com.ncatG);
   zero (com.fhK, com.npatt*com.ncatG);
   for (ig=0; ig<com.ngene; ig++) {
      SetPGene(ig, com.Mgene>1, com.Mgene>1, com.nalpha>1, x);
      for (ir=0; ir<com.ncatG; ir++) {
         EigenQ61(NULL, Root, U, V, com.kappa, com.rK[ir], PMat);
         PartialLikelihood (tree.origin, ig);
         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            if (com.fpatt[h]==0) continue;
            for (i=0,fh=0; i<com.ncode; i++)
               fh += com.pi[i]*nodes[tree.origin].lkl[h*com.ncode+i];
            com.fhK[ir*com.npatt+h] = com.freqK[ir]*fh;
         }
      }
   }
   for (h=0,lnL=0; h<com.npatt; h++) {
       for (ir=0,fh=0,it=0,t=0; ir<com.ncatG; ir++) {
          if (com.fhK[ir*com.npatt+h]>t) { t=com.fhK[ir*com.npatt+h]; it=ir; }
          fh+=com.fhK[ir*com.npatt+h];
       }
/*
fprintf (fout, "\n%4d%7.0f%9.1f: ", h+1, com.fpatt[h], com.rK[it]);
FOR (ir,com.ncatG)  fprintf (fout, " %9.4f", com.fhK[ir*com.npatt+h]/fh);
*/
       lnL-=com.fpatt[h]*log(fh);
       com.fhK[h]=com.rK[it];  com.fhK[com.npatt+h]=t/fh;
   }
   fprintf(fout,"\n\ndN/dS ratios along sequence:\n%4s%9s%9s\n",
         "Site","Rate","Prob");
   FOR (h, com.ls) {
      it=com.pose[h];
      fprintf(fout,"\n%4d%9.4f%9.4f",h+1,com.fhK[it],com.fhK[com.npatt+it]);
   }
   fprintf (fout,"\n\nlnL=%14.6f\n\n", -lnL);
   return (0);
}


/*

(1) Codon models for variable d_N/d_S ratios among lineages
    model=0: one d_N/d_S ratio for all lineages.
             fix_omega=1 will fix the single omega at the value provided
    model=1: each branch has its own omega, use fix_omega=0 only
    model=2: 2 to nbranch ratios.  Use consecutive branch marks 0,1,2,... 
             to specify omega.
             fix_omega=1 fixes the last omega at the value provided.
 
    Some computation can be saved under com.model==NS2 if only 2 or a few more 
    ratios are used instead of all b ratios.  Check whether it is worth the 
    effort.


(2) Amino acid models
    REVaa: The symmetrical part (S) of the rate matrix Q=S*PI are estimated, 
           making up 19*20/2-1=189 rate parameters for the matrix.  The aa 
           frequencies are estimated using the observed ones.  The Sij for 
           ijAAref=19*naa+9 (I-V) is set to one and others are relative rates;
    REVaa_0: AA1STEP[i*(i+1)+j] marks the aa pair i & j that are 
            interchangeable.  Sij for ijAAref=19*naa+9 (I-V) is set to one 
            and others are relative rates;


(3)
    Codon & amino acid models

    AAClasses: OmegaAA[i*(i-1)/2+j] marks the dN/dS ratio class for the pair 
            i & j.  Note kappa is before omega's in x[].
            OmegaAA[i*(i-1)/2+j]=-1, if AAs i & j are not interchangeable
                       =0,  for the background ratio
                       =1,...,nclass for AAs i & j specified in OmegaAA.dat.
            The total number of classes (N_OMEGA) is one plus the number 
            specified in the file OmegaAAf.

   N_OMEGA is the number of different dN/dS ratios in the NSb, NS2 models
      and in AAClasses. 
   OmegaBranch[ibranch]=0,1,...,(N_OMEGA-1) marks the dN/dS ratio for ibranch 
      in the NSb NS2 models
   AA1STEP[i*(i-1)/2+j] =1 if AAs i & j differ at one codon position;
                        =0 otherwise.

(4)
   Geometric and linear relationships between amino acid distance and 
   substitution rate:
      wij = a*(1-b*dij/dmax)
      wij = a*exp(-b*dij/dmax)
   aaDist = 0:equal, +:geometric; -:linear, {1-5:Grantham,Miyata,c,p,v}
*/
