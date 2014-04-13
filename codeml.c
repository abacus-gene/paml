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
#define NGENE         50
#define LSPNAME       30
#define NCODE         64
#define NCATG         40
/*
#define NP            (NBRANCH*2+NGENE-1+2)
*/
#define NP            (NBRANCH+NGENE-1+189+2)

extern char BASEs[],AAs[],GenetCode[][64],Nsensecodon[];
extern int noisy, NFunCall, *ancestor;
extern double *SeqDistance;

int Forestry (FILE *fout, double space[]);
void DetailOutput (FILE *fout, double x[], double var[], double space[]);
int GetOptions (char *ctlf);
int testx (double x[], int np);
int SetxBound (int np, double xb[][2]);
int SetxInitials (double x[]);
int GetInitials (double x[]);
int SetParameters (double x[]);
int SetPGene (int igene, int _pi, int _UVRoot, int _alpha, double x[]);
int PMatJC69like (double P[], double t, int n);
int setmark_61_64 (void);
int printfcode (FILE *fout, double fb61[], double space[]);
int printsmaCodon (FILE *fout,char * z[],int ns,int ls,int lline,int simple);
int InitializeCodon (FILE *fout, double space[]);
int AA2Codonf (double faa[20], double fcodon[]);
int DistanceMatAA (FILE *fout);
int DistanceMatNG86 (FILE *fout, double alpha);
int GetDaa(FILE *fout, double daa[]);
int QcStatistics (double branchl, double *S, double *dS, double *dN,
    double kappa, double omega, double Q[64*64]);
int EigenQc (double Root[], double U[], double V[], 
    double kappa, double omega, double Q[64*64]);
int EigenQaa(FILE *fout, double Root[], double U[], double V[],double rate[]);
int Qcodon2aa (double Qc[], double pic[], double Qaa[], double piaa[]);
int SetAA1STEP (void);
int GetOmegaAA(int OmegaAA[]);
int TestModelQc(FILE *fout, double x[]);
double lfun2dSdN (double x[], int np);
int VariancedSdN(double t, double omega, double vtw[2*2], double vdSdN[2*2]);
int GetFreqPairCodon(void);
int PairwiseCodon (FILE *fout, double space[]);
int PairwiseAA (FILE *fout);
double lfunNSsites (double x[], int np);
int lfunNSsites_rate (FILE* fout, double x[], int np);
void SimulateData2s61(void);
void Ina(void);

struct common_info {
   char *z[NS], *spname[NS], seqf[96],outf[96],treef[96],daafile[96];
   int seqtype, ns, ls, ngene, posG[NGENE+1], lgene[NGENE], npatt,*pose;
   int runmode,clock,verbose,print, codonf,aaDist,model,NSsites, cleandata;
   int icode, ncode, Mgene, ndata;
   int fix_rgene, fix_kappa, fix_omega, fix_alpha, fix_rho, getSE;
   int np, ntime, nrgene, nrate, nalpha, ncatG, maxnp;
   double *fpatt;
   double pi[NCODE],fb61[64],piG[NGENE][64],kappa,omega,alpha,rho,rgene[NGENE];
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
double PMat[NCODE*NCODE],U[NCODE*NCODE],V[NCODE*NCODE],Root[NCODE];
double *XCOPY;
  /* communication between SetParameters & PartialLikelihood & EigenQc.  
     Remove it? */
int LASTROUND=0;

double OMEGA_FIX=-1;  
/* fix the last dN/dS in the NSb, NS2 models with variable dN/dS ratios 
   for lineages.  Useful for testing whether w>1 for particular lineages. */
int N_OMEGA=-99, *OmegaBranch=NULL, OmegaAA[190], AA1STEP[190];

int *rateBranch=NULL, N_rateBranch=-1;

char *_nScale;     /* nScale[ns-1] for interior nodes */
double *_nScaleF=NULL;  /* nScaleF[npatt] for scale factors */

FILE *frub, *flfh, *frst, *frst1;
char *ratef="rates";
enum {Fequal, F1x4, F3x4, Fcodon} CodonFreqs;
enum {NS1, NSb, NS2} CodonModels;   /* , AAClasses=7 */
enum {Poisson, EqualInput, Empirical, Empirical_F,
     FromCodon=6, AAClasses=7, REVaa_0=8, REVaa=9} AAModel;

char *codonfreqs[]={"Fequal", "F1x4", "F3x4", "Fcodon"};
char *codonmodels[]={"One dN/dS ratio for branches", 
     "free dN/dS Ratios for branches", "several dN/dS ratios for branches",
     "", "", "", "","AAClasses"};
char *NSsitesmodels[]={"one-ratio","neutral model","selection model",
    "deleterious mutation", "deleterious mutation&selection"};
char *aamodels[]={"Poisson", "EqualInput", "Empirical", "Empirical_F", "",
     "", "FromCodon", "AAClasses", "REVaa_0", "REVaa"};

#define CODEML 1

#include "treesub.c"
#include "treespace.c"

int main(int argc, char *argv[])
{
   FILE *fout, *fseq;
   char ctlf[32]="codeml.ctl", *pmodel;
   char *seqtypestr[3]={"CODONML", "AAML", "CODON2AAML"};
   char *Mgenestr[]={"diff. rate", "separate data", "diff. rate & pi", 
                     "diff. rate & kappa", "diff. rate & pi & kappa"};
   int i, s1=0, s2=0, s3=0, idata;
   double *space=NULL;
   clock_t clock_start, clock_end;

/*
int j;
double d, dmax=0;
double aaprop[2][20]={{31, 124, 56, 54, 55, 85, 83, 3, 96, 111, 111,
119, 105, 132, 32.5, 32, 61, 170, 136, 84}, 
			{ -0.11, 0.079, -0.136, -0.285, -0.184, -0.067, -0.246, -0.073, 0.32, 0.001, -0.008, 0.049, -0.041, 0.438, -0.016, -0.153, -0.208, 0.493, 0.381, -0.155}};

for(i=0;i<20;i++,FPN(F0)) FOR(j,i) {
   d=fabs(aaprop[0][i]-aaprop[0][j]);
   if(d>dmax) dmax=d;
   printf(" %5.3f",d);
}
printf("dmax: %5.3f\n",dmax);
exit(0);
*/


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
   com.fix_rho=0;     com.rho=0.4; 
   com.getSE=0;       com.print=0;    com.verbose=1; 

   frub=fopen("rub","w");   frst=fopen("rst","w");  frst1=fopen("rst1","w");
   flfh=fopen("lfh", "w");

/*
printf ("\nrandom number seed? ");
scanf ("%d", &i);  SetSeed(i);
Ina();
SimulateData2s61();
*/
   if(argc>1) strcpy(ctlf,argv[1]);
   GetOptions (ctlf);
   if ((fout=fopen (com.outf, "w"))==NULL) error("outfile creation err.");
   if((fseq=fopen (com.seqf,"r"))==NULL)  {
      printf ("\n\nSequence file %s not found!\n", com.seqf);
      exit (-1);
   }
   if(noisy && com.seqtype==CODONseq) 
    { printcu(F0,NULL,com.icode); puts("The code is nice, uuh?"); }

   for (idata=0; idata<com.ndata; idata++) {
      if (com.ndata>1) {
         printf ("\nData set %d\n", idata+1);
         fprintf (fout, "\n\nData set %d\n", idata+1);
      }
      if (idata)  GetOptions (ctlf); /* Is this necessary? */

      ReadSeq((com.verbose?fout:NULL),fseq,com.seqtype); /*may change seqtype*/

      if(com.ngene>1 && com.Mgene==1)  printSeqsMgenes ();
      if(com.Mgene==1 && com.print) puts("Mgene & RateAncestor not tested");

      if (com.seqtype==CODONseq && (com.model==NSb||com.model==NS2)) {
        if(OmegaBranch) free(OmegaBranch);
        if((OmegaBranch=(int*)malloc(NBRANCH*sizeof(int)))==NULL) error("oom");
      }
     pmodel=(com.seqtype==CODONseq?codonmodels[com.model]:aamodels[com.model]);
    fprintf(fout,"%s %s Model: %s ",seqtypestr[com.seqtype-1],com.seqf,pmodel);
      if(com.seqtype==CODONseq||com.model==FromCodon) {
         if (com.fix_kappa)  fprintf(fout, " kappa = %.3f\n", com.kappa);
         if (com.fix_omega) fprintf(fout, " omega? = %.3f fixed\n", com.omega);
      }
      if (com.seqtype==AAseq && (com.model==Empirical||com.model==Empirical_F))
         fprintf (fout, "(%s) ", com.daafile);
      /* if (com.nrate) fprintf(fout, "(nrate:%d)  ", com.nrate); */
      if(com.clock) fprintf(fout," %s clock  ",com.clock==1?"Global":"Local");
      if (com.alpha && com.rho) fprintf (fout, "Auto-");
      if (com.alpha) fprintf (fout,"dGamma (ncatG=%d) ", com.ncatG);
      if (com.ngene>1) 
         fprintf (fout, " (%d genes: %s)  ", com.ngene, Mgenestr[com.Mgene]);
      if (com.alpha==0) com.nalpha=0;
      else              com.nalpha=(com.nalpha?com.ngene:!com.fix_alpha);
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

      if (com.clock==2) {
        if(rateBranch) free(rateBranch);
        if((rateBranch=(int*)malloc(com.ns*2*sizeof(int)))==NULL) error("oom");
      }
      com.maxnp = max((2*com.ns-2)*2+com.nrate+com.ngene-1+2, 50);
/*
printf ("\nmaxnp = %d\n", com.maxnp);
getchar ();
*/
      s3 = com.maxnp*(com.maxnp*2+2+12)*sizeof(double);
      s3 = max(s3, com.ls*(com.ns*3*sizeof(char)+sizeof(int)));
      i = com.ns*(com.ns-1)/2;
      s3 = max(s3, sizeof(double)*((com.ns*2-2)*(com.ns*2-2 + 4 + i) + i));
      if(com.seqtype==CODONseq) s3=max(s3,com.ns*com.ngene*sizeof(double));
      if (space) free(space);
      if ((space=(double*)malloc(s3))==NULL)  error("oom space");

      if (SeqDistance) free(SeqDistance);      
      if (ancestor) free(ancestor);
      SeqDistance=(double*)malloc(com.ns*(com.ns-1)/2*sizeof(double));
      ancestor=(int*)malloc(com.ns*(com.ns-1)/2*sizeof(int));
      if (SeqDistance==NULL||ancestor==NULL) error("oom");

      if(com.seqtype==AAseq) { 
         Initialize (fout, space, com.seqtype);
         if (com.model==FromCodon||com.model==AAClasses) 
            AA2Codonf(com.pi, com.fb61);  /* get codon freqs from aa freqs */ 
      }
      else
         if(InitializeCodon(fout,space)) error("giving up on stop codons");


      if (com.seqtype==CODONseq) DistanceMatNG86 (fout, 0);
      else                       DistanceMatAA (fout);
      fflush(fout);

      if (com.Mgene==3) FOR(i,com.ngene) xtoy(com.pi,com.piG[i],com.ncode);

      if (com.seqtype==AAseq && com.model==Poisson && !com.print) 
         PatternJC69like (NULL);

      if (com.alpha || com.NSsites) {
         if(com.fhK) free(com.fhK);
         com.fhK=(double*)malloc(s2=com.npatt*com.ncatG*sizeof(double));
         if(com.fhK==NULL) error("oom");
      }
      if (com.runmode==-2 && com.Mgene!=1) {
         if (com.seqtype==CODONseq) PairwiseCodon(fout, space);  
         else                       PairwiseAA(fout);  
      }
      else {
         printf("\n%9ld bytes for distance",com.ns*(com.ns-1)/2*sizeof(double));

         s1=(com.ns*(1+!com.cleandata)-1)*com.ncode*com.npatt*sizeof(double);
         if(com.chunk) free(com.chunk);
         if((com.chunk=(double*)malloc(s1))==NULL) error("oom chunk");
         printf("\n%9d bytes for partial likelihood\n", s1);
         printf ("%9d bytes for fhK\n%9d bytes for working space\n", s2,s3);

         if(com.ns>100) {
            if(_nScaleF) free(_nScaleF);
            _nScaleF=(double*)malloc(com.npatt*sizeof(double)+2*com.ns);
            if(_nScaleF==NULL) error("oom");
            _nScale=(char*)(_nScaleF+com.npatt);
         }

         if (com.Mgene==1)         MultipleGenes (fout, space);
         else if (com.runmode==0)  Forestry (fout, space);
         else if (com.runmode==3)  StepwiseAddition (fout, space);
         else if (com.runmode>=4)  Perturbation(fout, (com.runmode==4), space);
         else                      StarDecomposition (fout, space);
         FPN(F0);
      }
   }  /* for (idata) */
   if (noisy) putchar ('\a');

   clock_end=clock();
   printf("\n%.2f seconds used.\n",(((double)clock_end-clock_start))/CLOCKS_PER_SEC);

   return (0);
}


/* x[]: t[ntime]; rgene[ngene-1]; kappa; p[](NSsites); omega[]; 
        { alpha(for NSsites) !! alpha, rho || rK[], fK[] || rK[], MK[] }
*/

int Forestry (FILE *fout, double space[])
{
   FILE *ftree, *fin=(FILE*)fopen("in.codeml","r"), *frate=NULL;
   int  i,j,k, itree, ntree, np, iteration=1;
   int pauptree=0;
   double x[NP],xcom[NP-NBRANCH],xb[NP][2], lnL,lnL0=0,lnLm=0, e=1e-7,tl=0;
   double *var=space+NP, nchange=-1;

   if ((ftree=fopen (com.treef,"r"))==NULL) error ("treefile not found");
   GetTreeFileType(ftree, &ntree, &pauptree, 0);
   if (com.alpha) {
      if ((frate=(FILE*)fopen(ratef,"w"))==NULL) error("file err");
      fprintf(frate,"Rates for sites, from CODEML, %d trees\n", ntree);
   }
   if (ntree>10 && com.print) puts("\nlarge lfh file");
   fprintf (flfh, "%6d%6d%6d\n", ntree, com.ls, com.npatt);
   FOR(i,com.maxnp) xcom[i]=0.1;

   if(!com.cleandata) InitPartialLikelihood ();
   for(itree=0; pauptree||itree<ntree; itree++,iteration=1) {
      if(pauptree && PaupTreeRubbish(ftree)) 
         { puts("err or end of tree file."); break; }
      if(ReadaTreeN(ftree, &i, 1)) error("err or end of tree file.");

      printf ("\nTREE # %2d\n", itree+1);
      fprintf (fout,"\nTREE # %2d:  ", itree+1);
      fprintf (flfh,"\n\n%2d\n", itree+1);
      fprintf (frub,"\n\nTREE #%2d", itree+1);

      LASTROUND=0;

      if (i && !com.clock) {
         puts("\ntree has branch lengths, break to collect into in.codeml?\n");
         FOR(i,tree.nnode) 
            if(i!=tree.origin) x[nodes[i].ibranch]=nodes[i].branch; 
         matout(F0,x,1,tree.nbranch);
         getchar();
      }

      if(com.cleandata) nchange=MPScore(space);
      OutaTreeN(F0,0,0);   printf("   MP score: %.0f",nchange);
      OutaTreeN(fout,0,0); fprintf(fout,"   MP score: %.0f",nchange);
      if(!com.clock && nodes[tree.origin].nson<=2) {
         puts("\nThis is a rooted tree, without clock.  Check.\n");
         if(com.verbose) fputs("\nThis is a rooted tree.  Please check!",fout);
      }
      fflush (fout),  fflush (flfh);  fflush(frst);

      GetInitials(x);
      if ((np=com.np)>NP || np-com.ntime>NP-NBRANCH) error("raise NP");
      printf ("\nntime & nrate & np: %d  %d  %d ", com.ntime,com.nrate,com.np);
      if(itree)  for (i=0; i<np-com.ntime; i++) x[com.ntime+i]=xcom[i];

      if (com.seqtype==CODONseq && (com.model==NSb || com.model==NS2)) {
         printf ("\n%d dN/dS ratios for branches assumed:\n", N_OMEGA);
         FOR (i,tree.nbranch) printf("%4d", OmegaBranch[i]); FPN(F0);
         fprintf (fout, "\n%d dN/dS ratios for branches assumed:\n", N_OMEGA);
         FOR (i,tree.nbranch) fprintf(fout, "%4d", OmegaBranch[i]); FPN(fout);
      }
      if (com.clock==2) {
         printf("\n%d rates for branches assumed:\n",N_rateBranch);
         FOR (i,tree.nbranch) printf("%3d", rateBranch[i]); FPN(F0);
         FPN(F0); FOR(i,tree.nbranch)fprintf(fout,"%3d",rateBranch[i]); FPN(F0);
      }

      if (fin) {
         puts("\nInitials from in.codeml. Break if not correct");
/*         getchar ();
*/
         fscanf(fin,"%lf",&x[i=0]);
         if (x[0]==-1) iteration=0;  else i++;
         for( ;i<np;i++) if(fscanf(fin,"%lf",&x[i])!=1) break;
         if (i<np)  {
            printf("err at #%d in in.codeml. Edit or remove it.\n",i+1);
            exit(-1);
         }
      }
      if(iteration) SetxInitials (x);
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
/*
         gradient (np, x, lnL, space, com.plfun, space+np, 1);
         matout (F0, space, 1, np);
*/
      }
      if(iteration) {
         SetxBound (np, xb);
         j=ming2(noisy>2?frub:NULL,&lnL,com.plfun,NULL,x,xb, space,e,np);
         if(j) fprintf(fout,"\ncheck convergence..");
      }
      printf ("\nOut..\nlnL  = %12.6f\n", -lnL);
      if (itree==0) lnL0=lnLm=lnL;
      else if (lnL<lnLm) lnLm=lnL;
      if (itree==0)
         for (i=0; i<np-com.ntime; i++) xcom[i]=x[com.ntime+i];
      else if (!j)  
         for (i=0; i<np-com.ntime; i++) xcom[i]=xcom[i]*.2+x[com.ntime+i]*0.8;

      if (com.seqtype==CODONseq && com.NSsites==2) {
         /* if (com.getSE>0) */ LASTROUND=1;
         k=com.ntime+com.nrgene+!com.fix_kappa;  /* p0, p1, , w */
         tl=(x[k+0]=exp(x[k+0]))+(x[k+1]=exp(x[k+1]))+1;
         x[k+0]/=tl;  x[k+1]/=tl;
         /* for(i=1; i<com.ncatG-1;i++)   com.freqK[i]=tl/(com.ncatG-2); */
      }
      fprintf (fout,"\nlnL(ntime:%3d  np:%3d):%14.6f%+14.6f\n",
         com.ntime, np, -lnL, -lnL+lnL0);
      OutaTreeB (fout);  FPN (fout);
      FOR (i,np) fprintf(fout," %8.5f",x[i]); FPN(fout); fflush(fout);

fprintf(frst1," %10.3f", -lnL);
for(i=com.ntime;i<com.np;i++) fprintf(frst1," %8.4f", x[i]);  
FPN(frst1); fflush(frst1);

/*
fprintf(frst," %10.3f%7.3f", -lnL,x[com.ntime]);

for(i=com.ntime;i<com.np;i++) fprintf(frst," %8.4f", x[i]);  
fprintf(frst,"%12.3f", -lnL); 
*/

      if (com.getSE) {
         Hessian (np, x, lnL, space, var, com.plfun, var+np*np);
         matinv (var, np, np, var+np*np);
         fprintf (fout,"SEs for parameters:\n");
         FOR(i,np) fprintf(fout," %8.5f",var[i*np+i]>0.?sqrt(var[i*np+i]):-0);
         FPN (fout);
         if (com.getSE==2) matout2(fout, var, np, np, 12,7);
      }
      if (com.clock) {
         SetBranch (x);
         for (i=0,tl=0;i<tree.nnode;i++) if(i-tree.origin) tl+=nodes[i].branch;
      }
      else   tl=sum(x, com.ntime);       
      fprintf (fout,"\ntree length %s = %9.5f\n", 
         ((com.ngene>1)?"(for the 1st gene)":""),tl);
      OutaTreeN(fout,1,1);   fputs(";\n", fout);
      DetailOutput (fout, x, var, var+np*np);
      if (com.seqtype==AAseq && com.model>=2 /* AAClasses */)
         EigenQaa(fout, Root, U, V, x+com.ntime+com.nrgene); /* S & PAM */

      if (com.NSsites)       lfunNSsites_rate(frst,x,np);
      else if (com.print) {
         if(com.rho==0)      AncestralSeqs(frst,x,space);
         if(com.plfun!=lfun) lfunRates(frate,x,np);
      }
      com.print-=9;
      com.plfun (x, np);
      com.print+=9;
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

void DetailOutput (FILE *fout, double x[], double var[], double space[])
{
/* var[] is used for codon models if com.getSE=1 to calculate the variances 
   of dS and dN.
*/
   int i,j,k=com.ntime, np=com.np;
   double om,S,N,dS,dN, vtw[4],vSN[4], tnode;

   if(com.clock>=1) { /* SetBranch() should have been called before this. */
      for(i=0; i<tree.nnode;i++) 
         printf("Node %3d divtime %9.5f\n",i+1,nodes[i].divtime);

      printf("\nreference node & node time? ");
      scanf("%d%lf", &j,&tnode);
      if(--j>=com.ns && tnode>0)
         for(i=com.ns, tnode/=nodes[j].divtime; i<tree.nnode;i++) 
            printf("Node %3d Time %6.1f\n",i+1,nodes[i].divtime*tnode);
   } 
   fprintf(fout,"\nDetailed output identifying parameters\n");
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
      if(com.aaDist) {
         k=com.ntime+com.nrgene+!com.fix_kappa;
         fprintf (fout,"\nb = %9.5f", x[k++]);
         if (com.seqtype==CODONseq)  fprintf (fout,"\na = %9.5f\n", x[k++]);
      }
   }
/*
   k=com.ntime+com.nrate;
*/
   if (!com.fix_alpha)  fprintf (fout, "alpha (gamma) = %8.5f\n", x[k++]);
   if (!com.fix_rho)   fprintf (fout, "rho (correlation) = %8.5f\n", x[k]);
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
         fprintf(fout," %8.5f",com.MK[i*com.ncatG+j]);
   }

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
   if((com.seqtype==CODONseq /*||com.model==FromCodon||com.model==AAClasses */)
      && com.NSsites==0) {
      fprintf(fout,"\ndN & dS for each branch\n\n%7s   ","branch");
      fprintf(fout,"%9s%9s%9s%9s%9s%9s\n\n", "t","S","N","dN/dS","dN","dS");
      XCOPY=x;
      FOR (i,tree.nbranch) {
         fprintf(fout,"%4d..%-3d ",tree.branches[i][0]+1,tree.branches[i][1]+1);
         j=com.ntime+com.nrgene+!com.fix_kappa;
         if (com.model==AAClasses) om=-1;
         else if (com.model==NS1 ||com.model==FromCodon)
                                   om=(com.fix_omega?com.omega:x[j]);
         else if (com.model==NSb)  om=x[j+i];
         else if (!com.fix_omega || OmegaBranch[i]<N_OMEGA-1)
                              om=x[j+=OmegaBranch[i]];
         else                 om=OMEGA_FIX;
         QcStatistics(x[i],&S,&dS,&dN, com.kappa,om,space); 
         N=com.ls*3-S;
                         /* om not used in AAClasses model */

         if (com.getSE && !com.clock && com.model!=AAClasses) {
            puts ("calculate the SEs for dN & dS for each branch");
            vtw[0]=var[i*np+i];  vtw[3]=var[j*np+j]; 
            vtw[1]=vtw[2]=var[i*np+j]; 
            VariancedSdN(x[i], com.omega, vtw, vSN);
            fprintf(fout,"dN = %7.5f +- %.5f   dS = %7.5f +- %.5f",
               dN,(vSN[3]>0?sqrt(vSN[3]):-0),dS,(vSN[0]>0?sqrt(vSN[0]):-0));
            fprintf(fout," (by method 2)\n");
         }
         fprintf(fout,"%9.3f%9.1f%9.1f%9.3f%9.4f%9.4f\n",x[i],S,N,dN/dS,dN,dS);
         fprintf(frst,"%8.1f%8.1f %9.3f%9.4f%9.4f",N,S,dN/dS,dN,dS);
      }
   }
   FPN (fout);
}



extern double SmallDiff;

int GetOptions (char *ctlf)
{
   int i, nopt=30, lline=255;
   char line[255], *pline, opt[99], comment='*';
   char *optstr[] = {"seqfile", "outfile", "treefile", "seqtype", "noisy", 
        "cleandata", 
        "runmode", "clock", "getSE", "RateAncestor", "CodonFreq", "verbose",
        "model","aaDist","aaRatefile",
        "NSsites", "icode", "Mgene", "fix_kappa", "kappa",
        "fix_omega", "omega", "fix_alpha", "alpha","Malpha", "ncatG", 
        "fix_rho", "rho", "ndata", "Small_Diff"};
   double t;
   FILE  *fctl=fopen (ctlf, "r");
   char *daafiles[]={"", "grantham.dat", "miyata.dat", 
                     "g1974c.dat","g1974p.dat","g1974v.dat","g1974a.dat"};

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
                  case ( 5): com.cleandata=(int)t;  break;
                  case ( 6): com.runmode=(int)t;     break;
                  case ( 7): com.clock=(int)t;       break;
                  case ( 8): com.getSE=(int)t;       break;
                  case ( 9): com.print=(int)t;       break;
                  case (10): com.codonf=(int)t;      break;
                  case (11): com.verbose=(int)t;     break;
                  case (12): com.model=(int)t;       break;
                  case (13): com.aaDist=(int)t;      break;
                  case (14): sscanf(pline+2,"%s",com.daafile); break;
                  case (15): com.NSsites=(int)t;     break;
                  case (16): com.icode=(int)t;       break;
                  case (17): com.Mgene=(int)t;       break;
                  case (18): com.fix_kappa=(int)t;   break;
                  case (19): com.kappa=t;            break;
                  case (20): com.fix_omega=(int)t;   break;
                  case (21): com.omega=t;            break;
                  case (22): com.fix_alpha=(int)t;   break;
                  case (23): com.alpha=t;            break;
                  case (24): com.nalpha=(int)t;      break;
                  case (25): com.ncatG=(int)t;       break;
                  case (26): com.fix_rho=(int)t;     break;
                  case (27): com.rho=t;              break;
                  case (28): com.ndata=(int)t;       break;
                  case (29): SmallDiff=t;            break;
               }
               break;
            }
         }
         if (i==nopt)
            { printf ("\noption %s in %s\n", opt, ctlf);  exit(-1); }
      }
      fclose (fctl);
   }
   else
      if (noisy) puts ("\nno ctl file..\n");

if (com.NSsites>2) error("NSsites>2");

   if (noisy) FPN(F0);
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
         if (com.model>NS1 && com.NSsites) 
            error("dN/dS ratios among branches & sites not implemented");
         if (com.model>NS1 && com.alpha) 
            error("dN/dS ratios among branches not implemented for gamma");
         if (com.model>NS1 || com.NSsites) {
            if(com.fix_kappa>1) error("fix_kappa>1, not implemented.");
            if(com.Mgene>1) error("model: Diff NS+Mgene");
         }
         if (com.fix_omega) {
            OMEGA_FIX=com.omega;
            if (com.model==NSb || com.NSsites==1) error("fix_omega?");
         }
         if (com.model>NS2) error ("seqtype or model.");
         if (com.kappa<0)  error("kappa..");
         if (com.Mgene>=3 && com.nrate==0)  error("Mgene");

         com.nrate=!com.fix_kappa+!com.fix_omega;
         if (com.aaDist)  com.nrate++;
         if(com.runmode==-2 && (com.NSsites||com.alpha||com.aaDist))
            error("err: models for pairwise comparison.");
         if (com.NSsites) {
            com.nrate=!com.fix_kappa;
            if(com.NSsites<=2)  com.ncatG=com.NSsites+1;
            else { /* deleterious (& positive selection) models */
               if(com.alpha<=0) { com.alpha=.5; puts("NSsites: alpha reset.");}
               com.ncatG=1+5+(com.NSsites%2==0);  /* w3>0 for even numbers */
               com.nrate+=(com.NSsites%2==0&&!com.fix_omega);
          }
         }
      }
   }
   else
      error ("seqtype..");
   if(com.runmode==-2) com.cleandata=1; /* for both codon & aa seqs */
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

if(com.NSsites) error("NSsites in testx.  Check this is needed");

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
   double tb[]={4e-6,9}, rgeneb[]={0.1,99}, rateb[]={1e-4,99};
   double alphab[]={0.005,99}, rhob[]={0.01,0.99}, omegab[]={.008,89}; 
   double pb[]={.0001,.9999};

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
   if (com.NSsites) {
      if(com.NSsites==1) { xb[k][0]=1e-4; xb[k++][1]=.9999;} /* prob for w0*/
      else if(com.NSsites==2)
         { xb[k][0]=xb[k+1][0]=-99; xb[k][1]=xb[k+1][1]=99; k+=2; }
     xb[k][0]=omegab[0];  xb[k][1]=omegab[1];
   }
   else if((com.seqtype==CODONseq||com.model==FromCodon)&&com.model!=AAClasses)
     { if(!com.fix_omega) { xb[k][0]=omegab[0]; xb[k][1]=omegab[1]; } }
   else
      FOR(j,N_OMEGA) { xb[k+j][0]=omegab[0]; xb[k+j][1]=omegab[1]; }

   if (com.aaDist<0 && (com.seqtype==1||com.model==FromCodon)) {
      /* linear relationship between d_ij and w_ij */
      if(com.nrate != !com.fix_kappa+1+(com.seqtype==1)) error("in Setxbound");
      xb[com.ntime+com.nrgene+!com.fix_kappa][1]=1; /* 0<b<1 */
   }

   k=com.ntime+com.nrgene+com.nrate;
   for (i=0;i<com.nalpha;i++,k++)  FOR (j,2) xb[k][j]=alphab[j];
   if (!com.fix_rho)   FOR (j,2) xb[np-1][j]=rhob[j];

   if(noisy>2) {
      puts("\nBounds:\n");
      FOR(i,np) printf(" %9.6f", xb[i][0]);  FPN(F0);
      FOR(i,np) printf(" %9.6f", xb[i][1]);  FPN(F0);
   }

   return(0);
}


int SetxInitials (double x[])
{
/* This forces initial values into the boundary of the space
*/
   int i, k;
   double tb[]={.0002,3}, rgeneb[]={0.1,9}, rateb[]={.0001,89};
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

int GetInitials (double x[])
{
/* 
   perhaps try to restruct the code and make two sections for amino acids 
   and codons?
*/
   static int times=0;
   int i, j,k=0, nc=com.ncode, naa=20;
   int selection=(com.NSsites && com.NSsites%2==0);
   double t;

   if(_nScaleF) {
      zero(_nScaleF,com.npatt);  FOR(i,tree.nnode) _nScale[i]=0;
      SetnScale(tree.origin);
      if(noisy) { 
         puts ("\nNodes used for scaling: ");
         FOR(i,tree.nnode) if(_nScale[i]) printf(" %2d", i+1);
      }
   }
   com.plfun = (com.alpha==0 ? lfun : (com.rho==0?lfundG:lfunAdG));
   if (com.NSsites)  com.plfun=lfunNSsites; 

   if (times++==0) {
      if ((com.aaDist && (com.seqtype==CODONseq||com.model==FromCodon)) ||
          (com.seqtype==AAseq &&
          (com.model==Empirical||com.model==Empirical_F||com.model>=REVaa_0))){
         GetDaa(NULL, com.daa);
      }
   }

   if(com.clock==0)  com.ntime = tree.nbranch;
   else if (com.clock==1) com.ntime = tree.nnode-com.ns+(tree.origin<com.ns);
   else {  /* if(com.clock==2) */
      OutaTreeB(F0); FPN(F0);
      for (i=0,N_rateBranch=0; i<tree.nbranch; i++) {
            printf("Branch %2d: %2d..%-2d? ", i+1,
               tree.branches[i][0]+1, tree.branches[i][1]+1);
            scanf("%d", &rateBranch[i]);
            if (rateBranch[i]+1>N_rateBranch) N_rateBranch=rateBranch[i]+1;
         }
         printf ("\n\n%d rates for branches. ", N_rateBranch);
         com.ntime=tree.nnode-com.ns + N_rateBranch-1;
         if(com.ntime>tree.nbranch) 
            printf("\a\nntime=%d, too many rates??\n",com.ntime);
   }

   if (com.seqtype==CODONseq && (com.model==NSb||com.model==NS2)) {
      if (com.model==NS2) {
         printf("\nSpecify dN/dS ratios for all %d branches (0,1,2,etc.).\n",
                tree.nbranch);
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

   if (com.NSsites) {   /*  count p0,p1,p2, w */
      com.np++;         /* p0 */
      if (com.NSsites%2==0)  com.np+=1+(!com.fix_omega); /* plus p1 & w */
      if (com.NSsites>2)     com.np+=!com.fix_alpha;     /* plus alpha */
   }

   if (com.clock) {
      for(j=1,x[0]=1; j<tree.nnode-com.ns; j++) x[j]=0.8;
      for(; j<com.ntime; j++) x[j]=1;   /* rates for branches */
   }
   else {
      FOR (j,com.ntime) x[j]=.1;
      LSDistance (&t, x, testx);
   }

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
            if(j,com.nrate-!com.fix_kappa)  x[k++]=com.omega; 
         }
      }
      else {                                        /* CODONseq */
         if (com.nrate==0 && com.NSsites==0) 
            EigenQc(Root,U,V, com.kappa,com.omega,PMat);
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
            if(com.NSsites) {
               k=com.ntime+com.nrgene+!com.fix_kappa;
               if (com.NSsites==1) x[k++]=.7; /* p for w=0 */
               else {                         /* com.NSsites>=2 */
                  x[k++]=.7; x[k++]=.1;       /* p0 & p1 */
                  if(selection && !com.fix_omega) x[k++]=com.omega; /*w2*/
                  if(com.NSsites>2&&!com.fix_alpha) {
                     x[k++]=com.alpha;
                     if (com.fix_alpha) {
                        DiscreteGamma(com.freqK+1,com.rK+1,
                              com.alpha,com.alpha, com.ncatG-1-selection ,0);
                        for(j=1,com.rK[0]=0;j<com.ncatG-1-selection;j++) 
                           com.rK[j]*=com.omega;
                 }
              }
               }
          }
         }
      }
   }
   if (!com.fix_alpha)
      x[com.np-1-!com.fix_rho]=com.alpha;
   if (!com.fix_rho) x[com.np-1]=com.rho;
   if (com.rho)
      AutodGamma (com.MK, com.freqK, com.rK, &t, com.alpha, com.rho,com.ncatG);
   else if (com.alpha && com.fix_alpha && com.NSsites<=2)
      DiscreteGamma(com.freqK,com.rK,com.alpha,com.alpha,com.ncatG,0);

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
         if (!com.NSsites) 
            EigenQc(Root,U,V, com.kappa,com.omega,PMat);
      else 
         EigenQaa(NULL, Root, U, V, x+k);
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
   int K=com.ncatG, i,k, naa=20, selection=(com.NSsites>=2&&com.NSsites%2==0);
   double t;

   XCOPY=x;
   SetBranch (x);
   FOR(i,com.nrgene) com.rgene[i+1]=x[com.ntime+i];

   k=com.ntime+com.nrgene;
   if (com.nrate) {
      if (!com.fix_kappa) com.kappa=x[k]; 
      if (com.model!=AAClasses && !com.fix_omega) com.omega=x[k];
      if (com.seqtype==AAseq)
         EigenQaa(NULL, Root, U, V, x+com.ntime+com.nrgene);
      else if ((com.model==0 && com.NSsites==0) || com.model==AAClasses)
         /* CODONs, same N/S for branch & site */ 
         EigenQc(Root,U,V, com.kappa,com.omega,PMat);
      k+=com.nrate;
   }

   if(com.NSsites) {
      k=com.ntime+com.nrgene+!com.fix_kappa;
      if (com.NSsites>2 && !com.fix_alpha) {   /* continuous */
         com.alpha=x[com.ntime+com.nrgene+com.nrate+1+selection];
         DiscreteGamma(com.freqK+1,com.rK+1, com.alpha, com.alpha,K-2,0); 
         error("need to scale rK, not ready yet");
      }

      com.rK[0]=0; com.rK[1]=1; 
      if(selection) com.rK[K-1]=(com.fix_omega?OMEGA_FIX:x[k+2]);

      if (com.NSsites%2==1) {   /* p0 for w0=0, p1=1-p0 */
         com.freqK[0]=x[k];   com.freqK[1]=t=1-x[k];
         com.rK[0]=0;         com.rK[1]=1;   /* omega's */
         if(com.NSsites>2) {
            FOR(i,K-1) com.freqK[1+i]=t/(K-1);
         }
      }
      else if(selection) { /* selection: p0, p1 & p2 */
         if(LASTROUND) {
            com.freqK[0]=x[k++]; com.freqK[1]=t=x[k++];
            com.freqK[K-1]=1-com.freqK[0]-com.freqK[1];
            for(i=1; i<K-1;i++)  com.freqK[1]=t/(K-2);
         }
         else {
            t=(com.freqK[0]=exp(x[k+0]))+(com.freqK[1]=exp(x[k+1]))+1;
            com.freqK[K-1]=1/t;  com.freqK[0]/=t;  t=com.freqK[1]/t;
            for(i=1; i<K-1;i++)  com.freqK[i]=t/(K-2);
            k+=2;          /* p0 & p1 */
         }
      }
   }     /* if(NSsites) */

   if (!com.fix_alpha && com.NSsites==0) {
      com.alpha=x[k++];
      if (com.fix_rho) {
         DiscreteGamma (com.freqK, com.rK, com.alpha, com.alpha, com.ncatG, 0);
         if (com.NSsites==3)  FOR (i,com.ncatG) com.rK[i]*=com.omega;
      }
   }
   if (!com.fix_rho) {
      com.rho=x[k++];
      AutodGamma(com.MK, com.freqK, com.rK, &t, com.alpha, com.rho, com.ncatG);
   }
   return (0);
}

int QcStatistics (double branchl, double *S, double *dS, double *dN,
    double kappa, double omega, double Q[64*64])
{
/* This calculates the number of synonymous sites (S) and
   the numbers of synonymous and nonsynonymous substitutions (dS & dN).
   c[0-2] are rates for the 3 codon positions, useful for calculating ratios
   like c[1]/c[0].
   k3[0,1,2] are the kappa's for the three codon positions expected under
   the codon model.
*/
   int n=Nsensecodon[com.icode], i,j,k, ic1,ic2,ndiff,pos=0,from[3],to[3];
   double rs0,ra0,rs,ra, c0[3],c[3],ts[3],tv[3], t,space[64*3], *ri=space;
   double *pi=(com.seqtype==AAseq?com.fb61:com.pi);
   double *pomega=XCOPY+com.ntime+com.nrgene+!com.fix_kappa, w;

   if(branchl && (S==NULL||dS==NULL||dN==NULL)) error("QcStatistics");
   if (com.seqtype==AAseq) error("strange check");
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

      if ((from[pos]+to[pos]-1)*(from[pos]+to[pos]-5)==0) ts[pos]+=t*Q[i*n+j];
      else                                                tv[pos]+=t*Q[i*n+j];

      c[pos]+=t*Q[i*n+j];
      Q[j*n+i]=Q[i*n+j];
   }  /* for (i,j) */
/*
   t=rs+ra;
   if (com.NSsites==0)  FOR(i,n) FOR(j,n) Q[i*n+j]*=pi[j]/t;
   else                 FOR(i,n) FOR(j,n) Q[i*n+j]*=pi[j];
   FOR (i,n) Q[i*n+i]=-sum(Q+i*n,n);

  if(noisy>2) FOR(i,3) printf("\nts/tv[pos %d] = %.6f", i+1, ts[i]/tv[i]);
*/

   rs0=rs;
   t=rs0+ra0; rs0/=t;  ra0/=t;   *S=rs0*3*com.ls; 
   if (branchl) {  /* calculates dS & dN */
      t=rs+ra;   rs/=t;   ra/=t;
      *dS=branchl*rs/3/rs0;  *dN=branchl*ra/3/ra0;
   }
   else  { *dS=*dN=0; } 
   return (0);
}


int EigenQc (double Root[], double U[], double V[], 
    double kappa, double omega, double Q[64*64])
{
/* This contructs the rate matrix Q for codon substitution, using parameters 
   kappa and omega and com.pi (or com.fb61).
*/
   int n=Nsensecodon[com.icode], i,j,k, ic1,ic2,ndiff,pos=0,from[3],to[3];
   double t,space[64*3], *ri=space;
   double *pi=(com.seqtype==AAseq?com.fb61:com.pi);
   double *pomega=XCOPY+com.ntime+com.nrgene+!com.fix_kappa, w;

   FOR (i,n*n) Q[i]=0;
   for (i=0; i<n; i++) FOR (j,i) {  /* set up symmetrical part of Q first */
      ic1=FROM61[i]; from[0]=ic1/16; from[1]=(ic1/4)%4; from[2]=ic1%4;
      ic2=FROM61[j];   to[0]=ic2/16;   to[1]=(ic2/4)%4;   to[2]=ic2%4;
      for (k=0,ndiff=0; k<3; k++)  if (from[k]!=to[k]) { ndiff++; pos=k; }
      if (ndiff!=1)  continue; 
      Q[i*n+j]=1;
      if ((from[pos]+to[pos]-1)*(from[pos]+to[pos]-5)==0)  Q[i*n+j]=kappa;
      if((ic1=GenetCode[com.icode][ic1]-1)!=(ic2=GenetCode[com.icode][ic2]-1)){
         if (com.model==AAClasses) {
            if (ic1<ic2)  { k=ic2; ic2=ic1; ic1=k; }
            k=ic1*(ic1-1)/2+ic2;
            if (pomega[OmegaAA[k]]<0) {
               if (noisy)  printf ("ic1 & ic2 & iw & w: %d %d %d %.5f\n", 
                              ic1,ic2,OmegaAA[k],pomega[OmegaAA[k]]);
               pomega[OmegaAA[k]]=0;
            }
            if (com.seqtype==AAseq && com.nrate>65 && ic1*20+ic2==ijAAref)
                ;     /* if estimating grantham's matrix with aa sequences */
            else  Q[i*n+j]*=pomega[OmegaAA[k]];
         }
         else if (com.aaDist==0)  Q[i*n+j] *= omega;
         else {        /*  amino acid properties */
            /* w = pomega[0]*com.daa[ic1*20+ic2]*com.daa[ic1*20+ic2]; */
            w = pomega[0]*com.daa[ic1*20+ic2];
            if(com.aaDist>0)           Q[i*n+j] *= exp(-w);  /* geometric */
            else                       Q[i*n+j] *= 1-w;      /* linear */
            if (com.seqtype==CODONseq) Q[i*n+j] *= pomega[1];
/*
            Q[i*n+j]*=pomega[1]*(1-fabs(com.daa[ic1*20+ic2]-pomega[0]));
*/
         }
      }
      Q[j*n+i]=Q[i*n+j];
   } /* for(i) for(j) */
   if (com.seqtype==AAseq) return (0);
   FOR(i,n) FOR(j,n) Q[i*n+j]*=pi[j];
   for (i=0,t=0; i<n; i++) 
      { Q[i*n+i]=-sum(Q+i*n,n); t-=pi[i]*Q[i*n+i]; }
   if (com.NSsites==0)  FOR(i,n*n) Q[i]/=t;

   if (eigen (1,Q,n,Root,ri,U,V,space+n)) error("eigenQc err.");
   xtoy (U, V, n*n);
   matinv (V, n, n, space);

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
      EigenQc(Root,U,V, com.kappa,com.omega,PMat);
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
         fprintf (fout, "\n%-5s", getAAstr(aa3,i));
         FOR (j,i) fprintf (fout, " %4.0f", Q[i*naa+j]/t/com.pi[j]*100);
      }
      fputs("\n     ",fout);  FOR(i,naa) fprintf(fout,"%5s",getAAstr(aa3,i));
      FPN(fout);  fflush(fout);
   }
/*
   if (fout && frst) {
      fprintf(frst, "\nRate matrix (symmetrical part, Sij) for bubble plot\n");
      FOR (i,naa)  FOR (j,i) 
         fprintf(frst, "\t%d\t%d\t%.2f\n", i+1,j+1,Q[i*naa+j]/t/com.pi[j]*100);
   }
*/
   if (eigen(1,Q,naa,Root,space,U,V,space+com.ncode))  error("err: eigenQaa");
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
   if (com.cleandata && inode<com.ns) 
      for(h=pos0;h<pos1;h++) nodes[inode].lkl[h*n+com.z[inode][h]]=1;

   FOR (i, nodes[inode].nson) {
      ison=nodes[inode].sons[i];
      if (com.seqtype==CODONseq && (com.model==NSb||com.model==NS2)) {
         j=com.ntime+com.nrgene+!com.fix_kappa;
         om=XCOPY[j+OmegaBranch[k=nodes[ison].ibranch]];
         if (com.fix_omega && OmegaBranch[k]==N_OMEGA-1)  om=OMEGA_FIX;
         EigenQc(Root,U,V, com.kappa,om,PMat);
      }
      t=nodes[ison].branch*com.rgene[igene];
      if (com.seqtype==AAseq && com.model==Poisson)
         PMatJC69like (PMat, t, n);
      else 
         PMatUVRoot (PMat, t, n, U, V, Root);

      if (com.cleandata && nodes[ison].nson<1)       /* end node */
         for(h=pos0; h<pos1; h++) 
            FOR(j,n) nodes[inode].lkl[h*n+j]*=PMat[j*n+com.z[ison][h]];
      else
         for(h=pos0; h<pos1; h++) 
            FOR(j,n) {
               for(k=0,t=0; k<n; k++)    /* t is used as temp */
                  t+=PMat[j*n+k]*nodes[ison].lkl[h*n+k];
               nodes[inode].lkl[h*n+j]*=t;
            }
   }        /*  for (ison)  */

   if(_nScaleF && _nScale[inode])   /* scaling to avoid underflow */
      for(h=pos0; h<pos1; h++) {
         for(j=0,t=0;j<n;j++) 
            if(nodes[inode].lkl[h*n+j]>t) t=nodes[inode].lkl[h*n+j];
         for(j=0;j<n;j++)  nodes[inode].lkl[h*n+j]/=t;
         _nScaleF[h]+=log(t);
      }
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
   int i,n;

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
      FOR (is,ns)     {
       lt=0; 
       for(il=ig*lline,pz=z[is]+il; lt<lline && il<ls; il++,lt++,pz++) {
            b=(int) *pz;
            b=FROM61[b]; c[0]=b/16; c[1]=(b%16)/4; c[2]=b%4; c[3]=0;
            FOR(i,3) c[i]=BASEs[c[i]];
            if (is && simple)  {
               b=(int)z[0][il];
               b=FROM61[b]; c0[0]=b/16; c0[1]=(b%16)/4; c0[2]=b%4; c0[3]=0;
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
   double *fbsg=space, *fbs, fb3x4[3*4],fb3x4t[3*4],fb4[4],fb4t[4], t;
   char str[4]="";

   if(com.seqtype!=CODONseq) error("not codon seq.");
   PatternWeight (fout, space);
   for(j=0,zero(fbsg,com.ngene*com.ns*nc); j<com.ns; j++) {
      for(h=0,ig=0,lt=0; h<com.npatt; h++) {
         for(k=0;k<3;k++)  NucListall(com.z[j][h*3+k], &nb[k], ib[k]);
         if((k=nb[0]*nb[1]*nb[2])>1)  miss=1;   /* k = no. compatible codons */
         if(k>1) continue;
         ic=ib[0][0]*16+ib[1][0]*4+ib[2][0];
         fbsg[j*com.ngene*nc+ig*nc+ic]+=com.fpatt[h];
/*       ambiguity codons are ignored.
         FOR(i0,nb[0])  FOR(i1,nb[1])  FOR(i2,nb[2]) {
            ic=ib[0][i0]*16+ib[1][i1]*4+ib[2][i2];
            fbsg[j*com.ngene*nc+ig*nc+ic] += com.fpatt[h]/(double)k;
            if(GenetCode[com.icode][ic]==0) {
               printf("stop codon in seq %3d: ",j+1);
               for(k=0;k<3;k++)  printf("%c",com.z[j][h*3+k]);  FPN(F0);
               status=-1; 
            }
         }
*/
         if ((lt+=(int)com.fpatt[h])==com.lgene[ig]) ig++;
      }
   }

   if((fbs=(double*)malloc((com.ns+1)*nc*sizeof(double)))==NULL) error("fbs?");
   zero(fbs,(com.ns+1)*nc);
   FOR(j,com.ns) FOR(ig,com.ngene) 
      FOR(k,nc) fbs[j*nc+k]+=fbsg[j*com.ngene*nc+ig*nc+k];
   FOR(j,com.ns) FOR(k,nc) fbs[com.ns*nc+k]+=fbs[j*nc+k];
   fputs("\nCodon usage in species and their sums\n",fout);
   printcums(fout, com.ns+1, fbs, com.icode);
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
   else if (com.codonf==Fcodon) {
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
   free(fbs);
   if(miss) fputs("\n(Ambiguity codons are not used to calculate freqs.)\n",fout);

   if(com.cleandata) {
      FOR(j,com.ns) {
         if(transform(com.z[j],com.npatt*3,1,0)) error("strange??");
         FOR(h,com.npatt) {
            b[0]=com.z[j][h*3]; b[1]=com.z[j][h*3+1]; b[2]=com.z[j][h*3+2];
            if(FROM64[k=b[0]*16+b[1]*4+b[2]]==-1) {
               printf("\nstop codon %s in seq. %d\n", getcodon(str,k),j+1);
               exit(-1);
            }
           com.z[j][h]=FROM64[k];
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
   FOR(ic,64) if((iaa=GenetCode[com.icode][ic])>0) NCsyn[iaa-1]++;
   zero (fcodon, 64);
   for (ic=0; ic<Nsensecodon[com.icode]; ic++) {
      iaa=GenetCode[com.icode][FROM61[ic]]-1;
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
   This works with both coded (com.cleandata==TRUE) and uncoded data.
   In the latter case, the method does pairwise delection.
   alpha for gamma rates is used for dN only.
   Missing data (com.cleandata==FALSE) are not used (pairwise deletion).
*/
   FILE *fds,*fdn;
   char codon[2][3], dsf[32]="DistanceNG.dS",dnf[32]="DistanceNG.dN";
   int is,js, k,i0,h, wname=20, status=0, ndiff,nsd[4];
   int nb[3],ib[3][4], missing;
   double ns, na, nst, nat, S,N, St,Nt, dS,dN,dN_dS,y, bigD=3;

   fprintf(fout,"\nNei & Gojobori 1986. dN/dS (dN, dS)");
   fprintf(fout,"\n(This matrix is not used in later m.l. analysis.)\n");

   fputs("\nNumber of codon sites with 0,1,2,3 position differences\n",frst);

   fds=(FILE*)fopen(dsf,"w");  fdn=(FILE*)fopen(dnf,"w"); 
   if(fds==NULL || fdn==NULL) error("err DistanceMatNG86: file error");
   fprintf(fds,"%6d\n", com.ns);  fprintf(fdn,"%6d\n", com.ns);

   FOR (is,com.ns) {
      fprintf(fout,"\n%-*s", wname,com.spname[is]);
      fprintf(fds,   "%-*s ",wname,com.spname[is]);
      fprintf(fdn,   "%-*s ",wname,com.spname[is]);
      FOR (js,is) {
         FOR(k,4) nsd[k]=0;
         for (h=0,nst=nat=S=N=0; h<com.npatt; h++)  {

		      if(com.cleandata)
               FOR(i0,2) getcodon(codon[i0],FROM61[com.z[i0==0?is:js][h]]);
            else {
               FOR(i0,2) FOR(k,3) codon[i0][k]=com.z[i0==0?is:js][h*3+k];
               for(i0=0,missing=0;i0<2;i0++) {
                  FOR(k,3) NucListall(codon[i0][k], &nb[k], ib[k]);
                  if(nb[0]*nb[1]*nb[2]!=1)  missing=1;
               }
               if(missing) continue;
            }

            ndiff=difcodonNG(codon[0],codon[1],&St,&Nt,&ns,&na,0,com.icode);
            nsd[ndiff]+=(int)com.fpatt[h];

            S+=St*com.fpatt[h];
            N+=Nt*com.fpatt[h];
            nst+=ns*com.fpatt[h];
            nat+=na*com.fpatt[h];
         }

         y=com.ls*3./(S+N); S*=y; N*=y;       /* rescale for stop codons */
         fprintf (frst, "%4d vs. %4d %6d %5d %5d %5d  ", 
            is+1,js+1,nsd[0],nsd[1],nsd[2],nsd[3]);
         dS=1-4./3*nst/S;  dN=1-4./3*nat/N;
         if(noisy>=9 && (dS<=0||dN<=0)) { puts("\nNG86 distance."); status=-1;}
         if(dS==1) dS=0;
         else      dS=(dS<0?-1:3./4*(-log(dS)));
         if(dN==1) dN=0;
         else dN=(dN<0?-1:3./4*(alpha==0?-log(dN):alpha*(pow(dN,-1/alpha)-1)));

         dN_dS=(dS>0?dN/dS:-1);
         fprintf(fout,"%7.4f(%5.4f %5.4f)", dN_dS, dN, dS);
         fprintf(frst,"%7.4f(%6.4f %6.4f)\n", dN_dS, dN, dS);

         if(dN<0) dN=bigD; if(dS<0) dS=bigD;
         SeqDistance [is*(is-1)/2+js] = (S*dS+N*dN)*3/(S+N);
         fprintf(fds," %7.4f", dS);   fprintf(fdn," %7.4f", dN);
      }
      FPN(fds); FPN(fdn);
   }    /* for(is) */
   if(!com.cleandata) fputs("\n(Ambiguity data not used)",fout);
   FPN(fout);
   if(status) fprintf (fout, "NOTE: -1 means the method is inapplicable.\n");
   fputs("\n\n",frst); fflush (frst);
   fclose(fds);  fclose(fdn);
   return (0);
}


int GetDaa (FILE* fout, double daa[])
{
/* Get the amino acid distance (or substitution rate) matrix 
   (grantham, dayhoff, jones, etc).
*/
   FILE * fdaa;
   char aa3[4]="";
   int i,j,k, naa=20;
   double dmax=0,dmin=1e40;

   if ((fdaa=fopen(com.daafile, "r"))==NULL) 
      { printf("\nAA dist file %s not found.", com.daafile); exit(-1); }
   printf("\n\nReading matrix from %s..\n", com.daafile);
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
      FOR(i,naa) if(fscanf(fdaa,"%lf",&com.pi[i])!=1) error("err aaRatefile");
      if (fabs(1-sum(com.pi,20))>1e-4) {
         printf("\nSum of freq. = %.6f != 1 in aaRateFile\n",sum(com.pi,20)); 
         exit(-1);
      }
   }
   fclose (fdaa);

   if (fout) {
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
      printf("%d dN/dS ratios estimated from data.  Stop if wrong.",N_OMEGA);
      getchar ();
   }
   else {
      FOR (iomega, N_OMEGA-1) {
         fscanf(fin, "%d", &j);
         if (j!=iomega+1) { printf("err data file %s.", OmegaAAf); exit(-1); } 
         printf ("\nClass #%d: ", j);
         j=fgetc (fin);  if (j!=':') error("err expecting :");
         fgets (line, nline, fin);
      
         printf ("%s\n", line);
         for (j=0,npair=0; j<nline-1&&line[j]&&line[j]!='\n'; j++) {
            iaa=line[j];
            if (!isalpha(iaa)) continue;
            jaa=line[++j];  if(!isalpha(jaa)) error("err jaa");
            npair++;

            printf ("\npair %d: |%c%c| ", npair, iaa,jaa);
            iaa=CodeChara((char)iaa,AAseq); jaa=CodeChara((char)jaa,AAseq);
            if(iaa<0||iaa>19||jaa<0||jaa>19) error("aa not found");
            if (iaa<jaa)  { k=jaa, jaa=iaa; iaa=k; }
      
            printf ("|%c%c (%2d,%2d)| ", AAs[iaa], AAs[jaa],iaa,jaa);
            if (iaa==jaa) printf ("This pair has no effect.");
            if (OmegaAA[iaa*(iaa-1)/2+jaa]) error("This pair specified?");
            if (AA1STEP[iaa*(iaa-1)/2+jaa]==0) error("This pair has rate 0!");
            OmegaAA[iaa*(iaa-1)/2+jaa]=iomega+1;
            printf (" in class %d ", iomega+1);
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
      error ("TestModelQc: ns>5 || alpha>0.");
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
         for (j=0,y=com.pi[nodeb[tree.origin]]; j<tree.nbranch; j++) 
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


int GetFreqPairCodon(void)
{
/* Get the expected codon frequencies (com.pi[]) for a pair of encoded codon seqs
   com.z[0,1], using com.codonf.  Data are clean and coded.
*/
   int n=com.ncode, i,j, ic,b[3];
   double *pi=com.pi, fb3x4[12], fb4[4], y;

   if(!com.cleandata) error("clean data in GetFreqPairCodon?");

   if (com.codonf==0) { fillxc(pi,1./n,n); return 0; }   /* Fequal */
   for(i=0,zero(pi,n); i<com.npatt; i++)
      { pi[com.z[0][i]]+=(y=com.fpatt[i]/(2.*com.ls)); pi[com.z[1][i]]+=y; }
   if (com.codonf==3) return 0;                          /* Fcodon */

   for (i=0,zero(fb3x4,12),zero(fb4,4); i<n; i++) {
      ic=FROM61[i];  b[0]=ic/16; b[1]=(ic/4)%4; b[2]=ic%4;
      FOR (j,3) { fb3x4[j*4+b[j]]+=pi[i];  fb4[b[j]]+=pi[i]/3.; }
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
      EigenQc(Root,U,V, com.kappa,com.omega,PMat);

   FOR(k,n) expt[k]=exp(x[0]*Root[k]);
   for (h=0; h<com.npatt; h++) {
      if(com.fpatt[h]<1e-20) continue;
      for(k=0,fh=0;k<n;k++) fh+=U[com.z[0][h]*n+k]*expt[k]*V[k*n+com.z[1][h]];
      fh*=com.pi[com.z[0][h]];
      if(fh<=0) { printf("lfundSdN: fh = %.9f\n",fh); getchar(); }
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
   QcStatistics(t,&S,&dS,&dN,com.kappa,omega,PMat);

   eh=(t+1)*SmallDiff;
   QcStatistics(t+eh,&S,&dS1,&dN1,com.kappa,omega,PMat);
   QcStatistics(t-eh,&S,&dS2,&dN2,com.kappa,omega,PMat);
   JacobiSN[0*np+0]=(dS1-dS2)/(2*eh);
   JacobiSN[1*np+0]=(dN1-dN2)/(2*eh);
  
   eh=(omega+1)*SmallDiff;
   QcStatistics(t,&S,&dS1,&dN1,com.kappa,omega+eh,PMat);
   QcStatistics(t,&S,&dS2,&dN2,com.kappa,omega-eh,PMat);
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
   char *pz0[NS];   /* pz0, npatt0, & fpatt0 hold the old information */
   int npatt0=com.npatt;
   double *fpatt0;
   float fp[NCODE*NCODE];
   FILE *fds,*fdn;
   char dsf[32]="DistanceML.dS",dnf[32]="DistanceML.dN";
   int n=com.ncode, is,js,j,k,h, np, miss=0, wname=20;
   double x[3]={1,1,1}, xb[3][2], lnL, e=1e-6, *var=space+NP, S,dS,dN;
   double JacobiSN[2*3],T1[2*3],T2[2*3],vSN[2*2], dS1,dN1,dS2,dN2,y[3],eh; 
          /* for calculating SEs of dS & dN */
   double tb[2]={1e-5,9}, kappab[2]={.1,999}, omegab[2]={.001,99};

   fpatt0=(double*)malloc(npatt0*3*sizeof(double));
   FOR(k,com.ns) pz0[k]=com.z[k];
   com.z[0]=(char*)(fpatt0+npatt0);   com.z[1]=com.z[0]+npatt0;
   FOR (k,npatt0) fpatt0[k]=(float)com.fpatt[k];

   if(!com.cleandata) error("\nPairwiseCodon: clean data only");
   if (com.ngene>1 && com.Mgene==1) puts("ngene>1 to be tested.");
   if (noisy) printf("\n\npairwise comparison (Goldman & Yang 1994)..\n");
   fprintf(fout,"\npairwise comparison, codon frequencies: %s.\n",
      codonfreqs[com.codonf]);
   fds=(FILE*)fopen(dsf,"w");  fdn=(FILE*)fopen(dnf,"w"); 
   if(fds==NULL || fdn==NULL) error("err PairwiseCodon: file creation error");

   xb[0][0]=tb[0]; xb[0][1]=tb[1];
   if(!com.fix_kappa)  { xb[1][0]=kappab[0]; xb[1][1]=kappab[1]; }
   if(!com.fix_omega)  { k=1+!com.fix_kappa; xb[k][0]=omegab[0]; xb[k][1]=omegab[1]; }

   fprintf(fds,"%6d\n", com.ns);  fprintf(fdn,"%6d\n", com.ns);
   fprintf(frst, "\n\npairwise comparison (Goldman & Yang 1994)");
   fprintf(frst,
      "\nseq seq        N       S       dN      dS   dN/dS    Paras.\n");

   FOR(is,com.ns) {
      fprintf(fds,"%-*s ", wname,com.spname[is]);
      fprintf(fdn,"%-*s ", wname,com.spname[is]);
      FOR(js,is) {
         printf ("\n\n%4d vs. %3d", is+1, js+1);
         fprintf(fout,"\n\n%d (%s) ... %d (%s)",
              is+1,com.spname[is], js+1,com.spname[js]);
         fprintf (frst, "%3d %3d ", is+1, js+1);
         if(noisy>2) fprintf(frub, "\n\n%d (%s) ... %d (%s)",
                  is+1,com.spname[is], js+1,com.spname[js]);
         FOR(k,n*n) fp[k]=0;
         for(h=0; h<npatt0; h++) {
            j=max(pz0[is][h],pz0[js][h]); k=min(pz0[is][h],pz0[js][h]);
            fp[j*n+k]+=(float)fpatt0[h];
         }
         for(j=0,com.npatt=0;j<n;j++) FOR(k,j+1)  if(fp[j*n+k]) {
            com.z[0][com.npatt]=j;  com.z[1][com.npatt]=k; 
            com.fpatt[com.npatt++]=fp[j*n+k];
         }
         if(noisy>=2) printf("\n  npatt=%d ",com.npatt);
         GetFreqPairCodon();

         np=com.np=(com.ntime=1)+!com.fix_kappa+!com.fix_omega;  NFunCall=0;
         x[0]=min(SeqDistance[is*(is-1)/2+js], 1);
         xb[0][0]=max(xb[0][0],x[0]/10); xb[0][1]=min(xb[0][1],x[0]*10);
         
         if(noisy>=9) {
            FPN(F0);  FOR(k,np) printf(" %9.6f",xb[k][0]); FPN(F0);
                      FOR(k,np) printf(" %9.6f",xb[k][1]); FPN(F0);
         }

         if(!com.fix_kappa) if(x[1]<=xb[1][0]*2 || x[1]>=xb[1][1]/2) x[1]=2;
         if(!com.fix_omega) 
            { k=1+!com.fix_kappa; if(x[k]<=xb[k][0]*2 || x[k]>=xb[k][1]/2) x[k]=1; }
         if (x[0]) 
            ming2(noisy>2?frub:NULL,&lnL,lfun2dSdN,NULL,x,xb, space,e,np);
         else {  x[1]=x[2]=com.kappa=com.omega=0; lnL=0; }
         fprintf(fout,"\nlnL =%12.6f\n",-lnL);
         FOR(k,np) fprintf(fout," %8.5f",x[k]);  FPN(fout);

         if (x[0]&&com.getSE) {
            Hessian(np, x, lnL, space, var, lfun2dSdN, var+np*np);
            matinv(var, np, np, var+np*np);
            fprintf(fout,"SEs for parameters:\n");
            FOR(k,np) fprintf(fout," %8.5f",(var[k*np+k]>0.?sqrt(var[k*np+k]):-0));
            FPN(fout);
         }
         FPN(fout);
         QcStatistics(x[0],&S,&dS,&dN, com.kappa,com.omega,PMat);
         fprintf(fds," %7.4f", dS);   fprintf(fdn," %7.4f", dN);

         fprintf (fout,
             "t=%7.4f  S=%8.1f  N=%8.1f  dN/dS=%7.4f  dN=%7.4f  dS=%7.4f\n",
              x[0],S,com.ls*3-S,com.omega,dN,dS);

         fprintf(frst,"%8.1f%8.1f %8.4f%8.4f%8.4f",com.ls*3-S,S,dN,dS,com.omega);
         FOR(k,np) fprintf(frst," %8.4f",x[k]);  
         k=np-1;
         if (com.getSE)
            fprintf(frst," +-%6.4f",(var[k*np+k]>0.?sqrt(var[k*np+k]):-1));
         fprintf(frst," %9.3f\n",-lnL);
         if(com.getSE && !com.fix_omega) {
            FOR(k,np) {
               FOR(j,np) y[j]=x[j];
               y[k] += (eh=(x[k]+1)*SmallDiff);
               if(!com.fix_kappa) com.kappa=y[1];
               com.omega=y[1+!com.fix_kappa];
               QcStatistics(y[0],&S,&dS1,&dN1,com.kappa,com.omega,PMat);
               y[k] -= 2*eh;
               if(!com.fix_kappa) com.kappa=y[1];
               com.omega=y[1+!com.fix_kappa];
               QcStatistics(y[0],&S,&dS2,&dN2,com.kappa,com.omega,PMat);

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
      FPN(fds); FPN(fdn);   fflush(fds); fflush(fdn); 
   }     /* for (is) */
   fclose(fds); fclose(fdn);

   FOR(k,com.ns) com.z[k]=pz0[k];  
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
   char distf[32]="DistanceML.aa", *pz0[NS];
   int n=com.ncode, i, ns0=com.ns,is,js;
   double x, xb[2]={0,19}, lnL, e=1e-6, step;

   if (com.ngene>1 && com.Mgene==1) error("ngene>1 to be tested.");
   if (noisy) printf("\npairwise ML distances of AA seqs.\n");
   if((fdist=(FILE*)fopen(distf,"w"))==NULL) error("PairwiseAA: file error");

   if(com.model>Empirical_F)  error("PairwiseAA: model wrong");
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



double TreeScore(double x[], double space[])
{
   double xb[NP][2], e=1e-6, lnL;

   PointLklnodes ();
   puts("\aCheck TreeScore");

   com.ntime = com.clock ? tree.nnode-com.ns : tree.nbranch;  /* clock=2 */
   GetInitials (x);
   SetxBound (com.np, xb);
   if(!com.cleandata) InitPartialLikelihood ();

   ming2(NULL,&lnL,com.plfun,NULL,x,xb, space,e,com.np);
   return (lnL);
}


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
         EigenQc(Root,U,V, com.kappa,com.rK[ir],PMat);
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
/*
matout2(F0,x,1,np,9,5);
printf("lnL = %12.6f\n", -lnL);
*/
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
         EigenQc(Root,U,V, com.kappa,com.rK[ir],PMat);
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
         "Site","dN/dS","Prob");
   FOR (h, com.ls) {
      it=com.pose[h];
      fprintf(fout,"\n%4d%9.4f%9.4f",h+1,com.fhK[it],com.fhK[com.npatt+it]);
   }
   fprintf (fout,"\n\nlnL=%14.6f\n\n", -lnL);
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
 
    Some computation can be saved under com.model==NS2 if only 2 or a few more 
    ratios are used instead of all b ratios.  Check whether it is worth the 
    effort.

(*) Codon models for variable dN/dS ratios among sites
    NSsites = 0: one ratio for all sites
    NSsites = 1: neutral model, w0=0, w1=1
    NSsites = 2: positive-selection model, w0=0, w1=1, w2 estimated
    NSsites = 3 (odd): continuous neutral, w0=0, w1(from gamma)
    NSsites = 4 (even): continuous selection, w0=0, w1(from gamma); w2 est.ed

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

   N_OMEGA is the number of different dN/dS ratios in the NSb, NS2 models
      and in AAClasses. 
   OmegaBranch[ibranch]=0,1,...,(N_OMEGA-1) marks the dN/dS ratio for ibranch 
      in the NSb NS2 models
   AA1STEP[i*(i-1)/2+j] =1 if AAs i & j differ at one codon position;
                        =0 otherwise.

(*) Codon and amino acid models
    Geometric and linear relationships between amino acid distance and 
    substitution rate:
       wij = a*(1-b*dij/dmax)
       wij = a*exp(-b*dij/dmax)
    aaDist = 0:equal, +:geometric; -:linear, {1-5:Grantham,Miyata,c,p,v}
*/
