/* CODEML.c  (AAML.c & CODONML.c)
   Maximum likelihood parameter estimation for protein-coding DNA (codon) 
        sequences (seqtype=1) or amino-acid sequences (seqtype=2)
                 Copyright, Ziheng Yang (March 1993 onwards)

                 gcc -o codeml -O2 codeml.c tools.o eigen.o
                        codeml <ControlFileName>
*/
#include "tools.h"
#ifdef SMALLMACHINE
#define NS            30
#else
#define NS            50
#endif
#define NBRANCH       (NS*2-2)
#define NNODE         (NS*2-1)
#define NGENE         13
#define LSPNAME       30
#define NCODE         64
#define NCATG         8
#define NP            (NBRANCH+NGENE-1+189+2)

extern char NUCs[], AAs[], GenetCode[][64], Nsensecodon[];
extern int noisy, NFunCall, *ancestor;
extern double *SeqDistance;

int Forestry (FILE *fout, double space[]);
int GetOptions (char *ctlf);
int testx (double x[], int np);
int SetxBound (int np, double xb[][2]);
int GetInitials (double x[], double space[]);
int SetParameters (double x[]);
int SetPGene (int igene, int _pi, int _UVRoot, int _alfa, double x[]);
double lfunAdG (double x[], int np);
double lfundG (double x[], int np);
double lfun (double x[], int np);
int PartialLikelihood (int inode, int igene);
int PMatJC69like (double P[], double t, int n);
int setmark_61_64 (void);
int printfcode (FILE *fout, double fb61[], double space[]);
int InitializeCoding (FILE *fout, double space[]);
int DistanceMatCode (FILE *fout, double alfa);
int DistanceMatNG1986 (FILE *fout, double alfa);
int GetDaa (FILE *fout, char daafile[], double daa[]);
int EigenQ20 (FILE *fout, double Root[], double U[], double V[],double rate[]);
int EigenQ61 (FILE *fout, double Root[], double U[], double V[], 
    double kapa, double omega, double Q[64*64]);
int Qcodon2aa (double Qc[], double pic[], double Qaa[], double piaa[]);
int Setmark1Step (void);
int TestModelQ61 (FILE *fout, double x[]);
int Pairwise (FILE *fout, double space[]);

struct common_info {
   char *z[NS], spname[NS][LSPNAME+1], seqf[32], outf[32];
   int seqtype, ns, ls, ngene, posG[NGENE+1], lgene[NGENE], npatt,*fpatt,*pose;
   int runmode, clock, print, model, icode, ncode, Mgene;
   int given_rgene, given_kapa, given_alfa, given_rho, getSE;
   int np, ntime, nrgene, nrate, nalfa, ncatG, maxnp;
   double lmax, pi[NCODE], fb61[64], piG[NGENE][64];
   double kapa, omega, alfa, rho, rgene[NGENE];
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

int from61[64], from64[64], mark1Step[190];
double PMat[NCODE*NCODE], U[NCODE*NCODE], V[NCODE*NCODE], Root[NCODE], TIME2S;
FILE *frub, *flfh, *frst;
enum {Fequal, F1x4, F3x4, Fcodon} CodonModels;
enum {Poisson, EqualInput, Dayhoff, Jones, Dayhoff_Pi, Jones_Pi,
     FromCodon=6, NGrantham, REVaa0=8, REVaa=9} AAModel;
char *codonmodels[]={"Fequal", "F1x4", "F3x4", "Fcodon","","","","NGrantham"};
char *aamodels[]={"Poisson", "EqualInput", "Dayhoff", "Jones", "Dayhoff_Pi",
     "Jones_Pi", "FromCodon", "NGrantham", "REVaa_0", "REVaa"};

#define CODEML
#include "treesub.c"

int main(int argc, char *argv[])
{
   FILE *fout, *fseq;
   char ctlf[32]="codeml.ctl", *pmodel;
   char *seqtypestr[3]={"CODONML", "AAML", "CODON2AAML"};
   char *Mgenestr[]={"diff. rate", "separate data", "diff. rate & pi", 
                     "diff. rate & kapa", "diff. rate & pi & kapa"};

   int i, s2=0, s3=0;
   double *space;

   noisy=2;
   com.runmode=0;
   com.clock=0;       /* 1: clock, rooted tree;  0: no clock, unrooted tree  */
   com.given_rgene=0; /* 0: estimate rate factors for genes */

   com.seqtype=AAseq;
   com.model=Jones;
   com.icode=0;
   com.nrate=0;
   com.given_kapa=0;   com.kapa=1;     com.omega=2.1;
   com.given_alfa=1;   com.alfa=0.;   com.ncatG=4;   /* alfa=0 := inf */
   com.given_rho=0;    com.rho=0.4; 

   com.getSE=0;
   com.print=0;

   frub=fopen ("rub","w");   frst=fopen ("rst","w");
   flfh=fopen ("lfh", "w");

   if (argc>1) { strcpy(ctlf, argv[1]); printf ("\nctlfile reset to %s.\n", ctlf); }
   GetOptions (ctlf);
   if ((fout=fopen (com.outf, "w"))==NULL) error("outfile creation err.");
   if((fseq=fopen (com.seqf,"r"))==NULL)  error("No sequence file!");
   ReadSeq (NULL, fseq, 0);

   pmodel=(com.seqtype==CODONseq?codonmodels[com.model]:aamodels[com.model]);
   fprintf (fout, "%s  %s   Model: %s (nrate:%d)  ",
      seqtypestr[com.seqtype-1], com.seqf, pmodel, com.nrate);
   if (com.clock) fprintf (fout, " Clock  ");
   if (com.alfa && com.rho) fprintf (fout, "Auto-");
   if (com.alfa) fprintf (fout,"dGamma (ncatG=%d) ", com.ncatG);
   if (com.ngene>1) 
      fprintf (fout, " (%d genes: %s)  ", com.ngene, Mgenestr[com.Mgene]);
   if (com.alfa==0) com.nalfa=0;
   else             com.nalfa=(com.nalfa?com.ngene:!com.given_alfa);
   if (com.nalfa>1 && com.rho) error ("Malfa or rho");
   if (com.nalfa>1) fprintf (fout,"(%d gamma)", com.nalfa);
   if (com.Mgene>=3) com.nrate*=com.ngene;
   if (com.Mgene && com.ngene==1) error ("Mgene");

   com.maxnp = max((2*com.ns-2)+com.nrate+com.ngene-1+2, 50);
   s3 = (unsigned)com.maxnp*(com.maxnp*2+2+12)*sizeof(double);
   s3 = max (s3, (unsigned)com.ls*(com.ns*sizeof(char)+sizeof(int)));
   i = com.ns*(com.ns-1)/2;
   s3 = max (s3, sizeof(double)*((com.ns*2-2)*(com.ns*2-2 + 4 + i) + i));
   if ((space=(double*)malloc(s3))==NULL)  error ("oom space");

   if (com.ns<50) {
      SeqDistance=(double*)malloc(com.ns*(com.ns-1)/2*sizeof(double));
      ancestor=(int*)malloc(com.ns*(com.ns-1)/2*sizeof(int));
      if (SeqDistance==NULL||ancestor==NULL) error("oom");
   }

   if(com.seqtype==AAseq) Initialize (fout, space, com.seqtype);
   else                   InitializeCoding (fout, space); /* seqtype change!*/
   if (com.seqtype==CODONseq) DistanceMatNG1986 (fout, 0);
   else                       DistanceMatCode (fout, com.alfa);

   if (com.Mgene==3)  FOR (i,com.ngene) xtoy(com.pi, com.piG[i], com.ncode);

   if (com.seqtype==AAseq && com.model==Poisson && !com.print) 
      PatternJC69like (fout);

   if (com.alfa)
     if((com.fhK=(double*)malloc(s2=com.npatt*com.ncatG*sizeof(double)))==NULL)
         error ("oom");

   if (com.runmode==-2)
      if (com.seqtype==CODONseq) { Pairwise (fout, space);  exit(0); }
      else error ("runmode");

   printf ("\n%10d bytes for lkl\n",
           (com.ns-1)*com.ncode*com.npatt*sizeof(double));
   com.chunk = 
      (double*) malloc ((com.ns-1)*com.ncode*com.npatt*sizeof(double));
   if (com.chunk==NULL) error ("oom chunk");
   printf ("%10d bytes for fhK\n%10d bytes for space\n", s2, s3);


   if (com.Mgene==1)         MultipleGenes (fout, space);
   else if (com.runmode==0)  Forestry (fout, space);
   else                      TreeSearch (fout, space);
/*
   for (i=0; i<100; i++) {
      com.nrate=0;
      com.kapa=0.5*pow(1.5,(double)(i/10));
      com.omega=10.0*(1+i%10);
      Forestry (fout, space);
      fprintf (flfh,"\n%4d%8.3f%8.3f", i+1, com.kapa, com.omega);
   }
*/
   return (0);
}

int Forestry (FILE *fout, double space[])
{
   char treef[12]="trees. s";
   int  i, j, itree, ntree, np, nchange;
   double x[NP],xcom[NP-NBRANCH],xb[NP][2], lnL,lnL0=0,lnLm=0, e=1e-7,tl=0;
   double *var=space+NP;
   FILE *ftree;

   sprintf (treef, "trees.%ds", com.ns);
   if ((ftree=fopen (treef,"r"))==NULL)  error ("no tree file.");
   fscanf (ftree, "%d%d", &i, &ntree);
   if (i!=com.ns) error ("err:ns"); 
   if (ntree>10 && com.print) puts("\nlarge lfh file");
   fprintf (flfh, "%6d%6d%6d\n", ntree, com.ls, com.npatt);
   FOR (itree, ntree) {
      printf ("\nTREE # %2d\n", itree+1);
      fprintf (fout,"\nTREE # %2d:  ", itree+1);
      fprintf (flfh,"\n\n%2d\n", itree+1);
      fprintf (frub,"\n\nTREE #%2d", itree+1);

      if (ReadaTreeN(ftree, &i, 1)) error ("err tree..");

      nchange=MPScore (space);
      OutaTreeN (F0, 0, 0);     printf ("   MP score: %d", nchange);
      OutaTreeN (fout, 0, 0);   fprintf (fout, "   MP score: %d", nchange);

      fflush (fout),  fflush (flfh);
      GetInitials (x, space);
      np=com.np;
      PointLklnodes ();
      NFunCall=0;
printf ("\nntime & nrate & np: %d  %d  %d ", com.ntime, com.nrate, com.np);
/*
      fprintf (fout,"\n\nLS branch lengths\n");
      FOR (i, com.ntime) fprintf (fout,"%9.5f", x[i]);
*/
      if (itree)  for (i=0; i<np-com.ntime; i++) x[com.ntime+i]=xcom[i];
/*
      printf ("\nInitials (np =%3d)\n", np);
      FOR (i,np) scanf ("%lf", &x[i]);
      FOR (i,np) x[i]=max(x[i], 0.00001);
*/
      if (noisy) {
         printf ("\nnp%6d", np);	 
         if (noisy>2) matout (F0, x, 1, np); 
         lnL = com.plfun (x, np);
         printf ("\nlnL0 = %12.6f\n", -lnL);

      }
      if (com.clock==0 /* && com.model>=REVaa0 */ ) {
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

      fprintf (fout,"\nlnL(ntime:%3d  np:%3d):%14.6f%+14.6f\n",
         com.ntime, np, -lnL, -lnL+lnL0);
      OutaTreeB (fout);  FPN (fout);
      FOR (i, np) fprintf (fout,"%9.5f", x[i]);  FPN (fout);
      if (com.getSE) {
         Hessian (np, x, lnL, space, var, com.plfun, var+np*np);
         matinv (var, np, np, var+np*np);
         FOR (i, np)
            fprintf (fout,"%9.5f", var[i*np+i]>0.?sqrt(var[i*np+i]):0.);
         FPN (fout);
      }
      if (com.clock) {
         SetBranch (x, 0);
         for (i=0,tl=0;i<tree.nnode;i++) if(i-tree.origin) tl+=nodes[i].branch;
      }
      else   tl=sum(x, com.ntime);       
      fprintf (fout,"\ntree length: %9.5f\n", tl);
      OutaTreeN (fout, 1, 1);  FPN(fout);

      fprintf (fout, "\n# of lfun calls:%10d\n", NFunCall);

      if (com.plfun==lfun && com.print) {
         if (com.seqtype==AAseq) AncestralSeqs (frst, x, space);
      }
      else if (com.print)        lfunAdG_rate (frst, x, np);

      if(com.seqtype==CODONseq||com.model==FromCodon||com.model==NGrantham) { 
          EigenQ61 (fout, Root, U, V, com.kapa, com.omega, space);
         /* TestModelQ61 (fout, x);  */
      }
      if (com.seqtype==AAseq && com.model>=REVaa0)
         EigenQ20 (fout, Root, U, V, x+com.ntime+com.nrgene);
      else if (com.model==NGrantham) 
         for (i=0; i<20; i++,FPN(F0))
           for (j=0,printf("%4s",getaa(0,i+1)); j<i+1; j++)
               printf ("%5.0f", 100*com.daa[i*20+j]);

      com.print-=9;
      com.plfun (x, np);
      com.print+=9;

   }         /* for (itree) */
   fclose (ftree);
   if (ntree>2) fprintf (fout, "\nBest likelihood:%12.6f\n", -lnLm);
   return (0);
}

int GetOptions (char *ctlf)
{
   int i, nopt=20, lline=255;
   char line[255], *pline, opt[20], comment='*';
   char *optstr[] = {"seqfile", "outfile", "seqtype", "noisy", 
        "runmode", "clock", "getSE", "RateAncestor", "model","icode",
        "Mgene", "given_kapa", "kapa", "omega", "given_alfa", "alfa","Malfa",
        "ncatG", "given_rho", "rho"};
   double t;
   FILE  *fctl=fopen (ctlf, "r");

   if (fctl) {
      if (noisy) printf ("\n\nReading options from %s..\n", ctlf);
      for (;;) {
         if (fgets (line, lline, fctl) == NULL) break;
         for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
            if (isalnum(line[i]))  { t=1; break; }
            else if (line[i]==comment) break;
         if (t==0) continue;
         sscanf (line, "%s%*s%lf", opt, &t);
         if ((pline=strstr(line, "= "))==NULL) error ("option file.");

         for (i=0; i<nopt; i++) {
            if (strncmp(opt, optstr[i], 8)==0)  {
               if (noisy>2)
                  printf ("\n%3d %15s | %-20s %6.2f", i+1,optstr[i],opt,t);
               switch (i) {
                  case ( 0): sscanf(pline+2, "%s", com.seqf);    break;
                  case ( 1): sscanf(pline+2, "%s", com.outf);    break;
                  case ( 2): com.seqtype=(int)t;     break;
                  case ( 3): noisy=(int)t;           break;
                  case ( 4): com.runmode=(int)t;     break;
                  case ( 5): com.clock=(int)t;       break;
                  case ( 6): com.getSE=(int)t;       break;
                  case ( 7): com.print=(int)t;       break;
                  case ( 8): com.model=(int)t;       break;
                  case ( 9): com.icode=(int)t;       break;
                  case (10): com.Mgene=(int)t;       break;
                  case (11): com.given_kapa=(int)t;  break;
                  case (12): com.kapa=t;             break;
                  case (13): com.omega=t;            break;
                  case (14): com.given_alfa=(int)t;  break;
                  case (15): com.alfa=t;             break;
                  case (16): com.nalfa=(int)t;       break;
                  case (17): com.ncatG=(int)t;       break;
                  case (18): com.given_rho=(int)t;   break;
                  case (19): com.rho=t;              break;
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
   if (com.seqtype==AAseq || com.seqtype==CODON2AAseq) {
      com.ncode=20;
      switch (com.model) {
      case (Poisson):     case (EqualInput): case (Dayhoff): case (Jones):
      case (Dayhoff_Pi):  case (Jones_Pi):
         com.given_kapa=1; com.kapa=0; com.nrate=0;   break;
      case (FromCodon): 
         if (com.given_kapa==0)                           com.nrate=2;
         else if (com.given_kapa==1 || com.given_kapa==2) com.nrate=1;
         else if (com.given_kapa==3)                      com.nrate=0;
         else                                error ("err given_kapa");
         if (com.kapa<=0) error("kapa.."); 
         fillxc (com.fb61, 1./Nsensecodon[com.icode], Nsensecodon[com.icode]);
         break;
      case (NGrantham):
         Setmark1Step();  com.nrate+=2; 
         fillxc (com.fb61, 1./Nsensecodon[com.icode], Nsensecodon[com.icode]);
         break;
      case (REVaa0): com.given_kapa=0; com.kapa=0; Setmark1Step(); break; 
      case (REVaa):  com.given_kapa=0; com.kapa=0; com.nrate=189; break;
      }
      if ((com.Mgene>=2 && com.model==Fequal) 
        || (com.Mgene>=3 && com.model!=FromCodon))   error ("Mgene");
   }
   else if (com.seqtype==CODONseq) {
      com.ncode=Nsensecodon[com.icode];
      if (com.model==NGrantham) { Setmark1Step(); com.nrate+=2; }
      else {
         if (com.model>Fcodon) error ("seqtype or model.");
         com.nrate = (com.given_kapa<3)+(com.given_kapa==0);
      }
      if(com.kapa<=0) error("kapa..");
      if (com.Mgene>=3 && (com.model>Jones_Pi || com.nrate==0)) error("Mgene");
   }
   else  
      error ("seqtype..");
   if (com.given_alfa==1 && com.alfa==0) {
      if (com.rho) puts("rho set to 0.");  com.given_rho=1; com.rho=0; 
   }
   if(!com.given_alfa && com.alfa<=0) { com.alfa=0.5; puts("init alfa reset");}
   if(!com.given_rho && com.rho==0) { com.rho=0.001;  puts("init rho reset");}
   if(com.alfa)  if (com.ncatG<2 || com.ncatG>NCATG) error ("ncatG");

   if(com.model==NGrantham && com.given_kapa>=1) error ("NGrantham");

   return (0);
}

int testx (double x[], int np)
{
   int i,k;
   double tb[]={1e-5, 99}, rgeneb[]={0.1,20}, rateb[]={1e-5, 20}, omegab=1e-4;
   double alfab[]={0.005,10}, rhob[]={0.01,0.99};

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
   for (i=0; i<com.nalfa; i++,k++)
      if (x[k]<alfab[0] || x[k]>alfab[1]) return (-1);
   if (!com.given_rho) if (x[np-1]<rhob[0] || x[np-1]>rhob[1]) return (-1);
   return(0);
}

int SetxBound (int np, double xb[][2])
{
   int i, j, k;
   double tb[]={1e-5, 99}, rgeneb[]={0.1,99}, rateb[]={1e-5,20}, omegab=1e-4;
   double alfab[]={0.005,10}, rhob[]={0.01,0.99};

   FOR (i,com.ntime)  FOR (j,2) xb[i][j]=tb[j];
   FOR (i,com.nrgene) FOR (j,2) xb[com.ntime+i][j]=rgeneb[j]; 
   FOR (i,com.nrate)  FOR (j,2) xb[com.ntime+com.nrgene+i][j]=rateb[j];
   if (com.seqtype==CODONseq && com.nrate==2) 
      xb[com.ntime+com.nrgene+1][0]=omegab;
   k=com.ntime+com.nrgene+com.nrate;
   for (i=0;i<com.nalfa;i++,k++)  FOR (j,2) xb[k][j]=alfab[j];
   if (!com.given_rho)   FOR (j,2) xb[np-1][j]=rhob[j];
   return(0);
}

int GetInitials (double x[], double space[])
{
   static int times=0;
   int i, j, k=-1, n=com.ncode, naa=20, nr1=com.nrate;
   double t;
   char daafile[][32]={"grantham.dat", "dayhoff.dat", "jones.dat"};

   if (com.seqtype==CODONseq || com.model==FromCodon || com.model==NGrantham)
      k=0;
   else if (com.model==Dayhoff || com.model==Dayhoff_Pi)  k=1;
   else if (com.model>Dayhoff) k=2;
   if (times++==0) {
      if (k>=0) GetDaa (NULL, daafile[k], com.daa);
      if (k==0) for (i=0; i<20*20; i++) com.daa[i]*=.01;
   }
/*
if (k==0 && times==1) {
   for (i=0; i<20; i++)  { 
      for (j=0; j<i; j++)  com.daa[j*20+i]=com.daa[i*20+j]=1;
      com.daa[i*20+i]=0; 
   }
   printf ("\a\na.a. dij set to be equal\n");
}
*/
   com.plfun = (com.alfa==0 ? lfun : (com.rho==0?lfundG:lfunAdG));

   com.nrgene=(!com.given_rgene)*(com.ngene-1);
   if (!com.given_kapa && com.Mgene>=3) com.nrate*=com.ngene;

   com.np = com.ntime+com.nrgene+com.nrate+(!com.given_alfa)+(!com.given_rho);
   FOR (j, com.ntime)  x[j]=0.05+0.01*(com.ntime-j);

   if (com.ns<50) LSDistance (&t, x, testx);
/*
   FOR (j, com.ntime) x[j]=0.1;
*/
   FOR (j, com.nrgene)  x[com.ntime+j]=1;

   if (com.model==NGrantham) { /* last daa is not 1 to effect initial values */
      for (i=0,k=0; i<naa-1; i++)  for (j=i+1; j<naa; j++)  
         if (k<com.nrate-1 && mark1Step[k]==i*naa+j)  
            x[com.ntime+com.nrgene+k++]=com.daa[i*naa+j];
         else
            com.daa[i*naa+j]=com.daa[j*naa+i]=0;
      x[k=com.ntime+com.nrgene+com.nrate-2]=com.kapa; x[++k]=com.omega;
   }
   else {
     k=com.ntime+com.nrgene;
     if (com.seqtype==AAseq) {         /* AAseq */
        if (com.nrate==0)  EigenQ20 (NULL, Root, U, V, &t); /* once for all */
        if (com.model==REVaa0) {
           k=mark1Step[com.nrate];  t=com.daa[(k/n)*n+k%n];
           for (i=0,k=0; i<n-1; i++)  for (j=i+1; j<n&&k<com.nrate; j++)  
              if (mark1Step[k]==i*n+j)
                 x[com.ntime+com.nrgene+k++]=com.daa[i*n+j]/t;
        }
        else if (com.model==REVaa) { 
           for (i=0; i<n-1; i++)   for (j=i+1; j<n; j++)  
            if(i*n+j<18*n+19) x[k++]=max(0.02,com.daa[i*n+j]/com.daa[18*n+19]);
        }
        else if (com.model==FromCodon) {
           if (com.given_kapa==0)      { x[k++]=com.kapa;  x[k++]=com.omega; }
           else if (com.given_kapa==1)   x[k++]=com.omega;  /* kappa given */
           else if (com.given_kapa==2)   x[k++]=com.kapa;   /* omega given */
        }
     }
     else {                            /* CODONseq */
        if (com.nrate==0) 
           EigenQ61 (NULL, Root, U, V, com.kapa, com.omega, PMat);
        else {
           FOR (i, nr1*(com.Mgene>=3?com.ngene:1)) {
            if (com.given_kapa==1)   x[k+i*nr1]=com.omega;    /* kappa given */
            else if (com.given_kapa==2) x[k+i*nr1]=com.kapa;  /* omega given */
            else {  x[k+i*nr1]=com.kapa; x[k+1+i*nr1]=com.omega; }
	   }
        }
     }
   }
   if (!com.given_alfa) x[com.ntime+com.nrgene+com.nrate]=com.alfa;
   if (!com.given_rho) x[com.np-1]=com.rho;
   if (com.rho)
      AutodGamma (com.MK, com.freqK, com.rK, &t, com.alfa, com.rho, com.ncatG);
   else if (com.alfa)
      DiscreteGamma (com.freqK, com.rK, com.alfa, com.alfa, com.ncatG, 0);

   FOR (j, com.np) x[j]=max(x[j], 0.001);
   return (0);
}

int SetPGene (int igene, int _pi, int _UVRoot, int _alfa, double x[])
{
   int k, nr1=com.nrate/com.ngene;

   if (_pi) xtoy (com.piG[igene], com.pi, com.ncode);
   if (_UVRoot) {
      k=com.ntime+com.nrgene+(com.Mgene>=3)*igene*nr1;
      if (com.given_kapa==0 || com.given_kapa==2)  com.kapa=x[k++];
      if (com.given_kapa==0 || com.given_kapa==1)  com.omega=x[k++];
      if (com.seqtype==CODONseq)
         EigenQ61 (NULL, Root, U, V, com.kapa, com.omega, PMat);
      else 
         EigenQ20 (NULL, Root, U, V, x+k);
   }
   if (_alfa) {
      com.alfa=x[com.ntime+com.nrgene+com.nrate+igene];
      DiscreteGamma (com.freqK, com.rK, com.alfa, com.alfa, com.ncatG, 0);
   }
   return (0);
}

int SetParameters (double x[])
{
/* set com. variables and initialize U, V, Root etc,
*/
   int i,j, k, naa=20;
   double t;

   if (SetBranch (x, 1)) puts ("\nbranch len err..");
   FOR (i, com.nrgene) com.rgene[i+1]=x[com.ntime+i];

   if (com.model==NGrantham) {
      for (i=0,k=0; i<naa-1; i++)  for (j=i+1; j<naa&&k<com.nrate-2; j++)  
         if (mark1Step[k]==i*naa+j)
            com.daa[i*naa+j]=com.daa[j*naa+i]=x[com.ntime+com.nrgene+k++];
      com.kapa=x[k=com.ntime+com.nrgene+com.nrate-2]; com.omega=x[++k];
   }
   k=com.ntime+com.nrgene;
   if (com.nrate) {
      if (com.model!=NGrantham) {
         if (com.given_kapa==0 || com.given_kapa==2)  com.kapa=x[k++];
         if (com.given_kapa==0 || com.given_kapa==1)  com.omega=x[k++];
      }
      if (com.seqtype==AAseq)
         EigenQ20 (NULL, Root, U, V, x+com.ntime+com.nrgene);
      else 
         EigenQ61 (NULL, Root, U, V, com.kapa, com.omega, PMat);
   }
   if (!com.given_alfa) {
      com.alfa=x[k++];
      if (com.given_rho)
         DiscreteGamma (com.freqK, com.rK, com.alfa, com.alfa, com.ncatG, 0);
   }
   if (!com.given_rho) {
      com.rho=x[k++];
      AutodGamma (com.MK, com.freqK, com.rK, &t, com.alfa, com.rho, com.ncatG);
   }

   return (0);
}

int PartialLikelihood (int inode, int igene)
{
   int n=com.ncode, i,j,k,h, ison, pos0=com.posG[igene],pos1=com.posG[igene+1];
   double t;

   FOR (i,nodes[inode].nson)
      if (nodes[nodes[inode].sons[i]].nson>0)
         PartialLikelihood (nodes[inode].sons[i], igene);
   fillxc (nodes[inode].lkl+pos0*n, (double)(inode>=com.ns), (pos1-pos0)*n);
   if (inode<com.ns) 
      for (h=pos0; h<pos1; h++) nodes[inode].lkl[h*n+com.z[inode][h]-1]=1;
   FOR (i, nodes[inode].nson) {
      ison=nodes[inode].sons[i];
      if (com.seqtype==AAseq && com.model==Poisson)
         PMatJC69like (PMat, nodes[ison].branch*com.rgene[igene], n);
      else 
         PMatUVRoot (PMat, nodes[ison].branch*com.rgene[igene], n, U, V, Root);
      if (nodes[ison].nson<1)  /* end node */
         for (h=pos0; h<pos1; h++) 
            FOR (j,n) nodes[inode].lkl[h*n+j]*=PMat[j*n+com.z[ison][h]-1];
      else
         for (h=pos0; h<pos1; h++) 
            FOR (j,n) {
               for (k=0,t=0; k<n; k++)
                  t += PMat[j*n+k] * nodes[ison].lkl[h*n+k];
                  nodes[inode].lkl[h*n+j] *= t;
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
   int i, n;

   for (i=0,n=0; i<64; i++) {
      if (GenetCode[com.icode][i]==0) 
         from64[i]=-1; 
      else 
         { from61[n]=i; from64[i]=n++; }
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
   FOR (i, n) space[from61[i]] = fb61[i]*10000;
   printcu (fout, space, com.icode, space+64);
   return (0);
}

int InitializeCoding (FILE *fout, double space[])
{
/* for genes, site patterns and fpatts 
*/
   int h, i,j, it, ic[3], nc=com.ncode, nconstp, ig, lt, wname=20;
   double t, fb3x4[3*4], fb64s[(NS>NGENE?(NS+1):(NGENE+1))*64];
   double fb4t[4], fb3x4t[12], fb[64], *fbg=fb64s;
   double *fpi[NGENE+1];

   if (com.ls%3) error ("Coding seq?");
   FOR (j, com.ns) 
      if (transform (com.z[j], com.ls, 1, 0)) printf("in seq #%d\n", j+1);

   fprintf (fout, "\n\nns =%4d\tls =%4d", com.ns, com.ls);
   FOR (j,com.ns) fprintf (fout,"\n%s", com.spname[j]);
/*
   FOR (j,com.ns) prints (fout, com.z[j], com.ls, 60, 3, 1);
   printsma (fout, com.z, com.ns, com.ls, 60, 3, 1, 1, CODONseq);
*/
   zero (fb3x4t, 3*4);   zero (fb64s, (com.ns+1)*64);   zero (com.fb61, 64);
   fprintf (fout,"\n\ncodon position x base (3x4) table for each sequence");
   for (j=0; j<com.ns; j++) {
      for (h=0,zero(fb3x4, 3*4); h<com.ls; h+=3) {
         FOR (i,3) {
            ic[i]=(char)(com.z[j][h+i]-1);
            fb3x4[i*4+ic[i]] += 3./(double)com.ls;
            fb3x4t[i*4+ic[i]] += 3./((double)com.ns*com.ls);
         }
         it=ic[0]*16+ic[1]*4+ic[2];
         if (GenetCode[com.icode][it]==0) {
            printf ("\n%-20s code %1c%1c%1c",
             com.spname[j],NUCs[(int)ic[0]],NUCs[(int)ic[1]],NUCs[(int)ic[2]]);
            error ("stop codon inside.");
         }
         fb64s[j*64+it] ++;
         com.z[j][h/3]=(char)(from64[it]+1); /* recode codons into 1,2,...61 */
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
   FOR (i, 64) com.fb61[from64[i]]=fb64s[com.ns*64+i]/(com.ls*com.ns);
   fprintf (fout, "\n\nCodon usage for each species and their sums\n");
   printcums (fout, com.ns+1, fb64s, com.icode);
/*
   fprintf (fout, "\nAverage sense codon frequencies, see the table below.");
   matout (fout, com.fb61, 16, 4);
   printfcode (fout, com.fb61, space);
*/
   if (com.seqtype==CODON2AAseq) {
      FOR (j,com.ns)  FOR (h,com.ls) 
         com.z[j][h]=GenetCode[com.icode][from61[com.z[j][h]-1]];
      com.seqtype=AAseq;
      fprintf (fout, "\nTranslating into amino acid sequences..\n");
      printsma (fout, com.z, com.ns, com.ls, 60, 10, 1, 0, AAseq);
   }
   PatternWeight (fout, space);

   t=1./((double)com.lgene[0]*com.ns);  
   zero (com.pi, nc);  zero (fb, nc);  zero (fbg, (com.ngene+1)*64);
   for (h=0,ig=0,lt=0,nconstp=0; h<com.npatt; ) {
      for (j=1; j<com.ns; j++)  if(com.z[j][h]!=com.z[0][h]) break;
      if (j==com.ns) nconstp += com.fpatt[h];
      FOR (j,com.ns) fb[com.z[j][h]-1]+=(double)com.fpatt[h]*t;
      if (com.seqtype==CODONseq) 
         FOR (j,com.ns) fbg[ig*64+from61[com.z[j][h]-1]]+=com.fpatt[h];
      else
         FOR (j,com.ns) fbg[ig*nc+com.z[j][h]-1]+=com.fpatt[h];

      lt+=com.fpatt[h++];
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
   /* edit com.pi & com.piG for the model */
   if (com.model==Fequal) {
      fillxc (com.pi, 1./nc, nc);
      FOR (j,com.ngene) fillxc (com.piG[j], 1./nc, nc);
   }
   if (com.seqtype==CODONseq && (com.model==F1x4 || com.model==F3x4)) {
      fpi[0]=com.pi;
      FOR (it, com.ngene) fpi[it+1]=com.piG[it];
      FOR (it, com.ngene+1) {
         zero (fb4t, 4); zero(fb3x4t, 12);         
         FOR (i, nc) {
            j=from61[i]; ic[0]=j/16; ic[1]=(j/4)%4; ic[2]=j%4; 
            FOR (j,3) fb4t[ic[j]]+=fpi[it][i];
            FOR (j,3) fb3x4t[j*4+ic[j]]+=fpi[it][i];
         }
         if (com.model==F1x4) 
            FOR (i,nc) {
               j=from61[i];
               fpi[it][i]=fb4t[j/16]*fb4t[(j/4)%4]*fb4t[j%4];
            }
         else if (com.model==F3x4) 
            FOR (i,nc) {
               j=from61[i];
               fpi[it][i]=fb3x4t[j/16]*fb3x4t[4+(j/4)%4]*fb3x4t[2*4+j%4];
            }
         abyx (1/sum(fpi[it],nc), fpi[it], nc);
      }
   }
   fprintf (fout,"\n\n# constant sites: %6d (%6.2f%%)",
      nconstp, (double)nconstp*100./com.ls);
   for (h=0,com.lmax=-(double)com.ls*log((double)com.ls); h<com.npatt; h++)
      if (com.fpatt[h]>1)
         com.lmax+=(double)com.fpatt[h]*log((double)com.fpatt[h]);
   fprintf (fout,"\nmax{ln(L)}:%16.6f\n", com.lmax);
   return(0);
}

int DistanceMatCode (FILE *fout, double alfa)
{
   int i,j, h;
   double p;

   if (fout)  fprintf(fout,"\ndistances (alpha set at %.2f)\n", alfa);
   FOR (i, com.ns) {
      if (fout) fprintf (fout, "\n%-15s", com.spname[i]);
      FOR (j, i) {
         for (h=0,p=0; h<com.npatt; h++)  
            if (com.z[i][h]!=com.z[j][h]) p+=com.fpatt[h];
         p = p/com.ls;
         SeqDistance[i*(i-1)/2+j]=p= (alfa==0 ? p : alfa*(pow(1-p,-1/alfa)-1));
         if (fout) fprintf (fout, "%7.3f", p);
      }
   }
   if (fout) FPN(fout);
   return (0);
}

int DistanceMatNG1986 (FILE *fout, double alfa)
{
/* Nei & Gojobori (1986)
*/
   int i,j, h, wname=15;
   char codon1[4], codon2[4];
   double ns, na, nst, nat, S, N, St, Nt, largeD=5;

   if (fout) 
      fprintf(fout,"\nNei & Gojobori 1986. Ks/Ka (Ks, Ka)  alpha=%.2f", alfa);
   FOR (i, com.ns) {
      if (fout)  fprintf (fout, "\n%-*s", wname, com.spname[i]);
      FOR (j, i) {
         for (h=0,nst=nat=St=Nt=0; h<com.npatt; h++)  {
            strcpy (codon1, getcode(from61[com.z[i][h]-1]));
            strcpy (codon2, getcode(from61[com.z[j][h]-1]));
            difcodon (codon1, codon2, &S, &N, &ns, &na, 0, com.icode);
            St+=S*com.fpatt[h];
            Nt+=N*com.fpatt[h];
            nst+=ns*com.fpatt[h];
            nat+=na*com.fpatt[h];
         }
         if (noisy>=2)
            printf ("\n%4d%4d:%9.2f +%9.2f =%9.2f ..%9.2f +%9.2f =%9.2f",
               i+1,j+1,St,Nt,St+Nt,nst,nat, nst+nat);

         S = 1-4./3*nst/St;
         N = 1-4./3*nat/Nt;
         if (S<=0 || N<=0) puts ("\nvery different..");
         nst=(S<0?largeD:3./4*(alfa==0?-log(S):alfa*(pow(S,-1/alfa)-1)));
         nat=(N<0?largeD:3./4*(alfa==0?-log(N):alfa*(pow(N,-1/alfa)-1)));

         SeqDistance [i*(i-1)/2+j] = (0.21*nst+0.79*nat)*3;

         if (fout) fprintf (fout, "%7.3f(%6.3f%6.3f)", nst/nat, nst, nat);
/*
         fprintf (frst, "%4d%4d%7.3f%6.3f%6.3f\n", i+1,j+1, nst/nat,nst,nat);
*/
      }
   }
   if (fout) FPN (fout);
   return (0);
}

int GetDaa (FILE* fout, char daafile[], double daa[])
{ 
   FILE * fdaa;
   int i,j, n=20;

   if ((fdaa=fopen(daafile, "r"))==NULL) error ("\ndaa file not found");
   for (i=0; i<n; i++)  for (j=0; j<i+1; j++)  {
      fscanf(fdaa, "%lf", &daa[i*n+j]);
      daa[j*n+i]=daa[i*n+j];
   }
   fclose (fdaa);
   if (fout) {
      fprintf (fout, "\n%s\n", daafile);
      FOR (i,n) {
         fprintf (fout, "\n%4s", getaa (0, i+1));
         FOR (j,i+1)  fprintf (fout, "%5.0f", daa[i*n+j]);
      }
      FPN (fout);
   }
   return (0);
}

int EigenQ20 (FILE *fout, double Root[], double U[], double V[], double rate[])
{
   int n=20, i,j,k;
   double Q[20*20], mr, space[20*3];

   /* construct Q[20*20] */
   switch (com.model) {
   case (Poisson) :  case (EqualInput) : 
      fillxc (Q, 1.0, n*n);  break;
   case (Dayhoff)  : case (Jones):
         FOR (i, n) com.pi[i]=com.daa[i*n+i];
   case (Dayhoff_Pi):  case (Jones_Pi):
      FOR (i, n)
         FOR (j, i+1) Q[i*n+j]=Q[j*n+i]=com.daa[i*n+j]/100;
      break;
   case (FromCodon): case (NGrantham): 
      EigenQ61 (NULL, Root, U, V, com.kapa, com.omega, PMat);
      Qcodon2aa (PMat, com.fb61, Q, space);
      break;
   case (REVaa0)  : 
      for (i=0,k=0,zero(Q,n*n); i<n-1; i++) for (j=i+1; j<n&&k<com.nrate; j++)
         if (mark1Step[k]==i*n+j)  Q[i*n+j]=Q[j*n+i]=rate[k++];
      k=mark1Step[com.nrate];
      Q[(k/n)*n+k%n]=Q[(k%n)*n+k/n]=1;
      break;
   case (REVaa)  : 
      for (i=0,k=0; i<n-1; i++) for (j=i+1; j<n; j++)
         if (i*n+j!=18*n+19) Q[i*n+j]=Q[j*n+i]=rate[k++];
      Q[18*n+19]=Q[19*n+18]=1; 
      break;
   }
   FOR (i,n) FOR (j,n) Q[i*n+j]*=com.pi[j];
   for (i=0,mr=0; i<n; i++) {
      Q[i*n+i]=0;  Q[i*n+i]=-sum(Q+i*n,n);  mr-=com.pi[i]*Q[i*n+i]; 
   }
   if (fout) {
      matout (fout, Root, 1, 20);
      matout (fout, com.pi, 1, 20);
      FOR (i,n) {
         fprintf (fout, "\n%-5s", getaa(0, i+1));
         FOR (j,i) fprintf (fout, "%4.0f", Q[i*n+j]/mr/com.pi[j]*100);
         fprintf (fout, "%4.0f", -Q[i*n+i]/mr*100);
      }
      fprintf(fout, "\n     ");  FOR(i,n) fprintf(fout,"%4s",getaa(0,i+1));
      fprintf(fout, "\n     ");  FOR(i,n) fprintf(fout,"%4s",getaa(1,i+1));
      FPN (fout);
   }
   if (eigen (1, Q, n, Root, space, U, V, space+com.ncode))
      error ("eigenQ20 err");
   xtoy (U, V, n*n);
   matinv (V, n, n, space);
   FOR (i,n)  Root[i]=Root[i]/mr;

   if (fout) {
      fprintf (fout, "\nPAM matrix");
      PMatUVRoot (PMat, 0.01, n, U, V, Root);
      FOR (i,n) {
         fprintf (fout, "\n%5s", getaa(0, i+1));
         FOR (j,n)  fprintf (fout, "%6.0f", PMat[i*n+j]*100000);
      }
      FPN (fout);
   }
   return (0);
}

int EigenQ61 (FILE* fout, double Root[], double U[], double V[], 
    double kapa, double omega, double Q[64*64])
{
   int n=Nsensecodon[com.icode], i,j,k, ic1,ic2, ndiff, pos=0, from[3],to[3];
   double rs0, ra0, rs, ra, c0[3], c[3], t, space[64*3];
   double *pi=(com.seqtype==AAseq?com.fb61:com.pi);

   FOR (i,3) c[i]=c0[i]=0;  zero(Q,n*n);
   for (i=0, rs0=ra0=rs=ra=0; i<n; i++) FOR (j,i) {
      ic1=from61[i]; from[0]=ic1/16; from[1]=(ic1/4)%4; from[2]=ic1%4;
      ic2=from61[j];   to[0]=ic2/16;   to[1]=(ic2/4)%4;   to[2]=ic2%4;
      for (k=0,ndiff=0; k<3; k++)  if (from[k]!=to[k]) { ndiff++; pos=k; }
      if (ndiff!=1)  continue; 
      t=2*pi[i]*pi[j];
      if ((from[pos]+to[pos]-1)*(from[pos]+to[pos]-5)==0)  Q[i*n+j]=kapa;
      else   Q[i*n+j]=1;
      c0[pos]+=t*Q[i*n+j];
      if((ic1=GenetCode[com.icode][ic1]-1)!=(ic2=GenetCode[com.icode][ic2]-1)){
         ra0+=t*Q[i*n+j];

if (com.daa[ic1*20+ic2]==0) { printf ("\n%4d%4d", ic1+1,ic2+1); error (".."); }

         Q[i*n+j] *= exp(-com.daa[ic1*20+ic2]*omega);
/*
         Q[i*n+j] *= exp(-sqrt(com.daa[ic1*20+ic2])*omega);
         Q[i*n+j] *= exp(-fabs(com.daa[ic1*20+ic2]-0.05)*omega);
         Q[i*n+j] *= 2.15/com.daa[ic1*20+ic2]*omega;
         Q[i*n+j] *= exp(-com.daa[ic1*20+ic2]*omega);
*/
         ra+=t*Q[i*n+j];
      }
      else 
         rs+=t*Q[i*n+j];
      c[pos]+=t*Q[i*n+j];
      Q[j*n+i]=Q[i*n+j];
   }
   if (com.seqtype==AAseq && !fout) return (0);

   FOR (i,n) FOR (j,n) Q[i*n+j]*=pi[j]/(rs+ra);
   FOR (i,n) { Q[i*n+i]=0; Q[i*n+i]=-sum(Q+i*n, n); }

   if (fout) {
      rs0=rs;  t=rs0+ra0; rs0/=t;  ra0/=t;
      t=rs+ra; rs/=t;     ra/=t;

      fprintf(fout,"\nrs/ra[0]=%.3f (%.3f  %.3f)\n     [*]=%.3f (%.3f  %.3f)",
           rs0/ra0, rs0, ra0, rs/ra, rs, ra);
      fprintf (fout, "\nKs/Ka (Ks/t, Ka/t) %.3f (%.3f  %.3f)\n",
           (rs/ra)/(rs0/ra0), rs/3/rs0, ra/3/ra0);
      if (com.ns==2)  {
         printf ("%10.5f", (rs/ra)/(rs0/ra0));
         fprintf (frst, "%8.4f %8.4f%8.4f",
            (rs/ra)/(rs0/ra0), rs/3/rs0*TIME2S, ra/3/ra0*TIME2S);
      }
/*
      fprintf(fout,"\nc1:2:3 [0] = 1 :%6.3f :%6.3f  [*] = 1 :%6.3f :%6.3f\n",
           c0[1]/c0[0], c0[2]/c0[0], c[1]/c[0], c[2]/c[0]);
*/
   }
   if (com.seqtype==AAseq) return (0);
   if (eigen (1, Q, n, Root, space, U, V, space+n))  error("eigenQ61 err.");
   xtoy (U, V, n*n);
   matinv (V, n, n, space);

   return (0);
}

int Qcodon2aa (double Qc[], double pic[], double Qaa[], double piaa[])
{
/* Q61 -> Q20, only the symmetrical components of the rate matrices are used.
   Q20(aai,aaj) = SUMi SUMj (piC[i]*piC[j]]*Q61[i][j]) / (piAAi*piAAj)
*/
   int i, j, aai, aaj, nc=Nsensecodon[com.icode], naa=20;
   double ti, tij;

   zero(piaa,naa);  zero(Qaa,naa*naa);
   FOR (i,nc) piaa[GenetCode[com.icode][from61[i]]-1] += pic[i];

   FOR (i,nc) {
      aai=GenetCode[com.icode][from61[i]]-1;
      ti=pic[i]/piaa[aai];
      FOR (j, i) {
         aaj=GenetCode[com.icode][from61[j]]-1;
         if (Qc[i*nc+j]==0 || aai==aaj) continue;
         tij=ti*pic[j]*Qc[i*nc+j]/piaa[aaj];
         Qaa[aai*naa+aaj] += tij;    Qaa[aaj*naa+aai] += tij; 
      }
   }
   return (0);
}

int Setmark1Step (void)
{
/* mark1step[k] marks the k_th aa pair that differs by one mutation, 
   Q[i*naa+j] is the k_th nonzero element if mark1Step[k]=i*naa+j;
*/
   int n=Nsensecodon[com.icode], naa=20, i,j,k, ic1,ic2, ndiff, from[3], to[3];
   int *Q=(int*)PMat;

   FOR (i, naa*naa) Q[i]=0;
   for (i=0; i<n; i++) FOR (j,i) {
      ic1=from61[i]; from[0]=ic1/16; from[1]=(ic1/4)%4; from[2]=ic1%4;
      ic2=from61[j];   to[0]=ic2/16;   to[1]=(ic2/4)%4;   to[2]=ic2%4;
      for (k=0,ndiff=0; k<3; k++)  if (from[k]!=to[k]) ndiff++; 
      if (ndiff!=1)  continue; 
      ic1=GenetCode[com.icode][ic1]-1;  ic2=GenetCode[com.icode][ic2]-1;
      Q[ic1*naa+ic2]=Q[ic2*naa+ic1]=1;
   }
   for (i=0,k=0; i<naa-1; i++) for (j=i+1; j<naa; j++)
      if (Q[i*naa+j]) mark1Step[k++]=i*naa+j;
   com.nrate=k-1;
   return(0);
}

#if 0

int GetCategoryQ61 (char z[NS])
{
/* the category ID for a codon site with z[NS], transformed
   classified into 19 categories 
*/
   int i,j, icat, ncat=19, it, b[NS][3], nbase[3], markb[4];

   for (j=0; j<com.ns; j++) {
      it=from61[(int)z[j]];  b[j][0]=it/16; b[j][1]=(it/4)%4; b[j][2]=it%4;
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

   if (com.ns>5 || com.alfa || com.ngene>1)
      error ("TestModelQ61: ns>5 || alfa>0.");
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

int Pairwise (FILE *fout, double space[])
{
/* calculate Ks & Ka for all pairwise codon sequence comparisons
   use different npatt for different pairs
*/
   char *pz0[NS];
   int  n=com.ncode, ii, jj,k,h, np;
   int  ns0=com.ns, npatt0=com.npatt, *fpatt0, fp[NCODE*NCODE];
   double x[NP], xb[NP][2], lnL, e=1e-6, *var=space+NP;

   if (com.ngene>1) error("ngene>1 not supported");
   if (noisy) printf("\npairwise comparison (Goldman & Yang 1994)..\n");
   fprintf(fout,"\npairwise comparison");
   fprintf(frst, "pairwise comparison by the method of Goldman & Yang (1994)");
   fprintf(frst, "\nseq. seq.%8s %8s%8s%10s\tparas.\n", 
      "Ks/Ka", "Ks", "Ka", "lnL");
   tree.branches[0][0]=0;  tree.branches[0][1]=1;
   tree.nnode = (com.ntime=tree.nbranch=1) + 1;
   tree.origin=0;
   com.ns=2;
   com.clock=0;
   BranchToNode ();

   fpatt0=(int*) malloc(npatt0*3*sizeof(int));
   FOR (k, ns0) pz0[k]=com.z[k];
   com.z[0]=(char*)(fpatt0+npatt0);   com.z[1]=com.z[0]+npatt0;

   com.chunk =(double*) malloc ((com.ns-1)*n*com.npatt*sizeof(double));
   if (com.chunk==NULL||fpatt0==NULL) error ("oom");
   FOR (h, npatt0) fpatt0[h]=com.fpatt[h];
   PointLklnodes ();

   FOR (ii, ns0) {
      FOR (jj, ii) {
         printf ("%4d vs. %3d", ii+1, jj+1);
         fprintf (fout, "\n\n%4d vs. %4d", ii+1, jj+1);
         fprintf (frst, "%4d %4d ", ii+1, jj+1);

         FOR (k,n*n) fp[k]=0;
         for(h=0,zero(com.pi,n); h<npatt0; h++) {
            fp[(pz0[ii][h]-1)*n+pz0[jj][h]-1]+=fpatt0[h];
            fp[(pz0[jj][h]-1)*n+pz0[ii][h]-1]+=fpatt0[h];
            com.pi[pz0[ii][h]-1]+=fpatt0[h]/(2.*com.ls);  
            com.pi[pz0[jj][h]-1]+=fpatt0[h]/(2.*com.ls);
         }

         for (k=0,com.npatt=0; k<n; k++) for (h=0; h<=k; h++) {
            if (fp[k*n+h]) { 
               com.z[0][com.npatt]=k+1; com.z[1][com.npatt]=h+1;
               if (h==k) com.fpatt[com.npatt++]=fp[k*n+h]/2;
               else      com.fpatt[com.npatt++]=fp[k*n+h];
	    }
	 }
         com.posG[1]=com.npatt;

         GetInitials (x, space);   np=com.np;  NFunCall=0;
/*
         lnL = com.plfun (x, np);
         printf ("npatt:%4d%14.6f\n", com.npatt, lnL);

         ming1(noisy>2?frub:NULL,&lnL,com.plfun,NULL,testx,x,space,e,np);
*/
         SetxBound (np, xb);
         ming2(noisy>2?frub:NULL,&lnL,com.plfun,NULL,x,xb, space,e,np);

         fprintf (fout,"\nlnL%14.6f\n", -lnL);
         FOR (k, np)  fprintf (fout,"%9.5f", x[k]);

         if (com.getSE) {
            Hessian (np, x, lnL, space, var, com.plfun, var+np*np);
            matinv (var, np, np, var+np*np);
            FOR (k, np)
               fprintf (fout,"%9.5f", var[k*np+k]>0.?sqrt(var[k*np+k]):0.);
            FPN (fout);
	 }
TIME2S=x[0];
         EigenQ61 (fout, Root, U, V, com.kapa, com.omega, space);
         fprintf (frst, "%10.2f", -lnL);
         FOR (k,com.np) fprintf (frst, "%10.5f", x[k]); FPN (frst);
         fflush(fout);   fflush(frst);

         printf ("%14.6f%6d%6d", -lnL, com.npatt, NFunCall);
         FOR (k,com.np) printf ("%10.5f", x[k]);  FPN(F0);
      }  /* for (jj) */
   }     /* for (ii) */
   free(fpatt0);  free(com.chunk);
   return (0);
}
