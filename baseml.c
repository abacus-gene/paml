/* baseml.c 
      Maximum likelihood parameter estimation for aligned DNA (RNA) sequences,
                  combined with phylogenetic tree estimation.
                     Copyright, Ziheng YANG, July 1992 onwards

                           cc -c -fast tools.c eigen.c
                   cc -o baseml -fast baseml.c tools.o eigen.o -lm
                           baseml <ControlFileName>
*/
#include "tools.h"
#define NS            100
#define NBRANCH       (NS*2-2)
#define NNODE         (NS*2-1)
#define NGENE         15
#define LSPNAME       30
#define NCODE         4
#define NCATG         40
#define NP            (NBRANCH+NGENE+11)
/*
#define NP            (NBRANCH+3*NNODE+NGENE+11+NCATG*NCATG-1)
*/
extern int noisy, NFunCall, *ancestor;
extern double *SeqDistance; 

int Forestry (FILE *fout, double space[]);
void DetailOutput (FILE *fout, double x[]);
int GetOptions (char *ctlf);
int GetInitials (double x[]);
int SetxInitials (double x[]);
int SetParameters (double x[]);
int SetPGene (int igene, int _pi, int _UVRoot, int _alpha, double x[]);
int testx (double x[], int np);
double lfunAdG (double x[], int np);
double lfundG (double x[], int np);
double lfun (double x[], int np);
int PartialLikelihood (int inode, int igene);

void ReadPatternFreq (char* pexpf);
int BootstrapSeq (FILE *fout, double space[]);
int TestModel (FILE *fout, double x[], int nsep, double space[]);
int OldDistributions (int inode, double x[]);

struct CommonInfo {
   char *z[NS], spname[NS][LSPNAME+1], seqf[64],outf[64],treef[64];
   int seqtype, ns, ls, ngene, posG[NGENE+1], lgene[NGENE], *pose,npatt;
   int np, ntime, nrgene, nrate, nalpha, npi, nhomo, ncatG, maxnp, ncode,Mgene;
   int fix_kappa, fix_rgene, fix_alpha, fix_rho, nparK;
   int clock, model, getSE, runmode, print,verbose, ndata;
   float *fpatt;
   double pi[6], kappa, alpha, rho, rgene[NGENE], piG[NGENE][6];
   double freqK[NCATG], rK[NCATG], MK[NCATG*NCATG];
   double (*plfun)(double x[],int np), *chunk, *fhK;
}  com;
struct TREEB {
   int  nbranch, nnode, origin, branches[NBRANCH][2];
   double lnL;
}  tree;
struct TREEN {
   int father, nson, sons[NS], ibranch;
   double branch, divtime, kappa, pi[6], *lkl;
}  nodes[2*NS-1];

int nR=4;
static double PMat[16], Cijk[64], Root[4];
int PMatCijk (double PMat[], double t);

FILE *frub, *flfh, *frst;
char *models[]={"JC69","K80","F81","F84","HKY85","TN93","REV","UNREST"};
enum {JC69, K80, F81, F84, HKY85, TN93, REV, UNREST} MODELS;

#define BASEML 1
#include "treesub.c"
#include "treespace.c"

int main(int argc, char *argv[])
{
   FILE *fout, *fseq;
   int i, s1=0, s2=0, s3, READPATTERNFREQ=0, idata;
   char rstf[]="rst", ctlf[]="baseml.ctl";
   char *Mgenestr[]={"diff. rate", "separate data", "diff. rate & pi", 
                     "diff. rate & kappa", "diff. rate & pi & kappa"};
   double *space=NULL;

   com.ndata=1;
   noisy=0;
   com.runmode=2;     /* 0: user tree;  1: semi-automatic;  2: automatic  */

   com.clock=0;       /* 1: clock, rooted tree;  0: no clock, unrooted tree  */
   com.fix_rgene=0; /* 0: estimate rate factors for genes */
   com.nhomo=0;     
   com.getSE=0;       /* 1: want S.E.s of estimates;  0: don't want them */

   com.model=0;
   com.fix_kappa=0; com.kappa=5;
   com.fix_alpha=1; com.alpha=0.;  com.ncatG=4;    /* alpha=0 := inf */
   com.fix_rho=1;  com.rho=0.;   com.nparK=0;
   com.ncode=4;
   com.print=0;      com.verbose=1; 

   frub=fopen ("rub", "w");   flfh=fopen ("lfh", "w");
   frst=fopen (rstf, "w");

   if (argc>1)  strcpy(ctlf, argv[1]); 
   GetOptions (ctlf);
   if ((fout=fopen (com.outf, "w"))==NULL) error("outfile creation err.");
   if((fseq=fopen (com.seqf,"r"))==NULL)  {
      printf ("Sequence file %s not found!\n", com.seqf);
      exit (-1);
   }

   for (idata=0; idata<com.ndata; idata++) {
      if (com.ndata>1) {
         printf ("\nData set %d\n", idata+1);
         fprintf (fout, "\n\nData set %d\n", idata+1);
      }
      if (idata)  GetOptions (ctlf);
      if (READPATTERNFREQ) {
         if (com.rho || com.Mgene) error("model");
         ReadPatternFreq ("pexp.outt");
      }
      else
         ReadSeq ((com.verbose?fout:NULL), fseq, 0);
      if (com.ngene>1 && com.Mgene==1)  OutSeqsMgenes ();
      if (com.Mgene && com.ngene==1)  error ("option Mgene for 1 gene?");
      if(com.nhomo && com.ngene>1) error("nhomo not used with multiple genes");
      if(com.print && com.ngene>1) error("RateAncestor unavailable");
      if(com.nalpha && (!com.alpha || com.ngene==1 || com.fix_alpha))
         error("Malpha");
      if (com.alpha==0)  com.nalpha=0;
      else               com.nalpha=(com.nalpha?com.ngene:!com.fix_alpha);

      if (com.Mgene==1) com.nalpha=!com.fix_alpha;
      if (com.nalpha>1 && com.rho) error ("Malpha or rho");
      if (com.ngene==1) com.Mgene=0;
      if((com.Mgene>=2 && com.model==JC69) || (com.Mgene>=3 && com.model==F81) 
        || ((com.Mgene==2 || com.Mgene==4) && com.model==K80)
        || (com.Mgene>1 && com.nhomo>1) || (com.Mgene>=2 && com.model==UNREST))
         error ("model || Mgene");
      fprintf (fout,"BASEML %15s %8s ", com.seqf, models[com.model]);
      if (com.clock) fprintf (fout, " Clock  ");
      if (!com.nparK && com.alpha && com.rho) fprintf (fout, "  Auto-");
      if (com.alpha) fprintf (fout,"dGamma (ncatG=%d)", com.ncatG);
      if (com.nalpha>1) fprintf (fout,"(%d gamma)", com.nalpha);
      if (com.ngene>1) 
         fprintf (fout, " (%d genes: %s)  ", com.ngene, Mgenestr[com.Mgene]);
      if (com.nhomo)
       fprintf(fout,"\nNonhomo:%2d  fix_kappa%2d\n",com.nhomo,com.fix_kappa);
      if (com.nparK && com.ncatG>3)
         printf("\a\n%d rate categories for nonparametric?", com.ncatG);
      if (com.nparK) fprintf(fout,"\nnparK:%4d  K:%4d\n",com.nparK,com.ncatG);
   
      com.maxnp=(com.ns*2-2)*(1+(com.nhomo==2)) + (com.nhomo>2)*3*(com.ns*2-1) 
               +com.ngene + 12 + (com.nparK>=1)*com.ncatG*(com.ncatG-1);
      com.maxnp=max(50, com.maxnp);
      s3 = com.maxnp*(com.maxnp*2+2+12)*sizeof(double);
      s3 = max (s3, com.ls*(com.ns*sizeof(char)+sizeof(int)));
      if (space) free(space);
      if ((space=(double*)malloc(s3))==NULL)  error ("oom space");

      if (!READPATTERNFREQ) {
         if (SeqDistance) free(SeqDistance);   if (ancestor) free(ancestor);
         SeqDistance=(double*)malloc(com.ns*(com.ns-1)/2*sizeof(double));
         ancestor=(int*)malloc(com.ns*(com.ns-1)/2*sizeof(int));
         if (SeqDistance==NULL||ancestor==NULL) error("oom");
      }

      if (!READPATTERNFREQ) Initialize (fout, space, 0);
   
      if (com.Mgene==3) FOR (i,com.ngene) xtoy(com.pi, com.piG[i], 6);
      if (com.model==JC69 && !com.print && !READPATTERNFREQ) 
         PatternJC69like (fout);
   
      if (com.chunk) free(com.chunk);
      com.chunk=(double*)malloc((s1=(com.ns-1)*4*com.npatt)*sizeof(double));
      if (com.alpha || com.nparK) {
         s2 = com.npatt*com.ncatG;
         if (com.fhK) free(com.fhK);
         if ((com.fhK=(double *)malloc(s2*sizeof(double)))==NULL) error("oom");
      }
      if (noisy) printf("\nMem:\n%8d B lkl\n%8d B fhK\n%8d B space\n",
                         s1*sizeof(double), s2*sizeof(double), s3);
      if (com.chunk==NULL)  error ("oom chunk");

      FOR (i, com.ns*2-1) xtoy (com.pi, nodes[i].pi, 6);
      if(!READPATTERNFREQ) 
         DistanceMatNuc(fout,com.model,com.alpha);
/*
      BootstrapSeq (fout, space);
*/
      if (com.Mgene==1)        MultipleGenes (fout, space);
      else if (com.runmode==0) Forestry (fout, space);
      else if (com.runmode==3) StepwiseAddition (fout, space);
      else if (com.runmode>=4) Perturbation(fout, (com.runmode==4), space);
      else                  fprintf(frst,"%2d",StarDecomposition(fout,space));
   }
   fclose (fseq);
   if (noisy) putchar ('\a');
   return (0);
}

/* x[]: t[ntime], rgene[ngene-1], kappa[nbranch], pi[nnode*3], 
        { alpha, rho || rK[], fK[] || rK[], MK[] }
*/

int Forestry (FILE *fout, double space[])
{
   int  status=0, inbasemlg=0, i,j, itree, ntree, np, np0;
   double x[NP], xcom[NP-NBRANCH], lnL,lnL0=0,lnLm=0,e=1e-7, nchange;
   double xb[NP][2], tl=0,tlls, *var=space+com.maxnp;
   FILE *ftree, *finbasemlg=NULL, *fin=(FILE*)fopen("in.baseml","r"), *frates;

   if (com.alpha) 
      if ((frates=(FILE*)fopen("Rates.out","w"))==NULL) error("file err");
   if (com.alpha && com.rho==0 && com.nhomo==0 && com.nparK==0)
       { inbasemlg=1;  finbasemlg=fopen ("in.basemlg", "w"); }
   if ((ftree=fopen (com.treef,"r"))==NULL) error ("treefile error");
   fscanf (ftree, "%d%d", &i, &ntree);
   if (i!=com.ns) error ("ns in the tree file");
   fprintf (flfh,"\n%6d%6d%6d\n", ntree, com.ls, com.npatt);

   FOR (itree, ntree) {
      if (noisy) printf ("\nTREE # %2d: ", itree+1);
      fprintf (fout,"\nTREE # %2d:  ", itree+1);
      fprintf (frub,"\nTREE # %2d", itree+1);
      fprintf (flfh,"\n\n%2d\n", itree+1);

      if (ReadaTreeN (ftree, &i, 1)) error ("err tree..");

      nchange=MPScore (space);
      if (noisy) { OutaTreeN(F0,0,0);  printf("  MP score: %.2f\n",nchange); }
      OutaTreeN (fout, 0, 0);   fprintf (fout, "  MP score: %.2f", nchange);

      fflush (fout);   fflush (flfh);

      GetInitials (x);  
      if ((np=com.np)>NP) error("raise MP");
      if (itree) { np0=np; }
      if (itree && com.nhomo<=2)
         for (i=0; i<np-com.ntime; i++) x[com.ntime+i]=max(xcom[i],0.001);
      if (fin) {
         if (itree==0) {
            puts("\nInitials from in.baseml. Break if not correct");
            getchar ();
         }
         FOR (i,np) if (fscanf(fin,"%lf",&x[i])!=1) break;
         if (i<np)  error("in.baseml unusable. Edit or remove it.");
      }

      PointLklnodes ();
      NFunCall=0;
      SetxInitials (x);  /* pull the initial point into the feasible region */

      tlls = (com.clock ? 0 : sum(x,com.ntime));
      if(noisy)
         printf("ntime & nrate & np:%4d%4d%4d\n",com.ntime,com.nrate,com.np);
/*
      fprintf (fout,"\nLS branch lengths \n");
      FOR (i, com.ntime) fprintf (fout," %8.5f", x[i]);
*/
      if (noisy) {
         if (noisy>2) matout (F0, x, 1, np); 
         lnL=com.plfun(x,np);
         printf ("\nlnL0 = %12.6f\n", -lnL);
/*
         gradient (np, x, lnL, space, com.plfun, space+np, 1);
         matout (F0, space, 1, np);
*/
      }
      i=(com.clock==0 && !com.npi && (!com.nparK || com.ncatG==2) );
      if (i) {
         SetxBound (np, xb);
         j=ming2(noisy>2?frub:NULL,&lnL,com.plfun,NULL,x,xb, space,e,np);
      }
      else 
         j=ming1(noisy>2?frub:NULL,&lnL,com.plfun,NULL,testx,x,space,e,np);
      if (j || lnL<=0 || lnL>1e7) status=-1;  else status=0;
      if (itree==0) lnL0=lnLm=lnL;
      else if (lnL<lnLm) lnLm=lnL;
      if (noisy>1) printf ("\nOut...\nlnL  = %12.6f\n", -lnL);
      fprintf(fout,"\nlnL(ntime:%3d  np:%3d):%14.6f%+12.6f\n",
          com.ntime,np,-lnL,-lnL+lnL0);
      OutaTreeB (fout);  FPN (fout);
      FOR (i,np) fprintf(fout," %8.5f",x[i]);  FPN(fout);  fflush(fout);
      if (inbasemlg) matout (finbasemlg, x, 1, np);
/*
FOR(i,np) fprintf(frst," %5.3f", x[i]);  fprintf(frst,"%10.2f\n", lnL);
*/
      if (com.getSE) {
         Hessian (np, x, lnL, space, var, com.plfun, var+np*np);
         matinv (var, np, np, var+np*np);
         fprintf (fout,"SEs for parameters:\n");
         FOR (i, np)
            fprintf (fout," %8.5f", var[i*np+i]>0.?sqrt(var[i*np+i]):0.);
         FPN (fout);
         if (com.getSE==2) matout2(fout, var, np, np, 12,7);
      }
      if (com.clock) {
         SetBranch (x, 0);
         for(i=0,tl=0;i<tree.nnode;i++) if (i-tree.origin) tl+=nodes[i].branch;
      }
      else  tl=sum(x, com.ntime);
      fprintf (fout,"\ntree length: %9.5f", tl);
      if (com.ngene>1) fprintf (fout," (for 1st gene)");
      FPN(fout);  OutaTreeN (fout, 1, 1);  FPN(fout);
      if (com.verbose && com.np>com.ntime) DetailOutput (fout, x);
      if (status) {
         printf ("convergence?\n");  fprintf (fout,"check convergence..\n");
      }
      if (itree==0) 
         for (i=0; i<np-com.ntime; i++) xcom[i]=x[com.ntime+i];
      else if (!j)  
         for (i=0; i<np-com.ntime; i++) xcom[i]=xcom[i]*.8+x[com.ntime+i]*0.2;
/*
      TestModel (fout, x, 1, space);
      OldDistributions (tree.origin, x);
      fprintf (fout, "\n\n# of lfun calls:%10d\n", NFunCall);
*/
      if (com.model>=REV) 
         EigenREV(fout, x+com.ntime+com.nrgene, com.pi, &nR, Root, Cijk);
      
      com.print-=9;  com.plfun (x, np);  com.print+=9;
      if (com.print) {
         if (com.nhomo==0 && com.nparK==0 && com.model<=REV)
            AncestralSeqs (frst, x, space);
         if (com.plfun!=lfun)  
            lfunAdG_rate (frates, x, np);
      }
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
   if (inbasemlg) fclose (finbasemlg);
   fclose (ftree);  if(fin) fclose(fin);
   if (ntree>2)  fprintf (fout, "\nBest likelihood:%12.6f\n", -lnLm);
   return (0);
}


void DetailOutput (FILE *fout, double x[])
{
   int i,j,k=com.ntime, nr[]={0, 1, 0, 1, 1, 2, 5, 11};
   double Qfactor, *p=com.pi, t=0, k1,k2, S,V;

   fprintf (fout, "\nDetailed output identifying parameters\n");
   if (com.nrgene) {
      fprintf (fout, "\nrates for %d genes:%6.0f", com.ngene, 1.);
      FOR (i,com.nrgene) fprintf (fout, " %8.5f", x[k++]);
      FPN(fout);
   }
   if (com.nhomo==1 || com.nhomo>2) {
      fprintf (fout, "kappa under %s:", models[com.model]);
      FOR (i,com.nrate) fprintf (fout, " %8.5f", x[k++]);  FPN(fout);

      fprintf (fout, "base frequencies (%d sets):", com.npi);
      FOR (i, com.npi) {
         if (com.npi>1) fprintf (fout, "\n  Set #%d:", i+1);
         FOR (j,3) fprintf (fout, " %8.5f", x[k+3*i+j]);
         fprintf (fout, "%11.5f", 1-sum(x+k+3*i, 3));
      }
      FPN (fout);
      k+=3*com.npi;
   }
   else if (!com.fix_kappa) {
      fprintf (fout, "kappa or rate parameters in %s, ", models[com.model]);
      if (com.nhomo==2) {
         fprintf (fout, "\nbranch         t    kappa      TS     TV\n");
         FOR (i,tree.nbranch) {
            if (com.model==F84)  { k1=1+x[k+i]/p[4]; k2=1+x[k+i]/p[5]; }
            else                   k1=k2=x[k+i];
            S=2*p[0]*p[1]*k1+2*p[2]*p[3]*k2; V=2*p[4]*p[5];  Qfactor=1/(S+V);
            t=(com.clock ? nodes[tree.branches[i][1]].branch : x[i]);
            
            fprintf(fout,"%2d..%-2d %9.5f %8.5f %9.5f %8.5f\n",
               tree.branches[i][0]+1,tree.branches[i][1]+1, t,x[k+i],
               t*S/(S+V), t*V/(S+V));
	 }
      }
      else if (com.Mgene==3 || com.Mgene==4) {
         FOR (i,com.ngene) {
            fprintf (fout, "\n  gene #%d:", i+1);
            FOR (j,nr[com.model]) fprintf (fout, " %8.5f", x[k++]);
	 }
      }
      else  FOR (i,com.nrate) fprintf (fout, " %8.5f", x[k++]);
      FPN (fout);
   }
   if (!com.fix_alpha) {
      if (com.nalpha<=1) fprintf (fout, "gamma parameter alpha:");
      else  fprintf (fout, "gamma parameter alpha  for %d genes:", com.ngene);
      for (i=0; i<com.nalpha; i++) fprintf (fout, " %8.5f", x[k++]);
      FPN (fout);
   }
   if (!com.fix_rho) {
      fprintf (fout, "rho for the auto-discrete-gamma model: %9.5f",x[k]);
      FPN (fout);
   }
   if ((com.alpha || com.nparK>=1) && com.nalpha<=1) {
      fprintf (fout, "rates for categories:");
      FOR (i,com.ncatG) fprintf (fout, " %8.5f", com.rK[i]);
      fprintf (fout, "\nfreqs for categories:");
      FOR (i,com.ncatG) fprintf (fout, " %8.5f", com.freqK[i]);
      FPN (fout);
   }
   if (com.rho || com.nparK>=3 && com.nalpha<=1) {
      fprintf (fout, "transition probabilities between rate categories:\n");
      for(i=0;i<com.ncatG;i++,FPN(fout))  FOR(j,com.ncatG) 
         fprintf (fout, " %8.5f", com.MK[i*com.ncatG+j]);
   }
   FPN (fout);
}

extern double SmallDiff;

int GetOptions (char *ctlf)
{
   int i, nopt=24, lline=255; 
   char line[255], *pline, opt[20], comment='*';
   char *optstr[]={"seqfile","outfile","treefile","noisy","verbose","runmode",
        "clock","fix_rgene","Mgene","nhomo","getSE","RateAncestor",
        "model","fix_kappa","kappa","fix_alpha","alpha","Malpha","ncatG", 
        "fix_rho","rho","nparK", "ndata", "Small_Diff"};
   double t;
   FILE  *fctl=fopen (ctlf, "r");

   if (fctl) {
      printf ("\n\nReading options from %s..\n", ctlf);
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
                  case ( 2): sscanf(pline+2, "%s", com.treef);    break;
                  case ( 3): noisy=(int)t;           break;
                  case ( 4): com.verbose=(int)t;     break;
                  case ( 5): com.runmode=(int)t;     break;
                  case ( 6): com.clock=(int)t;       break;
                  case ( 7): com.fix_rgene=(int)t;   break;
                  case ( 8): com.Mgene=(int)t;       break;
                  case ( 9): com.nhomo=(int)t;       break;
                  case (10): com.getSE=(int)t;       break;
                  case (11): com.print=(int)t;       break;
                  case (12): com.model=(int)t;       break;
                  case (13): com.fix_kappa=(int)t;   break;
                  case (14): com.kappa=t;            break;
                  case (15): com.fix_alpha=(int)t;   break;
                  case (16): com.alpha=t;            break;
                  case (17): com.nalpha=(int)t;      break;
                  case (18): com.ncatG=(int)t;       break;
                  case (19): com.fix_rho=(int)t;     break;
                  case (20): com.rho=t;              break;
                  case (21): com.nparK=(int)t;       break;
                  case (22): com.ndata=(int)t;       break;
                  case (23): SmallDiff=t;            break;
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

   if (com.fix_alpha==1 && com.alpha==0) { 
      if (com.rho) puts ("rho reset to 0.");  com.fix_rho=1; com.rho=0;  
   }
   if (com.nparK>=1) { 
      com.fix_alpha=com.fix_rho=1; 
      if (com.alpha==0) com.alpha=0.5;   if (com.nparK>2) com.rho=0.4; 
   }
   if(com.model!=F84 && com.kappa<=0)  error("init kappa..");
   if(!com.fix_alpha && com.alpha<=0) { com.alpha=0.5; puts("init alpha reset");}
   if(!com.fix_rho && com.rho==0) { com.rho=0.001;  puts("init rho reset");}

   if (com.alpha)  if (com.ncatG<2 || com.ncatG>NCATG) error ("ncatG");

   if (com.nhomo==2) {
      if (com.model!=K80 && com.model!=F84 && com.model!=HKY85) error("nhomo");
   }
   else if (com.nhomo>2 && com.model!=F84 && com.model!=HKY85) error("nhomo");
   else
      if (com.nhomo==1 && com.model!=F81 && com.model!=F84 && com.model!=HKY85)
         error("nhomo");
   if (com.nhomo>1 && com.runmode>0)  error("nhomo incompatible with runmode");
   if (com.runmode==3 && (com.clock || com.npi || com.nparK))
      error("runmode+clock etc.");

   if (com.model==JC69 || com.model==K80 || com.model==UNREST)
      if (com.nhomo!=2)  com.nhomo=0;
   if (com.model==JC69 || com.model==F81) { com.fix_kappa=1; com.kappa=1; }
   if (com.model==TN93 || com.model==REV)   com.fix_kappa=0;

   return (0);
}

int GetInitials (double x[])
{
   int i,j,k, K=com.ncatG;
   double t;

   com.plfun=lfunAdG;
   if (com.alpha==0 && com.nparK==0)  com.plfun=lfun;
   else if ((com.alpha && com.rho==0) || com.nparK==1 || com.nparK==2)
      com.plfun=lfundG;

   com.nrgene = (!com.fix_rgene)*(com.ngene-1);
   com.nrate=0;
   if (!com.fix_kappa) {
      if (com.model<=HKY85)     com.nrate=1;
      else if (com.model==TN93) com.nrate=2;
      else                      com.nrate=(com.model==REV?5:11);
      if (com.Mgene>=3)         com.nrate*=com.ngene;
   }
   switch (com.nhomo) {
   case (0): com.npi=0;           break;   /* given 1 pi */
   case (1): com.npi=1;           break;   /* solve 1 pi */
   case (2): com.npi=0;  com.nrate=tree.nbranch;  break;
   case (3): com.npi=com.ns+(tree.origin>=com.ns)+(tree.nnode>com.ns+1);  
             com.nrate=(com.fix_kappa?1:tree.nbranch);
                                  break;   /* ns+2 pi    */
   case (4): com.npi=tree.nnode;  com.nrate=(com.fix_kappa?1:tree.nbranch);
                                  break;   /* nnode pi   */
   }

   if (com.clock) FOR (j,com.ntime) x[j]=0.01+0.0001*(double)(com.ntime-j);
   else           FOR (j,com.ntime) x[j]=0.015;

   LSDistance (&t, x, testx);
   FOR (i, com.ntime) x[i]=max(x[i], 0.0001);

   FOR (j,com.nrgene) x[com.ntime+j]=1;
   if (com.model<=TN93 && com.Mgene<=1)  
      EigenTN93 (com.model, com.kappa, com.kappa, com.pi, &nR, Root, Cijk);
   if (com.model==REV || com.model==UNREST)
      FOR (j, (com.Mgene>=3?com.ngene:1)) {
         k=com.ntime+com.nrgene+j*(com.model==REV?5:11);
         FOR (i,com.nrate) x[k+i]=0.2;
         if (com.model==REV)  x[k]=1;
         else x[k]=x[k+3]=x[k+8]=1;
      }
   else 
      FOR (i, com.nrate) x[com.ntime+com.nrgene+i]=com.kappa;
/*
   FOR (i,com.npi*3)  x[com.ntime+com.nrgene+com.nrate+i]=.25;
*/
   FOR (i,com.npi) xtoy(com.pi, x+com.ntime+com.nrgene+com.nrate+3*i, 3);

   com.np = k = com.ntime+com.nrgene+com.nrate+com.npi*3;
   if (com.alpha==0 && com.nparK==0) return(0);

   if (com.rho)
      AutodGamma (com.MK, com.freqK, com.rK, &t, com.alpha,com.rho,K);
   else 
      DiscreteGamma (com.freqK, com.rK, com.alpha, com.alpha, K, 0);

   for (i=0; i<com.nalpha; i++) x[k++]=com.alpha;
   if (!com.fix_rho)   x[k++]=com.rho;

   if (com.nparK) { xtoy (com.rK, x+k, K-1);  k+=K-1; }
   switch (com.nparK) {
   case (2):                            /* rK & fK */
      xtoy (com.freqK, x+k, K-1);       k+=K-1;         break;
   case (3):                            /* rK & MK (double stochastic) */
      fillxc (x+k, 1./K, (K-1)*(K-1));  k+=(K-1)*(K-1); break;
   case (4):                            /* rK & MK */
      fillxc (x+k, 1./K, K*(K-1));      k+=K*(K-1);     break;
   }
   com.np=k;
   return (0);
}

static complex cU[16], cV[16], cRoot[16];
int cPMat(double P[],double t,int n,complex cU[],complex cV[],complex cRoot[]);

int SetParameters (double x[])
{
/* set com. variables and initialize U, V, Root etc,
*/
   int i, j, k, NC=com.ncatG, status=0;
   double k1=com.kappa, k2=com.kappa, t, *p, space[NCATG*(NCATG+1)];

   k=SetBranch(x,1);
#if DEBUG
   if(k) puts ("\nbranch len err..");
#endif
   FOR (i, com.nrgene) com.rgene[i+1]=x[com.ntime+i];
   if (!com.fix_kappa && com.model<=TN93) {
       com.kappa=k1=k2=x[com.ntime+com.nrgene];
       if (com.model==TN93) k2=x[com.ntime+com.nrgene+1];
   }
   if (com.nhomo==1) {
      k=com.ntime+com.nrgene+com.nrate;
      xtoy (x+k, (p=com.pi), 3);
      if ( (p[3]=1-sum(p,3)) <= 0) status=-1;
      p[4]=p[0]+p[1];  p[5]=p[2]+p[3];
      if (com.model<=TN93) EigenTN93(com.model, k1, k2, p, &nR, Root, Cijk);
   }
   else if (com.nhomo==2) 
      for (i=0,k=com.ntime+com.nrgene; i<tree.nbranch; i++)
         nodes[tree.branches[i][1]].kappa=x[k+i];
   if (com.model<=TN93 && !com.fix_kappa && com.nhomo==0 && com.Mgene<=1)
      RootTN93 (com.model, k1, k2, com.pi, &t, Root);
   else if (com.nhomo>2) {
      for (i=0,k=com.ntime+com.nrgene; i<tree.nbranch; i++)
         nodes[tree.branches[i][1]].kappa=(com.fix_kappa?x[k]:x[k+i]);
      k+=com.nrate;
      if (com.nhomo==3) {
         FOR (i,com.ns) xtoy (x+k+i*3, nodes[i].pi, 3);
         for (i=com.ns; i<tree.nnode; i++)
            xtoy (x+k+com.ns*3, nodes[i].pi, 3);
         if (tree.origin>=com.ns && tree.nnode>com.ns+1)
            xtoy(x+k+3*(com.ns+1), nodes[tree.origin].pi, 3);
      }
      else
         FOR (i,tree.nnode) xtoy (x+k+i*3, nodes[i].pi, 3);
      FOR (i, tree.nnode) {
         p=nodes[i].pi;
         if ( (p[3]=1-sum(p,3)) <= 0) status=-1;
         p[4]=p[0]+p[1];  p[5]=p[2]+p[3];
      }
      xtoy (nodes[tree.origin].pi, com.pi, 6);
   }
   else if (com.model==REV && com.Mgene<=1)
      EigenREV (NULL, x+com.ntime+com.nrgene, com.pi, &nR, Root, Cijk);
   else if (com.model==UNREST && com.Mgene<=1)
      EigenUNREST (NULL, x+com.ntime+com.nrgene, com.pi, &nR, cRoot, cU, cV);

   if (com.nparK==0 && (com.alpha==0 || com.fix_alpha*com.fix_rho==1))
      return(status);
   if (com.nalpha>1) return (status);
   k = com.ntime+com.nrate+com.nrgene+com.npi*3;
   if (!com.fix_alpha) {
      com.alpha=x[k++];
      if (com.fix_rho)
         DiscreteGamma (com.freqK,com.rK,com.alpha,com.alpha,NC,0);
   }
   if (!com.fix_rho) {
      com.rho=x[k++];
      AutodGamma (com.MK, com.freqK, com.rK, &t,com.alpha,com.rho,NC);
   }
   if (com.nparK==0) return(status);

   xtoy (x+k, com.rK, NC-1);
   if (com.nparK==2) {
      xtoy (x+k+NC-1, com.freqK, NC-1);
      com.freqK[NC-1]=1-sum(com.freqK, NC-1);
   }
   else if (com.nparK==3) {   /* rK & MK (double stochastic matrix) */
      for (i=0, k+=NC-1; i<NC-1; k+=NC-1, i++) {
         xtoy (x+k, com.MK+i*NC, NC-1);
         com.MK[i*NC+NC-1]=1-sum(com.MK+i*NC,NC-1);
      }
      FOR (j, NC)
         for(i=0,com.MK[(NC-1)*NC+j]=1;i<NC-1; i++)
            com.MK[(NC-1)*NC+j]-=com.MK[i*NC+j];
   }
   else if (com.nparK==4) { /* rK & MK */
      for (i=0, k+=NC-1; i<NC; k+=NC-1, i++) {
         xtoy (x+k, com.MK+i*NC, NC-1);
         com.MK[i*NC+NC-1]=1-sum(com.MK+i*NC,NC-1);
      }
      PtoPi(com.MK, com.freqK, NC, space);
   }
   com.rK[NC-1]=(1-innerp(com.freqK, com.rK, NC-1))/com.freqK[NC-1];

   return (status);
}

int SetPGene (int igene, int _pi, int _UVRoot, int _alpha, double x[])
{
   int nr[]={0, 1, 0, 1, 1, 2, 5, 11};
   int k=com.ntime+com.nrgene+(com.Mgene>=3)*igene*nr[com.model];

   if (_pi) xtoy (com.piG[igene], com.pi, 6);
   if (_UVRoot) {
      if (com.model==K80) com.kappa=x[k]; 
      else if (com.model<=TN93) 
      EigenTN93(com.model,x[k],com.model==TN93?x[k+1]:0,com.pi,&nR,Root,Cijk);
      else if (com.model==REV)
         EigenREV (NULL, x+k, com.pi, &nR, Root, Cijk);
   }
   if (_alpha) {
      com.alpha=x[com.ntime+com.nrgene+com.nrate+com.npi+igene];
      DiscreteGamma (com.freqK, com.rK, com.alpha, com.alpha, com.ncatG, 0);
   }
   return (0);
}

int SetxInitials (double x[])
{
/* This forces initial values to the boundary 
*/
   int i,k;
   double tb[]={.0005, 3}, rgeneb[]={.001,9}, rateb[]={.1,20};
   double alphab[]={0.08, 99}, rhob[]={-0.1, 0.9};

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
   if (com.nparK) {
      for (i=0,k=com.ntime+com.nrgene+com.nrate; i<com.ncatG-1; i++,k++) { 
         if (x[k]>rateb[1]) x[k]=rateb[1];
         if (x[k]<rateb[0]) x[k]=rateb[0];
      }
      if (com.nparK==2) { if(x[k]>.99) x[k]=.99; if(x[k]<.01) x[k]=.01; k++; }
      /* if (com.nparK>2) error ("err SetxInitial.. do this now"); */
   }

   return(0);
}

int SetxBound (int np, double xb[][2])
{
   int i,j, k=com.ntime+com.nrgene+com.nrate;
   double tb[]={.4e-5, 20}, rgeneb[]={1e-4,99}, rateb[]={1e-5,999};
   double alphab[]={0.005, 99}, rhob[]={-0.2, 0.99};

   if (com.npi || (com.nparK && com.ncatG>2 && com.nparK%2==0))
       error ("npi || nparK in SetxBound");

   FOR (i,com.ntime)  FOR (j,2) xb[i][j]=tb[j];
   FOR (i,com.nrgene) FOR (j,2) xb[com.ntime+i][j]=rgeneb[j]; 
   FOR (i,com.nrate)  FOR (j,2) xb[com.ntime+com.nrgene+i][j]=rateb[j];
   for (i=0;i<com.nalpha;i++,k++)  FOR (j,2) xb[k][j]=alphab[j];
   if (!com.fix_rho)   FOR (j,2) xb[np-1][j]=rhob[j];
   if (com.nparK) {
      FOR (i,com.ncatG-1) { xb[k][0]=rateb[0]; xb[k++][1]=rateb[1]; }
      if (com.nparK==2) { xb[k][0]=1e-5; xb[k++][1]=1-1e-5; }
      if (com.nparK>2) error ("err SetxBound.. do this now");
   }

   return(0);
}

int testx (double x[], int np)
{
   int i,j, k, K=com.ncatG;
   double t, tb[]={.4e-5,20}, rgeneb[]={1e-4,99}, rateb[]={1e-5,99},pib=1e-5;
   double alphab[]={0.005, 40}, rhob[]={-0.2, 0.99}, Mb[]={1e-5, 1-1e-5};

   if (com.clock && SetBranch (x, 0))   return (-1);
   FOR (i,com.ntime)   if (x[i]<tb[0] || x[i]>tb[1]) return (-1);
   if (np==com.ntime) return (0); 
   FOR (i, com.nrgene)
      if (x[com.ntime+i]<rgeneb[0] || x[com.ntime+i]>rgeneb[1])
         return (-1);

   for (i=0,k=com.ntime+com.nrgene; i<com.nrate; i++,k++) 
      if (x[k]<rateb[0] || x[k]>rateb[1])  return (-1);
   FOR (i, com.npi)  {
      FOR (j,3) if (x[k+i*3+j] < pib) return (-1);
      if (sum(x+k+3*i,3) > 1-pib)     return(-1);
   }

   if (com.alpha==0 && com.nparK==0) return(0);
   k = com.ntime+com.nrate+com.nrgene+com.npi*3;
   for (i=0; i<com.nalpha; i++,k++)
      if (x[k]<alphab[0] || x[k]>alphab[1]) return (-1);
   if (!com.fix_rho) {
      if (x[k]<rhob[0] || x[k]>rhob[1]) return (-1);
      k++;
   }
   if (com.nparK==0) return (0);
   for (i=0; i<com.ncatG-1; i++,k++)        /* rK */
      if (x[k]<rateb[0] || x[k]>rateb[1])  return (-1);
   if (com.nparK==2) {                      /* fK */
      FOR (j,com.ncatG-1) if (x[k+j]<pib) return (-1);
      if (sum(x+k, com.ncatG-1) > 1-pib)  return(-1);
   }
   if (com.nparK<=2) return (0);

   for (i=0; i<(com.nparK==3?K-1:K); k+=K-1, i++) {
      xtoy (x+k, com.MK+i*K, K-1);
      com.MK[i*K+K-1]=1-sum(com.MK+i*K,K-1);
   }
   if (com.nparK==3) 
      FOR (j, K) {
         for (i=0,t=0; i<K-1; i++) t += com.MK[i*K+j];
         com.MK[(K-1)*K+j] = 1-t;
      }
   FOR (i, K*K) if (com.MK[i]<Mb[0]) return(-1);
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
      t=nodes[ison].branch*com.rgene[igene];
      if (com.nhomo>=2)
         EigenTN93(com.model,nodes[ison].kappa,1,nodes[ison].pi,&nR,Root,Cijk);
      if (com.model<=K80)
	      PMatK80 (PMat, t, (com.nhomo==2?nodes[ison].kappa:com.kappa));
      else if (com.model<=REV)  PMatCijk (PMat, t);
      else                      cPMat (PMat, t, n, cU, cV, cRoot);

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

int PMatCijk (double P[], double t)
{
/* P(t)ij = SUM Cijk * exp{Root*t}
*/
   int i,j,k, n=4, nr=nR;
   double expt[4];

   if (t<-1e-3) printf ("\nt = %.5f in PMatCijk", t);
   if (t<1e-7) { identity (P, n); return(0); }

   for (k=1; k<nr; k++) expt[k]=exp(t*Root[k]);
   FOR (i,n) FOR (j,n) {
      for (k=1,t=Cijk[i*n*nr+j*nr+0]; k<nr; k++)
         t+=Cijk[i*n*nr+j*nr+k]*expt[k];
      P[i*n+j]=t;
   }
   return (0);
}

#if 0

int CollapsSite (FILE *fout, int nsep, int ns, int *ncat, int SiteCat[])
{
   int j,k, it, h, b[NS], ndiff, n1, bit1;
/* n1: # of 1's   ...  bit1: the 1st nonzero bit */
   
   *ncat = 5 + (1<<(ns-1))-1 + nsep*11;
   if (fout) fprintf (fout, "\n# cat:%5d  # sep:%5d\n\n", *ncat, nsep);
 
   FOR (h, 1<<2*ns) {
      for (j=0,it=h; j<ns; b[ns-1-j]=it%4,it/=4,j++) ;
      for (j=1,ndiff=0; j<ns; j++)  {
         FOR (k,j) if (b[j]==b[k]) break;
         if (k==j) ndiff++;
      }
      switch (ndiff) {
      default  : SiteCat[h]=0;      break;
      case (0) : SiteCat[h]=b[0]+1; break;
      case (1) :
         for (j=1,it=0,n1=0,bit1=0; j<ns; j++) {
            k = (b[j]!=b[0]);
            it = it*2 + k;
            n1 += k;
            if (bit1==0 && k) bit1=ns-1-j;
         }
         it = 5 + it-1;
         if (nsep==0) { SiteCat[h]=it; break; }

         SiteCat[h]=it+min(bit1+1,nsep)*11;
         if (n1==1 && bit1<nsep) {
            SiteCat[h]-=11;
            SiteCat[h]+=(b[0]*4+b[ns-1-bit1]-b[0]-(b[0]<=b[ns-1-bit1]));
         }
         break;
      }
      if (fout) {
         FOR (j, ns) fprintf (fout, "%1c", NUCs[b[j]]);
         fprintf (fout, "%5d    ", SiteCat[h]);
         if (h%4==3) FPN (fout);
      }
   }
   return (0);
}

int GetPexpML (double x[], int ncat, int SiteCat[], double pexp[])
{
   int  j, it, h, nodeb[NNODE]; 
   int  isum, nsum=1<<(2*(tree.nbranch-com.ns+1));
   double fh, y, Pt[NBRANCH][16];

   if (com.ngene>1 || com.nhomo || com.alpha || com.nparK || com.model>REV)
      error ("Pexp()");
   SetParameters (x);
   FOR (j, tree.nbranch) 
      PMatCijk (Pt[j], nodes[tree.branches[j][1]].branch);

   for (h=0; h<(1<<2*com.ns); h++) {
      if (SiteCat[h] == 0) continue;
      for (j=0,it=h; j<com.ns; nodeb[com.ns-1-j]=it&3,it>>=2,j++) ;
      for (isum=0,fh=0; isum<nsum; isum++) {
         for (j=0,it=isum; j<tree.nbranch-com.ns+1; j++)
            { nodeb[com.ns+j]=it%4; it/=4; }
         for (j=0,y=com.pi[nodeb[tree.origin]]; j<tree.nbranch; j++) 
            y*=Pt[j][nodeb[tree.branches[j][0]]*4+nodeb[tree.branches[j][1]]];
         fh += y;
      }
      pexp[SiteCat[h]] += fh;
   }    
   pexp[0] = 1-sum(pexp+1,ncat-1);
   return (0);
}


int TestModel (FILE *fout, double x[], int nsep, double space[])
{
/* test of models, using com.
*/
   int j,h, it, ls=com.ls, ncat, *SiteCat;
   double *pexp=space, *nobs, lmax0, lnL0, X2, ef, de;

   SiteCat = (int *) malloc ((1<<2*com.ns)* sizeof (int));
   if (SiteCat == NULL)  error ("oom");
   CollapsSite (F0, nsep, com.ns, &ncat, SiteCat);
   fprintf (fout, "\n\nAppr. test of model.. ncat%6d  nsep%6d\n", ncat,nsep);

   nobs = pexp+ncat;
   zero (pexp, 2*ncat);
    /* nobs */
   FOR (h, com.npatt) {
      for (j=0,it=0; j<com.ns; j++) it = it*4+(com.z[j][h]-1);
      nobs[SiteCat[it]] += (double)com.fpatt[h];
   }
   GetPexpML (x, ncat, SiteCat, pexp);

   for (h=0,lnL0=0,X2=0,lmax0=-(double)ls*log((double)ls); h<ncat; h++) {
      if (nobs[h]>1) {
         lmax0 += nobs[h]*log(nobs[h]);
         lnL0 += nobs[h]*log(pexp[h]);
      }
      ef = com.ls*pexp[h];
      de = square(nobs[h]-ef)/ef;
      X2 += de;
      fprintf (fout, "\nCat #%3d%9.0f%9.2f%9.2f", h, nobs[h], ef,de);
   }
   fprintf (fout, "\n\nlmax0:%12.4f  D2:%12.4f   X2:%12.4f\n",
       lmax0, 2*(lmax0-lnL0), X2);
   free (SiteCat);

   return (0);
}

double freq[NNODE][4];

int OldDistributions (int inode, double x[])
{
/* reconstruct nucleotide frequencies at and down inode
   for models with com.nhomo>2
*/
   int i;
   double kappa=com.kappa;

   if (com.model>=REV)  error ("model");
   if (inode==tree.origin) {
      if (com.alpha || com.model>REV) puts ("\ninaccurate");
      SetParameters (x);
      xtoy (nodes[inode].pi, freq[inode], 4);
   }
   else {
      EigenTN93 (com.model, kappa, kappa, nodes[inode].pi,
          &nR, Root, Cijk);
      PMatCijk (PMat, nodes[inode].branch);

      matby (freq[nodes[inode].father], PMat, freq[inode], 1, 4, 4);
   }

   printf ("\nnode%2d b %8.5f k %8.5f sum PI =%12.6f =%12.6f\n",
      inode+1, nodes[inode].branch, nodes[inode].kappa,
      sum(freq[inode],4), sum(freq[inode]+4,2));

   fprintf (flfh, "\nnode %4d", inode+1);
   FOR (i, 4) fprintf (flfh, " %8.5f", freq[inode][i]);

   FOR (i, nodes[inode].nson)  
      OldDistributions (nodes[inode].sons[i], x);
   return (0);
}

#endif


#if 0

int BootstrapSeq (FILE *fout, double space[])
{
/* Bootstrap sequences, uses com.pose[] and destroys com.fpatt[]
*/
   int ib,h,it, nboot=100, npatt;

   printf ("\nBe patient: bootstrapping\n");
   for (ib=0; ib<nboot; ib++) {
      if (com.ngene>1) error ("ngene>1 in BootstrapSeq");
      for (h=0; h<com.npatt; h++)  com.fpatt[h]=0;
      for (h=0; h<com.ls; h++) 
         com.fpatt[com.pose[(int)(com.ls*rndu())]] ++;
      for (h=0,npatt=0; h<com.npatt; h++)  npatt+=(com.fpatt[h]>0);
      fprintf (fout, "\n%6d", npatt);

      it = StarDecomposition (fout, space);
      OutaTreeN (F0, 0, 0);   FPN(F0);
   }
   exit (0);
}

#endif


void ReadPatternFreq (char* pexpf)
{
   FILE *fpexp;
   int h, j, ch;
   float sum;

   printf ("\n\nRead site-pattern frequencies from %s\n", pexpf);
   if ((fpexp=fopen(pexpf, "r")) == NULL) error ("file err");
   fscanf (fpexp, "%d%d", &com.ns, &com.npatt);
   printf ("%3d species, %d site patterns\n", com.ns, com.npatt);

   if (com.model>1) {
      error("Should calculate base frequencies\a\n");
   }
   else {
      fillxc(com.pi, .25, 4);  com.pi[4]=com.pi[5]=.5;
   }
   xtoy (com.pi, com.piG[0], 6);

   if(com.npatt>(1<<(2*com.ns))) error("npatt");
   com.ls=com.npatt;
   com.fpatt=(float*) malloc (com.npatt*sizeof(float));
   FOR (j, com.ns) com.z[j]=(char*) malloc ((com.ls+1)*sizeof(char));
   if (com.fpatt==NULL || com.z[com.ns-1]==NULL) error ("oom");
   FOR (j,com.ns) com.z[j][com.ls]='\0';

   com.ngene=1; com.lgene[0]=com.ls; com.posG[0]=0; com.posG[1]=com.ls;
   com.rgene[0]=1;
   for (h=0; h<com.npatt; h++) {
      for (j=0; j<com.ns; j++) {
         ch=fgetc(fpexp);
         while (!isalpha(ch)) ch=fgetc(fpexp);
         com.z[j][h]=CodeChara((char)ch, 0);
      }
      if (fscanf(fpexp,"%f",&com.fpatt[h]) != 1) error ("ReadPatternFreq");
   }
   for (j=0,sum=0; j<com.npatt; j++) sum+=com.fpatt[j];
   printf ("\nSUM freq = %.6f = 1?\n\n", sum);
   fclose (fpexp);
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


/* problems and notes
   (1) Discrete-Gamma model: waste of memory; 
       com.fhK[com.npatt*sizeof(double)] is enough and so much is
       used in lfundG(), while com.fhK[npatt*sizeeof(double)*com.ncatG]
       is allocated, as this leads to one function lfunAdG_pr() for both 
       AdG and dG.

   (2) AdG model: generation of com.MK[K*K] is not independent
       of com.alpha or DiscreteGamma().

__________________________________________________________________________
non-homogeneous process models:
  nhomo            fix_kappa         models
 
  0 (1 pi given)   0 (solve 1 kappa)   JC69, K80, F81, F84, HKY85,
                                       REV(5), UNREST(11)
                   1 (1 kappa given)   K80, F84, HKY85
 
  1 (solve 1 pi)   0 (as above)       F84, HKY85, REV(5)
                   1                  F81(0), F84, HKY85        
 
  2 (ns+2 pi)      0 (solve 1 kappa)   F84 & HKY85  
                   1 (nbranch kappa)   F84 & HKY85  
 
  3 (nnode pi)     0,1  (as above)

space-time process models:
  nparK     fix_alpha       fix_rho        parameters 
  0         (0,1)          (0,1)           alpha & rho for AdG

  1         set to 1       set to 1        rK[]        (freqK=1/K) (K-1)
  2         set to 1       set to 1        rK[] & freqK[]          2(K-1)
  3         set to 1       set to 1        rK[] & MK[] (freqK=1/K) K(K-1)
  4         set to 1       set to 1        rK[] & MK[]             K*K-1
__________________________________________________________________________

*/
