/* baseml.c
     Maximum likelihood parameter estimation for aligned DNA (RNA) sequences,
                 combined with phylogenetic tree estimation.
                    Copyright, Ziheng YANG, July 1992 onwards

                         cc -c -fast tools.c eigen.c
                  cc -o baseml -fast baseml.c tools.o eigen.o -lm
                          baseml <ControlFileName>
*/

#ifdef __MWERKS__
/* Added by Andrew Rambaut to accommodate Macs -
   Brings up dialog box to allow command line parameters.
*/
#include <console.h>
#endif

#include "paml.h"

#define NS            1000
#define NBRANCH       (NS*2-2)
#define NNODE         (NS*2-1)
#define NGENE         1000
#define LSPNAME       30
#define NCODE         4
#define NCATG         8

#define NP            (NBRANCH+NGENE+11)
/*
#define NP            (NBRANCH+3*NNODE+NGENE+11+NCATG*NCATG-1)
*/
extern int noisy, NFunCall, NEigenQ, NPMatUVRoot, *ancestor, GeneticCode[][64];
extern char AAs[];
extern double *SeqDistance; 


int Forestry (FILE *fout);
void DetailOutput (FILE *fout, double x[], double var[]);
int GetOptions (char *ctlf);
int GetInitials(double x[], int *fromfile);
int SetxInitials (double x[]);
int SetParameters (double x[]);
int SetPGene (int igene, int _pi, int _UVRoot, int _alpha, double xcom[]);
int SetPSiteClass(int iclass, double xcom[]);
int testx (double x[], int np);
int PartialLikelihood (int inode, int igene);
int cPMat(double P[],double t,int n,complex cU[],complex cV[],complex cRoot[]);

int TransformxBack(double x[]);
void ReadPatternFreq (FILE* fout, char* pexpf);

int TestModel (FILE *fout, double x[], int nsep, double space[]);
int OldDistributions (int inode, double oldfreq[]);
int SubData(int argc, char *argv[]);
int GroupDistances();

struct CommonInfo {
   char *z[NS], *spname[NS], seqf[96],outf[96],treef[96], cleandata;
   int seqtype, ns,ls,ngene,posG[NGENE+1],lgene[NGENE], *pose,npatt;
   int np,ntime,nrgene,nrate,nalpha,npi,nhomo,ncatG,ncode,Mgene,sspace,slkl1;
   int fix_kappa, fix_rgene, fix_alpha, fix_rho, nparK, fix_blength;
   int clock, model, getSE, runmode, print,verbose, ndata, icode,coding;
   int method, npi0, readfpatt;
   double *fpatt, kappa, alpha, rho, rgene[NGENE], pi[4],piG[NGENE][4];
   double freqK[NCATG], rK[NCATG], MK[NCATG*NCATG];
   double pi_sqrt[NCODE];
   double (*plfun)(double x[],int np), *lkl, *lkl0, *fhK, *space;
}  com;
struct TREEB {
   int  nbranch, nnode, root, branches[NBRANCH][2];
   double lnL;
}  tree;
struct TREEN {
   int father, nson, sons[NS], ibranch;
   double branch, divtime, kappa, pi[4], *lkl, label;
}  *nodes;

int nR=4, LASTROUND;
double PMat[16], Cijk[64], Root[4];
complex cU[16], cV[16], cRoot[16];
char  StepMatrix[16];

double xcom[NP-NBRANCH];

int lklSiteClass=0;  /* is lkl memory allocated for each class? */
char _oldlkl[NNODE];     /* update lkl for nodes to save computation (0 yes; 1 no) */
int _nnodeScale=0;
char _nodeScale[NNODE];   /* nScale[ns-1] for interior nodes */
double *_nodeScaleF=NULL;  /* nScaleF[npatt] for scale factors */

FILE *frub, *flnf, *frst, *frst1, *fin=NULL;
char *ratef="rates";
char *models[]={"JC69","K80","F81","F84","HKY85","T92","TN93","REV","UNREST", "REVu","UNRESTu"};
enum {JC69, K80, F81, F84, HKY85, T92, TN93, REV, UNREST, REVu, UNRESTu} MODELS;

double _rateSite=1;
int *rateBranch=NULL, N_rateBranch=-1;
int N_PMatUVRoot=0;


#define BASEML 1
#include "treesub.c"
#include "treespace.c"

int main(int argc, char *argv[])
{
   FILE *fout, *fseq=NULL, *fpair[6];
   char pairfs[1][32]={"2base.t"};
   int i, slkl0=0,s2=0, idata;
   char rstf[96]="rst", ctlf[96]="baseml.ctl", timestr[64];
   char *Mgenestr[]={"diff. rate", "separate data", "diff. rate & pi", 
                     "diff. rate & kappa", "diff. rate & pi & kappa"};
   char *clockstr[]={"", "global clock", "local clock", "clock with dated seqs"};

#ifdef __MWERKS__
/* Added by Andrew Rambaut to accommodate Macs -
   Brings up dialog box to allow command line parameters.
*/
	argc=ccommand(&argv);
#endif

   starttime();
   com.ndata=1;
   com.cleandata=0;  noisy=0;  com.runmode=0;

   com.clock=0;
   com.fix_rgene=0;  /* 0: estimate rate factors for genes */
   com.nhomo=0;
   com.getSE=0;      /* 1: want S.E.s of estimates;  0: don't want them */

   com.seqtype=0;   com.model=0;
   com.fix_kappa=0; com.kappa=5;
   com.fix_alpha=1; com.alpha=0.;  com.ncatG=4;    /* alpha=0 := inf */
   com.fix_rho=1;   com.rho=0;     com.nparK=0;
   com.ncode=4;     com.icode=0;
   com.print=0;     com.verbose=1;  com.fix_blength=0;
   com.method=0;    com.space=NULL;
   SetSeed ((int)time(NULL));

   frub=gfopen("rub","w");  frst=gfopen(rstf,"w"); frst1=gfopen("rst1","w");

/*
com.pi[0]=0.45; com.pi[1]=0.45; com.pi[2]=0.05; com.pi[3]=0.05; 
com.model=4;
EigenTN93 (com.model,1,1,com.pi,&nR,Root,Cijk);
PMatCijk (PMat, 3./3.);
matout (F0, PMat, 4,4);
FOR(i,4) abyx(com.pi[i], PMat+i*4, 4);
matout (F0, PMat, 4,4);
exit(0);
*/

   if (argc>1)  strcpy(ctlf, argv[1]); 
   GetOptions (ctlf);

   fout=gfopen(com.outf, "w"); 
   if((fseq=fopen (com.seqf,"r"))==NULL)  {
      printf ("\n\nSequence file %s not found!\n", com.seqf);
      exit (-1);
   }
   fpair[0]=(FILE*)fopen(pairfs[0],"w");

   for (idata=0; idata<com.ndata; idata++) {
      if (com.ndata>1) {
         printf("\nData set %4d\n", idata+1);
         fprintf(fout, "\n\nData set %4d\n", idata+1);
         fprintf(frst1, "\t%4d", idata+1);
      }
      if (idata)  GetOptions (ctlf);
      if (com.readfpatt) {
         if (com.rho || com.Mgene) error2("rho or Mgene for readfpatt.");
         com.print=0;  com.cleandata=1;
         /* strcpy(com.seqf,"pexp.in"); */
         ReadPatternFreq (fout, com.seqf);
      }
      else
         ReadSeq ((com.verbose?fout:NULL), fseq);
      if(com.ndata==1) fclose(fseq);
      i=(com.ns*2-1)*sizeof(struct TREEN);
      if((nodes=(struct TREEN*)malloc(i))==NULL) error2("oom");

      /* BootstrapSeq("boot");  exit(0); */

      if(com.coding) {
         if(com.ls%3!=0 || (com.ngene!=1 && com.ngene!=3))
            error2("this is not a coding sequence.  Remove icode?");
      }

      if(com.ngene>1 && com.Mgene==1)  printSeqsMgenes (); 
      if(com.Mgene && com.ngene==1)  error2 ("option Mgene for 1 gene?");
      if(com.ngene>1 && com.nhomo) error2("nhomo for mutliple genes?");
      if(com.nalpha && (!com.alpha || com.ngene==1 || com.fix_alpha))
         error2("Malpha");
      if(com.nalpha>1 && com.rho) error2("Malpha or rho");
      if(com.alpha==0)  com.nalpha=0;
      else              com.nalpha=(com.nalpha?com.ngene:!com.fix_alpha);
      if(com.Mgene==1)  com.nalpha=!com.fix_alpha;

      if(com.ngene==1) com.Mgene=0;
      if((com.nhomo==1 && com.ngene>1) || (com.Mgene>1 && com.nhomo>=1))
         error2("nhomo does not work with Mgene options");

      if((com.Mgene>=2 && com.model==JC69) || (com.Mgene>=3 && com.model==F81) 
        || ((com.Mgene==2 || com.Mgene==4) && com.model==K80)
        || (com.Mgene>1 && com.nhomo>1) || (com.Mgene>=2 && com.model==UNREST))
         error2 ("model || Mgene");
      fprintf(fout,"BASEML (in %s)  %s  %s ", VerStr,com.seqf,models[com.model]);
      if(com.clock) fprintf(fout," %s ",clockstr[com.clock]);
      if (!com.nparK && com.alpha && com.rho) fprintf (fout, "  Auto-");
      if (com.alpha) fprintf (fout,"dGamma (ncatG=%d)", com.ncatG);
      if (com.nalpha>1) fprintf (fout,"(%d gamma)", com.nalpha);
      if (com.ngene>1) 
         fprintf (fout, " (%d genes: %s)  ", com.ngene, Mgenestr[com.Mgene]);
      if (com.nhomo>1)
         fprintf(fout,"\nNonhomo:%2d  fix_kappa%2d\n",com.nhomo,com.fix_kappa);
      if (com.nparK && com.ncatG>3)
         printf("\a\n%d rate categories for nonparametric model?", com.ncatG);
      if (com.nparK) fprintf(fout,"\nnparK:%4d  K:%4d\n",com.nparK,com.ncatG);

      if (com.clock==2)
        if((rateBranch=(int*)realloc(rateBranch,com.ns*2*sizeof(int)))==NULL) 
           error2("oom");
      com.sspace=800000;
      if((com.space=(double*)realloc(com.space,com.sspace))==NULL) 
         error2("oom space");

      i=com.ns*(com.ns-1)/2;
      SeqDistance=(double*)realloc(SeqDistance,i*sizeof(double));
      ancestor=(int*)realloc(ancestor, i*sizeof(int));
      if(SeqDistance==NULL||ancestor==NULL) error2("oom distance&ancestor");

      if (!com.readfpatt) Initialize (fout);
      if (com.Mgene==3) FOR (i,com.ngene) xtoy(com.pi,com.piG[i],4);

      if (com.model==JC69 && !com.print) 
         PatternJC69like(fout);

      if(!com.cleandata) {
         slkl0=com.ns*com.ncode*com.npatt*sizeof(double);
         if((com.lkl0=(double*)malloc(slkl0))==NULL) error2("oom lkl0");
      }
      com.slkl1=(com.ns-1)*com.ncode*com.npatt*sizeof(double);
      if((com.lkl=(double*)realloc(com.lkl,com.slkl1))==NULL) error2("oom lkl");
      if (com.alpha || com.nparK) {
         s2 = com.npatt*com.ncatG*sizeof(double);
         if((com.fhK=(double*)realloc(com.fhK,s2))==NULL) error2("oom");
      }

      printf ("\n%9ld bytes for distance ",
         com.ns*(com.ns-1)/2*(sizeof(double)+sizeof(int)));
      printf("\n%9d bytes for lkl0\n%9d bytes for lkl1\n",slkl0,com.slkl1);
      printf("%9d bytes for fhK\n%9d bytes for space\n",s2,com.sspace);

      /* FOR(i,com.ns*2-1) xtoy(com.pi,nodes[i].pi, 4);  check this? */

      DistanceMatNuc(fout,fpair[0],com.model,com.alpha);

      if (com.Mgene==1)        MultipleGenes (fout, fpair, com.space);
      else if (com.runmode==0) Forestry (fout);
      else if (com.runmode==2) fprintf(frst,"%2d",StarDecomposition(fout,com.space));
      else if (com.runmode==3) StepwiseAddition (fout, com.space);
      else if (com.runmode>=4) Perturbation(fout, (com.runmode==4), com.space);

      FPN(frst);  if((idata+1)%10==0) fflush(frst);
      free(nodes);
      printf("\nTime used: %s\n", printtime(timestr));
   }   /* for(idata) */
   if(com.ndata>1 && fseq) fclose(fseq);  
   fclose(fout); fclose(frst);   fclose(fpair[0]);
   if (noisy) putchar ('\a');

   return (0);
}


/* x[]: t[ntime], rgene[ngene-1], kappa[nbranch], pi[nnode*3], 
        { alpha, rho || rK[], fK[] || rK[], MK[] }
*/


static int ncorrect=0;

int Forestry (FILE *fout)
{
   static int times=0;
   int status=0, inbasemlg=0, i,j=0,k=0, itree,ntree, np,np0,iteration=1;
   int pauptree=0, btree=0, haslength;
   double x[NP], lnL,lnL0=0,lnLbest=0,e=1e-8, nchange=-1;
   double xb[NP][2], tl=0, *var=NULL;
   FILE *ftree, *finbasemlg=NULL, *frate=NULL;

   if(com.clock==3) GetSeqTimes();

   ftree=gfopen(com.treef,"r");
   GetTreeFileType(ftree,&ntree,&pauptree,0);
   if(com.alpha) frate=gfopen(ratef,"w");
   if (com.alpha && com.rho==0 && com.nhomo==0 && com.nparK==0) {
      inbasemlg=1;  finbasemlg=gfopen("in.basemlg","w");
   }
   flnf=gfopen("lnf","w+");
   fprintf(flnf,"%6d %6d %6d\n", ntree, com.ls, com.npatt);

   if(!com.cleandata) InitPartialLikelihood ();
   for(itree=0; ntree==-1||itree<ntree; itree++,iteration=1) {
      if((pauptree && PaupTreeRubbish(ftree)) || 
         ReadaTreeN(ftree,&haslength,&i,1))
            { printf("\nerr or end of tree file.\n"); break; }
      if(noisy) printf ("\nTREE # %2d: ", itree+1);
      fprintf (fout,"\nTREE # %2d:  ", itree+1);
      fprintf (frub,"\n\nTREE # %2d\n", itree+1);
      if(com.print) fprintf (frst,"\n\nTREE # %2d\n", itree+1);
      fprintf (flnf,"\n\n%2d\n", itree+1);

      LASTROUND=0;
      if (com.fix_blength==2 && !haslength) error2("no blengths in tree");
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

      if(com.cleandata) nchange=MPScore (com.space);
      if(noisy&&com.ns<99) 
         { OutaTreeN(F0,0,0); printf(" MP score: %.2f\n",nchange);}
      OutaTreeN(fout,0,0);  fprintf(fout,"  MP score: %.2f",nchange);
      if(!com.clock && com.model<=REV && com.nhomo<=2 
         && nodes[tree.root].nson<=2 && com.ns>2){
         puts("\nThis is a rooted tree, without clock.  Check.");
         if(com.verbose) fputs("\nThis is a rooted tree.  Please check!",fout);
      }
      fflush(fout);  fflush(flnf);

      GetInitials(x, &i);
      if(i==-1) iteration=0;
      if(iteration) SetxInitials (x); /* start within the feasible region */

      if((np=com.np)>NP) error2("raise NP and recompile");

      if((i=spaceming2(np))>com.sspace)
         if((com.space=(double*)realloc(com.space,com.sspace=i))==NULL)
            error2("oom space");
      if(itree) { np0=np; }
      if(itree && !fin && (com.nhomo==0 || com.nhomo==2))
         for (i=0; i<np-com.ntime; i++) x[com.ntime+i]=max2(xcom[i],0.001);
      if (com.clock==2) {
         printf("\n%d substitution rates for branches assumed:\n",N_rateBranch);
         FOR(i,tree.nbranch) printf("%3d", rateBranch[i]); FPN(F0);
         fprintf(fout,"\n%d substitution rates for branches assumed:\n",
            N_rateBranch);
         FOR(i,tree.nbranch) fprintf(fout,"%3d",rateBranch[i]); FPN(fout);
      }

      if(noisy>=2&&com.npi0) printf("\n%d zero frequencies.\n", com.npi0);
      PointLklnodes ();
      lnL = com.plfun (x,np);

      if(noisy) {
         printf("\nntime & nrate & np:%6d%6d%6d\n",com.ntime,com.nrate,com.np);
         if(noisy>2 && com.ns<50) { OutaTreeB(F0); FPN(F0); matout(F0,x,1,np); }
         printf("\nlnL0 = %12.6f\n",-lnL);

/*
         printf ("\nGradient:\n");
         gradient(np,x,lnL,com.space,com.plfun,com.space+np,1);
         matout(F0,com.space,1,np);
         printf ("\nHessian:\n");
         var=com.space+2*np;
         Hessian (np, x, lnL, com.space+np, var, com.plfun, var+np*np);
         FOR(i,np) printf(" %11.6f",var[i*np+i]);
         exit(0);
*/
      }

      if (iteration && np) {
         SetxBound (np, xb);
         if(com.method)
            j=minB(noisy>2?frub:NULL, &lnL,x,xb, e*100, com.space);
         else
            j=ming2(noisy>2?frub:NULL,&lnL,com.plfun,NULL,x,xb, com.space,e,np);

         if (j==-1 || lnL<=0 || lnL>1e7) status=-1;
         else status=0;
      }

      if (itree==0) { lnL0=lnLbest=lnL; btree=0; }
      else if (lnL<lnLbest) { lnLbest=lnL;  btree=itree; }
      if (noisy) printf ("Out...\nlnL  = %12.6f\n", -lnL);
      fprintf(fout,"\nlnL(ntime:%3d  np:%3d):%14.6f%+12.6f\n",
          com.ntime,np,-lnL,-lnL+lnL0);
      if(com.fix_blength<2) { OutaTreeB (fout);  FPN (fout); }
      if(LASTROUND==0) {
         LASTROUND=1;
         if((com.npi && com.model!=T92) || com.nparK>=2) TransformxBack(x);
      }
      if(com.clock) {  /* this applies to all clock models (clock=1,2,3) */
         for(i=com.ns; i<tree.nnode; i++) 
            if(i!=tree.root) x[i-com.ns]=nodes[i].divtime;
      }

for(i=0;i<np;i++) fprintf(frst1,"\t%.5f",x[i]); 
fprintf(frst1,"\t%.4f\t%d",-lnL,com.npatt);  FPN(frst1);
fflush(frst1);

      FOR(i,np) fprintf(fout," %8.5f",x[i]);  FPN(fout); fflush(fout);
      if(inbasemlg) matout(finbasemlg,x,1,np);

      if (com.getSE) {
         if(np>100) puts("Calculating SE's");
         if(com.sspace<np*(np+1)*(int)sizeof(double)) {
            com.sspace=np*(np+1)*sizeof(double);
            if((com.space=(double*)realloc(com.space,com.sspace))==NULL)
               error2("oom space for SE");
         }
         var=com.space+np;
         Hessian (np,x,lnL,com.space,var,com.plfun,var+np*np);
         matinv(var,np,np,var+np*np);
         fprintf(fout,"SEs for parameters:\n");
         FOR(i,np) fprintf(fout," %8.5f",(var[i*np+i]>0.?sqrt(var[i*np+i]):-1));
         FPN(fout);
         if (com.getSE==2) matout2(fout, var, np, np, 15, 10);
      }

      /* if(com.clock) SetBranch(x); */

      /* GroupDistances(); */

      for(i=0,tl=0;i<tree.nnode;i++) 
         if(i!=tree.root) tl+=nodes[i].branch;
      fprintf(fout,"\ntree length = %9.5f%s\n",tl,com.ngene>1?" (1st gene)":"");
      if(com.readfpatt)  { FPN(fout); OutaTreeN(fout,0,1); FPN(fout); }
      else {
         FPN(fout); OutaTreeN(fout,1,0); FPN(fout);
         FPN(fout); OutaTreeN(fout,1,1); FPN(fout);
      }

      if(com.np-com.ntime||com.clock)  DetailOutput(fout,x,var);
      if (status) {
         printf ("convergence?\n");  fprintf (fout,"check convergence..\n");
      }
      if (itree==0) 
         for (i=0; i<np-com.ntime; i++) xcom[i]=x[com.ntime+i];
      else if (!j)  
         for (i=0; i<np-com.ntime; i++) xcom[i]=xcom[i]*.8+x[com.ntime+i]*0.2;
/*
      TestModel(fout,x,1,com.space);
      fprintf(fout,"\n\n# of lfun calls:%10d\n", NFunCall);
*/

      com.print-=9;  com.plfun(x, np);  com.print+=9;
      if (com.print) {
         if (com.plfun!=lfun)  lfunRates (frate, x, np);
         /** think about more-general models **/
         if (com.nhomo==0 && com.nparK==0 && com.model<=REV && com.rho==0)
            AncestralSeqs(frst, x);
      }

      if(com.coding) { 
         fputs("\nTree with branch lengths for codon models:\n",fout);
         FOR(i,tree.nnode) nodes[i].branch*=(1+com.rgene[1]+com.rgene[2]);
         FPN(fout); OutaTreeN(fout,1,1); FPN(fout);
      }

   }         /* for (itree) */
   if(finbasemlg) fclose(finbasemlg);   if (frate) fclose(frate);
   fclose(ftree);  if(fin) { fclose(fin);  fin=NULL; }
   
   if(ntree>2)  { 
      fprintf(fout,"Best tree (#%d):%12.6f\n",btree+1,-lnLbest);
      if(btree==0) ncorrect++;
      fprintf(frst,"\t%d\t%d",btree+1,ncorrect);
   }
   if(ntree==-1) {
      ntree=itree;
      rewind(flnf);
      fprintf(flnf,"%6d", ntree);
   }

   if(noisy && ntree>1) { rewind(flnf);  rell(flnf,fout,ntree); }
   fclose(flnf);

   return(0);
}


int TransformxBack(double x[])
{
/* transform variables x[] back to their original definition after iteration,
   for output and for calculating SEs.
*/ 
   int i,k, K=com.ncatG;

   k=com.ntime+com.nrgene+com.nrate;
   for (i=0; i<com.npi; i++)  
      f_and_x(x+k+3*i,x+k+3*i,4,0,0);
   
   k+=com.npi*3 + K-1;        /* K-1 for rK */
   if (com.nparK==2)          /* rK & fK */
      f_and_x(x+k,x+k,K,0,0);
   else if (com.nparK==3)     /* rK & MK (double stochastic matrix) */
      for (i=0; i<K-1; k+=K-1,i++)  f_and_x(x+k,x+k,K,0,0);
   else if (com.nparK==4)     /* rK & MK */
      for (i=0; i<K;   k+=K-1,i++)  f_and_x(x+k,x+k,K,0,0);
   return(0);
}


void DetailOutput (FILE *fout, double x[], double var[])
{
   int i,j,k=com.ntime, nr[]={0, 1, 0, 1, 1, 1, 2, 5, 11};
   int n31pi=(com.model==T92?1:3);
   double Qfactor,*p=com.pi, t=0, k1,k2, S,V, Y=p[0]+p[1],R=p[2]+p[3],tnode;
   double *Qrate=x+com.ntime+com.nrgene, *oldfreq;

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
         fprintf(fout,"Node %3d Time %6.2f ",
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

   if(com.nhomo==1) {
      if(com.nrate) fprintf (fout, "kappa under %s:", models[com.model]);
      FOR(i,com.nrate) fprintf (fout, " %8.5f", x[k++]);  FPN(fout);
      fprintf (fout, "base frequency parameters\n");
      FOR (j,4) fprintf (fout, " %8.5f", com.pi[j]);
      k += n31pi;
   }
   else if(com.nhomo>=3) {
      fprintf (fout, "kappa under %s (in order of branches):", models[com.model]);
      FOR(i,com.nrate) fprintf(fout," %8.5f", x[k++]);  FPN(fout);
      SetParameters(x);
      if(com.alpha==0) {
         if((oldfreq=(double*)malloc(tree.nnode*4*sizeof(double)))==NULL) 
            error2("out of memory for OldDistributions()");
         OldDistributions (tree.root, oldfreq);
      }
      fputs("\n(frequency parameters for branches)  [frequencies at nodes] (see Yang & Roberts 1995 fig 1)\n\n",fout);
      for(i=0;i<tree.nnode;i++,FPN(fout)) {
         fprintf (fout, "Node #%d  (", i+1);
         FOR (j,4) fprintf (fout, " %7.5f ", nodes[i].pi[j]);
         if(com.alpha==0) {
            fprintf(fout,")  [");
            FOR(j,4) fprintf(fout," %7.5f",oldfreq[i*4+j]);
            fprintf(fout," ]");
         }
      }
      fprintf(fout,"\nNote: node %d is root.\n",tree.root+1);
      k += com.npi*n31pi;
      if(com.alpha==0) free(oldfreq);
   }
   else if (!com.fix_kappa) {
      fprintf(fout,"\nParameters in the rate matrix (%s) (Yang 1994 J Mol Evol 39:105-111):\n",
         models[com.model]);

      if (com.nhomo==2) {
         fprintf (fout, "\nbranch         t    kappa      TS     TV\n");
         FOR (i,tree.nbranch) {
            if (com.model==F84)  { k1=1+x[k+i]/R; k2=1+x[k+i]/R; }
            else                   k1=k2=x[k+i];
            S=2*p[0]*p[1]*k1+2*p[2]*p[3]*k2; V=2*Y*R;  Qfactor=1/(S+V);
            /* t=(com.clock ? nodes[tree.branches[i][1]].branch : x[i]); */
            t=nodes[tree.branches[i][1]].branch; 
            fprintf(fout,"%2d..%-2d %9.5f %8.5f %9.5f %8.5f\n",
               tree.branches[i][0]+1,tree.branches[i][1]+1, t,x[k+i],
               t*S/(S+V), t*V/(S+V));
         }
      }
      /* Mgene = 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff*/
      else if (com.Mgene>=2) {
         FOR (i,com.ngene) {
            fprintf (fout, "\n\nGene #%d\n", i+1);
            p = (com.Mgene==3 ? com.pi : com.piG[i]);
            Qrate = (com.Mgene==2 ? x+k : x+k+i*nr[com.model]);
            if (com.model<=TN93)
               FOR (j,nr[com.model]) fprintf(fout," %8.5f", Qrate[j]);
            else if (com.model==REV || com.model==REVu) 
               /* output Q matrix, no eigen calculation */
               EigenQREVbase(fout, Qrate, p, NULL, -1, &nR, Root, Cijk);
            else if (com.model==UNREST || com.model==UNRESTu) 
               EigenQunrest (fout, Qrate, p,&nR,cRoot,cU,cV);
         }
         if (com.Mgene>=3) k+=com.ngene*nr[com.model];
         else              k+=nr[com.model];
      }
      else {
         if (com.model<REV) FOR (i,com.nrate) fprintf (fout, " %8.5f", x[k++]);
         else k+=com.nrate;
      }
      FPN (fout);
   }

   if (com.Mgene<2) {
      if (com.model==REV || com.model==REVu) /* output Q, no eigen calculation */
         EigenQREVbase(fout, Qrate, com.pi, NULL, -1, &nR, Root, Cijk);
      else if (com.model==UNREST || com.model==UNRESTu) 
         EigenQunrest (fout, Qrate, com.pi, &nR, cRoot, cU, cV);
   }

   FOR(j,com.nalpha) {
      if (!com.fix_alpha)  
         fprintf(fout,"\nalpha (gamma, K=%d) = %8.5f", com.ncatG,(com.alpha=x[k++]));
      if(com.nalpha>1) 
         DiscreteGamma(com.freqK,com.rK,com.alpha,com.alpha,com.ncatG,0);
      fprintf(fout, "\nr: "); FOR(i,com.ncatG) fprintf(fout," %8.5f",com.rK[i]);
      fprintf(fout, "\nf: "); FOR(i,com.ncatG) fprintf(fout," %8.5f",com.freqK[i]);
      FPN(fout);
   }
   if (!com.fix_rho) {
      fprintf (fout, "rho for the auto-discrete-gamma model: %9.5f",x[k]);
      FPN (fout);
   }
   if (com.nparK>=1 && com.nalpha<=1) {
      fprintf (fout, "rates for categories:");
      FOR(i,com.ncatG) fprintf(fout," %8.5f",com.rK[i]);
      fprintf (fout, "\nfreqs for categories:");
      FOR(i,com.ncatG) fprintf(fout," %8.5f",com.freqK[i]);
      FPN(fout);
   }
   if (com.rho || com.nparK>=3 && com.nalpha<=1) {
      fprintf (fout, "transition probabilities between rate categories:\n");
      for(i=0;i<com.ncatG;i++,FPN(fout))  FOR(j,com.ncatG) 
         fprintf (fout, " %8.5f", com.MK[i*com.ncatG+j]);
   }
   FPN (fout);
}

extern double Small_Diff;
int GetStepMatrix(char*line)
{
/* read user definitions of the REV and UNREST models
   StepMatrix[4*4]:
      -1 at diagonals, 0 for default rate, and positive values for rates
*/
   char *p,*errstr="StepMatrix specification in the control file";
   int i,k, b1=-1,b2;

   p=strchr(line,'[');
   sscanf(++p,"%d", &com.nrate);
   if(com.nrate<0 || (com.model==REVu&&com.nrate>5) 
                  || (com.model==UNRESTu&&com.nrate>11)) error2(errstr);
   FOR(i,4) FOR(k,4) StepMatrix[i*4+k] = (i==k?-1:0);
   FOR(i,com.nrate) { 
      while(*p && *p!='(') p++;
      if(*p++ !='(') error2( "expecting (" );
      FOR(k,12) {
         while (isspace(*p)) p++;
         if(*p==')') break;
         b1=CodeChara(*p++,0);  b2=CodeChara(*p++,0);
         if(b1<0||b1>3||b2<0||b2>3) error2("bases out of range.");
         if(b1==b2||StepMatrix[b1*4+b2]>0) error2("pair already specified");
         if(com.model==REVu) StepMatrix[b1*4+b2]=StepMatrix[b2*4+b1]=i+1;
         else                StepMatrix[b1*4+b2]=i+1;
      }
      printf("rate %d: %d pairs\n", i+1,k);
   }

FOR(i,16) { printf("%3d", StepMatrix[i]); if((i+1)%4==0) FPN(F0); }

   return(0);
}

int GetOptions (char *ctlf)
{
   int i, nopt=29, lline=2048; 
   char line[2048], *pline, opt[20], comment='*';
   char *optstr[]={"seqfile","outfile","treefile","noisy",
        "cleandata", "verbose","runmode", "method",
        "clock","fix_rgene","Mgene","nhomo","getSE","RateAncestor",
        "model","fix_kappa","kappa","fix_alpha","alpha","Malpha","ncatG", 
        "fix_rho","rho","nparK", "ndata", "Small_Diff","icode", "readfpatt",
        "fix_blength"};
   double t;
   FILE *fctl;

   fctl=gfopen(ctlf,"r");
   printf ("Reading options from %s..\n", ctlf);
   for (;;) {
      if (fgets (line, lline, fctl) == NULL) break;
      for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
         if (isalnum(line[i]))  { t=1; break; }
         else if (line[i]==comment) break;
      if (t==0) continue;
      sscanf (line, "%s%*s%lf", opt, &t);
      if ((pline=strstr(line, "="))==NULL)
         error2("option file. add space around the equal sign =");

      for (i=0; i<nopt; i++) {
         if (strncmp(opt, optstr[i], 8)==0)  {
            if (noisy>=9)
               printf ("\n%3d %15s | %-20s %6.2f", i+1,optstr[i],opt,t);
            switch (i) {
               case ( 0): sscanf(pline+1, "%s", com.seqf);    break;
               case ( 1): sscanf(pline+1, "%s", com.outf);    break;
               case ( 2): sscanf(pline+1, "%s", com.treef);    break;
               case ( 3): noisy=(int)t;           break;
               case ( 4): com.cleandata=(char)t;  break;
               case ( 5): com.verbose=(int)t;     break;
               case ( 6): com.runmode=(int)t;     break;
               case ( 7): com.method=(int)t;      break;
               case ( 8): com.clock=(int)t;       break;
               case ( 9): com.fix_rgene=(int)t;   break;
               case (10): com.Mgene=(int)t;       break;
               case (11): com.nhomo=(int)t;       break;
               case (12): com.getSE=(int)t;       break;
               case (13): com.print=(int)t;       break;
               case (14): com.model=(int)t; 
                  if(com.model>UNREST) GetStepMatrix(line);  break;
               case (15): com.fix_kappa=(int)t;   break;
               case (16): com.kappa=t;            break;
               case (17): com.fix_alpha=(int)t;   break;
               case (18): com.alpha=t;            break;
               case (19): com.nalpha=(int)t;      break;
               case (20): com.ncatG=(int)t;       break;
               case (21): com.fix_rho=(int)t;     break;
               case (22): com.rho=t;              break;
               case (23): com.nparK=(int)t;       break;
               case (24): com.ndata=(int)t;       break;
               case (25): Small_Diff=t;           break;
               case (26): com.icode=(int)t; com.coding=1; break;
               case (27): com.readfpatt=(int)t; break;
               case (28): com.fix_blength=(int)t; break;
            }
            break;
         }
      }
      if (i==nopt)
        { printf ("\noption %s in %s not recognised\n", opt,ctlf); exit(-1); }
   }
   fclose (fctl);

   if (com.fix_alpha==1 && com.alpha==0) {
      if (!com.fix_rho || com.rho) error2("fix rho to 0 if alpha=0.");
   }
   if (com.nparK>=1) { 
      com.fix_alpha=com.fix_rho=1; 
      if(com.alpha==0) com.alpha=0.5; 
      if(com.nparK<=2) com.rho=0; else com.rho=0.4;
      if(com.nhomo>=1) error2("nhomo & nparK");
   }
   if(com.model!=F84 && com.kappa<=0)  error2("init kappa..");
   if(!com.fix_alpha && com.alpha<=0)  error2("init alpha..");
   if(!com.fix_rho && com.rho==0) { com.rho=0.001;  puts("init rho reset"); }

   if (com.alpha) 
      { if (com.ncatG<2 || com.ncatG>NCATG) error2 ("ncatG"); }
   else if (com.ncatG>1) com.ncatG=1;

   if(com.method &&(com.clock||com.rho)) 
      { com.method=0; puts("\aiteration method reset"); }

   if (com.nhomo==2) {
      if (com.model!=K80 && com.model!=F84 && com.model!=HKY85) error2("nhomo");
   }
   else if (com.nhomo>2 && !(com.model>=F84 && com.model<=T92)) error2("nhomo");
   else if (com.nhomo>2 && com.method) error2("nhomo & method.");
   else
      if (com.nhomo==1 && !(com.model>=F81 && com.model<=REV) && com.model!=REVu)
         error2("nhomo=1 and model");

   if (com.nhomo>1 && com.runmode>0)  error2("nhomo incompatible with runmode");
   if (com.clock && com.runmode>2)  error2("runmode & clock?");
   if (com.runmode==3 && (com.npi || com.nparK))
      error2("runmode incompatible with nparK or nhomo.");

   if (com.model==JC69 || com.model==K80 || com.model==UNREST)
      if (com.nhomo!=2)  com.nhomo=0;
   if (com.model==JC69 || com.model==F81) { com.fix_kappa=1; com.kappa=1; }
   if (com.model==TN93 || com.model==REV || com.model==REVu)  com.fix_kappa=0;
   if (com.nparK==3) {
      puts("\n\nnparK==3, double stochastic, may not work.  Use nparK=4?\n");
      getchar();
   }
   if(com.runmode!=-2) fin=fopen("in.baseml","r");
   return (0);
}



int GetInitials (double x[], int *fromfile)
{
   int i,j,k, K=com.ncatG, n31pi=(com.model==T92?1:3);
   int slkl1_new=(tree.nnode-com.ns)*com.ncode*com.npatt*com.ncatG*sizeof(double);
   double t=-1;

   NFunCall=NPMatUVRoot=NEigenQ=0;
   com.plfun=lfunAdG;
   if (com.alpha==0 && com.nparK==0)  com.plfun=lfun;
   else if ((com.alpha && com.rho==0) || com.nparK==1 || com.nparK==2)
      com.plfun=lfundG;

   if(com.method && com.fix_blength!=2 && com.plfun==lfundG) {
      lklSiteClass=1;
      if(com.slkl1!=slkl1_new) {
         com.slkl1=slkl1_new;
         printf("\n%9d bytes for lkl1, adjusted\n",com.slkl1);
         if((com.lkl=(double*)realloc(com.lkl,com.slkl1))==NULL)
            error2("oom lkl1");
      }
   }
   InitializeNodeScale();

   com.nrgene = (!com.fix_rgene)*(com.ngene-1);
   if (com.fix_kappa && (com.Mgene==3 || com.Mgene==4)) error2("Mgene options");

   if(com.model<=UNREST) {
      com.nrate=0;
      if (!com.fix_kappa) {
         if (com.model<=T92)        com.nrate=1;
         else if (com.model==TN93)  com.nrate=2;
         else                       com.nrate=(com.model==REV?5:11);
         if (com.Mgene>=3)          com.nrate*=com.ngene;
      }
   }
   switch (com.nhomo) {
   case (0): com.npi=0;           break;   /* given 1 pi */
   case (1): com.npi=1;           break;   /* solve 1 pi */
   case (2): com.npi=0;  com.nrate=tree.nbranch;  break;  /* b kappa's */
   case (3): com.npi=com.ns+(tree.root>=com.ns)+(tree.nnode>com.ns+1);
             com.nrate=(com.fix_kappa?1:tree.nbranch);  
             for(i=0; i<tree.nnode; i++)  nodes[i].label = (i<com.ns?i:com.ns);
             if(tree.root>=com.ns) nodes[tree.root].label = com.ns+1;
             break;   /* ns+2 pi */
   case (4): com.npi=tree.nnode;  com.nrate=(com.fix_kappa?1:tree.nbranch);
                                  break;   /* nnode pi   */
   case (5): com.nrate=(com.fix_kappa?1:tree.nbranch);
      for(i=0,com.npi=0; i<tree.nnode; i++) {
         j=(int)nodes[i].label;
         if(j+1>com.npi)  com.npi=j+1;
         if(j<0||j>tree.nnode-1) error2("node label in tree.");
      }
      printf("%d sets of frequency parameters\n",com.npi);
      break;   /* user-specified pi   */
   }

   if(com.fix_blength==2)  { com.ntime=0; com.method=0; }
   else if(com.clock==0)  com.ntime=tree.nbranch;
   else if(com.clock==1)  com.ntime=tree.nnode-com.ns+(tree.root<com.ns);
   else if(com.clock==3)  com.ntime=tree.nnode-com.ns+(tree.root<com.ns)+1;
   else {  /* if(com.clock==2) */
      if(com.nhomo>=2) error2("clock=2 & nhomo do not work together");
      for(i=0,k=0,N_rateBranch=0; i<tree.nbranch; i++) { 
         rateBranch[i]=j=(int)nodes[tree.branches[i][1]].label;
         if(j+1>N_rateBranch)  N_rateBranch=j+1;
         if(j<0||j>tree.nbranch-1) error2("branch label in the tree.");
         else if (j)  k=1;
      }
      if(k) puts("\nUsing branch marks in the tree file. Stop if wrong"); 
      else {
         OutaTreeB(F0); FPN(F0);
         for (i=0,N_rateBranch=0; i<tree.nbranch; i++) {
            printf("Branch %2d: %2d..%-2d? ", i+1,
               tree.branches[i][0]+1, tree.branches[i][1]+1);
            scanf("%d",&rateBranch[i]);
            if (rateBranch[i]+1>N_rateBranch) N_rateBranch=rateBranch[i]+1;
         }
      }
      printf("\n\n%d rates for branches. ",N_rateBranch);
      com.ntime=tree.nnode-com.ns + N_rateBranch-1;
      if(com.ntime>tree.nbranch) 
         printf("\a\nntime=%d, too many rates??\n",com.ntime);
   }

   if (com.clock) {
      for(j=1,x[0]=1; j<tree.nnode-com.ns+(tree.root<com.ns); j++)
         x[j]=0.2+0.8*rndu();
      for(; j<com.ntime; j++) x[j]=1;   /* rates for branches */
      if(com.clock==3) SetDivTime_OldSon(tree.root);
   }
   else if(com.fix_blength==-1) 
      FOR (j,com.ntime) x[j]=rndu();
   else if(com.fix_blength==0) {
      for(j=0,t=0,k=com.ns*(com.ns-1)/2; j<k; j++)  t+=SeqDistance[j]/k;
      t=max2(0.01,t);  t=min2(t,1.5);
      FOR (j,com.ntime) x[j]=t*rndu();
      if(com.ns<50) LSDistance (&t, x, testx);
   }

   FOR (j,com.nrgene) x[com.ntime+j]=1;
   if (com.model<=TN93 && com.Mgene<=1)  
      EigenTN93 (com.model,com.kappa,com.kappa,com.pi,&nR,Root,Cijk);
   if (com.model==REV || com.model==UNREST)
      FOR (j, (com.Mgene>=3?com.ngene:1)) {
         k=com.ntime+com.nrgene+j*(com.model==REV?5:11);
         FOR (i,com.nrate) x[k+i]=0.2+0.1*rndu();
         if (com.model==REV)  x[k]=1;
         else x[k]=x[k+3]=x[k+8]=1;
      }
   else 
      FOR(i,com.nrate) x[com.ntime+com.nrgene+i]=com.kappa;

   FOR(i,com.npi*n31pi) x[com.ntime+com.nrgene+com.nrate+i]=rndu()*.2;
   com.np = k = com.ntime+com.nrgene+com.nrate+com.npi*n31pi;

   if (com.alpha || com.nparK) {
      if (com.rho)
         AutodGamma(com.MK, com.freqK, com.rK, &t, com.alpha,com.rho,K);
      else 
         DiscreteGamma (com.freqK, com.rK, com.alpha, com.alpha, K, 0);

      for (i=0; i<com.nalpha; i++) x[k++]=com.alpha;
      if (!com.fix_rho)   x[k++]=com.rho;

      if (com.nparK) { xtoy(com.rK, x+k, K-1);  k+=K-1; }
      switch (com.nparK) {
      case (2):                            /* rK & fK */
         zero(x+k, K-1);       k+=K-1;         break;
      case (3):                            /* rK & MK (double stochastic) */
         zero(x+k, (K-1)*(K-1));  k+=(K-1)*(K-1); break;
      case (4):                            /* rK & MK */
         zero(x+k, K*(K-1));      k+=K*(K-1);     break;
      }
      com.np=k;
   }
   if(fin) readx(x,fromfile);
   else    *fromfile=0;

   return (0);
}




int SetParameters(double x[])
{
/* This sets parameters in com., nodes[], etc and is called before lfun() etc.
   Iinitialize U, V, Root etc, if necessary.
   For nhomo models (nhomo=1,3,4) 
      x[] has frequencies if (LASTROUND==1) or exp(pi)/(1+SUM(exp[pi])) if otherwise
*/
   int i, j, k, K=com.ncatG, status=0, n31pi=(com.model==T92?1:3);
   double k1=com.kappa, k2=com.kappa, t, space[NCATG*(NCATG+1)];

   if(com.fix_blength<2) SetBranch(x);
   if(com.np<=com.ntime) return(0);
   FOR (i, com.nrgene) com.rgene[i+1]=x[com.ntime+i];
   if (!com.fix_kappa && com.model<=TN93) {
       com.kappa=k1=k2=x[com.ntime+com.nrgene];
       if (com.model==TN93) k2=x[com.ntime+com.nrgene+1];
   }
   if (com.nhomo==1) {
      k=com.ntime+com.nrgene+com.nrate;
      if(com.model==T92)
         { com.pi[0]=com.pi[2]=(1-x[k])/2;  com.pi[1]=com.pi[3]=x[k]/2; }
      else {
         if (!LASTROUND) f_and_x(x+k,com.pi,4,0,0);
         else            xtoy (x+k, com.pi, 3);
         com.pi[3]=1-sum(com.pi,3);
      }
      if (com.model<=TN93)
         EigenTN93(com.model,k1,k2,com.pi,&nR,Root,Cijk);
      if (com.model==REV || com.model==REVu)
         getpi_sqrt (com.pi, com.pi_sqrt, com.ncode, &com.npi0);
   }
   else if (com.nhomo==2)
      for (i=0,k=com.ntime+com.nrgene; i<tree.nbranch; i++)
         nodes[tree.branches[i][1]].kappa=x[k+i];
 
   if (com.model<=TN93 && com.nhomo==0 && com.Mgene<=1)
      RootTN93 (com.model, k1, k2, com.pi, &t, Root);
   else if (com.nhomo>=3) {
      for (i=0,k=com.ntime+com.nrgene; i<tree.nbranch; i++)
         nodes[tree.branches[i][1]].kappa=(com.fix_kappa?x[k]:x[k+i]); /* ?? */
      k+=com.nrate;

      FOR (i,tree.nnode) {
         j = (com.nhomo==4 ? i : (int)nodes[i].label);
         if(com.model==T92) {
            nodes[i].pi[0]=nodes[i].pi[2] = (1-x[k+j])/2;
            nodes[i].pi[1]=nodes[i].pi[3] = x[k+j]/2;
         }
         else {
            if (!LASTROUND) f_and_x(x+k+j*3, nodes[i].pi, 4,0,0);
            else            xtoy   (x+k+j*3, nodes[i].pi, 3);
            nodes[i].pi[3]=1-sum(nodes[i].pi,3);
         }
      }
      xtoy(nodes[tree.root].pi, com.pi, 4);

/*
FOR (i,tree.nnode) {
printf("node %d (%2.0f ): ", i,nodes[i].label);
FOR(j,4) printf("%8.5f", nodes[i].pi[j]);
FPN(F0);
}
getchar();
*/
   }
   else if ((com.model==REV || com.model==REVu) && com.Mgene<=1)
      EigenQREVbase (NULL, x+com.ntime+com.nrgene, 
                     com.pi, com.pi_sqrt, com.npi0, &nR, Root, Cijk);
   else if ((com.model==UNREST || com.model==UNRESTu) && com.Mgene<=1)
      EigenQunrest (NULL, x+com.ntime+com.nrgene,com.pi,&nR,cRoot,cU,cV);

   if (com.nparK==0 && (com.alpha==0 || com.fix_alpha*com.fix_rho==1))
      return(status);
   if (com.nalpha>1) return (status);
   k = com.ntime+com.nrate+com.nrgene+com.npi*n31pi;
   if (!com.fix_alpha) {
      com.alpha=x[k++];
      if (com.fix_rho)
         DiscreteGamma (com.freqK,com.rK,com.alpha,com.alpha,K,0);
   }
   if (!com.fix_rho) {
      com.rho=x[k++];
      AutodGamma (com.MK, com.freqK, com.rK, &t,com.alpha,com.rho,K);
   }
   if (com.nparK==0) return(status);

   /* nparK models */
   xtoy (x+k, com.rK, K-1);

   if (com.nparK==2) {
      if (!LASTROUND)  f_and_x(x+k+K-1, com.freqK, K,0,0);
      else             xtoy  (x+k+K-1, com.freqK, K-1);
      com.freqK[K-1]=1-sum(com.freqK, K-1);
   }
   else if (com.nparK==3) {   /* rK & MK (double stochastic matrix) */
      for (i=0,k+=K-1; i<K-1; k+=K-1,i++) {
         if (!LASTROUND) f_and_x(x+k, com.MK+i*K, K,0,0);
         else            xtoy  (x+k, com.MK+i*K, K-1);
         com.MK[i*K+K-1]=1-sum(com.MK+i*K,K-1);
      }
      FOR(j, K) {
         for(i=0,com.MK[(K-1)*K+j]=1;i<K-1; i++)
            com.MK[(K-1)*K+j]-=com.MK[i*K+j];
         if (com.MK[(K-1)*K+j]<0)
            printf("SetPar: MK[K-1][j]=%.5f<0\n",com.MK[(K-1)*K+j]);
      }
   }
   else if (com.nparK==4) { /* rK & MK */
      for (i=0, k+=K-1; i<K; k+=K-1, i++) {
         if (!LASTROUND) f_and_x(x+k, com.MK+i*K, K,0,0);
         else            xtoy  (x+k, com.MK+i*K, K-1);
         com.MK[i*K+K-1]=1-sum(com.MK+i*K,K-1);
      }
      PtoPi(com.MK, com.freqK, K, space);
   }
   com.rK[K-1]=(1-innerp(com.freqK, com.rK, K-1))/com.freqK[K-1];
   return (status);
}


int SetPGene(int igene, int _pi, int _UVRoot, int _alpha, double xcom[])
{
/* xcom[] does not contain time parameters
   Note that com.piG[][] have been homogeneized if (com.Mgene==3)
*/
   int nr[]={0, 1, 0, 1, 1, 2, 5, 11};
   int k=com.nrgene+(com.Mgene>=3)*igene*nr[com.model];
   double ka1=xcom[k], ka2=(com.model==TN93?xcom[k+1]:-1);

   if(com.Mgene==2 && com.fix_kappa) ka1=ka2=com.kappa;

   if (_pi) {
      xtoy(com.piG[igene], com.pi, 4);
      if(com.model==REV || com.model==REVu)
         getpi_sqrt (com.pi, com.pi_sqrt, com.ncode, &com.npi0);
   }
   if (_UVRoot) {
      if (com.model==K80) com.kappa=ka1;
      else if (com.model<=TN93) 
         EigenTN93(com.model,ka1,ka2,com.pi,&nR,Root,Cijk);
      else if (com.model==REV || com.model==REVu)
         EigenQREVbase(NULL,xcom+k,com.pi, com.pi_sqrt, com.npi0, &nR, Root, Cijk);
   }
   if (_alpha) {
      com.alpha=xcom[com.nrgene+com.nrate+com.npi+igene]; /* check?? */
      DiscreteGamma(com.freqK,com.rK,com.alpha,com.alpha,com.ncatG,0);
   }
   return(0);
}


int SetxInitials (double x[])
{
/* This forces initial values to the boundary 
*/
   int i,k;
   double tb[]={.00006, 5}, rgeneb[]={.001,9}, rateb[]={.05,50};
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
   }
   return(0);
}

int SetxBound (int np, double xb[][2])
{
/* sets lower and upper bounds for variables during iteration
*/
   int i,j, k=com.ntime+com.nrgene+com.nrate, nf=0, n31pi=(com.model==T92?1:3);
   double tb[]={.4e-6, 99}, rgeneb[]={1e-4,999}, rateb[]={1e-5,999};
   double alphab[]={.04, 999}, rhob[]={-0.2, 0.99}, pb[]={.00001,.999999};
   double fb[]={-19,9}; /* transformed freqs.*/

   if(com.clock) {
      xb[0][0]=.00005;  xb[0][1]=tb[1];
      for(i=1;i<tree.nnode-com.ns;i++)  /* for node times, clock=1 or 2 */
         FOR(j,2) { xb[i][0]=pb[0]; xb[i][1]=pb[1]; }
      for(; i<com.ntime; i++)           /* for rates for branches, clock=2 */
         FOR(j,2) { xb[i][0]=tb[0]; xb[i][1]=tb[1]; }
   }
   else 
      FOR (i,com.ntime)  FOR (j,2) xb[i][j]=tb[j];
   FOR (i,com.nrgene) FOR (j,2) xb[com.ntime+i][j]=rgeneb[j];
   FOR (i,com.nrate)  FOR (j,2) xb[com.ntime+com.nrgene+i][j]=rateb[j];
   k=com.ntime+com.nrgene+com.nrate;
   for(i=0; i<com.npi*n31pi; i++) {
      xb[k][0]  =(com.model==T92?pb[0]:fb[0]);
      xb[k++][1]=(com.model==T92?pb[1]:fb[1]);
   }
   for (i=0;i<com.nalpha;i++,k++)  FOR (j,2) xb[k][j]=alphab[j];
   if (!com.fix_rho)   FOR (j,2) xb[np-1][j]=rhob[j];
   if (com.nparK) {
      FOR (i,com.ncatG-1) { xb[k][0]=rateb[0]; xb[k++][1]=rateb[1]; }
      if     (com.nparK==2) nf=com.ncatG-1;
      else if(com.nparK==3) nf=(com.ncatG-1)*(com.ncatG-1);
      else if(com.nparK==4) nf=(com.ncatG-1)*com.ncatG;
      FOR(i,nf) { xb[k][0]=fb[0]; xb[k++][1]=fb[1]; }
   }
   if(noisy>2 && np<20) {
      printf("\nBounds (np=%d):\n",np);
      FOR(i,np) printf(" %10.6f", xb[i][0]);  FPN(F0);
      FOR(i,np) printf(" %10.6f", xb[i][1]);  FPN(F0);
   }
   return(0);
}

int testx (double x[], int np)
{
/* This is now used for LS branch lengths estimation by nls2 only.
   To be removed.
*/
   int i;
   double tb[]={.4e-6,20};

   if (com.clock) error2 ("testx: clock not supported"); 
   FOR (i,com.ntime)  if (x[i]<tb[0] || x[i]>tb[1]) return (-1);
   return (0);
}

int PartialLikelihood (int inode, int igene)
{
   int n=com.ncode, i,j,k,h, ison, pos0=com.posG[igene],pos1=com.posG[igene+1];
   double t;

   FOR (i,nodes[inode].nson)
      if (nodes[nodes[inode].sons[i]].nson>0 && !_oldlkl[nodes[inode].sons[i]])
         PartialLikelihood (nodes[inode].sons[i], igene);
   fillxc (nodes[inode].lkl+pos0*n, (double)(inode>=com.ns), (pos1-pos0)*n);

   /* what are those two lines for? */
   if (com.cleandata && inode<com.ns) 
      for(h=pos0;h<pos1;h++) nodes[inode].lkl[h*n+com.z[inode][h]]=1;
   FOR (i, nodes[inode].nson) {
      t=nodes[ison=nodes[inode].sons[i]].branch*com.rgene[igene]*_rateSite;
      GetPMatBranch(PMat, NULL, t, ison);

      if (com.cleandata && nodes[ison].nson<1)       /* end node */
         for(h=pos0; h<pos1; h++)
            FOR(j,n) 
               nodes[inode].lkl[h*n+j]*=PMat[j*n+com.z[ison][h]];
      else {
         for(h=pos0; h<pos1; h++) 
            FOR(j,n) {
               for(k=0,t=0; k<n; k++)    /* t is used as temp */
                  t+=PMat[j*n+k]*nodes[ison].lkl[h*n+k];
               nodes[inode].lkl[h*n+j]*=t;
            }
      }
   }        /*  for (ison)  */
   if(_nnodeScale && _nodeScale[inode])  NodeScale(inode, pos0, pos1);
   return (0);
}


int PMatCijk (double P[], double t)
{
/* P(t)ij = SUM Cijk * exp{Root*t}
*/
   int i,j,k, n=4, nr=nR;
   double expt[4], pij;

   if (t<-.001 && noisy>3) 
      printf ("\nt = %.5f in PMatCijk", t);
   if (t<1e-200) { identity (P, n); return(0); }

   for (k=1; k<nr; k++) expt[k]=exp(t*Root[k]);
   FOR (i,n) FOR (j,n) {
      for (k=1,pij=Cijk[i*n*nr+j*nr+0]; k<nr; k++)
         pij+=Cijk[i*n*nr+j*nr+k]*expt[k];
      P[i*n+j] = (pij>0?pij:0);
   }
   return (0);
}


#ifdef UNDEFINED

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

         SiteCat[h]=it+min2(bit1+1,nsep)*11;
         if (n1==1 && bit1<nsep) {
            SiteCat[h]-=11;
            SiteCat[h]+=(b[0]*4+b[ns-1-bit1]-b[0]-(b[0]<=b[ns-1-bit1]));
         }
         break;
      }
      if (fout) {
         FOR (j, ns) fprintf (fout, "%1c", BASEs[b[j]]);
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
      error2 ("Pexp()");
   SetParameters (x);
   FOR (j, tree.nbranch) 
      PMatCijk (Pt[j], nodes[tree.branches[j][1]].branch);

   for (h=0; h<(1<<2*com.ns); h++) {
      if (SiteCat[h] == 0) continue;
      for (j=0,it=h; j<com.ns; nodeb[com.ns-1-j]=it&3,it>>=2,j++) ;
      for (isum=0,fh=0; isum<nsum; isum++) {
         for (j=0,it=isum; j<tree.nbranch-com.ns+1; j++)
            { nodeb[com.ns+j]=it%4; it/=4; }
         for (j=0,y=com.pi[nodeb[tree.root]]; j<tree.nbranch; j++) 
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

   SiteCat=(int*)malloc((1<<2*com.ns)* sizeof(int));
   if (SiteCat == NULL)  error2 ("oom");
   CollapsSite (F0, nsep, com.ns, &ncat, SiteCat);
   fprintf (fout, "\n\nAppr. test of model.. ncat%6d  nsep%6d\n", ncat,nsep);

   nobs = pexp+ncat;
   zero (pexp, 2*ncat);
    /* nobs */
   FOR (h, com.npatt) {
      for (j=0,it=0; j<com.ns; j++) it = it*4+(com.z[j][h]-1);
      nobs[SiteCat[it]] += com.fpatt[h];
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


#endif


int OldDistributions (int inode, double oldfreq[])
{
/* reconstruct nucleotide frequencies at and down inode
   for nonhomogeneous models com.nhomo==3 or 4.
   oldfreq[tree.nnode*4]
*/
   int i, n=4;
   double kappa=com.kappa;

   if (com.alpha || com.model>REV) {
      puts("OldDistributions() does not run when alpha > 0 or model >= TN93");
      return(-1);
   }
   if (inode==tree.root) {
      xtoy (nodes[inode].pi, oldfreq+inode*n, n);
   }
   else {
      if(!com.fix_kappa) kappa=nodes[inode].kappa;
      EigenTN93 (com.model, kappa, kappa, nodes[inode].pi, &nR, Root, Cijk);
      PMatCijk (PMat, nodes[inode].branch);
      matby (oldfreq+nodes[inode].father*n, PMat, oldfreq+inode*n, 1, n,n);
   }
   FOR (i,nodes[inode].nson)
      OldDistributions (nodes[inode].sons[i], oldfreq);
   return (0);
}


void ReadPatternFreq (FILE* fout, char* pexpf)
{
   FILE *fpexp;
   int h, j, ch;
   double t;

   printf ("\n\nRead site-pattern frequencies from %s\n", pexpf);
   fprintf (fout, "\nSite-pattern frequencies are read from %s\n", pexpf);
   fpexp=gfopen(pexpf,"r");
   fscanf (fpexp, "%d%d", &com.ns, &com.npatt);
   printf ("%3d species, %d site patterns\n", com.ns, com.npatt);

   com.ls=com.npatt;
   com.fpatt=(double*) malloc (com.npatt*sizeof(double));
   FOR (j, com.ns) com.z[j]=(char*) malloc (com.ls*sizeof(char));
   if (com.fpatt==NULL || com.z[com.ns-1]==NULL) error2 ("oom in ReadPatternFreq");

   com.ngene=1; com.lgene[0]=com.ls; com.posG[0]=0; com.posG[1]=com.ls;
   com.rgene[0]=1;
   zero(com.pi,4);   zero(com.fpatt,com.npatt);
   for (h=0; h<com.npatt; h++) {
      for (j=0; j<com.ns; j++) {
         ch=fgetc(fpexp);
         while (!isalpha(ch)) ch=fgetc(fpexp);
         com.z[j][h]=(char)CodeChara((char)ch, 0);
         if(com.z[j][h]<0) 
            error2("does not deal with ambiguity characters in patterns");
      }
      if (fscanf(fpexp,"%lf",&t) != 1) error2 ("ReadPatternFreq");
      com.fpatt[h] += t;
      for (j=0; j<com.ns; j++) com.pi[com.z[j][h]]+=t;
   }
   for (j=0,t=0; j<com.npatt; j++) t+=com.fpatt[j];
   printf ("\nSUM freq = %.6f = 1?\n\n", t);
   for (j=0; j<4; j++) com.pi[j]/=(t*(double)com.ns);
   fprintf(fout, "\nBase frequences"); 
   for (j=0; j<4; j++) fprintf(fout, "  %c: %9.5f", BASEs[j],com.pi[j]);
   FPN(fout);   FPN(fout);
   if(com.model==JC69 || com.model==K80)
     for(j=0;j<4;j++) com.pi[j]=.25;
   else if (com.model==REV || com.model==REVu) 
      getpi_sqrt (com.pi, com.pi_sqrt, com.ncode, &com.npi0);

   xtoy (com.pi, com.piG[0], 4);
   fclose (fpexp);
}



/* problems and notes

   (1) AdG model: generation of com.MK[K*K] is not independent
       of com.alpha or DiscreteGamma().

non-homogeneous process models:
  nhomo            fix_kappa         models
 
  0 (1 pi given)   0 (solve 1 kappa)   JC69, K80, F81, F84, HKY85,
                                       REV(5), UNREST(11)
                   1 (1 kappa given)   K80, F84, HKY85
 
  1 (solve 1 pi)   0 (as above)        F84, HKY85, REV(5)
                   1                   F81(0), F84, HKY85        
 
  2 (b kappa's)    ?                   K80, F84, HKY85

  3 (ns+2 pi)      0 (solve 1 kappa)   F84 & HKY85  
                   1 (nbranch kappa)   F84 & HKY85  
 
  4 (nnode pi)     0,1  (as above)


space-time process models:
  nparK     fix_alpha       fix_rho        parameters 
  0         (0,1)          (0,1)           alpha & rho for AdG

  1         set to 1       set to 1        rK[]        (freqK=1/K) (K-1)
  2         set to 1       set to 1        rK[] & freqK[]          2(K-1)
  3         set to 1       set to 1        rK[] & MK[] (freqK=1/K) K(K-1)
  4         set to 1       set to 1        rK[] & MK[]             K*K-1


Local clock models
parameters under local clock (com.clock=2)
   com.ntime = (#ancestral nodes) - 1 + (#rates) - 1
Parameters include (ns-1) node times t[] and rates for branches.
  x[0]: t0*r0
  x[1..(nid-1)]: ti/t0 ancestral node times expressed as ratios
  x[ns-1 .. ntime-1]: rates for branches

*/



#ifdef SUBDATA

int GetSubTreeN (int hasbranch, char keep[], int rennodes)
{
/* This extracts the sub tree when some species are kept (keep[i]=1) while 
   others are removed.  Branch lengths are preserved by adding them up 
   when nodes are removed.
   The algorithm modifies tree.branches[] to change the tree topology and 
   nodes[].branch and then updates nodes[].  If (rennodes), the nodes and 
   species are renumbered consecutively, in which case OutaTreeN(?,1,?) should 
   not be used as the pointers to species names are not updated.  
   If (!rennodes), the old nodes and species numberes are used, in which case 
   you can use both OutaTreeN(?,1,?) (which gives the correct species) or 
   OutaTreeN(?,0,?) (which uses the old species numbers).
   If the root has 2 sons (&& !com.clock), the tree is derooted.
   The program checks that species to be removed are removed and species 
   to be kept are kept.  Branch lengths are checked using small trees manually.

   Returns -1 if error (e.g., if <2 seqs are to be kept).
   Does not work if a current seq is ancestral and is removed.
   keep[] keeps species, and bkeep[] keeps branches.
   December 1998, Ziheng Yang
*/
   char bkeep[NBRANCH];
   int ib,j,k,father,inode,sib,nson0=nodes[tree.root].nson,newn,nmark[NNODE];

   if(nson0<2) puts("GetSubTreeN: strange root.");
   FOR(ib,tree.nbranch) bkeep[ib]=1;
   for(j=0,k=0;j<com.ns;j++) if(keep[j]) k++;
   if(k==com.ns) return(0); else if(k<2) return(-1);

   FOR(ib,tree.nbranch) {
      if(!bkeep[ib]) continue;
      inode=tree.branches[ib][1]; father=tree.branches[ib][0]; 
      if(inode>=com.ns || keep[inode]) continue;
      bkeep[ib]=0;

      nodes[father].nson--;  /* remove inode as a son */
      for(j=0,k=0;j<nodes[father].nson+1;j++) {
         if(k) nodes[father].sons[j-1]=nodes[father].sons[j];
         else if(nodes[father].sons[j]==inode)  k=1;
      }
      if(nodes[father].nson>=2)  continue;
      if(nodes[father].nson==0) {
         printf("\nib %d father %d\n", ib+1,father+1);
         error2("This should not be possible.");
      }
      else if(nodes[father].nson==1) {
         sib=nodes[father].sons[0];
         if(sib==inode) error2("This should not be possible 1");

         if(father==tree.root) {  /* remove old root */
            tree.root=sib; bkeep[nodes[sib].ibranch]=0;
         }
         else {  /* remove father node and father branch & reset sib branch */
            bkeep[nodes[father].ibranch]=0;
            j=nodes[sib].father=nodes[father].father;
            tree.branches[nodes[sib].ibranch][0]=j;
            for(k=0; k<nodes[j].nson; k++) 
               if(nodes[j].sons[k]==father) { nodes[j].sons[k]=sib; break; }
            if(hasbranch) nodes[sib].branch+=nodes[father].branch;
         }
      }
   }  /* for (ib) */
   if(nodes[tree.root].nson<2) error2("strange");
   else if(nodes[tree.root].nson==2 && !com.clock && 
      nodes[tree.root].nson<nson0) {
      /* deroot, k is new root */
      j=nodes[tree.root].sons[0]; k=nodes[tree.root].sons[1];
      if(j>k) {k=j; j=nodes[tree.root].sons[1]; }
      if(k<com.ns) if(noisy) puts("only two seqs left.");

      tree.root=k; 
      bkeep[nodes[k].ibranch]=0; 
      tree.branches[nodes[j].ibranch][0]=k;
      tree.branches[nodes[j].ibranch][1]=j;
      nodes[j].branch+=nodes[k].branch;
   }

   /* error2 checking, k=err */
   for(ib=0,k=0; ib<tree.nbranch; ib++) /* deleted species j? */
      if(bkeep[ib] && (j=tree.branches[ib][1])<com.ns && !keep[j]) k=1;
   if(k) { printf("\n\aDeleted species %d still in tree\n",j+1); getchar(); }
   FOR(k,com.ns) nmark[k]=0;    /* kept species? nmark[] for species in tree */
   FOR(ib,tree.nbranch)  if(bkeep[ib]) 
      FOR(j,2) if((k=tree.branches[ib][j])<com.ns) nmark[k]=1;
   FOR(j,com.ns)  if((int)keep[j]!=nmark[j])
      { printf("\a\nsp. %d should be kept?\n", j+1); getchar(); }

   if(rennodes) {   /* renumber nodes, using nmark[] to mark nodes */
      FOR(j,tree.nnode) nmark[j]=-1; /* absent */
      FOR(ib,tree.nbranch) if(bkeep[ib]) FOR(j,2) nmark[tree.branches[ib][j]]=1;
      for(j=0,newn=0; j<tree.nnode; j++) if(nmark[j]==1) nmark[j]=newn++;
      for(ib=0,newn=0; ib<tree.nbranch; ib++)  if(bkeep[ib]) {
         FOR(j,2) tree.branches[newn][j]=nmark[tree.branches[ib][j]];
         inode=tree.branches[ib][1];
         newn++;   /* new nbranch */
      }
      tree.nbranch=newn;  tree.root=nmark[tree.root];
      FOR(j,tree.nnode)  nodes[nmark[j]].branch=nodes[j].branch;
      BranchToNode();
   }
   else {
      FOR (j,tree.nnode) ClearNode (j);
      FOR (ib,tree.nbranch) {  /* BranchToNode(); */
         if(!bkeep[ib]) continue; 
         j=tree.branches[ib][0]; k=tree.branches[ib][1];
         nodes[j].sons[nodes[j].nson++]=k;
         nodes[k].father=j;  nodes[k].ibranch=ib;
      }
   }
   return (0);
}



struct TREEB tree0;
struct TREEN nodes0[2*NS-1];

int SubData(int argc, char *argv[])
{
   FILE *fin, *fout;
   char keep[NS], outfile[32]="s", infile[32]="mc.paup";
   int i,j, ntree, pauptree=0, format=1;
   int nkept=4, isample,nsample=5;

   puts("\n\nUsage: SubData <infile> <outfile> <#seqs> <#sample>\n");
   strcpy(com.seqf,infile);
   if(argc>1)  strcpy(com.seqf,argv[1]);
   if(argc>2)  strcpy(outfile, argv[2]);
   if(argc>3)  sscanf(argv[3],"%d",&nkept);
   if(argc>4)  sscanf(argv[4],"%d",&nsample);
   printf("\nreading %s\twriting %s\n", com.seqf,outfile);
   printf("%d samples of size %d\n", nsample,nkept);

   fin=gfopen(com.seqf,"r");
   fout=fopen(outfile,"w");
   ReadSeq (NULL, fin, 0);
   if(nkept>com.ns || nkept<2) error2("nkept");
   GetTreeFileType(fin, &ntree, &pauptree, 1);
   if (!pauptree) error2("Is this a paup data & tree file?");
   PaupTreeRubbish(fin);
   if (ReadaTreeN (fin, &i, 1)) error2 ("err reading tree.");

   fprintf(fout,"#NEXUS\n[%d small samples of size %d, from\n",nsample,nkept);
   OutaTreeN(fout,1,1); fputs("\n]\n\n",fout);
   appendfile(fout,"paupstart");
   tree0=tree;  memcpy(nodes0,nodes,sizeof(nodes));

   for (isample=0; isample<nsample; isample++) {
      printf("\ndata # %d: ", isample+1); 
      tree=tree0;  memcpy(nodes,nodes0,sizeof(nodes));
      FOR(i,com.ns) keep[i]=0;
      for(i=0; i<nkept; ) 
         { j=(int)(com.ns*rndu()); if(!keep[j]) { keep[j]=1; i++; } }

      fprintf(fout,"\n\n[data # %d: ", isample+1);  
      fputs("Species kept: ",fout);
      FOR(i,com.ns)  if(keep[i]) fprintf(fout,"%2d ",i+1);
      fputs("]\n",fout);
      printSeqs(fout,NULL,keep,1);

      /* do not renumber nodes, use GetSubTreeN(1,keep,0) &  OutaTreeN(?,1,?) */
      GetSubTreeN(1,keep,0);

      fputs("\nbegin trees;\n  tree smalltree=[&U] ", fout);   
      OutaTreeN(fout,1,1); fputs("\nend;\n", fout);

/*      OutaTreeN(F0,0,0);
*/
      appendfile(fout, "paupblock");

   }
   FPN(F0);
   return(0);
}
#endif
