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
#define NGENE         100
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
extern clock_t clock_start, clock_cur;

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
void ReadPatternFreq (char* pexpf);

struct CommonInfo {
   char *z[NS], *spname[NS], seqf[96],outf[96],treef[96];
   int seqtype, ns,ls,ngene,posG[NGENE+1],lgene[NGENE], *pose,npatt;
   int cleandata;
   int np,ntime,nrgene,nrate,nalpha,npi,nhomo,ncatG,ncode,Mgene,sspace,slkl1;
   int fix_kappa, fix_rgene, fix_alpha, fix_rho, nparK, fix_branch;
   int clock, model, getSE, runmode, print,verbose, ndata, icode,coding;
   int method, npi0;
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
   int father, nson, sons[NS], ibranch, label;
   double branch, divtime, kappa, pi[4], *lkl;
}  *nodes;

int nR=4, LASTROUND;
double PMat[16], Cijk[64], Root[4];
complex cU[16], cV[16], cRoot[16];
char  StepMatrix[16];

double xcom[NP-NBRANCH];

int lklSiteClass=0;  /* is lkl memory allocated for each class? */
char _oldlkl[NNODE];     /* update lkl for nodes? to save computation (0 yes; 1 no) */
int _nnodeScale=0;
char _nodeScale[NNODE];   /* nScale[ns-1] for interior nodes */
double *_nodeScaleF=NULL;  /* nScaleF[npatt] for scale factors */

FILE *frub, *flnf, *frst, *frst1, *fin=NULL;
char *ratef="rates";
char *models[]={"JC69","K80","F81","F84","HKY85","TN93","REV","UNREST", "REVu","UNRESTu"};
enum {JC69, K80, F81, F84, HKY85, TN93, REV, UNREST, REVu, UNRESTu} MODELS;

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
   int i, slkl0=0,s2=0, READPATTERNFREQ=0, idata;
   char rstf[96]="rst", ctlf[96]="baseml.ctl", timestr[64];
   char *Mgenestr[]={"diff. rate", "separate data", "diff. rate & pi", 
                     "diff. rate & kappa", "diff. rate & pi & kappa"};
   char *clockstr[]={"", "global clock", "local clock", "clock with dated seqs"};
   time_t tm;

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
   com.print=0;     com.verbose=1;
   com.method=0;    com.space=NULL;
   SetSeed ((int)time(&tm));

   frub=fopen("rub","w");  frst=fopen(rstf,"w"); frst1=fopen("rst1","w");
   if(frub==NULL || frst==NULL || frst1==NULL) error2("file open error");

   if (argc>1)  strcpy(ctlf, argv[1]); 
   GetOptions (ctlf);

   if ((fout=fopen (com.outf, "w"))==NULL) 
      error2("can't open the control file .ctl\n");
   if((fseq=fopen (com.seqf,"r"))==NULL)  {
      printf ("\n\nSequence file %s not found!\n", com.seqf);
      exit (-1);
   }
   if((fpair[0]=(FILE*)fopen(pairfs[0],"w"))==NULL) error2("2base.t file open error");

   for (idata=0; idata<com.ndata; idata++) {
      if (com.ndata>1) {
         printf("\nData set %4d\n", idata+1);
         fprintf(fout, "\n\nData set %4d\n", idata+1);
         fprintf(frst, "\t%4d", idata+1);
      }
      if (idata)  GetOptions (ctlf);
      if (READPATTERNFREQ) {
         if (com.rho || com.Mgene) error2("model");
         ReadPatternFreq ("pexp.in");
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

      /* if(com.ngene>1 && com.Mgene==1)  printSeqsMgenes (); */
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
         error2("nhomo does not work well with option G");

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
      com.sspace=max2(800000,com.ls*(int)(com.ns*sizeof(char)+sizeof(int)));
      if((com.space=(double*)realloc(com.space,com.sspace))==NULL) 
         error2("oom space");

      if (!READPATTERNFREQ) {
         i=com.ns*(com.ns-1)/2;
         SeqDistance=(double*)realloc(SeqDistance,i*sizeof(double));
         ancestor=(int*)realloc(ancestor, i*sizeof(int));
         if(SeqDistance==NULL||ancestor==NULL) error2("oom distance&ancestor");
      }

      if (!READPATTERNFREQ) Initialize (fout, com.space);
      if (com.Mgene==3) FOR (i,com.ngene) xtoy(com.pi,com.piG[i],4);
/*
      if (com.model==JC69 && !com.print && !READPATTERNFREQ) 
         PatternJC69like(fout);
*/
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

      if(!READPATTERNFREQ) 
         DistanceMatNuc(fout,fpair[0],com.model,com.alpha);

      if (com.Mgene==1)        MultipleGenes (fout, fpair, com.space);
      else if (com.runmode==0) Forestry (fout);
      else if (com.runmode==3) StepwiseAddition (fout, com.space);
      else if (com.runmode>=4) Perturbation(fout, (com.runmode==4), com.space);
      else                  fprintf(frst,"%2d",StarDecomposition(fout,com.space));
      
      FPN(frst);  if((idata+1)%10==0) fflush(frst);
      free(nodes);
      printf("\nTime used: %s\n", printtime(timestr));
      starttime();
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
   int pauptree=0, btree=0, length_label;
   double x[NP], lnL,lnL0=0,lnLbest=0,e=1e-8, nchange=-1;
   double xb[NP][2], tl=0, *var=NULL;
   FILE *ftree, *finbasemlg=NULL, *frate=NULL;

   if(com.clock==3) GetSeqTimes();

   if ((ftree=fopen(com.treef,"r"))==NULL) {
      printf("\ntree file %s not found.\n", com.treef);
      exit(-1);
   }
   GetTreeFileType(ftree,&ntree,&pauptree,0);
   if(com.alpha)
      if((frate=(FILE*)fopen(ratef,"w"))==NULL) error2("rate file open err");
   if (com.alpha && com.rho==0 && com.nhomo==0 && com.nparK==0) {
      inbasemlg=1;  if((finbasemlg=fopen("in.basemlg","w"))==NULL) error2("file err"); 
   }
   if((flnf=fopen("lnf","w+"))==NULL) error2("lnf file open error2");
   fprintf(flnf,"%6d %6d %6d\n", ntree, com.ls, com.npatt);

   if(!com.cleandata) InitPartialLikelihood ();
   for(itree=0; ntree==-1||itree<ntree; itree++,iteration=1) {
      if((pauptree && PaupTreeRubbish(ftree)) || 
         ReadaTreeN(ftree,&length_label,1))
            { printf("\nerr or end of tree file.\n"); break; }
      if(noisy) printf ("\nTREE # %2d: ", itree+1);
      fprintf (fout,"\nTREE # %2d:  ", itree+1);
      fprintf (frub,"\n\nTREE # %2d\n", itree+1);
      if(com.print) fprintf (frst,"\n\nTREE # %2d\n", itree+1);
      fprintf (flnf,"\n\n%2d\n", itree+1);

      LASTROUND=0;
      if (times++==0 && (length_label==1 || length_label==3)) {
         if(com.clock) puts("\nBranch lengths in tree are ignored");
         else {
            puts("\ntree has branch lengths.\n0:ignore?\n1:initials?\n2:fixed?");
            com.fix_branch=1;
            scanf("%d",&com.fix_branch);
			
            if(com.fix_branch==1) {
               FOR(i,tree.nnode) 
                  if((x[nodes[i].ibranch]=nodes[i].branch)<0)
                     error2("blength<0 in tree");
            }
         }
      }

      if(com.cleandata) nchange=MPScore (com.space);
      if(noisy&&com.ns<99) 
         { OutaTreeN(F0,0,0); printf(" MP score: %.2f\n",nchange);}
      OutaTreeN(fout,0,0);  fprintf(fout,"  MP score: %.2f",nchange);
      if(!com.clock && com.model<=REV && com.nhomo<=2 && nodes[tree.root].nson<=2 && com.ns>2){
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
      if(com.fix_branch<2) { OutaTreeB (fout);  FPN (fout); }
      LASTROUND=1;
      if(com.npi || com.nparK>=2) TransformxBack(x);
      if(com.clock) {  /* this applies to all clock models (clock=1,2,3) */
         for(i=com.ns; i<tree.nnode; i++) 
            if(i!=tree.root) x[i-com.ns]=nodes[i].divtime;
      }
/*
for(i=com.ntime;i<np;i++) fprintf(frst,"%9.5f",x[i]); 
fprintf(frst,"%12.4f %4d ",-lnL,com.npatt); 
*/

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
         FOR (i,np) var[i*np+i]=(var[i*np+i]>0.?sqrt(var[i*np+i]):-1);
         FOR (i,np) fprintf(fout," %8.5f", var[i*np+i]);  FPN (fout);
         /* if (com.getSE==2) matout2(fout,var,np,np,12,7); */
      }

/*    if(com.clock) SetBranch(x); */

      for(i=0,tl=0;i<tree.nnode;i++) 
         if(i!=tree.root) tl+=nodes[i].branch;
      fprintf(fout,"\ntree length = %9.5f%s\n",tl,com.ngene>1?" (1st gene)":"");
      /* FPN(fout); OutaTreeN(fout,0,1); FPN(fout); */
      FPN(fout); OutaTreeN(fout,1,0); FPN(fout);
      FPN(fout); OutaTreeN(fout,1,1); FPN(fout);

      if(com.np-com.ntime||com.clock)  DetailOutput(fout,x,var);
      if (status) {
         printf ("convergence?\n");  fprintf (fout,"check convergence..\n");
      }
      if (itree==0) 
         for (i=0; i<np-com.ntime; i++) xcom[i]=x[com.ntime+i];
      else if (!j)  
         for (i=0; i<np-com.ntime; i++) xcom[i]=xcom[i]*.8+x[com.ntime+i]*0.2;

      com.print-=9;  com.plfun(x, np);  com.print+=9;
      if (com.print) {
         if (com.plfun!=lfun)  lfunRates (frate, x, np);
         if (com.nhomo==0 && com.nparK==0 && com.model<=REV && com.rho==0)
            AncestralSeqs (frst, x);
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


void DetailOutput(FILE *fout, double x[], double var[])
{
   int i,j,k=com.ntime, nr[]={0, 1, 0, 1, 1, 2, 5, 11};
   double Qfactor,*p=com.pi, t=0, k1,k2, S,V, Y=p[0]+p[1],R=p[2]+p[3],tnode;
   double *Qrate;

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
   if (com.nhomo==1 || com.nhomo>2) {
      fprintf (fout, "kappa under %s:", models[com.model]);
      FOR (i,com.nrate) fprintf (fout, " %8.5f", x[k++]);  FPN(fout);

      fprintf (fout, "base frequencies (%d set):", com.npi);
      FOR (i, com.npi) {
         if (com.npi>1) fprintf (fout, "\n  Set #%d:", i+1);
         FOR (j,3) fprintf (fout, " %8.5f", x[k+3*i+j]);
         fprintf (fout, "%11.5f", 1-sum(x+k+3*i, 3));
      }
      FPN (fout);
      k+=3*com.npi;
   }
   else if (!com.fix_kappa) {
      fprintf(fout,"\nParameters in the rate matrix (%s) (for definition, see Yang 1994 J Mol Evol 39:105-111):\n",
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
               EigenQREVbase(fout, Qrate, p, &nR, Root, Cijk);
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
      if (com.model==REV || com.model==REVu) 
         EigenQREVbase(fout, x+com.ntime+com.nrgene, com.pi, &nR, Root, Cijk);
      else if (com.model==UNREST || com.model==UNRESTu) 
         EigenQunrest (fout, x+com.ntime+com.nrgene,com.pi,&nR,cRoot,cU,cV);
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
   int i, nopt=27, lline=2048; 
   char line[2048], *pline, opt[20], comment='*';
   char *optstr[]={"seqfile","outfile","treefile","noisy",
        "cleandata", "verbose","runmode", "method",
        "clock","fix_rgene","Mgene","nhomo","getSE","RateAncestor",
        "model","fix_kappa","kappa","fix_alpha","alpha","Malpha","ncatG", 
        "fix_rho","rho","nparK", "ndata", "Small_Diff","icode"};
   double t;
   FILE *fctl;

   if((fctl=fopen(ctlf,"r"))==NULL) error2("\nctl file open error.\n");
   printf ("\n\nReading options from %s..\n", ctlf);
   for (;;) {
      if (fgets (line, lline, fctl) == NULL) break;
      for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
         if (isalnum(line[i]))  { t=1; break; }
         else if (line[i]==comment) break;
      if (t==0) continue;
      sscanf (line, "%s%*s%lf", opt, &t);
      if ((pline=strstr(line, "="))==NULL)
         error2("err: option file. add space around the equal sign =");

      for (i=0; i<nopt; i++) {
         if (strncmp(opt, optstr[i], 8)==0)  {
            if (noisy>=9)
               printf ("\n%3d %15s | %-20s %6.2f", i+1,optstr[i],opt,t);
            switch (i) {
               case ( 0): sscanf(pline+1, "%s", com.seqf);    break;
               case ( 1): sscanf(pline+1, "%s", com.outf);    break;
               case ( 2): sscanf(pline+1, "%s", com.treef);    break;
               case ( 3): noisy=(int)t;           break;
               case ( 4): com.cleandata=(int)t;   break;
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
                  if(com.model>REV) GetStepMatrix(line);  break;
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
   else if (com.nhomo>2 && com.model!=F84 && com.model!=HKY85) error2("nhomo");
   else if (com.nhomo>2 && com.method) error2("nhomo & method.");
   else
      if (com.nhomo==1 && com.model!=F81 && com.model!=F84 && com.model!=HKY85 
                       && com.model!=REV && com.model!=REVu)
         error2("nhomo=1 and model");

   if (com.nhomo>1 && com.runmode>0)  error2("nhomo incompatible with runmode");
   if (com.clock && com.runmode>2)  error2("runmode & clock?");
   if (com.runmode==3 && (com.npi || com.nparK))
      error2("runmode etc.");

   if (com.model==JC69 || com.model==K80 || com.model==UNREST)
      if (com.nhomo!=2)  com.nhomo=0;
   if (com.model==JC69 || com.model==F81) { com.fix_kappa=1; com.kappa=1; }
   if (com.model==TN93 || com.model==REV)   com.fix_kappa=0;
   if (com.nparK==3) {
      puts("\n\nnparK==3, double stochastic, may not work.  Use nparK=4?\n");
      getchar();
   }
   if(com.runmode!=-2) fin=fopen("in.baseml","r");

   return (0);
}



int GetInitials (double x[], int *fromfile)
{
   int i,j,k, K=com.ncatG;
   int slkl1_new=(tree.nnode-com.ns)*com.ncode*com.npatt*com.ncatG*sizeof(double);
   double t=-1;

   NFunCall=NPMatUVRoot=NEigenQ=0;
   com.plfun=lfunAdG;
   if (com.alpha==0 && com.nparK==0)  com.plfun=lfun;
   else if ((com.alpha && com.rho==0) || com.nparK==1 || com.nparK==2)
      com.plfun=lfundG;

   if(com.method && com.fix_branch!=2 && com.plfun==lfundG) {
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
         if (com.model<=HKY85)     com.nrate=1;
         else if (com.model==TN93) com.nrate=2;
         else                      com.nrate=(com.model==REV?5:11);
         if (com.Mgene>=3)         com.nrate*=com.ngene;
      }
   }
   switch (com.nhomo) {
   case (0): com.npi=0;           break;   /* given 1 pi */
   case (1): com.npi=1;           break;   /* solve 1 pi */
   case (2): com.npi=0;  com.nrate=tree.nbranch;  break;  /* b kappa's */
   case (3): com.npi=com.ns+(tree.root>=com.ns)+(tree.nnode>com.ns+1);  
             com.nrate=(com.fix_kappa?1:tree.nbranch);  break;   /* ns+2 pi */
   case (4): com.npi=tree.nnode;  com.nrate=(com.fix_kappa?1:tree.nbranch);
                                  break;   /* nnode pi   */
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
   else if(com.fix_branch==0) {
      for(j=0,t=0,k=com.ns*(com.ns-1)/2; j<k; j++)  t+=SeqDistance[j]/k;
      FOR (j,com.ntime) x[j]=min2(t,5)*rndu();
      if(com.ns<50) LSDistance (&t, x, testx);
/*
matout (F0, x, 1, com.ntime);
printf("SS=%.8f\n", t);
*/
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

   FOR (i,com.npi) fillxc(x+com.ntime+com.nrgene+com.nrate+3*i,0, 3);

   com.np = k = com.ntime+com.nrgene+com.nrate+com.npi*3;

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
/* This is called before lfun etc to set parameters in com. etc and
   initialize U, V, Root etc,
   For nhomo models (nhomo=1,3,4) 
      x[] has freqs pi's if (LASTROUND)
              exp(pi)/(1+SUM(exp[pi])) if otherwise
   This function does not change x[].

*/
   int i, j, k, K=com.ncatG, status=0;
   double k1=com.kappa, k2=com.kappa, t, space[NCATG*(NCATG+1)];

   if(com.fix_branch<2) SetBranch(x);
   if(com.np<=com.ntime) return(0);
   FOR (i, com.nrgene) com.rgene[i+1]=x[com.ntime+i];
   if (!com.fix_kappa && com.model<=TN93) {
       com.kappa=k1=k2=x[com.ntime+com.nrgene];
       if (com.model==TN93) k2=x[com.ntime+com.nrgene+1];
   }
   if (com.nhomo==1) {
      k=com.ntime+com.nrgene+com.nrate;

      if (!LASTROUND) f_and_x(x+k,com.pi,4,0,0);
      else            xtoy (x+k, com.pi, 3);
      com.pi[3]=1-sum(com.pi,3);
      if (com.model<=TN93) EigenTN93(com.model,k1,k2,com.pi,&nR,Root,Cijk);
   }
   else if (com.nhomo==2) 
      for (i=0,k=com.ntime+com.nrgene; i<tree.nbranch; i++)
         nodes[tree.branches[i][1]].kappa=x[k+i];
 
   if (com.model<=TN93 && com.nhomo<=1 && com.Mgene<=1)
      RootTN93 (com.model, k1, k2, com.pi, &t, Root);
   else if (com.nhomo>2) {
      for (i=0,k=com.ntime+com.nrgene; i<tree.nbranch; i++)
         nodes[tree.branches[i][1]].kappa=(com.fix_kappa?x[k]:x[k+i]);
      k+=com.nrate;
      if (com.nhomo==3) {
         FOR (i,com.ns) 
            if (!LASTROUND) f_and_x(x+k+i*3,nodes[i].pi,4,0,0);
            else            xtoy  (x+k+i*3, nodes[i].pi, 3);
         for (i=com.ns; i<tree.nnode; i++)
            if (!LASTROUND) f_and_x(x+k+com.ns*3,nodes[i].pi,4,0,0);
            else            xtoy  (x+k+com.ns*3, nodes[i].pi, 3);

         if (tree.root>=com.ns && tree.nnode>com.ns+1)
            if (!LASTROUND) f_and_x(x+k+(com.ns+1)*3,nodes[tree.root].pi,4,0,0);
            else            xtoy   (x+k+(com.ns+1)*3, nodes[tree.root].pi,3);
      }
      else { /*  if(com.nhomo==4) */
         FOR (i,tree.nnode)
            if (!LASTROUND) f_and_x(x+k+i*3,nodes[i].pi,4,0,0);
            else            xtoy   (x+k+i*3, nodes[i].pi, 3);
      }

      FOR(i,tree.nnode) nodes[i].pi[3]=1-sum(nodes[i].pi,3);
      xtoy(nodes[tree.root].pi, com.pi, 4);
   }
   else if ((com.model==REV || com.model==REVu) && com.Mgene<=1)
      EigenQREVbase (NULL, x+com.ntime+com.nrgene, com.pi, &nR, Root, Cijk);
   else if ((com.model==UNREST || com.model==UNRESTu) && com.Mgene<=1)
      EigenQunrest (NULL, x+com.ntime+com.nrgene,com.pi,&nR,cRoot,cU,cV);

   if (com.nparK==0 && (com.alpha==0 || com.fix_alpha*com.fix_rho==1))
      return(status);
   if (com.nalpha>1) return (status);
   k = com.ntime+com.nrate+com.nrgene+com.npi*3;
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
   int j,k=com.nrgene+(com.Mgene>=3)*igene*nr[com.model];
   double ka1=xcom[k], ka2=(com.model==TN93?xcom[k+1]:-1);

   if(com.Mgene==2 && com.fix_kappa) ka1=ka2=com.kappa;

   if (_pi) {
      xtoy(com.piG[igene], com.pi, 4);
      if(com.model==REV || com.model==REVu)
      for(j=0,com.npi0=0; j<com.ncode; j++)
         if(com.pi[j]) com.pi_sqrt[com.npi0++]=sqrt(com.pi[j]);
      com.npi0=com.ncode-com.npi0;
   }
   if (_UVRoot) {
      if (com.model==K80) com.kappa=ka1;
      else if (com.model<=TN93) 
         EigenTN93(com.model,ka1,ka2,com.pi,&nR,Root,Cijk);
      else if (com.model==REV)
         EigenQREVbase(NULL,xcom+k,com.pi,&nR,Root,Cijk);
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
   int i,j, k=com.ntime+com.nrgene+com.nrate, nf=0;
   double tb[]={.4e-6, 99}, rgeneb[]={1e-4,999}, rateb[]={1e-5,999};
   double alphab[]={.04, 99}, rhob[]={-0.2, 0.99}, pb[]={.00001,.999999};
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
   for(i=0,k=com.ntime+com.nrgene+com.nrate; i<com.npi*3; i++) {
      xb[k][0]=fb[0];  xb[k++][1]=fb[1];
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
   if (com.cleandata && inode<com.ns) 
      for(h=pos0;h<pos1;h++) nodes[inode].lkl[h*n+com.z[inode][h]]=1;

   FOR (i, nodes[inode].nson) {
      t=nodes[ison=nodes[inode].sons[i]].branch*com.rgene[igene]*_rateSite;
      GetPMatBranch(PMat, NULL, t, ison);

      if (com.cleandata && nodes[ison].nson<1)       /* end node */
         for(h=pos0; h<pos1; h++)
            FOR(j,n) nodes[inode].lkl[h*n+j]*=PMat[j*n+com.z[ison][h]];
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

   if (t<-2e-4 && noisy>3) 
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

void ReadPatternFreq (char* pexpf)
{
   FILE *fpexp;
   int h, j, ch;
   double sum;

   printf ("\n\nRead site-pattern frequencies from %s\n", pexpf);
   if ((fpexp=fopen(pexpf,"r")) == NULL) error2 ("pexp file open err");
   fscanf (fpexp, "%d%d", &com.ns, &com.npatt);
   printf ("%3d species, %d site patterns\n", com.ns, com.npatt);

   if (com.model>1) {
      error2("Should calculate base frequencies\a\n");
   }
   else 
      fillxc(com.pi, .25, 4);  
   xtoy (com.pi, com.piG[0], 4);

   if(com.npatt>(1<<(2*com.ns))) error2("npatt");
   com.ls=com.npatt;
   com.fpatt=(double*) malloc (com.npatt*sizeof(double));
   FOR (j, com.ns) com.z[j]=(char*) malloc ((com.ls+1)*sizeof(char));
   if (com.fpatt==NULL || com.z[com.ns-1]==NULL) error2 ("oom");
   FOR (j,com.ns) com.z[j][com.ls]='\0';

   com.ngene=1; com.lgene[0]=com.ls; com.posG[0]=0; com.posG[1]=com.ls;
   com.rgene[0]=1;
   for (h=0; h<com.npatt; h++) {
      for (j=0; j<com.ns; j++) {
         ch=fgetc(fpexp);
         while (!isalpha(ch)) ch=fgetc(fpexp);
         com.z[j][h]=(char)CodeChara((char)ch, 0);
      }
      if (fscanf(fpexp,"%lf",&com.fpatt[h]) != 1) error2 ("ReadPatternFreq");
   }
   for (j=0,sum=0; j<com.npatt; j++) sum+=com.fpatt[j];
   printf ("\nSUM freq = %.6f = 1?\n\n", sum);
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
