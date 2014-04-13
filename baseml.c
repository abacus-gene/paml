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

#include "tools.h"
#define NS            1000
#define NBRANCH       (NS*2-2)
#define NNODE         (NS*2-1)
#define NGENE         50
#define LSPNAME       30
#define NCODE         4
#define NCATG         40

#define NP            (NBRANCH+NGENE+11)
/*
#define NP            (NBRANCH+3*NNODE+NGENE+11+NCATG*NCATG-1)
*/
extern int noisy, NFunCall, *ancestor;
extern char GenetCode[][64], AAs[];

extern double *SeqDistance; 

int Forestry (FILE *fout, double space[]);
void DetailOutput (FILE *fout, double x[]);
int GetOptions (char *ctlf);
int GetInitials (double x[]);
int SetxInitials (double x[]);
int SetParameters (double x[]);
int SetPGene (int igene, int _pi, int _UVRoot, int _alpha, double xcom[]);
int SetPSiteClass(int iclass, double xcom[]);
int testx (double x[], int np);
int PartialLikelihood (int inode, int igene);

int TransformxBack(double x[]);
void ReadPatternFreq (char* pexpf);

int TestModel (FILE *fout, double x[], int nsep, double space[]);
int OldDistributions (int inode, double x[]);
int SubData(int argc, char *argv[]);

struct CommonInfo {
   char *z[NS], *spname[NS], seqf[96],outf[96],treef[96];
   int seqtype, ns,ls,ngene,posG[NGENE+1],lgene[NGENE], *pose,npatt;
   int cleandata;
   int np,ntime,nrgene,nrate,nalpha,npi,nhomo,ncatG,ncode,Mgene,sspace;
   int fix_kappa, fix_rgene, fix_alpha, fix_rho, nparK, fix_branch;
   int clock, model, getSE, runmode, print,verbose, ndata, icode,coding;
   int method;
   double *fpatt, kappa, alpha, rho, rgene[NGENE], pi[4],piG[NGENE][4];
   double freqK[NCATG], rK[NCATG], MK[NCATG*NCATG];
   double (*plfun)(double x[],int np), *lkl, *lkl0, *fhK;
}  com;
struct TREEB {
   int  nbranch, nnode, origin, branches[NBRANCH][2];
   double lnL;
}  tree;
struct TREEN {
   int father, nson, sons[NS], ibranch, label;
   double branch, divtime, kappa, pi[4], *lkl;
}  nodes[2*NS-1];

int nR=4, LASTROUND;
double PMat[16], Cijk[64], Root[4];


int lklSiteClass=0;  /* is lkl memory allocated for each class? */
char _oldlkl[NNODE];     /* update lkl for nodes? to save computation */
int _nnodeScale=0;
char _nodeScale[NNODE];   /* nScale[ns-1] for interior nodes */
double *_nodeScaleF=NULL;  /* nScaleF[npatt] for scale factors */

FILE *frub, *flnf, *frst;
char *ratef="rates";
char *models[]={"JC69","K80","F81","F84","HKY85","TN93","REV","UNREST"};
enum {JC69, K80, F81, F84, HKY85, TN93, REV, UNREST} MODELS;

double _rateSite=1;
int *rateBranch=NULL, N_rateBranch=-1;

#define BASEML 1
#include "treesub.c"
#include "treespace.c"

int main(int argc, char *argv[])
{
   FILE *fout, *fseq;
   int i, slkl0=0,slkl1=0,s2=0, READPATTERNFREQ=0, idata;
   char rstf[]="rst", ctlf[]="baseml.ctl";
   char *Mgenestr[]={"diff. rate", "separate data", "diff. rate & pi", 
                     "diff. rate & kappa", "diff. rate & pi & kappa"};
   double *space=NULL,t;
   clock_t clock_start, clock_end;

#ifdef __MWERKS__
/* Added by Andrew Rambaut to accommodate Macs -
   Brings up dialog box to allow command line parameters.
*/
	argc=ccommand(&argv);
#endif

   clock_start=clock();
   com.ndata=1;
   com.cleandata=0;  /* 0: delete; 1:treat missing data */
   noisy=0;
   com.runmode=2;    /* 0: user tree;  1: semi-automatic;  2: automatic  */

   com.clock=0;      /* 1: clock, rooted tree; 2:local clock; 0: no clock, unrooted tree */
   com.fix_rgene=0;  /* 0: estimate rate factors for genes */
   com.nhomo=0;
   com.getSE=0;      /* 1: want S.E.s of estimates;  0: don't want them */

   com.model=0;
   com.fix_kappa=0; com.kappa=5;
   com.fix_alpha=1; com.alpha=0.;  com.ncatG=4;    /* alpha=0 := inf */
   com.fix_rho=1;   com.rho=0.;    com.nparK=0;
   com.ncode=4;
   com.print=0;     com.verbose=1;  
   com.method=1;
   SetSeed(123);

   frub=fopen("rub","w");  frst=fopen(rstf,"w");

/*
printf ("\nrandom number seed? ");
scanf ("%d", &i);
SetSeed(i);
SubData(argc,argv);
exit(0);
*/

   if (argc>1)  strcpy(ctlf, argv[1]); 
   GetOptions (ctlf);

   if ((fout=fopen (com.outf, "w"))==NULL) error("outfile creation err.");

   if((fseq=fopen (com.seqf,"r"))==NULL)  {
      printf ("\n\nSequence file %s not found!\n", com.seqf);
      exit (-1);
   }
   for (idata=0; idata<com.ndata; idata++) {
      if (com.ndata>1) {
         printf("\nData set %d\n", idata+1);
         fprintf(fout, "\n\nData set %d\n", idata+1);
         fprintf(frst, "\t%d", idata+1);
      }
      if (idata)  GetOptions (ctlf);
      if (READPATTERNFREQ) {
         if (com.rho || com.Mgene) error("model");
         ReadPatternFreq ("pexp.in");
      }
      else
         ReadSeq ((com.verbose?fout:NULL), fseq, 0);

      /* BootstrapSeq("boot");  exit(0); */

      if(com.coding) {
         if(com.ls%3!=0 || (com.ngene!=1 && com.ngene!=3))
            error("this is not a coding sequence.  Remove icode?");
      }

      if (com.ngene>1 && com.Mgene==1)  printSeqsMgenes ();
      if (com.Mgene && com.ngene==1)  error ("option Mgene for 1 gene?");
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
      if(com.clock) fprintf(fout," %s clock  ",com.clock==1?"Global":"Local");
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
           error("oom");
      com.sspace=max2(300000,com.ls*(int)(com.ns*sizeof(char)+sizeof(int)));
      if((space=(double*)realloc(space,com.sspace))==NULL) error("oom space");

      if (!READPATTERNFREQ) {
         i=com.ns*(com.ns-1)/2;
         SeqDistance=(double*)realloc(SeqDistance,i*sizeof(double));
         ancestor=(int*)realloc(ancestor, i*sizeof(int));
         if(SeqDistance==NULL||ancestor==NULL) error("oom distance&ancestor");
      }

      if (!READPATTERNFREQ) Initialize (fout, space, 0);
      if (com.Mgene==3) FOR (i,com.ngene) xtoy(com.pi, com.piG[i],4);

      if (com.model==JC69 && !com.print && !READPATTERNFREQ) 
         PatternJC69like(fout);

      if(!com.cleandata) {
         slkl0=com.ns*com.ncode*com.npatt*sizeof(double);
         if((com.lkl0=(double*)malloc(slkl0))==NULL) error("oom lkl0");
      }
      slkl1=(com.ns-1)*com.ncode*com.npatt*sizeof(double);
      if((com.lkl=(double*)realloc(com.lkl,slkl1))==NULL) error("oom lkl1");

      if (com.alpha || com.nparK) {
         s2 = com.npatt*com.ncatG*sizeof(double);
         if((com.fhK=(double*)realloc(com.fhK,s2))==NULL) error("oom");
      }

      printf ("\n%9ld bytes for distance ",
         com.ns*(com.ns-1)/2*(sizeof(double)+sizeof(int)));
      printf("\n%9d bytes for lkl0\n%9d bytes for lkl1\n",slkl0,slkl1);
      printf("%9d bytes for fhK\n%9d bytes for space\n",s2,com.sspace);
      /* FOR(i,com.ns*2-1) xtoy(com.pi,nodes[i].pi, 4);  check this? */

      if(!READPATTERNFREQ) 
         DistanceMatNuc(fout,com.model,com.alpha);

      if (com.Mgene==1)        MultipleGenes (fout, space);
      else if (com.runmode==0) Forestry (fout, space);
      else if (com.runmode==3) StepwiseAddition (fout, space);
      else if (com.runmode>=4) Perturbation(fout, (com.runmode==4), space);
      else                  fprintf(frst,"%2d",StarDecomposition(fout,space));
      
      FPN(frst);  if((idata+1)%10==0) fflush(frst);
   }   /* for(idata) */
   fclose(fseq);  fclose(fout); fclose(frst); 
   if (noisy) putchar ('\a');

   clock_end=clock(); i=(int)(t=(clock_end*1.0-clock_start)/CLOCKS_PER_SEC);
   printf("\nTime used: %02d:%02d:%02.1f.\n", i/3600,(i%3600)/60, t-(i/60)*60);
   return (0);
}


/* x[]: t[ntime], rgene[ngene-1], kappa[nbranch], pi[nnode*3], 
        { alpha, rho || rK[], fK[] || rK[], MK[] }
*/
int ncorrect=0;

int Forestry (FILE *fout, double space[])
{
   static int times=0;
   int status=0, inbasemlg=0, i,j=0,k=0, itree,ntree, np,np0,iteration=1;
   int pauptree=0, btree=0, length_label;
   double x[NP], xcom[NP-NBRANCH], lnL,lnL0=0,lnLbest=0,e=1e-7, nchange=-1;
   double xb[NP][2], tl=0,tlls, *var=NULL;
   FILE *ftree,*finbasemlg=NULL,*fin=fopen("in.baseml","r"),*frate=NULL;

   if ((ftree=fopen(com.treef,"r"))==NULL) error("treefile not found");
   GetTreeFileType(ftree,&ntree,&pauptree,0);
   if(com.alpha) {
      if((frate=(FILE*)fopen(ratef,"w"))==NULL) error("rate file open err");
      fprintf(frate,"Rates for sites from BASEML, %d trees\n", ntree);
   }
   if (com.alpha && com.rho==0 && com.nhomo==0 && com.nparK==0) {
      inbasemlg=1;  if((finbasemlg=fopen("in.basemlg","w"))==NULL) error("file err"); 
   }
   if((flnf=fopen("lnf","w+"))==NULL) error("lnf file open error");

   if(!com.cleandata) InitPartialLikelihood ();
   for(itree=0; pauptree||itree<ntree; itree++,iteration=1) {
      if(pauptree && PaupTreeRubbish(ftree)) 
         { puts("err or end of tree file."); break; }
      if(ReadaTreeN(ftree,&length_label,1)) error("err or end of tree file.");
      if(noisy) printf ("\nTREE # %2d: ", itree+1);
      fprintf (fout,"\nTREE # %2d:  ", itree+1);
      fprintf (frub,"\n\nTREE # %2d\n", itree+1);
      fprintf (flnf,"\n\n%2d\n", itree+1);

      LASTROUND=0;
      if (times++==0 && (length_label==1 || length_label==3)) {
         if(com.clock) puts("\nBranch lengths in tree are ignored");
         else {
            puts("\ntree has branch lengths.\n0:ignore?\n1:initials?\n2:fixed?");
            scanf("%d",&com.fix_branch);
            /* com.fix_branch=2; */
            if(com.fix_branch==1) {
               FOR(i,tree.nnode) 
                  if((x[nodes[i].ibranch]=nodes[i].branch)<0)
                     error("blength<0 in tree");
            }
         }
      }

      if(com.cleandata) nchange=MPScore (space);
      if (noisy) { OutaTreeN(F0,0,0); printf("  MP score: %.2f\n",nchange);}
      OutaTreeN(fout,0,0);  fprintf(fout,"  MP score: %.2f",nchange);
      if(!com.clock && com.nhomo<=2 && nodes[tree.origin].nson<=2) {
         puts("\nThis is a rooted tree, without clock.  Check.\n");
         if(com.verbose) fputs("\nThis is a rooted tree.  Please check!",fout);
      }
      fflush(fout);   fflush(flnf);

      GetInitials (x);
      if((np=com.np)>NP) error("raise NP and recompile");

      if((i=spaceming2(np))>com.sspace)
         if((space=(double*)realloc(space,com.sspace=i))==NULL) 
            error("oom space");
      if(itree) { np0=np; }
      if(itree && (com.nhomo==0 || com.nhomo==2))
         for (i=0; i<np-com.ntime; i++) x[com.ntime+i]=max2(xcom[i],0.001);
      if(fin) {
         if(itree==0)
            puts("\nInitials from in.baseml.  Stop if not correct");
         fscanf(fin,"%lf",&x[i=0]);
         if (x[0]!=-1) i++;  else  iteration=0;
         for( ;i<np;i++) if(fscanf(fin,"%lf",&x[i])!=1) break;
         if (i<np)  {
            printf("err at #%d in in.baseml. Edit or remove it.\n",i+1);
            exit(-1);
         }
      }
      if (com.clock==2) {
         printf("\n%d substitution rates for branches assumed:\n",N_rateBranch);
         FOR(i,tree.nbranch) printf("%3d", rateBranch[i]); FPN(F0);
         fprintf(fout,"\n%d substitution rates for branches assumed:\n",
            N_rateBranch);
         FOR(i,tree.nbranch) fprintf(fout,"%3d",rateBranch[i]); FPN(fout);
      }

      PointLklnodes ();
      NFunCall=0;

      tlls = (com.clock?0:sum(x,com.ntime));
      if(noisy) {
         printf("\nntime & nrate & np:%4d%4d%4d\n",com.ntime,com.nrate,com.np);
         OutaTreeB(F0);  FPN(F0);
         if (noisy>2) matout(F0, x, 1, np);

         lnL = com.plfun (x,np);
         printf("\nlnL0 = %12.6f\n",-lnL);
/*
         printf ("\nGradient:\n");
         gradient(np,x,lnL,space,com.plfun,space+np,1);
         matout(F0,space,1,np);
         printf ("\nHessian:\n");
         var=space+2*np;
         Hessian (np, x, lnL, space+np, var, com.plfun, var+np*np);
         FOR(i,np) printf(" %11.6f",var[i*np+i]);
         exit(0);
*/
      }

      if (iteration && np) {
         SetxBound (np, xb);
         if(com.method==1 || com.ntime==0)
            j=ming2(noisy>2?frub:NULL,&lnL,com.plfun,NULL,x,xb, space,e,np);
         else
            minB(noisy>2?frub:NULL, &lnL,x,xb, space);
         if (j || lnL<=0 || lnL>1e7) status=-1;  else status=0;
      }

      if (itree==0) { lnL0=lnLbest=lnL; btree=0; }
      else if (lnL<lnLbest) { lnLbest=lnL;  btree=itree; }
      if (noisy) printf ("\nOut...\nlnL  = %12.6f\n", -lnL);
      fprintf(fout,"\nlnL(ntime:%3d  np:%3d):%14.6f%+12.6f\n",
          com.ntime,np,-lnL,-lnL+lnL0);
      if(com.fix_branch<2) { OutaTreeB (fout);  FPN (fout); }
      LASTROUND=1;
      if(com.npi || com.nparK>=2) TransformxBack(x);
/*
for(i=com.ntime;i<np;i++) fprintf(frst,"%9.5f",x[i]); 
fprintf(frst,"%12.4f %4d ",-lnL,com.npatt); 
*/

      FOR(i,np) fprintf(fout," %8.5f",x[i]);  FPN(fout); fflush(fout);
      if(inbasemlg) matout(finbasemlg,x,1,np);

      if (com.getSE) {
         var=space+np;
         Hessian (np,x,lnL,space,var,com.plfun,var+np*np);
         matinv(var,np,np,var+np*np);
         fprintf(fout,"SEs for parameters:\n");
         FOR (i, np)
            fprintf (fout," %8.6f", var[i*np+i]>0.?sqrt(var[i*np+i]):-0);
         FPN (fout);
         if (com.getSE==2) matout2(fout,var,np,np,12,7);
      }

      if(com.clock) SetBranch(x);
      for(i=0,tl=0;i<tree.nnode;i++) if(i!=tree.origin) tl+=nodes[i].branch;
      fprintf(fout,"\ntree length = %9.5f%s\n",tl,com.ngene>1?" (1st gene)":"");
      FPN(fout); OutaTreeN(fout,1,0); FPN(fout);
      FPN(fout); OutaTreeN(fout,1,1); FPN(fout);

      if(com.np-com.ntime||com.clock)  DetailOutput(fout,x);
      if (status) {
         printf ("convergence?\n");  fprintf (fout,"check convergence..\n");
      }
      if (itree==0) 
         for (i=0; i<np-com.ntime; i++) xcom[i]=x[com.ntime+i];
      else if (!j)  
         for (i=0; i<np-com.ntime; i++) xcom[i]=xcom[i]*.8+x[com.ntime+i]*0.2;
/*
      TestModel(fout,x,1,space);
      OldDistributions (tree.origin, x);
      fprintf(fout,"\n\n# of lfun calls:%10d\n", NFunCall);
*/
      if (com.model>=REV) 
         EigenREV(fout, x+com.ntime+com.nrgene, com.pi, &nR, Root, Cijk);
      
      com.print-=9;  com.plfun(x, np);  com.print+=9;
      if (com.print) {
         if (com.plfun!=lfun)  lfunRates (frate, x, np);
         if (com.nhomo==0 && com.nparK==0 && com.model<=REV && com.rho==0)
            AncestralSeqs (frst, x, space);
      }

      if(com.coding) { 
         fputs("\nExpected tree length for codon models:\n",fout);
         FOR(i,tree.nnode) nodes[i].branch*=(1+com.rgene[1]+com.rgene[2]);
         FPN(fout); OutaTreeN(fout,1,1); FPN(fout);
      }

   }         /* for (itree) */
   if(finbasemlg) fclose(finbasemlg);   if (frate) fclose(frate);
   if(fin) fclose(fin);  fclose(ftree);
   if(ntree>2)  { 
      fprintf(fout,"\nBest tree (#%d):%12.6f\n",btree+1,-lnLbest);
      if(btree==0) ncorrect++;
      fprintf(frst,"\t%d\t%d",btree+1,ncorrect);
   }

   if(noisy && ntree>1) { rewind(flnf);  rell(flnf,fout,ntree,space); }
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
   for (i=0; i<com.npi; i++)  xtoFreq(x+k+3*i, x+k+3*i, 4);
   k+=com.npi*3 + K-1;        /* K-1 for rK */
   if (com.nparK==2)          /* rK & fK */
      xtoFreq(x+k, x+k, K);
   else if (com.nparK==3)     /* rK & MK (double stochastic matrix) */
      for (i=0; i<K-1; k+=K-1,i++)  xtoFreq(x+k, x+k, K);
   else if (com.nparK==4)     /* rK & MK */
      for (i=0; i<K;   k+=K-1,i++)   xtoFreq(x+k, x+k, K);
   return(0);
}


void DetailOutput(FILE *fout, double x[])
{
   int i,j,k=com.ntime, nr[]={0, 1, 0, 1, 1, 2, 5, 11};
   double Qfactor,*p=com.pi, t=0, k1,k2, S,V, Y=p[0]+p[1],R=p[2]+p[3],tnode;

   /* date estimation under molecular clocks */
   if(com.clock>=1 && noisy>=9) { /* SetBranch() called before this. */
      fputs("\nNode & distance to present\n",fout);
      for(i=com.ns; i<tree.nnode;i++) 
         printf("Node %3d distance %9.5f\n",i+1,nodes[i].divtime);
      for (; ;) {
         printf("\nreference node & node time (-1 -1 to quit)? ");
         scanf("%d%lf", &j,&tnode);
         if(j<com.ns) break;
         FPN(F0);
         fprintf(fout,"\n\nNode %d Time %.3f\n\n",j,tnode);
         if(--j>=com.ns && tnode>0) {
            for(i=com.ns, tnode/=nodes[j].divtime; i<tree.nnode;i++) {
               printf("Node %3d Time %6.2f\n",i+1,nodes[i].divtime*tnode);
               fprintf(fout,"Node %3d Time %6.2f\n",i+1,nodes[i].divtime*tnode);
            }
         }
      }
   } 
   fprintf(fout,"\nDetailed output identifying parameters\n");
   if(com.clock>=2) {
      fprintf (fout,"rates for branches:    1");
      for(k=tree.nnode-com.ns; k<com.ntime; k++) fprintf(fout," %8.5f",x[k]);
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
      fprintf(fout,"kappa or rate parameter(s) in %s: ",models[com.model]);
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
      else  fprintf (fout, "gamma parameter alpha for %d genes:", com.ngene);
      for (i=0; i<com.nalpha; i++) fprintf (fout, " %8.5f", x[k++]);
      FPN (fout);
   }
   if (!com.fix_rho) {
      fprintf (fout, "rho for the auto-discrete-gamma model: %9.5f",x[k]);
      FPN (fout);
   }
   if ((com.alpha || com.nparK>=1) && com.nalpha<=1) {
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

int GetOptions (char *ctlf)
{
   int i, nopt=27, lline=255; 
   char line[255], *pline, opt[20], comment='*';
   char *optstr[]={"seqfile","outfile","treefile","noisy",
        "cleandata", "verbose","runmode", "method",
        "clock","fix_rgene","Mgene","nhomo","getSE","RateAncestor",
        "model","fix_kappa","kappa","fix_alpha","alpha","Malpha","ncatG", 
        "fix_rho","rho","nparK", "ndata", "Small_Diff","icode"};
   double t;
   FILE *fctl;

   if((fctl=fopen(ctlf,"r"))==NULL) error("\nctl file open error.\n");
   printf ("\n\nReading options from %s..\n", ctlf);
   for (;;) {
      if (fgets (line, lline, fctl) == NULL) break;
      for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
         if (isalnum(line[i]))  { t=1; break; }
         else if (line[i]==comment) break;
      if (t==0) continue;
      sscanf (line, "%s%*s%lf", opt, &t);
      if ((pline=strstr(line, "= "))==NULL) error("option file.");

      for (i=0; i<nopt; i++) {
         if (strncmp(opt, optstr[i], 8)==0)  {
            if (noisy>2)
               printf ("\n%3d %15s | %-20s %6.2f", i+1,optstr[i],opt,t);
            switch (i) {
               case ( 0): sscanf(pline+2, "%s", com.seqf);    break;
               case ( 1): sscanf(pline+2, "%s", com.outf);    break;
               case ( 2): sscanf(pline+2, "%s", com.treef);    break;
               case ( 3): noisy=(int)t;           break;
               case ( 4): com.cleandata=(int)t;  break;
               case ( 5): com.verbose=(int)t;     break;
               case ( 6): com.runmode=(int)t;     break;
               case ( 7): com.method=(int)t;      break;
               case ( 8): com.clock=(int)t;       break;
               case ( 9): com.fix_rgene=(int)t;   break;
               case (10): com.Mgene=(int)t;       break;
               case (11): com.nhomo=(int)t;       break;
               case (12): com.getSE=(int)t;       break;
               case (13): com.print=(int)t;       break;
               case (14): com.model=(int)t;       break;
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
         { printf ("\noption %s in %s\n", opt, ctlf);  exit (-1); }
   }
   fclose (fctl);

   if (com.fix_alpha==1 && com.alpha==0) {
      if (!com.fix_rho || com.rho) error("fix rho to 0 if alpha=0.");
   }
   if (com.nparK>=1) { 
      com.fix_alpha=com.fix_rho=1; 
      if(com.alpha==0) com.alpha=0.5; 
      if(com.nparK<=2) com.rho=0; else com.rho=0.4;
      if(com.nhomo>=1) error("nhomo & nparK");
   }
   if(com.model!=F84 && com.kappa<=0)  error("init kappa..");
   if(!com.fix_alpha && com.alpha<=0)  error("init alpha..");
   if(!com.fix_rho && com.rho==0) { com.rho=0.001;  puts("init rho reset"); }

   if (com.alpha)  if (com.ncatG<2 || com.ncatG>NCATG) error ("ncatG");
   if(com.method==0 &&(com.clock||com.rho)) 
      { com.method=1; puts("method reset"); }
   if (com.nhomo==2) {
      if (com.model!=K80 && com.model!=F84 && com.model!=HKY85) error("nhomo");
   }
   else if (com.nhomo>2 && com.model!=F84 && com.model!=HKY85) error("nhomo");
   else
      if (com.nhomo==1 && com.model!=F81 && com.model!=F84 && com.model!=HKY85)
         puts("nhomo: probably wrong.");

   if (com.nhomo>1 && com.runmode>0)  error("nhomo incompatible with runmode");
   if (com.clock && com.runmode>2)  error("runmode & clock?");
   if (com.runmode==3 && (com.npi || com.nparK))
      error("runmode etc.");

   if (com.model==JC69 || com.model==K80 || com.model==UNREST)
      if (com.nhomo!=2)  com.nhomo=0;
   if (com.model==JC69 || com.model==F81) { com.fix_kappa=1; com.kappa=1; }
   if (com.model==TN93 || com.model==REV)   com.fix_kappa=0;
   if (com.nparK==3) {
      puts("\n\nnparK==3, double stochastic, may not work.  Use nparK=4?\n");
      getchar();
   }

   return (0);
}

int GetInitials (double x[])
{
   static int times=0;
   int i,j,k, K=com.ncatG;
   int slkl1=(tree.nnode-com.ns)*com.ncode*com.npatt*sizeof(double);
   double t;

   com.plfun=lfunAdG;
   if (com.alpha==0 && com.nparK==0)  com.plfun=lfun;
   else if ((com.alpha && com.rho==0) || com.nparK==1 || com.nparK==2)
      com.plfun=lfundG;

   if(times++==0 && com.method==0 && com.fix_branch!=2 && com.plfun==lfundG) {
      lklSiteClass=1;
      printf("\n%9d bytes for lkl1, adjusted\n",slkl1*com.ncatG);
      if((com.lkl=(double*)realloc(com.lkl,slkl1*com.ncatG))==NULL)
         error("oom lkl1");
      sleep(1000);
   }
   InitializeNodeScale();

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
   case (2): com.npi=0;  com.nrate=tree.nbranch;  break;  /* b kappa's */
   case (3): com.npi=com.ns+(tree.origin>=com.ns)+(tree.nnode>com.ns+1);  
             com.nrate=(com.fix_kappa?1:tree.nbranch);  break;   /* ns+2 pi */
   case (4): com.npi=tree.nnode;  com.nrate=(com.fix_kappa?1:tree.nbranch);
                                  break;   /* nnode pi   */
   }

   if(com.fix_branch==2)   com.ntime=0;
   else if(com.clock==0)  com.ntime=tree.nbranch;
   else if(com.clock==1)  com.ntime=tree.nnode-com.ns+(tree.origin<com.ns);
   else {  /* if(com.clock==2) */
      for(i=0,k=0,N_rateBranch=0; i<tree.nbranch; i++) { 
         rateBranch[i]=j=nodes[tree.branches[i][1]].label;
         if(j+1>N_rateBranch)  N_rateBranch=j+1;
         if(j<0||j>tree.nbranch-1) error("branch label in the tree.");
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
      for(j=1,x[0]=1; j<tree.nnode-com.ns; j++) x[j]=0.7;
      for(; j<com.ntime; j++) x[j]=1;   /* rates for branches */
   }
   else if(com.fix_branch==0) {
      FOR (j,com.ntime) x[j]=.1;
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

   FOR (i,com.npi) fillxc(x+com.ntime+com.nrgene+com.nrate+3*i,0, 3);

   com.np = k = com.ntime+com.nrgene+com.nrate+com.npi*3;
   if (com.alpha==0 && com.nparK==0) return(0);

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
   return (0);
}

static complex cU[16], cV[16], cRoot[16];
int cPMat(double P[],double t,int n,complex cU[],complex cV[],complex cRoot[]);



void xtoFreq(double x[], double freq[], int n)
{
/* This transforms the iterative variable x[] to freq[0,1,n-2]:
      freq[k] = exp(x[k])/(1+SUM(exp(x[k]))), k=0,1,...,n-2
   x[] and freq[] may be the same vector, and the last element 
   freq[n-1] and x[n-1] is not used.  
*/
   int k;
   double t;

   for(k=0,t=1;k<n-1;k++) t+=exp(x[k]);
   for(k=0;k<n-1;k++) freq[k]=exp(x[k])/t;
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
   FOR (i, com.nrgene) com.rgene[i+1]=x[com.ntime+i];
   if (!com.fix_kappa && com.model<=TN93) {
       com.kappa=k1=k2=x[com.ntime+com.nrgene];
       if (com.model==TN93) k2=x[com.ntime+com.nrgene+1];
   }
   if (com.nhomo==1) {
      k=com.ntime+com.nrgene+com.nrate;
      if (!LASTROUND) xtoFreq(x+k, com.pi, 4);
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
            if (!LASTROUND) xtoFreq(x+k+i*3, nodes[i].pi, 4);
            else            xtoy   (x+k+i*3, nodes[i].pi, 3);
         for (i=com.ns; i<tree.nnode; i++)
            if (!LASTROUND) xtoFreq(x+k+com.ns*3, nodes[i].pi, 4);
            else            xtoy   (x+k+com.ns*3, nodes[i].pi, 3);

         if (tree.origin>=com.ns && tree.nnode>com.ns+1)
            if (!LASTROUND) xtoFreq(x+k+(com.ns+1)*3, nodes[tree.origin].pi,4);
            else            xtoy   (x+k+(com.ns+1)*3, nodes[tree.origin].pi,3);
      }
      else { /*  if(com.nhomo==4) */
         FOR (i,tree.nnode)
            if (!LASTROUND) xtoFreq(x+k+i*3, nodes[i].pi, 4);
            else            xtoy   (x+k+i*3, nodes[i].pi, 3);
      }

      FOR(i,tree.nnode) nodes[i].pi[3]=1-sum(nodes[i].pi,3);
      xtoy(nodes[tree.origin].pi, com.pi, 4);
   }
   else if (com.model==REV && com.Mgene<=1)
      EigenREV (NULL, x+com.ntime+com.nrgene, com.pi, &nR, Root, Cijk);
   else if (com.model==UNREST && com.Mgene<=1)
      EigenUNREST (NULL, x+com.ntime+com.nrgene,com.pi,&nR,cRoot,cU,cV);

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
      if (!LASTROUND) xtoFreq(x+k+K-1, com.freqK, K);
      else            xtoy   (x+k+K-1, com.freqK, K-1);
      com.freqK[K-1]=1-sum(com.freqK, K-1);
   }
   else if (com.nparK==3) {   /* rK & MK (double stochastic matrix) */
      for (i=0,k+=K-1; i<K-1; k+=K-1,i++) {
         if (!LASTROUND) xtoFreq(x+k, com.MK+i*K, K);
         else            xtoy   (x+k, com.MK+i*K, K-1);
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
         if (!LASTROUND) xtoFreq(x+k, com.MK+i*K, K);
         else            xtoy   (x+k, com.MK+i*K, K-1);
         com.MK[i*K+K-1]=1-sum(com.MK+i*K,K-1);
      }
      PtoPi(com.MK, com.freqK, K, space);
   }
   com.rK[K-1]=(1-innerp(com.freqK, com.rK, K-1))/com.freqK[K-1];

   return (status);
}

int SetPGene(int igene, int _pi, int _UVRoot, int _alpha, double xcom[])
{
/* xc[] does not contain time parameters
*/
   int nr[]={0, 1, 0, 1, 1, 2, 5, 11};
   int k=com.nrgene+(com.Mgene>=3)*igene*nr[com.model];
   double ka1=xcom[k], ka2=(com.model==TN93?xcom[k+1]:-1);

   if (_pi) xtoy(com.piG[igene], com.pi, 4);
   if (_UVRoot) {
      if (com.model==K80) com.kappa=ka1;
      else if (com.model<=TN93) 
         EigenTN93(com.model,ka1,ka2,com.pi,&nR,Root,Cijk);
      else if (com.model==REV)
         EigenREV(NULL,xcom+k,com.pi,&nR,Root,Cijk);
   }
   if (_alpha) {
      com.alpha=xcom[com.nrgene+com.nrate+com.npi+igene]; /* check?? */
      DiscreteGamma(com.freqK,com.rK,com.alpha,com.alpha,com.ncatG,0);
   }
   return(0);
}


int SetxBound (int np, double xb[][2])
{
/* sets lower and upper bounds for variables during iteration
*/
   int i,j, k=com.ntime+com.nrgene+com.nrate, nf=0;
   double tb[]={.4e-5, 99}, rgeneb[]={1e-4,99}, rateb[]={1e-5,999};
   double alphab[]={.04, 99}, rhob[]={-0.2, 0.99}, pb[]={.00001,.99999};
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
      xb[k][0]=fb[0]; xb[k++][1]=fb[1];
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
   double tb[]={.4e-5,20};

   if (com.clock) error ("testx: clock not supported"); 
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
      if (com.nhomo==2)
         EigenTN93(com.model,nodes[ison].kappa,1,com.pi,&nR,Root,Cijk);
      else if (com.nhomo>2)
         EigenTN93(com.model,nodes[ison].kappa,1,nodes[ison].pi,&nR,Root,Cijk);
      if (com.model<=K80)
         PMatK80(PMat, t, (com.nhomo==2?nodes[ison].kappa:com.kappa));
      else if (com.model<=REV)  PMatCijk(PMat, t);
      else                      cPMat (PMat, t, n, cU, cV, cRoot);

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

   if(_nnodeScale && _nodeScale[inode]) {  /* scaling to avoid underflow */
      for(i=0,k=0; i<tree.nnode; i++)   /* inode is k-th node for scaling */
         if(i==inode) break;
         else if(_nodeScale[i]) k++;

      for(h=pos0; h<pos1; h++) {
         for(j=0,t=0;j<n;j++)
            if(nodes[inode].lkl[h*n+j]>t) t=nodes[inode].lkl[h*n+j];
         
         if(t<1e-300)  _nodeScaleF[k*com.npatt+h]=-500;  /* w=0 for sites */
         else {  
            for(j=0;j<n;j++)  nodes[inode].lkl[h*n+j]/=t;
            _nodeScaleF[k*com.npatt+h]=log(t);
         }
      }
   }
   return (0);
}


int PMatCijk (double P[], double t)
{
/* P(t)ij = SUM Cijk * exp{Root*t}
*/
   int i,j,k, n=4, nr=nR;
   double expt[4];

   if (t<-1e-3) 
      printf ("\nt = %.5f in PMatCijk", t);
   if (t<1e-7) { identity (P, n); return(0); }

   for (k=1; k<nr; k++) expt[k]=exp(t*Root[k]);
   FOR (i,n) FOR (j,n) {
      for (k=1,t=Cijk[i*n*nr+j*nr+0]; k<nr; k++)
         t+=Cijk[i*n*nr+j*nr+k]*expt[k];
      P[i*n+j]=t;
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

   SiteCat=(int*)malloc((1<<2*com.ns)* sizeof(int));
   if (SiteCat == NULL)  error ("oom");
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

   fprintf (flnf, "\nnode %4d", inode+1);
   FOR (i, 4) fprintf (flnf, " %8.5f", freq[inode][i]);
   FOR (i, nodes[inode].nson)  
      OldDistributions (nodes[inode].sons[i], x);
   return (0);
}

#endif







void ReadPatternFreq (char* pexpf)
{
   FILE *fpexp;
   int h, j, ch;
   double sum;

   printf ("\n\nRead site-pattern frequencies from %s\n", pexpf);
   if ((fpexp=fopen(pexpf,"r")) == NULL) error ("pexp file open err");
   fscanf (fpexp, "%d%d", &com.ns, &com.npatt);
   printf ("%3d species, %d site patterns\n", com.ns, com.npatt);

   if (com.model>1) {
      error("Should calculate base frequencies\a\n");
   }
   else 
      fillxc(com.pi, .25, 4);  
   xtoy (com.pi, com.piG[0], 4);

   if(com.npatt>(1<<(2*com.ns))) error("npatt");
   com.ls=com.npatt;
   com.fpatt=(double*) malloc (com.npatt*sizeof(double));
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
      if (fscanf(fpexp,"%lf",&com.fpatt[h]) != 1) error ("ReadPatternFreq");
   }
   for (j=0,sum=0; j<com.npatt; j++) sum+=com.fpatt[j];
   printf ("\nSUM freq = %.6f = 1?\n\n", sum);
   fclose (fpexp);
}


double TreeScore(double x[], double space[])
{
   double xb[NP][2], e=1e-6, lnL;

   PointLklnodes ();
   if(com.clock==2) error("local clock in TreeScore");
   com.ntime = com.clock ? tree.nnode-com.ns : tree.nbranch;
   GetInitials (x);
   SetxBound (com.np, xb);
   if(!com.cleandata) InitPartialLikelihood ();

   ming2(NULL,&lnL,com.plfun,NULL,x,xb, space,e,com.np);
   return (lnL);
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


parameters under local clock (com.clock=2)
   com.ntime = (#ancestral nodes) - 1 + (#rates) - 1
Think of node times t[] (ns-1 of them) and rates for branches.
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
   int ib,j,k,father,inode,sib,nson0=nodes[tree.origin].nson,newn,nmark[NNODE];

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
         error("This should not be possible.");
      }
      else if(nodes[father].nson==1) {
         sib=nodes[father].sons[0];
         if(sib==inode) error("This should not be possible 1");

         if(father==tree.origin) {  /* remove old root */
            tree.origin=sib; bkeep[nodes[sib].ibranch]=0;
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
   if(nodes[tree.origin].nson<2) error("strange");
   else if(nodes[tree.origin].nson==2 && !com.clock && 
      nodes[tree.origin].nson<nson0) {
      /* deroot, k is new root */
      j=nodes[tree.origin].sons[0]; k=nodes[tree.origin].sons[1];
      if(j>k) {k=j; j=nodes[tree.origin].sons[1]; }
      if(k<com.ns) if(noisy) puts("only two seqs left.");

      tree.origin=k; 
      bkeep[nodes[k].ibranch]=0; 
      tree.branches[nodes[j].ibranch][0]=k;
      tree.branches[nodes[j].ibranch][1]=j;
      nodes[j].branch+=nodes[k].branch;
   }

   /* error checking, k=err */
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
      tree.nbranch=newn;  tree.origin=nmark[tree.origin];
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

   if ((fin=fopen(com.seqf,"r"))==NULL) error("infile error");
   if ((fout=fopen(outfile,"w"))==NULL) error("outfile creation error");
   ReadSeq (NULL, fin, 0);
   if(nkept>com.ns || nkept<2) error("err");
   GetTreeFileType(fin, &ntree, &pauptree, 1);
   if (!pauptree) error("Is this a paup data & tree file?");
   PaupTreeRubbish(fin);
   if (ReadaTreeN (fin, &i, 1)) error ("err reading tree.");

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
