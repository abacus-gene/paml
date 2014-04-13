/* mcmctree.c 

   Markov chain Monte Carlo calcualtioin of posterior probabilities of
   trees, with a prior (coalescent, Yule, and birth-death processes).

                         Ziheng YANG, September 1995 onwards

                cc -c -fast tools.c 
                cc -o mcmctree -fast mcmctree.c tools.o -lm
                     mcmctree <ControlFileName>

     x[birth, death, mut, kappa etc.]  
*/

#include "paml.h"

#ifdef __MWERKS__
/* Added by Andrew Rambaut to accommodate Macs -
   Brings up dialog box to allow command line parameters.
*/
#include <console.h>
#endif

#define NS            9
#define NBRANCH      (NS*2-2)
#define NNODE        (NS*2-1)
#define NGENE         1
#define LSPNAME       40
#define NCODE         4
#define NCATG         8

extern int noisy, NFunCall;
extern char BASEs[];

int GetOptions (char *ctlf);
int GetLHistoryI (int iLH);
int GetIofLHistory (void);
double OneTree(double delta);
int CountLHistory(char LHistories[], double space[]);
int ReorderNodes (char LHistory[]);
void MCMCtrees (FILE* fout, double space[]);
void MCMCunrootedtree (double space[]);
void EvaluateLHs(FILE*fout, FILE*fLH, int IofLHs[], double lnPLHs[],
    int ntreekept, double delta);
int InitPartialLikelihood (void);
int PartialLikelihood (int inode, int igene);
int PMatCijk (double P[], double t);
double lfun (void);
double lfundG (void);

struct CommonInfo {
   char *z[NS], *spname[NS], seqf[96],outf[96],treef[96],LHf[96];
   int seqtype, ns, ls, ngene, posG[NGENE+1],lgene[NGENE],*pose,npatt;
   int np, ncode, ntime, nrate, nrgene, nalpha, npi, ncatG, print, MCMC, seed;
   int cleandata;
   int priTt,hier, model, clock, fix_alpha, fix_kappa, fix_rgene, Mgene;
   double *fpatt;
   double pi[6],lmax, birth,death,sample, mut, kappa,alpha, beta,delta0,delta1;
   double rgene[NGENE],piG[NGENE][6], freqK[NCATG], rK[NCATG], *lkl, *fhK;
}  com;
struct TREEB {
   int  nbranch, nnode, root, branches[NBRANCH][2];
   double lnL;
}  tree;                
struct TREEN {
   int father, nson, sons[NS], ibranch, label;
   double branch, divtime, *lkl;
}  nodes[2*NS-1];

FILE *frub, *frst;
char *models[]={"JC69", "K80", "F81", "F84", "HKY85", "TN93", "REV"};
enum {JC69, K80, F81, F84, HKY85, TN93, REV} MODELS;

int nR=4, nreplicateMC;
double lnpBestTree, PMat[16], Cijk[64], Root[4];
double _rateSite=1;
int LASTROUND=0; /* no use for this */

#define REALSEQUENCE
#define PARSIMONY
#define NODESTRUCTURE
#define BIRTHDEATH
#include "treesub.c"
#include "treespace.c"

int main(int argc, char *argv[])
{
   char ctlf[]="mcmctree.ctl", rstf[96]="rst";
   char *priorTt[]={"Coalescent", "Yule", "Birth-death"};
   int idata, ndata=1;
   double *space;
   FILE  *fout, *fseq;

#ifdef __MWERKS__
/* Added by Andrew Rambaut to accommodate Macs -
   Brings up dialog box to allow command line parameters.
*/
	argc=ccommand(&argv);
#endif

   noisy=3;          
   com.seqtype=0;    com.model=F84;    com.kappa=1.63;
   com.alpha=0.;     com.ncatG=1;
   com.priTt=2;      /* 0: coalesent; 1: Yule; 2: birth-death */
   com.hier=0;   /* hierarchical, with priors for birth & death rates */ 
/*   com.birth=8.2;    com.death=4.1;   com.sample=9./185;   com.mut=0.24; */
   com.birth=6.7;    com.death=2.5;   com.sample=9./150;   com.mut=0.24; 
   com.ncode=4;      com.clock=1;
   com.print=1;

   frub=fopen ("rub", "w");   frst=fopen (rstf, "w");
/*
   printf ("\nNumber of data sets? ");
   scanf ("%d", &ndata);
*/
   if (argc>1) 
      { strcpy(ctlf, argv[1]); printf ("\nctlfile reset to %s.\n", ctlf); }
   GetOptions (ctlf);
   if ((fout=fopen(com.outf,"w"))==NULL) error2("outfile creation err.");
   if((fseq=fopen (com.seqf,"r"))==NULL)  {
      printf ("\n\nSequence file %s not found!\n", com.seqf);
      exit (-1);
   }

   if ((space=(double*)malloc(50000*sizeof(double)))==NULL) error2("oom");

   for (idata=0; idata<ndata; idata++) {
      if (ndata>1) {
         printf ("\nData set %d\n", idata+1);
         fprintf (fout, "\n\nData set %d\n", idata+1);
      }
      ReadSeq(NULL,fseq);
      
      fprintf(fout,"MCMCTREE %15s %8s (%s prior) ", 
                     com.seqf, models[com.model], priorTt[2]);
      if (com.clock) fprintf (fout, " Clock  ");
      if (com.alpha) fprintf (fout,"dGamma (ncatG=%d)", com.ncatG);
      if (com.nalpha>1) fprintf (fout,"(%d gamma)", com.nalpha);
      if (com.ngene>1) error2 ("multiple genes not supported");
      
      fprintf(fout,"\n%s Bayesian analysis:", 
          (com.hier?"Hierarchical":"Empirical"));
      fprintf (fout, "\n birth & death & sample & mut: %9.4f%9.4f%9.4f%9.4f",
          com.birth, com.death, com.sample, com.mut);
      printf("\n%s Bayesian analysis:", (com.hier?"Hierarchical":"Empirical"));
      printf ("\n birth & death & sample & mut: %9.4f%9.4f%9.4f%9.4f",
          com.birth, com.death, com.sample, com.mut);
      
      Initialize (fout, space);
      if (com.model==JC69) PatternJC69like (fout);
      if (com.lkl) free(com.lkl);
      com.lkl=(double*)malloc((com.ns-1)*com.ncode*com.npatt*sizeof(double));
      if(com.alpha){
         if (com.fhK) free(com.fhK);
         com.fhK=(double*)malloc(com.npatt*com.ncatG*sizeof(double));
         DiscreteGamma (com.freqK, com.rK, com.alpha, com.alpha, com.ncatG, 0);
      }
      if (com.lkl==NULL || (com.alpha && com.fhK==NULL)) error2 ("oom");
      
      if (com.model>1)
         EigenTN93(com.model,com.kappa,com.kappa,com.pi,&nR,Root,Cijk);
      
      printf ("\n%smolecular clock assumed.\n", (com.clock?"":"no "));
      fflush (fout);
      SetSeed (com.seed);
      if(!com.cleandata) InitPartialLikelihood ();
      MCMCtrees(fout, space);
   }
   fclose (fseq);
   putchar ('\a');
   return (0);
}


int PMatCijk (double P[], double t)
{
   int i,j,k, n=4, nr=nR;
   double expt[4];

   if (t<1e-7) { identity(P,n); return(0); }
   for (k=1; k<nr; k++) expt[k]=exp(t*Root[k]);
   FOR (i,n) FOR (j,n) {
      for (k=1,t=Cijk[i*n*nr+j*nr+0]; k<nr; k++)
         t+=Cijk[i*n*nr+j*nr+k]*expt[k];
      P[i*n+j]=t;
   }
   return (0);
}

int InitPartialLikelihood (void)
{
/* set partial likelihood at tips of the tree, considering missing data.
   This need to be modified when more proper algorithm for dealing with
   missing data is worked out.
   Need testing if sequences in the data are ancestors.
*/
   int n=com.ncode, is,j,k,h, b;
   char *pch=BASEs;

   FOR(is,com.ns) nodes[is].lkl=com.lkl+n*com.npatt*is;
   for (is=0;is<com.ns;is++) {
      zero(nodes[is].lkl, com.npatt*n);
      for (h=0;h<com.npatt;h++) {
            k=strchr(pch,com.z[is][h])-pch; 
            if(k<0) { printf("Character %d\n",com.z[is][h]); exit(-1); }
            if(k<n)
               nodes[is].lkl[h*n+k]=1;
            else
               FOR(j,nBASEs[k]) {
                  b=strchr(BASEs,EquateNUC[k][j])-BASEs;
                  nodes[is].lkl[h*n+b]=1.;  /* Joe's idea */
               }
      }  /* for (h) */
   }     /* for (is) */
   return(0);
}

int PartialLikelihood (int inode, int igene)
{
   int n=com.ncode, i,j,k,h, ison, pos0=com.posG[igene],pos1=com.posG[igene+1];
   double t;


   FOR (i,nodes[inode].nson)
      if (nodes[nodes[inode].sons[i]].nson>0)
         PartialLikelihood (nodes[inode].sons[i], igene);
   fillxc (nodes[inode].lkl+pos0*n, (double)(inode>=com.ns), (pos1-pos0)*n);
   if (com.cleandata && inode<com.ns) 
      for(h=pos0;h<pos1;h++) nodes[inode].lkl[h*n+com.z[inode][h]]=1;

   FOR (i, nodes[inode].nson) {
      ison=nodes[inode].sons[i];
      t=nodes[ison].branch*com.rgene[igene]*_rateSite;
      if (com.model<=K80)       PMatK80 (PMat, t, com.kappa);
      else if (com.model<=REV)  PMatCijk (PMat, t);

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
    return (0);
}


double lfundG (void)
{
/* discrete gamma rates for sites
*/
   int  h, ir, i, ig;
   double lnL, fh=0;

   zero (com.fhK, com.npatt);
   for (ig=0; ig<com.ngene; ig++) {
      for (ir=0; ir<com.ncatG; ir++) {
         _rateSite=com.rK[ir];
         PartialLikelihood (tree.root, ig);
         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            for (i=0,fh=0; i<com.ncode; i++)
               fh += com.pi[i]*nodes[tree.root].lkl[h*com.ncode+i];
            com.fhK[h] += com.freqK[ir]*fh;
         }
      }
   }
   for (h=0,lnL=0; h<com.npatt; h++) {
      if (com.fhK[h]<=0)
         printf ("\nlfundG: h=%4d  fhK=%9.4f", h, com.fhK[h]);
      lnL-=log(com.fhK[h])*com.fpatt[h];
   }
   return (lnL);
}

double lfun (void)
{
   int  h, i, ig;
   double lnL=0, fh;

   for (ig=0; ig<com.ngene; ig++) {
      PartialLikelihood (tree.root, ig);
      for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
         for (i=0,fh=0; i<com.ncode; i++)
            fh += com.pi[i]*nodes[tree.root].lkl[h*com.ncode+i];
         if (fh<=0) {
            printf ("  lfun: h=%4d fh=%9.4f: %8.4f%8.4f%8.4f\n",
                    h, fh, com.birth, com.death, com.mut);
            OutaTreeN(F0,0,1); FPN(F0);
	 }
         lnL-=log(fh)*com.fpatt[h];
      }
   }
   return (lnL);
}

int GetOptions (char *ctlf)
{
   int i, nopt=19, lline=255;
   char line[255],*pline, opt[20], comment='*';
   char *optstr[] = {"seqfile","outfile","treefile","LHfile", "MCMC", "beta",
        "seed", "delta0", "delta1", "model", "kappa", "alpha", "ncatG", 
        "hierarch", "birth", "death", "sample", "mutate", "cleandata"};
   double t=1;
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
         if ((pline=strstr(line, "="))==NULL) error2 ("option file.");

         for (i=0; i<nopt; i++) {
            if (strncmp(opt, optstr[i], 8)==0)  {
               if (noisy>2)
                  printf ("\n%3d %15s | %-20s %6.2f", i+1,optstr[i],opt,t);
               switch (i) {
                  case ( 0): sscanf(pline+2, "%s", com.seqf);    break;
                  case ( 1): sscanf(pline+2, "%s", com.outf);    break;
                  case ( 2): sscanf(pline+2, "%s", com.treef);   break;
                  case ( 3): sscanf(pline+2, "%s", com.LHf);     break;
                  case ( 4): com.MCMC=(int)t;       break;
                  case ( 5): com.beta=t;            break;
                  case ( 6): com.seed=(int)t;       break;
                  case ( 7): com.delta0=t;          break;
                  case ( 8): com.delta1=t;          break;
                  case ( 9): com.model=(int)t;      break;
                  case (10): com.kappa=t;           break;
                  case (11): com.alpha=t;           break;
                  case (12): com.ncatG=(int)t;      break;
                  case (13): com.hier=(int)t;       break;
                  case (14): com.birth=t;           break;
                  case (15): com.death=t;           break;
                  case (16): com.sample=t;          break;
                  case (17): com.mut=t;             break;
                  case (18): com.cleandata=(int)t;  break;
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

   if (com.model>HKY85)  error2("model.");
   if (com.model!=F84 && com.kappa<=0)  error2("init kappa.");
   if (com.alpha==0)  com.nalpha=0;
   else              com.nalpha=(com.nalpha?com.ngene:!com.fix_alpha);
   if (com.Mgene==1) error2 ("separate analyses not implemented.  DIY.");
   if (com.alpha)  if (com.ncatG<2 || com.ncatG>NCATG) error2 ("ncatG");
   if (com.model==JC69 || com.model==F81) { com.fix_kappa=1; com.kappa=1; }
   if (com.delta1>10 || com.delta1<0.1) puts ("\adelta1 not good..");

   return (0);
}

double OneTree (double delta)
{
/* returns log(P{data|tree}) for one tree (labeled history), calculated 
   by Monte Carlo integration over the node times or branch lengths.
   Stopping rule: slnp=s.e.(log(P))<delta
*/
   int i, ir, maxnr=1000000, minnr=1000, nrgap=100;
   double lnp, ScaleF,scalegap=20, t, p,mp,sp, mlnp,slnp, d=1;
   double birth=com.birth, death=com.death;

   PointLklnodes ();
   for (ir=0,ScaleF=mp=mlnp=sp=slnp=0; ir<maxnr; ir++) {
      if (com.hier)  { birth=2*com.birth*rndu();  death=2*com.death*rndu(); }
      if (com.clock) {
         /* check whether rooted or unrooted tree should be used here */
         BranchLengthBD (1, birth, death, com.sample, com.mut);
         FOR(i,com.ns) nodes[i].branch=max2(1e-7,nodes[i].branch);
      }
      else
         FOR(i,tree.nnode) if(i!=tree.root) nodes[i].branch=rndu()*com.mut;
      lnp=(com.alpha ? -lfundG() : -lfun());
      if (ir==0) ScaleF=lnp+scalegap;
      else if (lnp-ScaleF>scalegap) {
         mp*=(t=exp(ScaleF-(lnp+scalegap)));
         sp*=t*t;
         ScaleF=lnp+scalegap; 
      }
      p=exp(lnp-ScaleF);
      sp+=square(p-mp)*ir/(ir+1.);  mp=(mp*ir+p)/(ir+1.);
      slnp+=square(lnp-mlnp)*ir/(ir+1.);  mlnp=(mlnp*ir+lnp)/(ir+1.);

      if (ir+1>=minnr && (ir+1)%nrgap==0) {

         d=lnpBestTree-ScaleF-log(mp); 
         if (d>100) break;
         d=(d<10?1:1+d/50.);  d=min2(d,20);
         /* if (sqrt(slnp/ir/(ir+1.))<delta*d)  break; */
         if (sqrt(sp/ir/(ir+1.))/mp<delta*d)  break;
      }
   }
   nreplicateMC=min2(ir+1,maxnr);
   if ((mp=ScaleF+log(mp))>lnpBestTree) lnpBestTree=mp;

   /* printf ("mlnp: %14.4f%14.4f\n", mp, mlnp); */

   return (mp);
}


#define MAXTREEKEPT 500


void MCMCtrees (FILE* fout, double space[])
{
/* MCMC calculation in the tree space.  The chain moves to another labeled 
   history of the same rooted tree topology or to a random labeled history 
   of another direct neighbor by NNI.
*/
   struct TREE {struct TREEB tree; struct TREEN nodes[2*NS-1]; } tree0;
   int irun,burnin=200,nrun=100000, i,j;
   int ntreeNNI=2*(com.ns-2-!com.clock), iLH=-1,jLH=0, norder=500;
   int *IofLHs, *counts, ntreekept=0, nLHi=-1,nLHj;
   double *lnPLHs,lnPcur=0,lnP=0,alpha, nchange;
   char LH_NNI[(NS-1)*500], TL='T', wantL;
   FILE *ftree, *fLH=NULL;

   if(MAXTREEKEPT<burnin) puts("\nMAXTREEKEPT is too small??");
   lnPLHs=(double*)malloc(MAXTREEKEPT*3*sizeof(double));
   if(lnPLHs==NULL) error2("oom MCMCtrees");
   IofLHs=(int*)(lnPLHs+MAXTREEKEPT);
   counts=IofLHs+MAXTREEKEPT;

   for (i=com.ns,j=1; i>2; i--)  j*=i*(i-1)/2;
   printf ("\n%d labeled histories.\n", j);
   lnpBestTree=-1e40;
   if (com.MCMC)
      printf ("\nRun MCMC to generate candidate labeled histories.\n");
   else {
      printf("\nLabeled histories read from the file %s.\n", com.LHf);
      if ((fLH=fopen(com.LHf,"r"))==NULL) error2("LHs file open error2");
      EvaluateLHs(fout, fLH, IofLHs, lnPLHs, ntreekept, com.delta1);
      free(lnPLHs);
      return ;
   }
   if ((ftree=fopen(com.treef,"r"))==NULL) error2 ("no treefile");
   fscanf(ftree, "%d%d", &i, &irun);  /* irun (ntree) ignored */
   if (i!=com.ns) error2 ("ns in the tree file");
   if (ReadaTreeN (ftree, &i, 1)) error2 ("err tree..");  fclose(ftree);
   OutaTreeN (F0, 0, 0);  FPN(F0);   OutaTreeB (F0);  FPN(F0);
   if (com.ns*2-1!=tree.nnode) error2("Use a rooted tree to start.");
   nLHj=CountLHistory(LH_NNI, space);
   FOR (i,nLHj) {
      printf ("\nLH#%3d/%d:", i+1, nLHj);
      FOR (j,com.ns-1) printf("%3d", LH_NNI[i*(com.ns-1)+j]+1);
   }
   FPN (F0);
   if (nLHj>1) {
      printf ("\nchoose an ordering (1 -- %d)?\n", nLHj);
      scanf("%d", &jLH);  jLH--;
   }
   printf ("\nStart the chain with ordering #%d...\n", jLH+1);
   if (jLH<0 || jLH>nLHj-1) error2 ("not in range."); 

   FOR (i,MAXTREEKEPT) counts[i]=0;
   for (irun=-burnin,ntreekept=0,wantL=0; irun<nrun; irun++) {
      /* find next potential LH (jLH), given at the start of the chain */
      if (irun>-burnin) {
         tree=tree0.tree;  memcpy(nodes,tree0.nodes,sizeof(nodes));
         if ((nLHi==1 || rndu()>com.beta) && !wantL) { /* change topology*/
            TL='T';
            NeighborNNI((int)(ntreeNNI*rndu()));
         }
         else  TL='L';                                 /* change LH */
         nLHj=CountLHistory(LH_NNI, space);  jLH=(int)(nLHj*rndu());
         if (nLHj>norder) { printf("%9d>%9d",nLHj,norder); error2("norder"); }
      }
      ReorderNodes (LH_NNI+jLH*(com.ns-1));
      jLH=GetIofLHistory ();
      if (jLH==iLH) { ++wantL; irun--; continue; } 
      else            wantL=0;
      nchange=MPScore (space);

      printf ("%7d: %9d >%9d ", irun+1,iLH,jLH);
      OutaTreeN (F0, 0, 0);
      printf("%2c %3d %5.0f %9.3f > ", TL,nLHj,nchange,lnPcur);  fflush(F0);

      /* calculate lnP for jLH. */
      for (i=0,nreplicateMC=0; i<ntreekept; i++) 
         if (IofLHs[i]==jLH) { lnP=lnPLHs[i]; break; }
      if (i==ntreekept) {
         lnP=OneTree(com.delta0);
         IofLHs[ntreekept]=jLH; lnPLHs[ntreekept++]=lnP;
         if(ntreekept>=MAXTREEKEPT)  
           { puts("\n\nReached max2 number of trees kept.");  break; }
      }
      alpha=(irun==-burnin ? 1 : exp(lnP-lnPcur)*nLHj/(double)nLHi);
      if (alpha>1) alpha=1;

      printf("%9.3f %8.4f", lnP,alpha);

      /* accept or reject jLH? */
      if (rndu()<alpha) {  /* go to tree j */
         tree0.tree=tree;  memcpy (tree0.nodes,nodes,sizeof(nodes));
         iLH=jLH;  lnPcur=lnP;  nLHi=nLHj;
         printf (" Y");
      } 
      else   printf (" N");

      if (irun>=0) {
         for(i=0;i<ntreekept;i++) if (iLH==IofLHs[i]) break;
         counts[i]++;
      }
      printf("%8d%5d\n", nreplicateMC, ntreekept);
   }  /* for (irun) */

   printf ("\n\n%d LHs collected into the file %s\n", ntreekept,com.LHf);
   printf ("\n\nNo. of counts out of %d\n", irun);
   fprintf (fout, "\n%d LHs are collected from the MCMC run\n", ntreekept);
   fprintf (fout, "\n\nNo. of counts out of %d\n", irun);
   for (i=0,j=0; i<ntreekept; i++) {
      if (counts[i]==0) continue;
      IofLHs[j]=IofLHs[i];  lnPLHs[j]=lnPLHs[i];
      counts[j]=counts[i];  GetLHistoryI (IofLHs[j]);
      printf("\n%3d %9d %12.4f%9d %6.1f%% ",
           j+1,IofLHs[j],lnPLHs[j],counts[j],(double)counts[j]/irun);
      OutaTreeN(F0,0,0);
      fprintf(fout,"\n%3d %9d %12.4f%9d %6.1f%% ",
           j+1,IofLHs[j],lnPLHs[j],counts[j],(double)counts[j]/irun);
      OutaTreeN(fout,0,0);
      j++;
   }
   ntreekept=j;

   fLH=fopen(com.LHf, "w");
   fprintf (fLH, "%6d%6d\n", com.ns, ntreekept);
   for (j=0; j<ntreekept; j++) {
      GetLHistoryI (IofLHs[j]);
      fprintf(fLH,"%10d  ",IofLHs[j]);  OutaTreeN(fLH,0,0);  FPN(fLH);
   }
   fclose (fLH);  fflush(fout);
   lnpBestTree=-1e40;
   EvaluateLHs(fout, NULL, IofLHs, lnPLHs, ntreekept, com.delta1);
   free(lnPLHs);
   return ;
}

void EvaluateLHs(FILE*fout, FILE*fLH, int IofLHs[], double lnPLHs[],
    int ntreekept, double delta)
{
/* calculate relative posterior probabilities of LHs by more accurate
   calculation of their log{P}.
   if (fLH) the trees are read from a file.
*/
   char line[201];
   int j;
   double maxlnP=-1e99;
   
   printf("\n\nEvaluate LHs (delta=%.3f):\n\n", delta);
   fprintf(fout, "\n\nEvaluate LHs (delta=%.3f):\n\n", delta);
   if (fLH) {
      fscanf(fLH, "%d%d", &j, &ntreekept);
      if (j!=com.ns || ntreekept>MAXTREEKEPT)  error2 ("LHfile err.");
   }
   for (j=0; j<ntreekept; j++) {
      if (fLH) { fscanf(fLH,"%d",&IofLHs[j]); fgets(line,200,fLH); }
      GetLHistoryI (IofLHs[j]);
      nreplicateMC=0;
      lnPLHs[j]=OneTree(delta);

      printf ("#%3d/%3d %10d%12.4f  ", j+1,ntreekept,IofLHs[j],lnPLHs[j]);
      OutaTreeN(F0,0,0);  printf ("%14d\n", nreplicateMC);
      fprintf (fout, "#%3d%10d%12.4f  ", j+1,IofLHs[j],lnPLHs[j]);
      OutaTreeN(fout,0,0);  fprintf (fout, "%14d\n", nreplicateMC);
   }
   FOR (j,ntreekept) if (lnPLHs[j]>maxlnP) maxlnP=lnPLHs[j];
   FOR (j,ntreekept) lnPLHs[j]=exp(lnPLHs[j]-maxlnP);
   abyx (1/sum(lnPLHs,ntreekept), lnPLHs, ntreekept);

   puts ("\nprobabilities of LHs\n");
   fputs ("\nprobabilities of LHs\n",fout);
   for (j=0; j<ntreekept; j++) {
      GetLHistoryI (IofLHs[j]);  
      printf ("#%3d%10d%12.5f  ", j+1,IofLHs[j],lnPLHs[j]);
      OutaTreeN(F0,0,0);  OutaTreeB(F0);  FPN(F0);
      fprintf (fout, "#%3d%10d%12.5f  ", j+1,IofLHs[j],lnPLHs[j]);
      OutaTreeN(fout,0,0);  OutaTreeB(fout);  FPN(fout);
   }
}
