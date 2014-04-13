/* evolver.c
   Copyright, Ziheng Yang, April 1995.

     cc -fast -o evolver evolver.c tools.o -lm
     evolver
*/
/*
#define CODONLRT
*/
#include "tools.h"
#define NS            1000
#define NBRANCH       (NS*2-2)
#define NCODE         64

struct CommonInfo {
   char *z[2*NS-1], spname[NS][10], daafile[96];
   int ns,ls,npatt,*fpatt,np,ntime,ncode,clock,rooted,model,icode,cleandata;
   int seqtype, *pose;
   double kappa, omega, alpha, pi[64],*rates, *chunk, daa[20*20];
}  com;
struct TREEB {
   int nbranch, nnode, origin, branches[NBRANCH][2];
}  tree;
struct TREEN {
   int father, nson, sons[3], ibranch;
   double branch, divtime, omega, *lkl, daa[20*20];
}  nodes[2*NS-1];

extern char BASEs[],GenetCode[][64];

#define NODESTRUCTURE
#define BIRTHDEATH
#include "treesub.c"
#include "treespace.c"
void TreeDistances (void);
void Simulate (char*ctlf);
void MakeSeq(char*z, int ls);
int PMatExpQt (double Qt[], double space[]);
int GetPMatCodon(double t, double kappa, double omega, double pi[]);
int GetPMatAA (double t, double pi[]);

char *MCctlf[]={"MCbase.dat","MCcodon.dat","MCaa.dat"};

int main(int argc, char*argv[])
{
   FILE* ftree;
   char ctlf[96]="MCcodon.dat";
   int i, option=6, Ib[NS], rooted, wantbranch=0;
   double bfactor=1;

   if(argc>1) strcpy(ctlf,argv[1]);
   for (puts(""); ;) {
      printf("\n\t(1) Get a random UNROOTED tree?\n"); 
      printf("\t(2) Get a random ROOTED tree?\n"); 
      printf("\t(3) List all UNROOTED trees into file trees?\n");
      printf("\t(4) List all ROOTED trees into file trees?\n");
      printf("\t(5) Simulate nucleotide data sets (use %s)?\n",MCctlf[0]);
      printf("\t(6) Simulate codon data sets      (use %s)?\n",MCctlf[1]);
      printf("\t(7) Simulate amino acid data sets (use %s)?\n",MCctlf[2]);
      printf("\t(8) Calculate identical bi-partitions between trees?\n");
      printf("\t(0) Quit?\n");

#ifdef CODONLRT
      MCctlf[1]="MCcodonLRT.dat";
      option=6;
#else
      scanf("%d", &option);
#endif

      if (option==0) exit(0);
      if (option<5)  { printf ("No. of species: ");   scanf ("%d", &com.ns); }
      if (com.ns>NS) error ("No. of species too large.. raise NS.");

      rooted=!(option%2);
      if (option<3) {
         printf("\nrandom number seed? ");  scanf("%d",&i);  SetSeed(i);
         printf ("Want random (uniform) branch lengths (0/1): ");
         scanf ("%d", &wantbranch);
         if (wantbranch) { printf ("mean branch: "); scanf("%lf",&bfactor); }
      }
      switch (option) {
      case(1):
      case(2):
         for (i=0; i<com.ns-2-!rooted; i++) Ib[i]=(int)((3.+i)*rndu());
         MakeTreeIb (com.ns, Ib, rooted); 
         if (wantbranch) FOR (i,tree.nnode) nodes[i].branch=2.*bfactor*rndu();
         FPN(F0);  OutaTreeN(F0, 0, wantbranch);  FPN(F0);  FPN(F0);
         break;
      case(3):
      case(4): 
         if ((ftree=fopen("trees","w"))==NULL) error("file creation err");
         ListTrees(ftree, com.ns, rooted);
         fclose(ftree);
         break;
      case(5):  com.seqtype=0; Simulate(MCctlf[option-5]);   break;
      case(6):  com.seqtype=1; Simulate(MCctlf[option-5]);   break;
      case(7):  com.seqtype=2; Simulate(MCctlf[option-5]);   break;
      case(8):  TreeDistances ();                            break;
      default:  exit(0);
      }
   }
   return(0);
}


void TreeDistances (void)
{
   int i,j,ntree, k,*nib, parti2B[NS], nsame, IBsame[NS],nIBsame[NS], lparti=0;
   char treef[64], *partition, outf[32]="dT";
   FILE *ftree, *fout=(FILE*)fopen(outf,"w");
   double psame, mp, vp;

   puts("\nNumber of identical bi-partitions between trees.\nTree file name?");
   scanf ("%s", treef);
   if ((ftree=fopen (treef,"r"))==NULL) error ("Tree file not found.");
   if(fout==NULL) error("outfile creation error");
   fscanf (ftree, "%d%d", &com.ns, &ntree);

   if(ntree<2) error("ntree");
   printf ("\n%d species, %d trees\n", com.ns, ntree);
   puts("\n\t1: first vs. rest?\n\t2: all pairwise comparsons?\n");
   scanf("%d", &k);

   lparti=(com.ns-1)*com.ns*sizeof(char);
   i=(k==1?2:ntree)*lparti;
   printf("\n%d bytes of space requested.\n", i);
   partition=(char*)malloc(i);
   nib=(int*)malloc(ntree*sizeof(int));
   if (partition==NULL || nib==NULL) error("out of memory");

   if(k==2) {    /* pairwise comparisons */
     fputs("Number of identical bi-partitions in pairwise comparisons\n",fout);
      for (i=0; i<ntree; i++) {
         ReadaTreeN (ftree, &k, 1); 
         nib[i]=tree.nbranch-com.ns;
         BranchPartition(partition+i*lparti, parti2B);
      }
      for (i=0; i<ntree; i++,FPN(F0),FPN(fout)) {
         printf("%2d (%2d):", i+1,nib[i]);
         fprintf(fout,"%2d (%2d):", i+1,nib[i]);
         for (j=0; j<i; j++) {
            nsame=NSameBranch(partition+i*lparti,partition+j*lparti, 
                  nib[i],nib[j],IBsame);
            fprintf(fout," %2d", nsame);
         }
      }
   }
   else {  /* first vs. others */
      ReadaTreeN (ftree, &k, 1); 
      nib[0]=tree.nbranch-com.ns;
      if (nib[0]==0) error("1st tree is a star tree..");
      BranchPartition (partition, parti2B);
      fputs ("Comparing the first tree with the others\nFirst tree:\n",fout);
      OutaTreeN(fout,0,0);  FPN(fout);  OutaTreeB(fout);  FPN(fout); 
      fputs ("\nInternal branches in the first tree:\n",fout);
      FOR(i,nib[0]) { 
         k=parti2B[i];
         fprintf(fout,"%3d (%2d..%-2d): ( ",
            i+1,tree.branches[k][0]+1,tree.branches[k][1]+1);
         FOR(j,com.ns) if(partition[i*com.ns+j]) fprintf(fout,"%d ",j+1);
         fputs(")\n",fout);
      }
      if(nodes[tree.origin].nson<=2) 
         fputs("\nRooted tree, results may not be correct.\n",fout);
      fputs("\nCorrect internal branches compared with the 1st tree:\n",fout);
      FOR(k,nib[0]) nIBsame[k]=0;
      for (i=1,mp=vp=0; i<ntree; i++,FPN(fout)) {
         ReadaTreeN (ftree, &k, 1); 
         nib[1]=tree.nbranch-com.ns;
         BranchPartition (partition+lparti, parti2B);
         nsame=NSameBranch (partition,partition+lparti, nib[0],nib[1],IBsame);

         psame=nsame/(double)nib[0];
         FOR(k,nib[0]) nIBsame[k]+=IBsame[k];
         fprintf(fout,"1 vs. %3d: %4d: ", i+1,nsame);
         FOR(k,nib[0]) if(IBsame[k]) fprintf(fout," %2d", k+1);
         printf("1 vs. %5d: %6d/%d  %10.4f\n", i+1,nsame,nib[0],psame);
         vp += square(psame - mp)*(i-1.)/i;
         mp=(mp*(i-1.) + psame)/i;
      }
      vp=(ntree<=2 ? 0 : sqrt(vp/((ntree-1-1)*(ntree-1.))));
      fprintf(fout,"\nmean and S.E. of proportion of identical partitions\n");
      fprintf(fout,"between the 1st and all the other %d trees ", ntree-1);
      fprintf(fout,"(ignore these if not revelant):\n %.4f +- %.4f\n", mp, vp);
      fprintf(fout,"\nNumbers of times, out of %d, ", ntree-1);
      fprintf(fout,"interior branches of tree 1 are present");
      fputs("\n(This may be bootstrap support for nodes in tree 1)\n",fout);
      FOR(k,nib[0]) { 
         i=tree.branches[parti2B[k]][0]+1;  j=tree.branches[parti2B[k]][1]+1; 
         fprintf(fout,"%3d (%2d..%-2d): %6d (%5.1f%%)\n",
            k+1,i,j,nIBsame[k],nIBsame[k]*100./(ntree-1.));
      }
   }
   printf("\nResults in file %s\n", outf);
   free(partition);  free(nib);  fclose(ftree);  fclose(fout);
}





static double Pcode[64*64], spaceMC[64*64*3];

int GetPMatCodon(double t, double kappa, double omega, double pi[])
{
/* Construct the rate matrix Q[] and then get P(t), Pcode[], from PMatExpQt().
   64 codon are used, and stop codons have 0 freqs.
*/
   int n=com.ncode, i,j,k, c[2],ndiff,pos=0,from[3],to[3];
   double mr, *Q=spaceMC, *space=Q+n*n;
   
   FOR (i,n*n) Q[i]=0;
   for (i=0; i<n; i++) FOR (j,i) {
      from[0]=i/16; from[1]=(i/4)%4; from[2]=i%4;
      to[0]=j/16;   to[1]=(j/4)%4;   to[2]=j%4;
      c[0]=GenetCode[com.icode][i];   c[1]=GenetCode[com.icode][j];
      if (c[0]*c[1]==0)  continue;
      for (k=0,ndiff=0; k<3; k++)  if (from[k]!=to[k]) { ndiff++; pos=k; }
      if (ndiff!=1)  continue;
      Q[i*n+j]=1;
      if ((from[pos]+to[pos]-1)*(from[pos]+to[pos]-5)==0)  Q[i*n+j]*=kappa;
      if(c[0]!=c[1])  Q[i*n+j]*=omega;
      Q[j*n+i]=Q[i*n+j];
   }
   FOR(i,n) FOR(j,n) Q[i*n+j]*=com.pi[j];
   for (i=0,mr=0;i<n;i++) { 
      Q[i*n+i]=-sum(Q+i*n,n); mr-=pi[i]*Q[i*n+i]; 
   }
   FOR(i,n*n) Q[i]*=t/mr;

   PMatExpQt (Q, space);

   return (0);
}


int PMatExpQt (double Qt[], double space[])
{
/* This calculates the transition probability matrix P(t), in Pcode,
   using polynomial expansion. 
    P(t) = I + Qt + (Qt)^2/2 + (Qt)^3/3! + ...
*/
   int nterms=100, n=com.ncode, i,j;
   double *P=Pcode, *Q=Qt, *Q1=space, *Q2=Q1+n*n, mr,div,small=1e-6;

   xtoy(Q,P,n*n);  FOR(i,n) P[i*n+i]++;   /* I + Qt */
   xtoy(Q,Q1,n*n);
   for (i=2,div=2;i<nterms;i++,div*=i) {  /* k is divisor */
      matby(Q, Q1, Q2, n, n, n);
      for(j=0,mr=0;j<n*n;j++) { P[j]+=Q2[j]/div; mr+=fabs(Q2[j]); }
      mr/=div;
      if (mr<small) break;
      xtoy(Q2,Q1,n*n);
   }

   if(testTransP(P, n)) puts("err P..");

   return(0);
}



void MakeSeq(char*z, int ls)
{
/* makes a codon or AA sequence using com.pi[].  Codons are encoded as 
   0,1,2,...,63
*/
   int j,h, n=com.ncode;
   double p[64],r;

   FOR(j,n) p[j]=com.pi[j];
   for (j=1; j<n; j++) p[j]+=p[j-1];
   if (fabs(p[n-1]-1) > 1e-6) 
      { printf("\nsum pi = %.6f != 1!\n",p[n-1]); exit(-1); }
   for (h=0; h<com.ls; h++) {
      for(j=0,r=rndu();j<n;j++) if(r<p[j]) break;
      z[h]=(char)j;
   }
}

void Evolve1 (int inode)
{
/* evolve sequence com.z[tree.origin] along the tree to generate com.z[], 
   using nodes[].branch, nodes[].omega, & com.model
   Needs com.z[0,1,...,nnode-1], while com.z[0] -- com.z[ns-1] constitute
   the data.
   For codon or amino acid seqs.
*/
   int is, h,i,j, ison, from, n=com.ncode;
   double t, r;
   
   for (is=0; is<nodes[inode].nson; is++) {
      ison=nodes[inode].sons[is];
      FOR(h,com.ls) com.z[ison][h]=com.z[inode][h];
      if((t=nodes[ison].branch)<1e-7) continue;
      if(com.seqtype==BASEseq) 
         EvolveHKY85(com.z[ison],com.z[ison],com.ls,t,com.rates,
                     com.pi,com.kappa,1);
      else 
         FOR (h,com.ls) {
            if (h==0 || (com.rates && com.rates[h]!=com.rates[h-1])) {
               r=(com.rates?com.rates[h]:1);
               if(com.seqtype==CODONseq) 
                  GetPMatCodon(t*r,com.kappa,nodes[ison].omega,com.pi);
               else 
                  GetPMatAA(t*r,com.pi);
               FOR(i,n) for(j=1; j<n; j++)  Pcode[i*n+j]+=Pcode[i*n+j-1];
	    }
            for(j=0,from=com.z[ison][h],r=rndu(); j<n; j++)
               if(r<Pcode[from*n+j]) break;
            com.z[ison][h]=j;
         }
      if (nodes[ison].nson) Evolve1(ison); 
   }  /* for (is) */
}


enum {Poisson, EqualInput, Empirical, Empirical_F} AAModel;
char *aamodels[]={"Poisson", "EqualInput", "Empirical", "Empirical_F"};


int GetDaa (FILE* fout, double daa[])
{
/* Get the amino acid substitution rate matrix (grantham, dayhoff, jones, etc).
*/
   FILE * fdaa;
   char aa3[4]="";
   int i,j, naa=20;

   if ((fdaa=fopen(com.daafile, "r"))==NULL) 
      { printf("\nAA dist file %s not found.", com.daafile); exit(-1); }
   printf("\nReading rate matrix from %s\n", com.daafile);

   for (i=0; i<naa; i++)  for (j=0,daa[i*naa+i]=0; j<i; j++)  {
      fscanf(fdaa, "%lf", &daa[i*naa+j]);
      daa[j*naa+i]=daa[i*naa+j];
   }
   if (com.model==Empirical) {
      FOR(i,naa) if(fscanf(fdaa,"%lf",&com.pi[i])!=1) error("err aaRatefile");
      if (fabs(1-sum(com.pi,20))>1e-4) error("\nSum of aa freq. != 1\n");
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

   return (0);
}


int GetPMatAA (double t, double pi[])
{
/* Construct the rate matrix Q[] and then get P[] from PMatExpQt().
*/
   int naa=20, i,j;
   double mr, *Q=spaceMC, *space=Q+naa*naa;

   FOR (i,naa*naa) Q[i]=0;
   switch (com.model) {
   case (Poisson)   : case (EqualInput) : 
      fillxc (Q, 1., naa*naa);  break;
   case (Empirical)   : case (Empirical_F):
      FOR(i,naa) FOR(j,i) Q[i*naa+j]=Q[j*naa+i]=com.daa[i*naa+j]/100;
      break;
   }
   FOR (i,naa) FOR (j,naa) Q[i*naa+j]*=com.pi[j];
   for (i=0,mr=0; i<naa; i++) {
      Q[i*naa+i]=0; Q[i*naa+i]=-sum(Q+i*naa,naa);  mr-=com.pi[i]*Q[i*naa+i]; 
   }
   FOR(i,naa*naa) Q[i]*=t/mr;

   PMatExpQt (Q, space);

   return (0);
}


void Simulate (char*ctlf)
{
/* simulate nr data sets of nucleotide, codon, or AA sequences.
   ls: number of codons
   All 64 codons are used.
*/
   enum {PAML,PAUP} SEQDATAFORMAT;
   FILE *fin, *fseq, *ftree;
   char seqf[32]="mc.paml";
   char *paupstart="paupstart",*paupblock="paupblock",*paupend="paupend";
   int i,j,k, h,ir,n,nr,fixtree=1,size=10000,format=PAML,rooted=0,b[3], ncatG;
   double birth=0, death=0, sample=1, mut=1, t, *space;

   com.alpha=0; com.cleandata=1;

   printf("\nReadind options from data file %s.\n", ctlf);
   com.ncode=n=(com.seqtype==0 ? 4 : (com.seqtype==1?64:20));
   if((fin=(FILE*)fopen(ctlf,"r"))==NULL) 
      { printf("\ndata file %s does not exist?\n", ctlf);  exit(-1); }
   fscanf (fin, "%d", &format);
   sprintf(seqf, "mc.%s", format?"paup":"paml");
   printf("\nSimulated data will go into %s.\n", seqf);
   if(format==PAUP) printf("%s, %s, & %s will be appended if existent\n",
                       paupstart,paupblock,paupend);

   fscanf (fin, "%d%d%d%d", &i,&com.ns, &com.ls, &nr);
   SetSeed(i);
   if(com.ns>NS) error("too many seqs?");
   printf ("\n%d seqs, %d sites, %d replicates\n", com.ns,com.ls,nr);
   k=(com.ns*com.ls* (com.seqtype==CODONseq?4:1) *nr)/1024+1;
   printf ("Seq file will be about %dK bytes.",k);

   if (fixtree) {
      if(ReadaTreeN(fin,&i,1))  error("err tree..");
      if(!i) {
         FOR(i,tree.nbranch)
            fscanf(fin,"%lf",&nodes[tree.branches[i][1]].branch);
      }
      printf("\nModel tree & branch lengths:\n"); OutaTreeN(F0,0,1); FPN(F0);
      FPN(F0); OutaTreeB(F0); FPN(F0);
      FOR(i,tree.nbranch) printf("%9.5f",nodes[tree.branches[i][1]].branch); 
      if(com.seqtype==CODONseq) {
#ifdef CODONLRT
         FOR(i,tree.nbranch)
            fscanf(fin,"%lf",&nodes[tree.branches[i][1]].omega);
         printf("\nand dN/dS ratios:\n"); 
         FOR(i,tree.nbranch) printf("%9.5f",nodes[tree.branches[i][1]].omega);
         FPN(F0);
#endif
      }
   }
   else {   /* random trees */
      printf ("\nbirth rate, death rate, sampling fraction, and ");
      printf ("mutation rate (tree height)?\n");
      fscanf (fin, "%lf%lf%lf%lf", &birth, &death, &sample, &mut);
   }

   if(com.seqtype==BASEseq) {
      fscanf(fin,"%lf",&com.kappa);
      printf("\nkappa = %9.5f\n",com.kappa);
   }
   else if(com.seqtype==CODONseq) {
#ifndef CODONLRT
      fscanf(fin,"%lf",&com.omega);
      printf("\nomega = %9.5f\n",com.omega);
      FOR(i,tree.nbranch) nodes[tree.branches[i][1]].omega=com.omega;
#endif
      fscanf(fin,"%lf",&com.kappa);
      printf("\nkappa = %9.5f\n",com.kappa);
   }

   if(com.seqtype==BASEseq || com.seqtype==AAseq) {
      fscanf(fin,"%lf%d", &com.alpha, &ncatG);
      if (com.alpha) {
         if((com.rates=(double*)malloc(com.ls*sizeof(double)))==NULL) 
             error("oom for rates");
         printf("\nGamma rates at sites, alpha =%.4f (K=%d)\n",com.alpha,ncatG);
      }
   }
   if(!com.alpha)  puts("\n\nRates are constant over sites.");

   /* get aa substitution model and rate matrix */
   if(com.seqtype==AAseq) {
      fscanf(fin,"%d",&com.model);
      printf("\nmodel: %s",aamodels[com.model]); 
      if(com.model>=2)  { fscanf(fin,"%s",com.daafile); GetDaa(NULL,com.daa); }
   }
   /* get freqs com.pi[] */
   if(com.seqtype!=AAseq || com.model==1 || com.model==3) 
      FOR(k,com.ncode) fscanf(fin,"%lf",&com.pi[k]);
   else if(com.model==0) 
      fillxc(com.pi,1./com.ncode,com.ncode);
   printf("\nsum pi = 1 = %.6f:", sum(com.pi,com.ncode));
   matout2(F0,com.pi,1,com.ncode,7,4);


   puts("\nAll parameters are read, now ready to simulate.\n");
   FOR(i,tree.nnode)sprintf(com.spname[i],"seq.%d",i+1);
   FOR(j,com.ns*2-1) com.z[j]=(char*)malloc(com.ls*sizeof(char));
   size=max(size, com.ls*com.ns*sizeof(char));
   space=(double*)malloc(size);
   if (com.z[com.ns*2-1-1]==NULL || space==NULL) error("oom for space");
   if((fseq=fopen(seqf,"w"))==NULL) error("seq file creation error.");
   if(format==PAUP) appendfile(fseq,paupstart);
   for (ir=0; ir<nr; ir++) {
      if (!fixtree) {    /* right now tree is fixed */
         RandomLHistory (rooted, space);
         if (rooted && com.ns<10) j=GetIofLHistory ();
         BranchLengthBD (1, birth, death, sample, mut);
         if (rooted && com.ns<10) printf ("\ntree used (LH #%d):\n", j);
         else                     printf ("\ntree used: "); 
         OutaTreeN(F0,0,1);  puts(";");
      }
      MakeSeq(com.z[tree.origin],com.ls);

      if (com.alpha) 
         Rates4Sites (com.rates,com.alpha,ncatG,com.ls, 0,space);

      Evolve1(tree.origin);

      if (format==PAUP) fprintf(fseq,"\n\n[Replicate # %d]\n", ir+1);
      printSeqs(fseq, NULL, NULL, format);
      if(format==PAUP && !fixtree) {
         fprintf(fseq,"\nbegin tree;\n   tree true_tree = [&U] "); 
         OutaTreeN(fseq,1,1); fputs(";\n",fseq);
         fprintf(fseq,"end;\n\n");
      }
      if(format==PAUP) appendfile(fseq,paupblock);
      printf ("did data set %d.\n", ir+1);
   }   /* for (ir) */

   if(format==PAUP) appendfile(fseq,paupend);
   fclose(fseq); 
   FOR(j,com.ns*2-1) free(com.z[j]);
   free(space);
   exit (0);
}
