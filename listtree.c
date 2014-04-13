/* listtree.c
   Copyright, Ziheng Yang, April 1995.

     cc -O2 -o listtree listtree.c tools.o -lm
     listtree
*/

#include "tools.h"
#ifdef SMALLMACHINE
#define NS            200
#else
#define NS            500
#endif
#define NBRANCH       (NS*2-2)
#define NCODE         4

struct CommonInfo {
   char *z[2*NS-1], spname[NS][10];
   int ns, ls, npatt, *fpatt, np, ntime, ncode, clock, rooted, model;
   double kappa, alpha, pi[6],*rates, *chunk;
}  com;
struct TREEB {
   int nbranch, nnode, origin, branches[NBRANCH][2];
}  tree;
struct TREEN {
   int father, nson, sons[3], ibranch;
   double branch, divtime, *lkl;
}  nodes[2*NS-1];

char *models[]={"JC69","K80","F81","F84","HKY85","TN93","REV","UNREST"};
enum {JC69, K80, F81, F84, HKY85, TN93, REV, UNREST} MODELS;

#define NODESTRUCTURE
#define EVOLVE
#include "treesub.c"
#include "treespace.c"
int TreeDistances (void);
int SimulateData (void);

int main (void)
{
   FILE* ftree;
   int i, nopt=7, option, Ib[NS], rooted, wantbranch=0;
   double bfactor=1;

   for (; ;) {
      puts ("\n\n");
      puts ("\t\t(1) Get a random UNROOTED tree?\n"); 
      puts ("\t\t(2) Get a random ROOTED tree?\n"); 
      puts ("\t\t(3) List all UNROOTED trees into file trees?\n"); 
      puts ("\t\t(4) List all ROOTED trees into file trees?\n");
      puts ("\t\t(5) Simulate data sets?\n");
      puts ("\t\t(6) Calculate identical bi-partitions between trees?\n");
      puts ("\t\t(0) Quit?\n");

      scanf("%d", &option);
      if (option==0) exit(0);
      else if (option<1 || option>nopt) error ("option.");

      if (option<5)  { printf ("No. of species: ");   scanf ("%d", &com.ns); }
      if (com.ns>NS) error ("No. of species too large.. raise NS.");

      rooted=!(option%2);
      if (option==5) SimulateData();
      else if (option<3) {
         printf ("Want random (uniform) branch lengths (0/1): ");
         scanf ("%d", &wantbranch);
         if (wantbranch) { printf ("mean branch: "); scanf("%lf",&bfactor); }
      }
      switch (option) {
      case (1):   case (2):
         for (i=0; i<com.ns-2-!rooted; i++) Ib[i]=(int)((3.+i)*rndu());
         MakeTreeIb (com.ns, Ib, rooted); 
         if (wantbranch) FOR (i,tree.nnode) nodes[i].branch=2.*bfactor*rndu();
         FPN(F0);  OutaTreeN(F0, 0, wantbranch);  FPN(F0);  FPN(F0);
         break;
      case (3):   case (4): 
         if ((ftree=fopen("trees","w"))==NULL) error("file creation err");
         ListTrees (ftree, com.ns, rooted);
         fclose (ftree);
         break;
      case (6): TreeDistances ();  break;
      }
exit (0);
   }
   return(0);
}


int TreeDistances (void)
{
   int itree,ntree, nib[2], parti2B[NS], nsame;
   char treef[64], *partition[2];
   FILE *ftree;
   double psame, mp, vp;

   printf ("\nNumber and proportion of identical bi-partitions ");
   printf ("between the 1st tree and the others\n\nTreeFileName? ");
   scanf ("%s", treef);
   if ((ftree=fopen (treef,"r"))==NULL) error ("Tree file not found.");
   fscanf (ftree, "%d%d", &com.ns, &ntree);

   partition[0]=(char*)malloc(2*(com.ns-1)*com.ns*sizeof(char));
   partition[1]=partition[0]+(com.ns-1)*com.ns;
   if (partition[0]==NULL) error("oom");

   for (itree=0,mp=vp=0; itree<ntree; itree++) {
      ReadaTreeN (ftree, &nsame, 1);   /* nsame used for other purpose */
/*
      printf ("  %6d%6d  ", tree.nbranch, tree.nnode);
      OutaTreeN (F0, 0, 0);  FPN(F0); 
*/
      nib[itree>0]=tree.nbranch-com.ns;
      BranchPartitionLarge (partition[itree>0], parti2B);
      if (itree) {
         nsame=NSameBranchLarge (partition[0],partition[1], nib[0],nib[1]);
         psame=nsame/(double)nib[0];
         printf ("1 vs. %5d: %6d%10.4f\n", itree+1, nsame, psame);
         vp += square(psame - mp)*(itree-1.)/itree;
         mp=(mp*(itree-1.) + psame)/itree;
      }
      else if (nib[0]==0) puts("1st tree is a star tree..");
   }
   free (partition[0]);   fclose(ftree);
   vp=sqrt(vp/((ntree-1-1)*(ntree-1.)));
   printf ("\nmean and S.E. of proportion of identical partitions\n");
   printf ("between the 1st and all the other %d trees ", ntree-1);
   printf ("(ignore these if not revelant):\n");
   printf (" %.4f +- %.4f\n", mp, vp);

   return (0);
}



int SimulateData (void)
{
/* simulate nr data sets, using one fixed tree, or using a random tree from 
   the birth death process with species sampling.
*/
   FILE *fseq, *ftree, *fpaupb=NULL;
   char seqf[32], treef[32], *paupblock="paupblock";
   int i,j, ir, nr, fixtree=0, ncatG, size=10000, format=1;
   double birth=0, death=0, sample=1, mut=1, *space;

   printf ("\nrandom number generator seed? ");
   scanf ("%d", &i);  SetSeed (i);
   printf ("\n# species, # sites, # replicates? ");
   scanf ("%d%d%d", &com.ns, &com.ls, &nr);
   printf ("\nThe file will be about %d bytes.\n", com.ns*com.ls*nr);

   sprintf(treef, "trees.%ds", com.ns);
   printf ("\ntree from file %s (1) or birth-death process (0)? ", treef);
   scanf ("%d", &fixtree);
   if (fixtree) {
      if ((ftree=fopen(treef,"r"))==NULL) error ("tree file error.");
      fscanf (ftree, "%d%d", &i, &j); 
      if (i!=com.ns || ReadaTreeN(ftree, &i, 1))  error ("err tree..");
      printf("\nModel tree (stop if incorrect): "); OutaTreeN(F0,0,1); FPN(F0);
   }
   else {
      printf ("\nbirth rate, death rate, sampling fraction, and ");
      printf ("mutation rate (tree height)?\n");
      scanf ("%lf%lf%lf%lf", &birth, &death, &sample, &mut);
   }
   
   printf ("\nsubstitution model (0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85)? ");
   scanf ("%d", &com.model);
   com.kappa=1;
   if (com.model!=0 && com.model!=2) {
      printf ("\nkappa (transition/transversion rate ratio)? ");
      scanf ("%lf", &com.kappa);
   }
   fillxc(com.pi, .25, 4);
   if (com.model>=2) {
      printf ("\nbase frequencies? ");
      FOR (i,4) scanf ("%lf", &com.pi[i]);
   }
   com.pi[4]=com.pi[0]+com.pi[1];    com.pi[5]=com.pi[2]+com.pi[3]; 

   printf ("\nalpha (0 for one rate), # categories (0 for continuous)? ");
   scanf ("%lf%d", &com.alpha, &ncatG);

   printf ("\noutput format (0 for PAML, 1 for PAUP)? ");
   scanf ("%d", &format);
   sprintf (seqf, "simulateddata.%s", format?"paup":"paml");
   printf ("\n\nSimulated data will go into %s (Return to continue).\n", seqf);
   if(format && (fpaupb=fopen(paupblock,"r")))
      printf("File %s is appended after each data set\n",paupblock);

   fseq=fopen(seqf, "w");

   if (com.alpha)  { 
      com.rates=(double*)malloc(com.ls*sizeof(double));
      if (com.rates==NULL) error("oom for rates");
   }

   FOR (j, com.ns*2-1) com.z[j]=(char*)malloc((com.ls+1)*sizeof(char));
   size=max(size, (com.ls+1)*com.ns*sizeof(char));
   space=(double*)malloc(size);
   if (com.z[com.ns*2-1-1]==NULL || space==NULL) error("oom for space");
   for (i=0;i<com.ns;i++) sprintf(com.spname[i], "seq.%d", i+1);

   for (ir=0,com.ncode=4,com.clock=1; ir<nr; ir++) {
      if (!fixtree) {
         RandomLHistory (1, space);
         if (com.ns<10) j=GetIofLHistory ();
         BranchLengthBD (1, birth, death, sample, mut);
         if (com.ns<10) printf ("\ntree used (LH #%d):\n", j);
         else           printf ("\ntree used: "); 
         OutaTreeN (F0,0,1);  FPN(F0);
      }
      if (com.alpha) 
         Rates4Sites (com.rates,com.alpha,ncatG,com.ls, !com.model,space);

      GenerateSeq ();
      OutSeqs (fseq, 0, NULL, format);
      if(format && fpaupb) {
         for (;;) {
            if((j=fgetc(fpaupb))==EOF) break; 
            fputc (j, fseq); 
	 }
         rewind (fpaupb);
      }
      printf ("did data set %d\n", ir+1);
/*
 printsma (F0, com.z, com.ns*2-1, com.ls, 60, 10, 1, 0, 0);
 FOR (i,com.ls) com.z[tree.origin][i]+='1';
 printaSeq (F0, com.z[tree.origin], com.ls, 60, 10);
*/
   }
   fclose (fseq);  if (format && fpaupb) fclose (fpaupb);
   FOR (j, com.ns*2-1) free(com.z[j]);
   free (space);
   return (0);
}
