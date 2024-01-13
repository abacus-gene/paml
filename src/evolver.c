/* evolver.c
   Copyright, Ziheng Yang, April 1995.

     cl -Ox evolver.c tools.c
     cl -Ox -DCodonNSbranches    -FeevolverNSbranches.exe    evolver.c tools.c
     cl -Ox -DCodonNSsites       -FeevolverNSsites.exe       evolver.c tools.c
     cl -Ox -DCodonNSbranchsites -FeevolverNSbranchsites.exe evolver.c tools.c

     cc -O3 -o evolver evolver.c tools.c -lm
     cc -O4 -DCodonNSbranches -o evolverNSbranches evolver.c tools.c -lm
     cc -O4 -DCodonNSsites -o evolverNSsites evolver.c tools.c -lm
     cc -O4 -DCodonNSbranchsites -o evolverNSbranchsites evolver.c tools.c -lm

     evolver
     evolver 5 MCbase.dat
     evolver 6 MCcodon.dat
     evolver 7 MCaa.dat
     evolver 9 <TreesFile> <MainTreeFile>
*/

/*
#define CodonNSbranches
#define CodonNSsites
#define CodonNSbranchsites
*/

#include "paml.h"

#define NS            5000
#define NBRANCH       (NS*2-2)
#define MAXNSONS      20
#define LSPNAME       96
#define NCODE         64
#define NCATG         40


struct CommonInfo {
   char* z[2 * NS - 1];
   char* spname[NS], daafile[512], cleandata, readpattern;
   int ns, ls, npatt, np, ntime, ncode, clock, rooted, model, icode;
   int seqtype, * pose, ncatG, NSsites;
   int ngene, lgene[1], posG[1 + 1];  /* not used */
   double piG[1][4], rgene[1];  /* not used */
   double* fpatt, kappa, omega, alpha, pi[64], * conP, daa[20 * 20];
   double freqK[NCATG], rK[NCATG];
   char* siteID;        /* used if ncatG>1 */
   double* siterates;   /* rates for gamma or omega for site or branch-site models */
   double* omegaBS, * QfactorBS;     /* omega IDs for branch-site models */
}  com;


struct TREEB {
   int nbranch, nnode, root, branches[NBRANCH][2];
}  tree;

#if(1)
struct TREEN {
   int father, nson, sons[MAXNSONS], ibranch;
   double branch, age, omega, label, label2, * conP;
   char * name, * annotation, fossil;
}  *nodes;
#else
struct NODE {
   int father, nson, sons[2], ibranch;
   double branch, age, omega, label, label2, * conP;
   char* annotation, fossil;
};

struct TREE {
   struct NODE* nodes;
   int root;
};
#endif

extern char BASEs[];
extern int GeneticCode[][64], noisy;
int LASTROUND = 0; /* not used */

#define EVOLVER
#define NODESTRUCTURE
#define BIRTHDEATH
#include "treesub.c"
#include "treespace.c"

void TreeDistances(FILE* fout);
void Simulate(char* ctlf);
void MakeSeq(char* z, int ls);
int EigenQbase(double rates[], double pi[], double Root[], double U[], double V[], double Q[]);
int EigenQcodon(int getstats, double kappa, double omega, double pi[], double Root[], double U[], double V[], double Q[]);
int EigenQaa(double pi[], double Root[], double U[], double V[], double Q[]);
void CladeMrBayesProbabilities(char treefile[]);
int between_f_and_x(void);
void LabelClades(FILE* fout);

char* MCctlf0[] = { "MCbase.dat","MCcodon.dat","MCaa.dat" };
char* seqf[] = { "mc.txt", "mc.txt", "mc.nex", "mc.nex" };

enum { JC69, K80, F81, F84, HKY85, T92, TN93, REV } BaseModels;
char* basemodels[] = { "JC69","K80","F81","F84","HKY85","T92","TN93","REV" };
enum { Poisson, EqualInput, Empirical, Empirical_F } AAModels;
char* aamodels[] = { "Poisson", "EqualInput", "Empirical", "Empirical_F" };


double PMat[NCODE * NCODE], U[NCODE * NCODE], V[NCODE * NCODE], Root[NCODE];
static double Qfactor = -1, Qrates[5];  /* Qrates[] hold kappa's for nucleotides */


int main(int argc, char* argv[])
{
   char* MCctlf = NULL, outf[512] = "evolver.out", treefile[512] = "mcmc.txt", maintreefile[512] = "\0";
   int i, option = -1, ntree = 1, rooted, BD = 0, gotoption = 0, pick1tree = 0;
   double birth = -1, death = -1, sample = -1, mut = -1, * space;
   FILE* fout = zopen(outf, "w");

   printf("EVOLVER in %s\n", pamlVerStr);
   com.alpha = 0; com.cleandata = 1; com.model = 0; com.NSsites = 0;

   if (argc > 2 && !strcmp(argv[argc - 1], "--stdout-no-buf"))
      setvbuf(stdout, NULL, _IONBF, 0);
   if (argc > 1) {
      gotoption = 1;   sscanf(argv[1], "%d", &option);
   }
   if (argc == 1)
      printf("Results for options 1-4 & 8 go into %s\n", outf);
   else if (option != 5 && option != 6 && option != 7 && option != 9) {
      puts("Usage: \n\tevolver \n\tevolver option datafile"); exit(-1);
   }
   if (option >= 4 && option <= 6)
      MCctlf = argv[2];
   else if (option == 9) {
      strcpy(treefile, argv[2]);
      if (argc > 3) strcpy(maintreefile, argv[3]);
      if (argc > 4) sscanf(argv[4], "%d", &pick1tree);
   }

   for (i = 0; i < NS; i++) {
      com.spname[i] = (char*)malloc(LSPNAME * sizeof(char*));
      if (com.spname[i] == NULL) zerror("oom");
   }

#if defined (CodonNSbranches)
   option = 6;  com.model = 1;
   MCctlf = (argc == 3 ? argv[2] : "MCcodonNSbranches.dat");
   gotoption = 1;
#elif defined (CodonNSsites)
   option = 6;  com.NSsites = 3;
   MCctlf = (argc == 3 ? argv[2] : "MCcodonNSsites.dat");
   gotoption = 1;
#elif defined (CodonNSbranchsites)
   option = 6;  com.model = 1; com.NSsites = 3;
   MCctlf = (argc == 3 ? argv[2] : "MCcodonNSbranchsites.dat");
   gotoption = 1;
#endif

   if (!gotoption) {
      for (; ;) {
         fflush(fout);
         printf("\n\t(1) Get random UNROOTED trees?\n");
         printf("\t(2) Get random ROOTED trees?\n");
         printf("\t(3) List all UNROOTED trees?\n");
         printf("\t(4) List all ROOTED trees?\n");
         printf("\t(5) Simulate nucleotide data sets (use %s)?\n", MCctlf0[0]);
         printf("\t(6) Simulate codon data sets      (use %s)?\n", MCctlf0[1]);
         printf("\t(7) Simulate amino acid data sets (use %s)?\n", MCctlf0[2]);
         printf("\t(8) Calculate identical bi-partitions between trees?\n");
         printf("\t(9) Calculate clade support values (evolver 9 treefile maintreefile <pick1tree>)?\n");
         printf("\t(11) Label clades?\n");
         printf("\t(0) Quit?\n");

         option = 5;
         scanf("%d", &option);

         if (option == 0) exit(0);
         if (option >= 5 && option <= 7) break;
         if (option < 5) {
            printf("No. of species: ");
            scanf("%d", &com.ns);
         }
         if (com.ns > NS) zerror("Too many species.  Raise NS.");
         if ((space = (double*)malloc(10000 * sizeof(double))) == NULL) zerror("oom");
         rooted = !(option % 2);
         if (option < 3) {
            printf("\nnumber of trees & random number seed? ");
            scanf("%d%d", &ntree, &i);
            SetSeed(i, 1);
            printf("Want branch lengths from the birth-death process (0/1)? ");
            scanf("%d", &BD);
         }
         if (option <= 4) {
            if (com.ns < 3) zerror("no need to do this?");
            i = (com.ns * 2 - 1) * sizeof(struct TREEN);
            if ((nodes = (struct TREEN*)malloc(i)) == NULL)
               zerror("oom");
         }
         switch (option) {
         case(1):   /* random UNROOTED trees */
         case(2):   /* random ROOTED trees */
            /* default names */
            if (com.ns <= 52)
               for (i = 0; i < com.ns; i++)  sprintf(com.spname[i], "%c", (i < 26 ? 'A' + i : 'a' + i - 26));
            else
               for (i = 0; i < com.ns; i++)  sprintf(com.spname[i], "S%d", i + 1);

            if (BD) {
               printf("\nbirth rate, death rate, sampling fraction, and ");
               printf("mutation rate (tree height)?\n");
               scanf("%lf%lf%lf%lf", &birth, &death, &sample, &mut);
            }
            for (i = 0; i < ntree; i++) {
               RandomLHistory(rooted, space);
               if (BD)
                  BranchLengthBD(1, birth, death, sample, mut);
               if (com.ns < 20 && ntree < 10) { OutTreeN(F0, 0, BD); puts("\n"); }
               OutTreeN(fout, 1, BD);  fprintf(fout, "\n");
            }
            /*
            for (i=0; i<com.ns-2-!rooted; i++)
               Ib[i] = (int)((3.+i)*rndu());
            MakeTreeIb (com.ns, Ib, rooted);
            */
            break;
         case(3):
         case(4):
            ListTrees(fout, com.ns, rooted);
            break;
         case(8):  TreeDistances(fout);  break;
         case(9):
            printf("tree file names? ");
            scanf("%s%s", treefile, maintreefile);
            break;
         case(10): between_f_and_x();    break;
         case(11): LabelClades(fout);    break;
         default:  exit(0);
         }
      }
   }

   if (option >= 5 && option <= 7) {
      com.seqtype = option - 5;  /* 0, 1, 2 for bases, codons, & amino acids */
      Simulate(MCctlf ? MCctlf : MCctlf0[option - 5]);
   }
   else if (option == 9) {

      CladeSupport(fout, treefile, 1, maintreefile, pick1tree);
      /* CladeMrBayesProbabilities("\data\Lakner2008SB/M1510YY03.nex.t"); */
   }
   return(0);
}


int between_f_and_x(void)
{
   /* this helps with the exponential transform for frequency parameters */
   int i, n, fromf = 0;
   double x[100];

   for (;;) {
      printf("\ndirection (0:x=>f; 1:f=>x; -1:end)  &  #classes? ");
      scanf("%d", &fromf);
      if (fromf == -1) return(0);
      scanf("%d", &n);  if (n > 100) zerror("too many classes");
      printf("input the first %d values for %s? ", n - 1, (fromf ? "f" : "x"));
      for (i = 0; i < n - 1; i++) scanf("%lf", &x[i]);
      x[n - 1] = (fromf ? 1 - sum(x, n - 1) : 0);
      f_and_x(x, x, n, fromf, 1);
      matout(F0, x, 1, n);
   }
}


void LabelClades(FILE* fout)
{
   /* This reads in a tree and scan species names to check whether they form a
      paraphyletic group and then label the clade.
      It assumes that the tree is unrooted, and so goes through two rounds to check
      whether the remaining seqs form a monophyletic clade.
   */
   FILE* ftree;
   int unrooted = 1, iclade, sizeclade = -1, mrca = -1, paraphyl, is, imrca, i, j, k, lasts=-1, haslength;
   char key[96] = "A", treef[64] = "/A/F/flu/HA.all.prankcodon.tre", * p, chosen[NS], * endstr = "end";
   int* anc[NS - 1], loc, bitmask, SI = sizeof(int) * 8;
   int debug;

   printf("Tree file name? ");
   scanf("%s", treef);
   printf("Treat tree as unrooted (0 no, 1 yes)? ");
   scanf("%d", &unrooted);

   ftree = zopen(treef, "r");
   fscanf(ftree, "%d%d", &com.ns, &j);
   if (com.ns <= 0) zerror("need ns in tree file");
   debug = (com.ns < 20);

   i = (com.ns * 2 - 1) * sizeof(struct TREEN);
   if ((nodes = (struct TREEN*)malloc(i)) == NULL) zerror("oom");
   for (i = 0; i < com.ns * 2 - 1; i++)  nodes[i].annotation = NULL;
   for (i = 0; i < com.ns - 1; i++) {
      anc[i] = (int*)malloc((com.ns / SI + 1) * sizeof(int));
      if (anc[i] == NULL)  zerror("oom");
   }
   ReadTreeN(ftree, &haslength, 1, 0);
   fclose(ftree);
   if (debug) { OutTreeN(F0, 1, PrNodeNum);  printf("\n"); }

   for (iclade = 0; iclade < com.ns - 1; iclade++) {
      printf("\nString for selecting sequences (followed by non-digit) (end to end)? ");
      scanf("%s", key);
      if (strcmp(endstr, key) == 0)
         break;
      for (i = 0; i < com.ns; i++)
         chosen[i] = '\0';

      k = (int)strlen(key);
      for (i = 0; i < com.ns; i++) {
         if ((p = strstr(com.spname[i], key))
            && !isdigit(p[k]))
            chosen[i] = 1;
      }

      /*
      for(i=0; i<com.ns; i++)
         if(strstr(com.spname[i], key)) chosen[i] = 1;
      */

      /* look for MRCA, going through two rounds, assuming unrooted tree */
      for (imrca = 0; imrca < 1 + unrooted; imrca++) {
         if (imrca)
            for (i = 0; i < com.ns; i++) chosen[i] = 1 - chosen[i];

         for (i = 0, sizeclade = 0; i < com.ns; i++)
            if (chosen[i]) {
               sizeclade++;
               lasts = i;
            }

         if (sizeclade <= 1 || sizeclade >= com.ns - 1) {
            puts("unable to form a clade.  <2 seqs.");
            break;
         }
         for (i = 0; i < com.ns - 1; i++) for (j = 0; j < com.ns / SI + 1; j++)
            anc[i][j] = 0;
         for (is = 0; is < com.ns; is++) {
            if (chosen[is] == 0) continue;
            loc = is / SI;  bitmask = 1 << (is % SI);
            for (j = nodes[is].father; j != -1; j = nodes[j].father) {
               anc[j - com.ns][loc] |= bitmask;
               if (is == lasts) {
                  for (i = 0, k = 0; i < com.ns; i++)
                     if (anc[j - com.ns][i / SI] & (1 << (i % SI)))
                        k++;
                  if (k == sizeclade) {
                     mrca = j;  break;
                  }
               }
            }
         }
         if (imrca == 0 && mrca != tree.root) /* 1st round is enough */
            break;
      }

      if (sizeclade <= 1 || sizeclade >= com.ns - 1 || mrca == tree.root) {
         printf("Unable to label.  Ignored.");
         continue;
      }

      if (debug)
         for (is = 0; is < com.ns - 1; is++) {
            printf("\nnode %4d: ", is + com.ns);
            for (j = 0; j < com.ns; j++) {
               loc = j / SI;  bitmask = 1 << (j % SI);
               printf(" %d", (anc[is][loc] & bitmask) != 0);
            }
         }

      printf("\nClade #%d (%s): %d seqs selected, MRCA is %d\n", iclade + 1, key, sizeclade, mrca + 1);
      for (is = 0, paraphyl = 0; is < com.ns; is++) {
         if (chosen[is] == 0)
            for (j = nodes[is].father; j != -1; j = nodes[j].father)
               if (j == mrca) { paraphyl++;  break; }
      }
      if (paraphyl)
         printf("\nThis clade is paraphyletic, & includes %d other sequences\n", paraphyl);

      nodes[mrca].label = iclade + 1;
      if (debug) OutTreeN(F0, 1, haslength | PrLabel);
   }

   for (i = 0; i < com.ns - 1; i++)  free(anc[i]);
   OutTreeN(fout, 1, haslength | PrLabel);
   fprintf(fout, "\n");
   printf("Printed final tree with labels in evolver.out\n");
   exit(0);
}

void TreeDistanceDistribution(FILE* fout)
{
   /* This calculates figure 3.7 of Yang (2006).
      This reads the file of all trees (such as 7s.all.trees), and calculates the
      distribution of partition distance in all pairwise comparisons.
   */
   int i, j, ntree, k, * nib, nsame, IBsame[NS], lpart = 0;
   char treef[64] = "5s.all.trees", * partition;
   FILE* ftree;
   double mPD[NS], PD1[NS];  /* distribution of partition distances */

   puts("Tree file name?");
   scanf("%s", treef);

   ftree = zopen(treef, "r");
   fscanf(ftree, "%d%d", &com.ns, &ntree);
   printf("%2d sequences %2d trees.\n", com.ns, ntree);
   i = (com.ns * 2 - 1) * sizeof(struct TREEN);
   if ((nodes = (struct TREEN*)malloc(i)) == NULL) zerror("oom");

   lpart = (com.ns - 1) * com.ns * sizeof(char);
   i = ntree * lpart;
   printf("\n%d bytes of space requested.\n", i);
   partition = (char*)malloc(i);
   nib = (int*)malloc(ntree * sizeof(int));
   if (partition == NULL || nib == NULL) zerror("out of memory");

   puts("\ntree #: mean prop of tree pairs with 0 1 2 ... shared bipartitions\n");
   fputs("\ntree #: prop of tree pairs with 0 1 2 ... shared bipartitions\n", fout);
   for (i = 0; i < ntree; i++) {
      ReadTreeN(ftree, &j, 0, 1);
      nib[i] = tree.nbranch - com.ns;
      Tree2Partition(partition + i * lpart);
   }
   for (k = 0; k < com.ns - 3; k++) mPD[k] = 0;
   for (i = 0; i < ntree; i++) {
      for (k = 0; k < com.ns - 3; k++) PD1[k] = 0;
      for (j = 0; j < ntree; j++) {
         if (j == i) continue;
         nsame = NSameBranch(partition + i * lpart, partition + j * lpart, nib[i], nib[j], IBsame);
         PD1[nsame] ++;
      }
      for (k = 0; k < com.ns - 3; k++) PD1[k] /= (ntree - 1.);
      for (k = 0; k < com.ns - 3; k++) mPD[k] = (mPD[k] * i + PD1[k]) / (i + 1.);
      printf("%8d (%5.1f%%):", i + 1, (i + 1.) / ntree * 100);
      for (k = 0; k < com.ns - 3; k++) printf(" %7.4f", mPD[k]);
      fprintf(fout, "%8d:", i + 1);  for (k = 0; k < com.ns - 3; k++) fprintf(fout, " %7.4f", PD1[k]);
      printf("%s", (com.ns < 8 || (i + 1) % 100 == 0 ? "\n" : "\r"));
      fprintf(fout, "\n");
   }
   free(partition); free(nodes); free(nib); fclose(ftree);
   exit(0);
}


void TreeDistances(FILE* fout)
{
   /* I think this is broken after i changed the routine Tree2Partition().
   */
   int i, j, ntree, k, * nib, parti2B[NS] = { 0 }, nsame, IBsame[NS], nIBsame[NS], lpart = 0;
   char treef[64] = "5s.all.trees", * partition;
   FILE* ftree;
   double psame, mp, vp;

   /*
   TreeDistanceDistribution(fout);
   */

   puts("\nNumber of identical bi-partitions between trees.\nTree file name?");
   scanf("%s", treef);

   ftree = zopen(treef, "r");
   fscanf(ftree, "%d%d", &com.ns, &ntree);
   printf("%2d sequences %2d trees.\n", com.ns, ntree);
   i = (com.ns * 2 - 1) * sizeof(struct TREEN);
   if ((nodes = (struct TREEN*)malloc(i)) == NULL) zerror("oom");

   if (ntree < 2) zerror("ntree");
   printf("\n%d species, %d trees\n", com.ns, ntree);
   puts("\n\t1: first vs. rest?\n\t2: all pairwise comparisons?\n");
   k = 2;
   scanf("%d", &k);

   lpart = (com.ns - 1) * com.ns * sizeof(char);
   i = (k == 1 ? 2 : ntree) * lpart;
   printf("\n%d bytes of space requested.\n", i);
   partition = (char*)malloc(i);
   nib = (int*)malloc(ntree * sizeof(int));
   if (partition == NULL || nib == NULL) zerror("out of memory");

   if (k == 2) {    /* pairwise comparisons */
      fputs("Number of identical bi-partitions in pairwise comparisons\n", fout);
      for (i = 0; i < ntree; i++) {
         ReadTreeN(ftree, &j, 0, 1);
         nib[i] = tree.nbranch - com.ns;
         Tree2Partition(partition + i * lpart);
      }
      for (i = 0; i < ntree; i++) {
         printf("%2d (%2d):", i + 1, nib[i]);
         fprintf(fout, "%2d (%2d):", i + 1, nib[i]);
         for (j = 0; j < i; j++) {
            nsame = NSameBranch(partition + i * lpart, partition + j * lpart, nib[i], nib[j], IBsame);
            printf(" %2d", nsame);
            fprintf(fout, " %2d", nsame);
         }
         printf("\n");
         fprintf(fout, "\n");
      }
   }
   else {  /* first vs. others */
      ReadTreeN(ftree, &j, 0, 1);
      nib[0] = tree.nbranch - com.ns;
      if (nib[0] == 0) zerror("1st tree is a star tree..");
      Tree2Partition(partition);
      fputs("Comparing the first tree with the others\nFirst tree:\n", fout);
      OutTreeN(fout, 0, 0);  fprintf(fout, "\n");
      OutTreeB(fout);        fprintf(fout, "\n");
      fputs("\nInternal branches in the first tree:\n", fout);
      for (i = 0; i < nib[0]; i++) {
         k = parti2B[i];
         fprintf(fout, "%3d (%2d..%-2d): ( ",
            i + 1, tree.branches[k][0] + 1, tree.branches[k][1] + 1);
         for (j = 0; j < com.ns; j++)
            if (partition[i * com.ns + j]) fprintf(fout, "%d ", j + 1);
         fprintf(fout, ")\n");
      }
      if (nodes[tree.root].nson <= 2)
         fputs("\nRooted tree, results may not be correct.\n", fout);
      fputs("\nCorrect internal branches compared with the 1st tree:\n", fout);
      for (k = 0; k < nib[0]; k++) nIBsame[k] = 0;
      for (i = 1, mp = vp = 0; i < ntree; i++) {
         ReadTreeN(ftree, &j, 0, 1);
         nib[1] = tree.nbranch - com.ns;
         Tree2Partition(partition + lpart);
         nsame = NSameBranch(partition, partition + lpart, nib[0], nib[1], IBsame);

         psame = nsame / (double)nib[0];
         for (k = 0; k < nib[0]; k++) nIBsame[k] += IBsame[k];
         fprintf(fout, "1 vs. %3d: %4d: ", i + 1, nsame);
         for (k = 0; k < nib[0]; k++)
            if (IBsame[k]) fprintf(fout, " %2d", k + 1);
         printf("1 vs. %5d: %6d/%d  %10.4f\n", i + 1, nsame, nib[0], psame);
         vp += square(psame - mp) * (i - 1.) / i;
         mp = (mp * (i - 1.) + psame) / i;
         fprintf(fout, "\n");
      }
      vp = (ntree <= 2 ? 0 : sqrt(vp / ((ntree - 1 - 1) * (ntree - 1.))));
      fprintf(fout, "\nmean and S.E. of proportion of identical partitions\n");
      fprintf(fout, "between the 1st and all the other %d trees ", ntree - 1);
      fprintf(fout, "(ignore these if not revelant):\n %.4f +- %.4f\n", mp, vp);
      fprintf(fout, "\nNumbers of times, out of %d, ", ntree - 1);
      fprintf(fout, "interior branches of tree 1 are present");
      fputs("\n(This may be bootstrap support for nodes in tree 1)\n", fout);
      for (k = 0; k < nib[0]; k++) {
         i = tree.branches[parti2B[k]][0] + 1;  j = tree.branches[parti2B[k]][1] + 1;
         fprintf(fout, "%3d (%2d..%-2d): %6d (%5.1f%%)\n",
            k + 1, i, j, nIBsame[k], nIBsame[k] * 100. / (ntree - 1.));
      }
   }
   free(partition);  free(nodes); free(nib);  fclose(ftree);
}


int EigenQbase(double rates[], double pi[], double Root[], double U[], double V[], double Q[])
{
   /* Construct the rate matrix Q[] for nucleotide model REV.
   */
   int i, j, k;
   double mr, space[4];

   zero(Q, 16);
   for (i = 0, k = 0; i < 3; i++)
      for (j = i + 1; j < 4; j++)
         if (i * 4 + j != 11)
            Q[i * 4 + j] = Q[j * 4 + i] = rates[k++];
   Q[3 * 4 + 2] = Q[2 * 4 + 3] = 1;
   for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
         Q[i * 4 + j] *= pi[j];
   for (i = 0, mr = 0; i < 4; i++) {
      Q[i * 4 + i] = 0;
      Q[i * 4 + i] = -sum(Q + i * 4, 4);
      mr -= pi[i] * Q[i * 4 + i];
   }
   abyx(1 / mr, Q, 16);

   eigenQREV(Q, com.pi, 4, Root, U, V, space);
   return (0);
}


static double freqK_NS = -1;

int EigenQcodon(int getstats, double kappa, double omega, double pi[],
   double Root[], double U[], double V[], double Q[])
{
   /* Construct the rate matrix Q[].
      64 codons are used, and stop codons have 0 freqs.
   */
   int n = com.ncode, i, j, k, c[2], ndiff, pos = 0, from[3], to[3];
   double mr, space[64];

   for (i = 0; i < n * n; i++) Q[i] = 0;
   for (i = 0; i < n; i++)
      for (j = 0; j < i; j++) {
         from[0] = i / 16; from[1] = (i / 4) % 4; from[2] = i % 4;
         to[0] = j / 16;   to[1] = (j / 4) % 4;   to[2] = j % 4;
         c[0] = GeneticCode[com.icode][i];   c[1] = GeneticCode[com.icode][j];
         if (c[0] == -1 || c[1] == -1)  continue;
         for (k = 0, ndiff = 0; k < 3; k++)  if (from[k] != to[k]) { ndiff++; pos = k; }
         if (ndiff != 1)  continue;
         Q[i * n + j] = 1;
         if ((from[pos] + to[pos] - 1) * (from[pos] + to[pos] - 5) == 0)  Q[i * n + j] *= kappa;
         if (c[0] != c[1])  Q[i * n + j] *= omega;
         Q[j * n + i] = Q[i * n + j];
      }
   for (i = 0; i < n; i++) for (j = 0; j < n; j++)
      Q[i * n + j] *= com.pi[j];
   for (i = 0, mr = 0; i < n; i++) {
      Q[i * n + i] = -sum(Q + i * n, n);
      mr -= pi[i] * Q[i * n + i];
   }

   if (getstats)
      Qfactor += freqK_NS * mr;
   else {
      if (com.ncatG == 0)
         for (i = 0; i < n * n; i++)
            Q[i] *= 1 / mr;
      else
         for (i = 0; i < n * n; i++)
            Q[i] *= Qfactor;  /* NSsites models */
      eigenQREV(Q, com.pi, n, Root, U, V, space);
   }
   return (0);
}



int EigenQaa(double pi[], double Root[], double U[], double V[], double Q[])
{
   /* Construct the rate matrix Q[]
   */
   int n = 20, i, j;
   double mr, space[20];

   zero(Q, n * n);
   switch (com.model) {
   case (Poisson): case (EqualInput):
      fillxc(Q, 1., n * n);  break;
   case (Empirical): case (Empirical_F):
      for (i = 0; i < n; i++) for (j = 0; j < i; j++)
         Q[i * n + j] = Q[j * n + i] = com.daa[i * n + j] / 100;
      break;
   }
   for (i = 0; i < n; i++) for (j = 0; j < n; j++)
      Q[i * n + j] *= com.pi[j];
   for (i = 0, mr = 0; i < n; i++) {
      Q[i * n + i] = 0; Q[i * n + i] = -sum(Q + i * n, n);  mr -= com.pi[i] * Q[i * n + i];
   }

   eigenQREV(Q, com.pi, n, Root, U, V, space);
   for (i = 0; i < n; i++)  Root[i] /= mr;

   return (0);
}


int GetDaa(FILE* fout, double daa[])
{
   /* Get the amino acid substitution rate matrix (grantham, dayhoff, jones, etc).
   */
   FILE* fdaa;
   char aa3[4] = "";
   int i, j, n = 20;

   fdaa = zopen(com.daafile, "r");
   printf("\nReading rate matrix from %s\n", com.daafile);

   for (i = 0; i < n; i++)  for (j = 0, daa[i * n + i] = 0; j < i; j++) {
      fscanf(fdaa, "%lf", &daa[i * n + j]);
      daa[j * n + i] = daa[i * n + j];
   }
   if (com.model == Empirical) {
      for (i = 0; i < n; i++)
         if (fscanf(fdaa, "%lf", &com.pi[i]) != 1)
            zerror("err aaRatefile");
      if (fabs(1 - sum(com.pi, 20)) > 1e-4) zerror("\nSum of aa freq. != 1\n");
   }
   fclose(fdaa);

   if (fout) {
      fprintf(fout, "\n%s\n", com.daafile);
      for (i = 0; i < n; i++) {
         fprintf(fout, "\n%4s", getAAstr(aa3, i));
         for (j = 0; j < i; j++)  fprintf(fout, "%5.0f", daa[i * n + j]);
      }
      fprintf(fout, "\n");
   }

   return (0);
}




void MakeSeq(char* z, int ls)
{
   /* generate a random sequence of nucleotides, codons, or amino acids by
      sampling com.pi[], or read the ancestral sequence from the file RootSeq.txt
      if the file exists.
   */
   int i, j, h, n = com.ncode, ch, n31 = (com.seqtype == 1 ? 3 : 1), lst;
   double p[64], r, small = 5e-5;
   char rootseqf[] = "RootSeq.txt", codon[4] = "   ";
   FILE* frootseq = (FILE*)fopen(rootseqf, "r");
   static int times = 0;

   if (frootseq) {
      if (times++ == 0) printf("Reading sequence at the root from file.\n\n");
      if (com.siterates && com.ncatG > 1)
         zerror("sequence for root doesn't work for site-class models");

      for (lst = 0; ; ) {
         for (i = 0; i < n31; i++) {
            while ((ch = fgetc(frootseq)) != EOF && !isalpha(ch));
            if (ch == EOF) zerror("EOF when reading root sequence.");
            if (isalpha(ch))
               codon[i] = (char)(ch = CodeChara((char)ch, com.seqtype));
         }
         if (com.seqtype == 1) ch = codon[0] * 16 + codon[1] * 4 + codon[2];
         if (ch<0 || ch>n - 1)
            printf("error when reading site %d\n", lst + 1);
         if (com.seqtype == 1 && com.pi[ch] == 0)
            printf("you seem to have a stop codon in the root sequence\n");

         z[lst++] = (char)ch;
         if (lst == com.ls) break;
      }
      fclose(frootseq);
   }
   else {
      for (j = 0; j < n; j++)  p[j] = com.pi[j];
      for (j = 1; j < n; j++)  p[j] += p[j - 1];
      if (fabs(p[n - 1] - 1) > small) {
         printf("\nsum pi = %.6f != 1!\n", p[n - 1]);
         exit(-1);
      }
      for (h = 0; h < com.ls; h++) {
         for (j = 0, r = rndu(); j < n - 1; j++)
            if (r < p[j]) break;
         z[h] = (char)j;
      }
   }
}



void Evolve(int inode)
{
   /* evolve sequence com.z[tree.root] along the tree to generate com.z[],
      using nodes[].branch, nodes[].omega, & com.model
      Needs com.z[0,1,...,nnode-1], while com.z[0] -- com.z[ns-1] constitute
      the data.
      For codon sequences, com.siterates[] has w's for NSsites and NSbranchsite models.
   */
   int is, h, i, j, ison, from, n = com.ncode, longseq = 100000;
   double t, rw;

   for (is = 0; is < nodes[inode].nson; is++) {
      ison = nodes[inode].sons[is];
      memcpy(com.z[ison], com.z[inode], com.ls * sizeof(char));
      t = nodes[ison].branch;

      if (com.seqtype == 1 && com.model && com.NSsites) { /* branch-site models */
         Qfactor = com.QfactorBS[ison];
         for (h = 0; h < com.ls; h++)
            com.siterates[h] = com.omegaBS[ison * com.ncatG + com.siteID[h]];
      }

      for (h = 0; h < com.ls; h++) {
         /* decide whether to recalcualte PMat[]. */
         if (h == 0 || (com.siterates && com.siterates[h] != com.siterates[h - 1])) {
            rw = (com.siterates ? com.siterates[h] : 1);

            if (com.seqtype == BASEseq) {
               if (com.model <= TN93)
                  PMatTN93(PMat, t * Qfactor * rw * Qrates[0], t * Qfactor * rw * Qrates[1], t * Qfactor * rw, com.pi);
               else if (com.model == REV)
                  PMatUVRoot(PMat, t * rw, com.ncode, U, V, Root);
            }
            else if (com.seqtype == CODONseq) { /* Watch out for NSsites model */
               if (com.model || com.NSsites) { /* no need to update UVRoot if M0 */
                  if (com.model && com.NSsites == 0) /* branch */
                     rw = nodes[ison].omega;  /* should be equal to com.rK[nodes[].label] */

                  EigenQcodon(0, com.kappa, rw, com.pi, Root, U, V, PMat);
               }
               PMatUVRoot(PMat, t, com.ncode, U, V, Root);
            }
            else if (com.seqtype == AAseq) {
               PMatUVRoot(PMat, t * rw, com.ncode, U, V, Root);
            }

            for (i = 0; i < n; i++)
               for (j = 1; j < n; j++)
                  PMat[i * n + j] += PMat[i * n + j - 1];
         }
         for (j = 0, from = com.z[ison][h], rw = rndu(); j < n - 1; j++)
            if (rw < PMat[from * n + j]) break;
         com.z[ison][h] = j;
      }

      if (com.ls > longseq) printf("\r   nodes %2d -> %2d, evolving . .   ", inode + 1, ison + 1);

      if (nodes[ison].nson) Evolve(ison);
   }  /* for (is) */

   if (inode == tree.root && com.ls > longseq)  printf("\r%s", strc(50, ' '));
}



void Simulate(char* ctlf)
{
   /* simulate nr data sets of nucleotide, codon, or AA sequences.
      ls: number of nucleotides, codons, or AAs in each sequence.
      All 64 codons are used for codon sequences.
      When com.alpha or com.ncatG>1, sites are randomized after sequences are
      generated.
      space[com.ls] is used to hold site marks.
      format:  0: paml sites; 1: paml patterns; 2: paup nex; 3: paup JC69 format
   */
   char* ancf = "ancestral.txt", * siteIDf = "siterates.txt";
   char* aaf = "mc.aa.txt", * c12f = "mc.codon12.txt";
   FILE* fin, * fseq, * fanc = NULL, * faa = NULL, * fc12 = NULL, * fsiteID = NULL;
   char* paupstart = "paupstart", * paupblock = "paupblock", * paupend = "paupend";
   char line[32000];
   int lline = 32000, i, j, k, ir, n, nr, fixtree = 1, sspace = 10000, rooted = 1;
   int h = 0, format = 0, nrate = 1, counts[NCATG];
   int* siteorder = NULL;
   char* tmpseq = NULL, * pc;
   double birth = 0, death = 0, sample = 1, mut = 1, tlength, * space, * blengthBS;
   double T, C, A, G, Y, R, Falias[NCATG];
   int    Lalias[NCATG];

   noisy = 1;
   printf("\nReading options from data file %s\n", ctlf);
   com.ncode = n = (com.seqtype == 0 ? 4 : (com.seqtype == 1 ? 64 : 20));
   fin = (FILE*)zopen(ctlf, "r");
   fscanf(fin, "%d", &format);
   fgets(line, lline, fin);
   printf("\nSimulated data will go into %s.\n", seqf[format]);
   if (format == 2) printf("%s, %s, & %s will be appended if existent.\n", paupstart, paupblock, paupend);
   if (com.seqtype == 1) {
      printf("translated aa sequences will go into %s.\n", aaf);
      faa = (FILE*)zopen(aaf, "w");
      fc12 = (FILE*)zopen(c12f, "w");
   }
   fscanf(fin, "%d", &i);
   fgets(line, lline, fin);
   SetSeed(i, 1);
   fscanf(fin, "%d%d%d", &com.ns, &com.ls, &nr);
   fgets(line, lline, fin);
   i = (com.ns * 2 - 1) * sizeof(struct TREEN);
   if ((nodes = (struct TREEN*)malloc(i)) == NULL) zerror("oom");

   if (com.ns > NS) zerror("too many seqs?");
   printf("\n%d seqs, %d sites, %d replicate(s)\n", com.ns, com.ls, nr);
   k = (com.ns * com.ls * (com.seqtype == CODONseq ? 4 : 1) * nr) / 1000 + 1;
   printf("Seq file will be about %dK bytes.\n", k);
   for (i = 0; i < com.ns; i++)          /* default spname */
      sprintf(com.spname[i], "S%d", i + 1);

   if (fixtree) {
      fscanf(fin, "%lf", &tlength);   fgets(line, lline, fin);
      if (ReadTreeN(fin, &i, 1, 1))  /* might overwrite spname */
         zerror("err tree..");

      if (i == 0) zerror("use : to specify branch lengths in tree");
      for (i = 0, T = 0; i < tree.nnode; i++)
         if (i != tree.root) T += nodes[i].branch;
      if (tlength > 0) {
         for (i = 0; i < tree.nnode; i++)
            if (i != tree.root) nodes[i].branch *= tlength / T;
      }
      printf("tree length = %.3f\n", (tlength > 0 ? tlength : T));
      if (com.ns < 100) {
         printf("\nModel tree & branch lengths:\n");
         OutTreeN(F0, 1, 1);  printf("\n");
         OutTreeN(F0, 0, 1);  printf("\n");
      }
      if (com.seqtype == CODONseq && com.model && !com.NSsites) { /* branch model */
         for (i = 0; i < tree.nnode; i++)
            nodes[i].omega = nodes[i].label;
         printf("\n");
         OutTreeN(F0, 1, PrBranch | PrLabel);
         printf("\n");
      }
   }
   else {   /* random trees, broken or need testing? */
      printf("\nbirth rate, death rate, sampling fraction, mutation rate (tree height)?\n");
      fscanf(fin, "%lf%lf%lf%lf", &birth, &death, &sample, &mut);
      fgets(line, lline, fin);
      printf("%9.4f %9.4f %9.4f %9.4f\n", birth, death, sample, mut);
   }

   if (com.seqtype == BASEseq) {
      fscanf(fin, "%d", &com.model);
      fgets(line, lline, fin);
      if (com.model<0 || com.model>REV) zerror("model err");
      if (com.model == T92) zerror("T92: please use HKY85 with T=A and C=G.");

      printf("\nModel: %s\n", basemodels[com.model]);
      if (com.model == REV)        nrate = 5;
      else if (com.model == TN93)  nrate = 2;
      for (i = 0; i < nrate; i++)
         fscanf(fin, "%lf", &Qrates[i]);
      fgets(line, lline, fin);
      if (nrate <= 2)
         for (i = 0; i < nrate; i++) printf("kappa %9.5f\n", Qrates[i]);
      printf("\n");
      if (nrate == 5) {
         printf("a & b & c & d & e: ");
         for (i = 0; i < nrate; i++) printf("%9.5f", Qrates[i]);
         printf("\n");
      }
      if ((com.model == JC69 || com.model == F81) && Qrates[0] != 1)
         zerror("kappa should be 1 for this model");
   }
   else if (com.seqtype == CODONseq) {
      for (i = 0; i < 64; i++)
         getcodon(CODONs[i], i);
      if (com.model == 0 && com.NSsites) {  /* site model */
         fscanf(fin, "%d", &com.ncatG);   fgets(line, lline, fin);
         if (com.ncatG > NCATG) zerror("ncatG>NCATG");
         for (i = 0; i < com.ncatG; i++) fscanf(fin, "%lf", &com.freqK[i]);
         fgets(line, lline, fin);
         for (i = 0; i < com.ncatG; i++) fscanf(fin, "%lf", &com.rK[i]);
         fgets(line, lline, fin);
         printf("\n\ndN/dS (w) for site classes (K=%d)", com.ncatG);
         printf("\nf: ");  for (i = 0; i < com.ncatG; i++) printf("%9.5f", com.freqK[i]);
         printf("\nw: ");  for (i = 0; i < com.ncatG; i++) printf("%9.5f", com.rK[i]);
         printf("\n");
      }
      else if (com.model && com.NSsites) {  /* branchsite model */
         fscanf(fin, "%d", &com.ncatG);
         fgets(line, lline, fin);
         if (com.ncatG > min2(NCATG, 127))
            zerror("ncatG too large");
         for (i = 0; i < com.ncatG; i++)
            fscanf(fin, "%lf", &com.freqK[i]);
         fgets(line, lline, fin);
         printf("\n%d site classes.\nFreqs: ", com.ncatG);
         for (i = 0; i < com.ncatG; i++)
            printf("%9.5f", com.freqK[i]);

         if ((com.omegaBS = (double*)malloc((com.ncatG + 2) * tree.nnode * sizeof(double))) == NULL)
            zerror("oom");
         com.QfactorBS = com.omegaBS + com.ncatG * tree.nnode;
         blengthBS = com.QfactorBS + tree.nnode;

         for (i = 0; i < tree.nnode; i++)
            blengthBS[i] = nodes[i].branch;
         for (k = 0; k < com.ncatG; k++) {
            ReadTreeN(fin, &i, 0, 1);
            if (i) zerror("do not include branch lengths except in the first tree.");
            /* if (!j) zerror("Use # to specify omega's for branches"); */
            for (i = 0; i < tree.nnode; i++)
               com.omegaBS[i * com.ncatG + k] = nodes[i].label;
         }
         for (i = 0; i < tree.nnode; i++) {
            nodes[i].branch = blengthBS[i];
            nodes[i].label = nodes[i].omega = 0;
         }
         for (i = 0; i < tree.nnode; i++) {  /* print out omega as node labels. */
            nodes[i].annotation = pc = (char*)malloc(20 * com.ncatG * sizeof(char));
            sprintf(pc, "'[%.2f", com.omegaBS[i * com.ncatG + 0]);
            for (k = 1, pc += strlen(pc); k < com.ncatG; k++, pc += strlen(pc))
               sprintf(pc, ", %.2f", com.omegaBS[i * com.ncatG + k]);
            sprintf(pc, "]'");
         }
         printf("\n");
         OutTreeN(F0, 1, PrBranch | PrLabel);
         printf("\n");
      }
      else if (com.model == 0) {  /* M0 */
         fscanf(fin, "%lf", &com.omega);
         fgets(line, lline, fin);
         printf("omega = %9.5f\n", com.omega);
         for (i = 0; i < tree.nbranch; i++)
            nodes[tree.branches[i][1]].omega = com.omega;
      }

      fscanf(fin, "%lf", &com.kappa);   fgets(line, lline, fin);
      printf("kappa = %9.5f\n", com.kappa);
   }

   if (com.seqtype == BASEseq || com.seqtype == AAseq) {
      fscanf(fin, "%lf%d", &com.alpha, &com.ncatG);
      fgets(line, lline, fin);
      if (com.alpha)
         printf("Gamma rates, alpha =%.4f (K=%d)\n", com.alpha, com.ncatG);
      else {
         com.ncatG = 0;
         puts("Rates are constant over sites.");
      }
   }
   if (com.alpha || com.ncatG) { /* this is used for codon NSsites as well. */
      k = com.ls;
      if (com.seqtype == 1 && com.model && com.NSsites) k *= tree.nnode;
      if ((com.siterates = (double*)malloc(k * sizeof(double))) == NULL) zerror("oom1");
      if ((siteorder = (int*)malloc(com.ls * sizeof(int))) == NULL) zerror("oom2");
   }

   if (com.seqtype == AAseq) { /* get aa substitution model and rate matrix */
      fscanf(fin, "%d", &com.model);
      printf("\nmodel: %s", aamodels[com.model]);
      if (com.model >= 2) { fscanf(fin, "%s", com.daafile); GetDaa(NULL, com.daa); }
      fgets(line, lline, fin);
   }

   /* get freqs com.pi[] */
   if ((com.seqtype == BASEseq && com.model > K80) ||
      com.seqtype == CODONseq ||
      (com.seqtype == AAseq && (com.model == 1 || com.model == 3)))
      for (k = 0; k < com.ncode; k++) fscanf(fin, "%lf", &com.pi[k]);
   else if (com.model == 0 || (com.seqtype == BASEseq && com.model <= K80))
      fillxc(com.pi, 1. / com.ncode, com.ncode);

   printf("sum pi = 1 = %.6f:", sum(com.pi, com.ncode));
   matout2(F0, com.pi, com.ncode / 4, 4, 9, 6);
   if (com.seqtype == CODONseq) {
      fscanf(fin, "%d", &com.icode);   fgets(line, lline, fin);
      printf("genetic code: %d\n", com.icode);
      for (k = 0; k < com.ncode; k++)
         if (GeneticCode[com.icode][k] == -1 && com.pi[k])
            zerror("stop codons should have frequency 0?");
   }

   if (com.seqtype == BASEseq) {
      if (com.model < REV) {
         T = com.pi[0]; C = com.pi[1]; A = com.pi[2]; G = com.pi[3]; Y = T + C; R = A + G;
         if (com.model == F84) {
            Qrates[1] = 1 + Qrates[0] / R;   /* kappa2 */
            Qrates[0] = 1 + Qrates[0] / Y;   /* kappa1 */
         }
         else if (com.model <= HKY85) Qrates[1] = Qrates[0];
         Qfactor = 1 / (2 * T * C * Qrates[0] + 2 * A * G * Qrates[1] + 2 * Y * R);
      }
      else
         if (com.model == REV) EigenQbase(Qrates, com.pi, Root, U, V, PMat);
   }

   /* get Qfactor for NSsites & NSbranchsite models */
   if (com.seqtype == CODONseq && com.NSsites) {
      if (!com.model) {  /* site models */
         for (k = 0, Qfactor = 0; k < com.ncatG; k++) {
            freqK_NS = com.freqK[k];
            EigenQcodon(1, com.kappa, com.rK[k], com.pi, NULL, NULL, NULL, PMat);
         }
         Qfactor = 1 / Qfactor;
         printf("Qfactor for NSsites model = %9.5f\n", Qfactor);
      }
      else {            /* branch-site models */
         for (i = 0; i < tree.nnode; i++) {
            if (i == tree.root) { com.QfactorBS[i] = -1; continue; }
            for (k = 0, Qfactor = 0; k < com.ncatG; k++) {
               freqK_NS = com.freqK[k];
               EigenQcodon(1, com.kappa, com.omegaBS[i * com.ncatG + k], com.pi, NULL, NULL, NULL, PMat);
            }
            com.QfactorBS[i] = 1 / Qfactor;  Qfactor = 0;
            printf("node %2d: Qfactor = %9.5f\n", i + 1, com.QfactorBS[i]);
         }
      }
   }
   if (com.seqtype == CODONseq && com.ncatG <= 1 && com.model == 0)
      EigenQcodon(0, com.kappa, com.omega, com.pi, Root, U, V, PMat);
   else if (com.seqtype == AAseq)
      EigenQaa(com.pi, Root, U, V, PMat);

   puts("\nAll parameters are read.  Ready to simulate\n");
   sspace = max2(sspace, 8000000);
   sspace = max2(sspace, sizeof(double) * com.ls);
   space = (double*)malloc(sspace);
   if (com.alpha || com.ncatG) tmpseq = (char*)space;
   if (space == NULL) {
      printf("oom for space, %d bytes needed.", sspace);
      exit(-1);
   }

   fseq = zopen(seqf[format], "w");
   if (format == 2 || format == 3) appendfile(fseq, paupstart);

   fanc = (FILE*)zopen(ancf, "w");
   if (fixtree) {
      fputs("\nAncestral sequences generated during simulation ", fanc);
      fprintf(fanc, "(check against %s)\n", seqf[format]);
      OutTreeN(fanc, 0, 0);
      fprintf(fanc, "\n");
      OutTreeB(fanc); fprintf(fanc, "\n");
   }
   if (com.alpha || com.NSsites) {
      fsiteID = (FILE*)zopen(siteIDf, "w");
      if (com.seqtype == 1)
         fprintf(fsiteID, "\nSite class IDs\n");
      else
         fprintf(fsiteID, "\nRates for sites\n");
      if (com.seqtype == CODONseq && com.NSsites) {
         if (!com.model) matout(fsiteID, com.rK, 1, com.ncatG);
         if ((com.siteID = (char*)malloc(com.ls * sizeof(char))) == NULL)
            zerror("oom siteID");
      }
   }

   for (ir = 0; ir < nr; ir++) {
      for (j = 0; j < com.ns * 2 - 1; j++) {
         com.z[j] = (char*)realloc(com.z[j], com.ls * sizeof(char));
         if (com.z[j] == NULL) zerror("memory problem for z[]");
      }
      if (!fixtree) {    /* right now tree is fixed */
         RandomLHistory(rooted, space);
         if (rooted && com.ns < 10) j = GetIofLHistory();
         BranchLengthBD(1, birth, death, sample, mut);
         if (com.ns < 20) {
            printf("\ntree used: ");
            OutTreeN(F0, 1, 1);
            printf("\n");
         }
      }
      MakeSeq(com.z[tree.root], com.ls);

      if (com.alpha)
         Rates4Sites(com.siterates, com.alpha, com.ncatG, com.ls, 0, space);
      else if (com.seqtype == 1 && com.NSsites) { /* for NSsites */
         /* the table for the alias algorithm is the same, but ncatG is small. */
         MultiNomialAliasSetTable(com.ncatG, com.freqK, Falias, Lalias, space);
         MultiNomialAlias(com.ls, com.ncatG, Falias, Lalias, counts);

         for (i = 0, h = 0; i < com.ncatG; i++)
            for (j = 0; j < counts[i]; j++) {
               com.siteID[h] = (char)i;
               com.siterates[h++] = com.rK[i]; /* overwritten later for branchsite */
            }
      }

      Evolve(tree.root);

      /* randomize sites for site-class model */
      if (com.siterates && com.ncatG > 1) {
         if (format == 1 && ir == 0)
            puts("\nrequested site pattern counts as output for site-class model.\n");
         randorder(siteorder, com.ls, (int*)space);
         for (j = 0; j < tree.nnode; j++) {
            memcpy(tmpseq, com.z[j], com.ls * sizeof(char));
            for (h = 0; h < com.ls; h++) com.z[j][h] = tmpseq[siteorder[h]];
         }
         if (com.alpha || com.ncatG > 1) {
            memcpy(space, com.siterates, com.ls * sizeof(double));
            for (h = 0; h < com.ls; h++) com.siterates[h] = space[siteorder[h]];
         }
         if (com.siteID) {
            memcpy((char*)space, com.siteID, com.ls * sizeof(char));
            for (h = 0; h < com.ls; h++) com.siteID[h] = *((char*)space + siteorder[h]);
         }
      }

      /* print sequences*/
      if (format == 1 || format == 3) {
         for (i = 0; i < com.ns; i++) for (h = 0; h < com.ls; h++)    com.z[i][h] ++;  /* code as 1, 2, ... */
         PatternWeightSimple();
         for (i = 0; i < com.ns; i++) for (h = 0; h < com.npatt; h++) com.z[i][h] --;  /* code as 0, 1, ... */
         if (format == 3)
            PatternWeightJC69like();
      }
      if (format == 2 || format == 3) fprintf(fseq, "\n\n[Replicate # %d]\n", ir + 1);
      printSeqs(fseq, com.z, com.spname, com.ns, com.ls, com.npatt, com.fpatt, NULL, NULL, format); /* printsma not usable as it codes into 0,1,...,60. */
      if ((format == 2 || format == 3) && !fixtree) {
         fprintf(fseq, "\nbegin tree;\n   tree true_tree = [&U] ");
         OutTreeN(fseq, 1, 1); fputs(";\n", fseq);
         fprintf(fseq, "end;\n\n");
      }
      if (format == 2 || format == 3) appendfile(fseq, paupblock);
      /* print translated protein sequences, and codonpositions 1 & 2 */
      if (com.seqtype == 1 && faa) {
         fprintf(faa, "%6d %6d\n", com.ns, com.ls);
         for (j = 0; j < com.ns; j++) {
            fprintf(faa, "\n%-10s  ", com.spname[j]);
            for (h = 0; h < com.ls; h++) {
               int iaa = GeneticCode[com.icode][(int)com.z[j][h]];
               fprintf(faa, "%c", AAs[iaa]);
               if ((h + 1) % 10 == 0)  fprintf(faa, " ");
            }
         }
         fprintf(faa, "\n\n");

         fprintf(fc12, "%6d %6d\n", com.ns, com.ls * 2);
         for (j = 0; j < com.ns; j++) {
            fprintf(fc12, "\n%-10s  ", com.spname[j]);
            for (h = 0; h < com.ls; h++) {
               int ic = com.z[j][h];
               fprintf(fc12, "%c%c", CODONs[ic][0], CODONs[ic][1]);
               if ((h + 1) % 5 == 0)  fprintf(fc12, " ");
            }
         }
         fprintf(fc12, "\n\n");
      }

      /* print ancestral seqs, rates for sites. */
      if (format != 1 && format != 3) {  /* don't print ancestors if site patterns are printed. */
         j = (com.seqtype == CODONseq ? 3 * com.ls : com.ls);
         fprintf(fanc, "[replicate %d]\n", ir + 1);

         if (!fixtree) {
            if (format < 2) {
               OutTreeN(fanc, 1, 1);
               fprintf(fanc, "\n\n");
            }
         }
         else {
            fprintf(fanc, "%6d %6d\n", tree.nnode - com.ns, j);
            for (j = com.ns; j < tree.nnode; j++) {
               fprintf(fanc, "node%-26d  ", j + 1);
               print1seq(fanc, com.z[j], com.ls, NULL);
               fprintf(fanc, "\n");
            }
            fprintf(fanc, "\n");

            if (fsiteID) {
               if (com.seqtype == CODONseq && com.NSsites && com.model == 0) { /* site model */
                  k = 0;
                  if (com.rK[com.ncatG - 1] > 1)
                     for (h = 0; h < com.ls; h++)
                        if (com.rK[(int)com.siteID[h]] > 1) k++;
                  fprintf(fsiteID, "\n[replicate %d: %2d]\n", ir + 1, k);
                  if (k)  for (h = 0, k = 0; h < com.ls; h++) {
                     if (com.rK[(int)com.siteID[h]] > 1) {
                        fprintf(fsiteID, "%4d ", h + 1);
                        if (++k % 15 == 0) fprintf(fsiteID, "\n\n");
                     }
                  }
                  fprintf(fsiteID, "\n\n");
               }
               else if (com.seqtype == CODONseq && com.NSsites && com.model) { /* branchsite */
                  fprintf(fsiteID, "\n[replicate %d]\n", ir + 1);
                  for (h = 0; h < com.ls; h++) {
                     fprintf(fsiteID, " %4d ", com.siteID[h] + 1);
                     if (h == com.ls - 1 || (h + 1) % 15 == 0)
                        fprintf(fsiteID, "\n\n");
                  }
               }
               else {       /* gamma rates */
                  fprintf(fsiteID, "\n[replicate %d]\n", ir + 1);
                  for (h = 0; h < com.ls; h++) {
                     fprintf(fsiteID, "%7.4f ", com.siterates[h]);
                     if (h == com.ls - 1 || (h + 1) % 10 == 0)
                        fprintf(fsiteID, "\n\n");
                  }
               }
            }
         }
      }
      printf("\rdid data set %d %s", ir + 1, (com.ls > 100000 || nr < 100 ? "\n" : ""));
   }   /* for (ir) */
   if (format == 2 || format == 3) appendfile(fseq, paupend);

   fclose(fseq);
   if (!fixtree) fclose(fanc);
   if (faa) fclose(faa);
   if (fc12) fclose(fc12);
   if (com.alpha || com.NSsites) fclose(fsiteID);
   for (j = 0; j < com.ns * 2 - 1; j++) free(com.z[j]);
   free(space);
   if (com.model && com.NSsites) /* branch-site model */
      for (i = 0; i < tree.nnode; i++)  free(nodes[i].annotation);
   free(nodes);
   if (com.alpha || com.ncatG) {
      free(com.siterates);
      com.siterates = NULL;
      free(siteorder);
      if (com.siteID) free(com.siteID);
      com.siteID = NULL;
   }
   if (com.seqtype == 1 && com.model && com.NSsites) free(com.omegaBS);
   com.omegaBS = NULL;

   exit(0);
}


int GetSpnamesFromMB(FILE* fmb, char line[], int lline)
{
   /* This reads species names from MrBayes output file fmb, like the following.

         Taxon  1 -> 1_Arabidopsis_thaliana
         Taxon  2 -> 2_Taxus_baccata
   */
   int j, ispecies;
   char* p = NULL, * mbstr1 = "Taxon ", * mbstr2 = "->";

   puts("Reading species names from mb output file.\n");
   rewind(fmb);
   for (ispecies = 0; ; ) {
      if (fgets(line, lline, fmb) == NULL) return(-1);
      if (strstr(line, mbstr1) && strstr(line, mbstr2)) {
         p = strstr(line, mbstr1) + 5;
         sscanf(p, "%d", &ispecies);
         p = strstr(line, mbstr2) + 3;
         if (com.spname[ispecies - 1][0])
            zerror("species name already read?");

         for (j = 0; isgraph(*p) && j < lline; ) com.spname[ispecies - 1][j++] = *p++;
         com.spname[ispecies - 1][j] = 0;

         printf("\tTaxon %2d:  %s\n", ispecies, com.spname[ispecies - 1]);
      }
      else if (ispecies)
         break;
   }
   com.ns = ispecies;
   rewind(fmb);

   return(0);
}

char* GrepLine(FILE* fin, char* query, char* line, int lline)
{
   /* This greps infile to search for query[], and returns NULL or line[].
   */

   rewind(fin);
   for (; ; ) {
      if (fgets(line, lline, fin) == NULL)
         return(NULL);
      if (strstr(line, query))
         return(line);
   }
   return(NULL);
}


void CladeMrBayesProbabilities(char treefile[])
{
   /* This reads a tree from treefile and then scans a set of MrBayes output files
      (mbfiles) to retrieve posterior probabilities for every clade in that tree.
      It first scans the first mb output file to get the species names.

      Sample mb output:
      6 -- ...........................*************   8001 1.000 0.005 (0.000)
      7 -- ....................********************   8001 1.000 0.006 (0.000)

      Note 4 Jan 2014: This uses parti2B[], and is broken after i rewrote Tree2Partition().
   */
   int lline = 100000, i, j, k, nib, inode, parti2B[NS] = { 0 };
   char line[100000], * partition, * p;
   char symbol[2] = ".*", cladestr[NS + 1] = { 0 };
   FILE* ftree, * fmb[20];
   double* Pclade, t;
   int nmbfiles = 2;
   char* mbfiles[] = { "mb-1e-4.out", "mb-1e-1.out" };

   printf("tree file is %s\nmb output files:\n", treefile);
   ftree = zopen(treefile, "r");
   for (k = 0; k < nmbfiles; k++)
      fmb[k] = zopen(mbfiles[k], "r");
   for (k = 0; k < nmbfiles; k++) printf("\t%s\n", mbfiles[k]);

   GetSpnamesFromMB(fmb[0], line, lline);  /* read species names from mb output */

   fscanf(ftree, "%d%d", &i, &k);
   if (i && i != com.ns) zerror("do you mean to specify ns in the tree file?");
   i = (com.ns * 2 - 1) * sizeof(struct TREEN);
   if ((nodes = (struct TREEN*)malloc(i)) == NULL) zerror("oom");
   ReadTreeN(ftree, &i, 0, 1);

   printf("\n");
   OutTreeN(F0, 0, 0);
   printf("\n");
   nib = tree.nbranch - com.ns;
   for (i = 0; i < tree.nnode; i++) {
      nodes[i].annotation = NULL;
      if (i > com.ns) nodes[i].annotation = (char*)malloc(100 * sizeof(char));
   }

   partition = (char*)malloc(nib * com.ns * sizeof(char));
   if (partition == NULL) zerror("oom");
   if ((Pclade = (double*)malloc(nib * nmbfiles * sizeof(double))) == NULL)
      zerror("oom");
   for (i = 0; i < nib * nmbfiles; i++) Pclade[i] = 0;

   Tree2Partition(partition);

   for (i = 0; i < nib; i++) {
      inode = tree.branches[parti2B[i]][1];
      if (partition[i * com.ns + 0])
         for (j = 0; j < com.ns; j++) cladestr[j] = symbol[1 - partition[i * com.ns + j]];
      else
         for (j = 0; j < com.ns; j++) cladestr[j] = symbol[(int)partition[i * com.ns + j]];
      printf("#%2d branch %2d node %2d  %s", i + 1, parti2B[i], inode, cladestr);

      for (k = 0; k < nmbfiles; k++) {
         if (GrepLine(fmb[k], cladestr, line, lline)) {
            p = strstr(line, cladestr);
            sscanf(p + com.ns, "%lf%lf", &t, &Pclade[i * nmbfiles + k]);
         }
      }
      for (k = 0; k < nmbfiles; k++) printf("%6.2f", Pclade[i * nmbfiles + k]);
      printf("\n");
      for (k = 0, p = nodes[inode].annotation; k < nmbfiles; k++) {
         sprintf(p, "%3.0f%s", Pclade[i * nmbfiles + k] * 100, (k < nmbfiles - 1 ? "/" : ""));
         p += 4;
      }
   }
   printf("\n");
   OutTreeN(F0, 1, PrLabel);
   printf("\n");

   for (i = 0; i < tree.nnode; i++) free(nodes[i].annotation);
   free(nodes); free(partition);  free(Pclade);
   fclose(ftree);
   for (k = 0; k < nmbfiles; k++) fclose(fmb[k]);
   exit(0);
}
