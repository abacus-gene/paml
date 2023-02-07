/* yn00.c
   Pairwise estimation of dS and dN by the method of Yang & Nielsen
   (2000 Mol. Biol. Evol. 17:32-43)

   Copyright, 1998, Ziheng Yang

                 cc -o yn00 -fast yn00.c tools.o -lm
                 cl -O2 yn00.c tools.o
                 yn00 <SequenceFileName>

  Codon sequences are encoded as 0,1,...,61, as in codeml.c.
*/
#include "paml.h"
#define NS            1000
#define LSPNAME       30
#define NCODE         64
#define NGENE         2000

int GetOptions(char* ctlf);
int EncodeSeqCodon(void);
int Statistics(FILE* fout, double space[]);
int DistanceMatLWL85(FILE* fout);
int DistanceYN00(int is, int js, double* S, double* N, double* dS, double* dN,
   double* SEdS, double* SEdN, double* t, double space[]);
int GetKappa(void);
int GetFreqs(int is1, int is2, double f3x4[], double pi[]);
int CountSites(char z[], double pi[], double* Stot, double* Ntot,
   double fbS[], double fbN[]);
int GetPMatCodon(double P[], double t, double kappa, double omega, double space[]);
int CountDiffs(char z1[], char z2[],
   double* Sdts, double* Sdtv, double* Ndts, double* Ndtv, double PMat[]);
int DistanceF84(double n, double P, double Q, double pi[],
   double* k_HKY, double* t, double* SEt);
double dsdnREV(int is, int js, double space[]);

int ExpPattFreq(double t, double kappa, double omega, double pi[], double space[]);
int ConsistencyMC(void);
int InfiniteData(double t, double kappa, double omega, double f3x4_0[],
   double space[]);
void SimulateData2s64(FILE* fout, double f3x4_0[], double space[]);

struct common_info {
   char *z[NS], *spname[NS], seqf[512], outf[512];
   int ns, ls, npatt, codonf, icode, ncode, getSE, * pose, verbose, seqtype, readpattern;
   int cleandata, fcommon, kcommon, weighting, ndata, print;
   double* fpatt, pi[NCODE], f3x4s[NS][12], kappa, omega;
   int ngene, posG[NGENE + 1], lgene[NGENE], fix_rgene, model;
   double rgene[NGENE], piG[NGENE][NCODE], alpha;
}  com;


int FROM61[64], FROM64[64], FourFold[4][4];
double PMat[NCODE * NCODE];
char* codonfreqs[] = { "Fequal", "F1x4", "F3x4", "Fcodon" };
enum { Fequal, F1x4, F3x4, Fcodon } CodonFreqs;

FILE* frst, * frst1, * frub;
extern char BASEs[], AAs[];
extern int noisy, GeneticCode[][64];
int Nsensecodon;
enum { NODEBUG, KAPPA, SITES, DIFF } DebugFunctions;
int debug = 0;

double omega_NG, dN_NG, dS_NG;  /* what are these for? */


#define YN00
#define REALSEQUENCE
#include "treesub.c"


int main(int argc, char* argv[])
{
   char dsf[512] = "2YN.dS", dnf[512] = "2YN.dN", tf[512] = "2YN.t";
   FILE* fout, * fseq, * fds, * fdn, * ft;
   char ctlf[96] = "yn00.ctl", timestr[64];
   int    is, js, idata, wname = 20, sspace;
   double t = 0.4, dS = 0.1, dN = 0.1, S, N, SEdS, SEdN, f3x4[12], * space = NULL;

   /* ConsistencyMC(); */

   printf("YN00 in %s\n", pamlVerStr);
   starttimer();
   if (argc > 1)  strcpy(ctlf, argv[1]);
   com.seqtype = 1;  com.cleandata = 1;  /* works for clean data only? */
   com.ndata = 1;  com.print = 0;
   noisy = 1; com.icode = 0;  com.fcommon = 0;  com.kcommon = 1;
   GetOptions(ctlf);
   setmark_61_64();
   fout = fopen(com.outf, "w");
   frst = fopen("rst", "w");
   frst1 = fopen("rst1", "w");
   frub = fopen("rub", "w");
   if (fout == NULL || frst == NULL) zerror("outfile creation err.");
   fds = (FILE*)fopen(dsf, "w");
   fdn = (FILE*)fopen(dnf, "w");
   ft = (FILE*)fopen(tf, "w");
   if (fds == NULL || fdn == NULL || ft == NULL) zerror("file open error");

   if ((fseq = fopen(com.seqf, "r")) == NULL) {
      printf("\n\nSequence file %s not found!\n", com.seqf);
      exit(-1);
   }
   for (idata = 0; idata < com.ndata; idata++) {
      if (com.ndata > 1) {
         printf("\nData set %d\n", idata + 1);
         fprintf(fout, "\n\nData set %d\n", idata + 1);
         fprintf(frst, "\t%d", idata + 1);
      }

      ReadSeq((com.verbose ? fout : NULL), fseq, com.cleandata, 0, 0);
      SetMapAmbiguity(com.seqtype, 0);

      sspace = max2(200000, 64 * com.ns * sizeof(double));
      sspace = max2(sspace, 64 * 64 * 5 * sizeof(double));
      if ((space = (double*)realloc(space, sspace)) == NULL) zerror("oom space");

      com.kappa = 4.6;
      com.omega = 1;
      fprintf(fout, "YN00 %15s", com.seqf);
      Statistics(fout, space);

      if (noisy) printf("\n\n(A) Nei-Gojobori (1986) method\n");
      fprintf(fout, "\n\n(A) Nei-Gojobori (1986) method\n");
      DistanceMatNG86(fout, NULL, NULL, NULL, 0);
      fflush(fout);

      if (noisy) printf("\n\n(B) Yang & Nielsen (2000) method\n\n");
      fprintf(fout, "\n\n(B) Yang & Nielsen (2000) method\n\n");
      fprintf(fout, "Yang Z, Nielsen R (2000) Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. Mol. Biol. Evol. 17:32-43\n");
      if (!com.weighting) fputs("\n(equal weighting of pathways)\n", fout);

      if (com.fcommon)  GetFreqs(-1, -1, f3x4, com.pi);
      if (com.kcommon) {
         GetKappa();
         printf("kappa = %.2f\n\n", com.kappa);
         /* puts("kappa?"); scanf("%lf", &com.kappa); */
      }

      fputs("\nseq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE\n\n", fout);
      fprintf(fds, "%6d\n", com.ns);
      fprintf(fdn, "%6d\n", com.ns);
      fprintf(ft, "%6d\n", com.ns);
      for (is = 0; is < com.ns; is++) {
         fprintf(fds, "%-*s ", wname, com.spname[is]);
         fprintf(fdn, "%-*s ", wname, com.spname[is]);
         fprintf(ft, "%-*s ", wname, com.spname[is]);
         for (js = 0; js < is; js++) {
            if (noisy) printf("%3d vs. %3d\n", is + 1, js + 1);
            fprintf(fout, " %3d  %3d ", is + 1, js + 1);

            if (!com.fcommon) GetFreqs(is, js, f3x4, com.pi);
            if (!com.kcommon) GetKappa();
            DistanceYN00(is, js, &S, &N, &dS, &dN, &SEdS, &SEdN, &t, space);

            fprintf(fout, "%7.1f %7.1f %8.4f %7.4f %7.4f %6.4f +- %6.4f %7.4f +- %6.4f\n",
               S, N, t, com.kappa, com.omega, dN, SEdN, dS, SEdS);
            fprintf(frst, " YN: %8.4f%8.4f%8.4f %6.4f +- %6.4f %7.4f +- %6.4f\n",
               t, com.kappa, com.omega, dN, SEdN, dS, SEdS);

            fprintf(fds, " %7.4f", dS); fprintf(fdn, " %7.4f", dN); fprintf(ft, " %7.4f", t);
         }    /* for (js) */
         fprintf(fds, "\n"); fprintf(fdn, "\n"); fprintf(ft, "\n");
         fflush(fds); fflush(fdn); fflush(ft);
      }       /* for (is) */
      fprintf(fds, "\n"); fprintf(fdn, "\n"); fprintf(ft, "\n");

      if (noisy) printf("\n\n(C) LWL85, LPB93 & LWLm methods\n\n");
      fprintf(fout, "\n\n(C) LWL85, LPB93 & LWLm methods\n\n");
      fprintf(fout, "Li W.-H., C.-I. Wu, Luo (1985) A new method for estimating synonymous and nonsynonymous rates of nucleotide substitutions considering the relative likelihood of nucleotide and codon changes. Mol. Biol. Evol. 2: 150-174.\n");
      fprintf(fout, "Li W-H (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. J. Mol. Evol. 36:96-99\n");
      fprintf(fout, "Pamilo P, Bianchi NO (1993) Evolution of the Zfx and Zfy genes - rates and interdependence between the genes. Mol. Biol. Evol. 10:271-281\n");
      fprintf(fout, "Yang Z (2006) Computational Molecular Evolution. Oxford University Press, Oxford. Eqs. 2.12 & 2.13\n");

      DistanceMatLWL85(fout);

      fflush(frst);
      if (noisy) printf("\nTime used: %s\n", printtime(timestr));
   }
   return (0);
}



int GetOptions(char* ctlf)
{
   int i, nopt = 9, lline = 4096;
   char line[4096], * pline, opt[20], comment = '*';
   char* optstr[] = { "seqfile","outfile", "verbose", "noisy", "icode",
        "weighting","commonkappa", "commonf3x4", "ndata" };
   double t;
   FILE* fctl;

   if ((fctl = fopen(ctlf, "r")) == NULL) zerror("\nctl file open error.\n");
   printf("\nReading options from %s..\n", ctlf);
   for (;;) {
      if (fgets(line, lline, fctl) == NULL) break;
      for (i = 0, t = 0, pline = line; i < lline && line[i]; i++)
         if (isalnum(line[i])) { t = 1; break; }
         else if (line[i] == comment) break;
      if (t == 0) continue;
      sscanf(line, "%s%*s%lf", opt, &t);
      if ((pline = strstr(line, "=")) == NULL) zerror("option file.");

      for (i = 0; i < nopt; i++) {
         if (strncmp(opt, optstr[i], 8) == 0) {
            if (noisy > 2)
               printf("\n%3d %15s | %-20s %6.2f", i + 1, optstr[i], opt, t);
            switch (i) {
            case (0): sscanf(pline + 2, "%s", com.seqf);    break;
            case (1): sscanf(pline + 2, "%s", com.outf);    break;
            case (2): com.verbose = (int)t;     break;
            case (3): noisy = (int)t;           break;
            case (4): com.icode = (int)t;       break;
            case (5): com.weighting = (int)t;   break;
            case (6): com.kcommon = (int)t;     break;
            case (7): com.fcommon = (int)t;     break;
            case (8): com.ndata = (int)t;       break;
            }
            break;
         }
      }
      if (i == nopt)
      {
         printf("\noption %s in %s\n", opt, ctlf);  exit(-1);
      }
   }

   for (i = 0, Nsensecodon = 0; i < 64; i++)
      if (GeneticCode[com.icode][i] != -1) Nsensecodon++;
   com.ncode = Nsensecodon;
   fclose(fctl);
   printf("\n");
   return (0);
}

int DistanceYN00(int is, int js, double* S, double* N, double* dS, double* dN,
   double* SEdS, double* SEdN, double* t, double space[])
{
   /* calculates dS, dN, w, t by weighting.
      com.kappa & com.pi[] are calculated beforehand are not updated.
   */
   int j, k, ir, nround = 10, status = 0;
   double fbS[4], fbN[4], fbSt[4], fbNt[4], St, Nt, Sdts, Sdtv, Ndts, Ndtv, kappaS, kappaN;
   double w0 = 0, dS0 = 0, dN0 = 0, accu = 5e-4, minomega = 1e-5, maxomega = 99;

   if (*t == 0) *t = .5;
   if (com.omega <= 0) com.omega = 1;
   for (k = 0; k < 4; k++) fbS[k] = fbN[k] = 0;
   if (debug) printf("\nCountSites\n");
   if (noisy > 3) printf("\n");
   for (k = 0, *S = *N = 0; k < 2; k++) {
      CountSites(com.z[k == 0 ? is : js], com.pi, &St, &Nt, fbSt, fbNt);
      *S += St / 2;
      *N += Nt / 2;
      for (j = 0; j < 4; j++) {
         fbS[j] += fbSt[j] / 2;
         fbN[j] += fbNt[j] / 2;
      }
      if (noisy > 3) printf("Seq. %d: S = %9.3f N=%9.3f\n", k + 1, St, Nt);
   }
   if (noisy > 3) {
      printf("Ave.  : S = %9.3f N=%9.3f\n\n", *S, *N);
      printf("Base freqs at syn & nonsyn sites\n%10s%10s%10s%10s\n", "T", "C", "A", "G");
      for (k = 0; k < 4; k++) printf(" %9.6f", fbS[k]);
      printf("\n");
      for (k = 0; k < 4; k++) printf(" %9.6f", fbN[k]);
      printf("\n");
   }
   if (noisy > 3)
      printf(" #    Sdts   Sdtv   Ndts   Ndtv |       t   kappa       w      dN      dS |   kappaS  kappaN\n");

   /* initial values? */
   if (com.weighting) {
      if (*t < 0.001 || *t>5) *t = 0.5;
      if (com.omega < 0.01 || com.omega>5) com.omega = .5;
   }
   for (ir = 0; ir < (com.weighting ? nround : 1); ir++) {   /* weighting or iteration */
      if (com.weighting)
         GetPMatCodon(PMat, *t, com.kappa, com.omega, space);
      else
         for (j = 0; j < com.ncode * com.ncode; j++)
            PMat[j] = 1;

      CountDiffs(com.z[is], com.z[js], &Sdts, &Sdtv, &Ndts, &Ndtv, PMat);

      if (DistanceF84(*S, Sdts / *S, Sdtv / *S, fbS, &kappaS, dS, SEdS)) status = -1;
      if (DistanceF84(*N, Ndts / *N, Ndtv / *N, fbN, &kappaN, dN, SEdN)) status = -1;

      if (*dS < 1e-9) {
         status = -1;
         com.omega = maxomega;
      }
      else
         com.omega = max2(minomega, *dN / *dS);
      *t = *dS * 3 * *S / (*S + *N) + *dN * 3 * *N / (*S + *N);
      if (noisy > 3) {
         printf("%2d %7.2f%7.2f%7.2f%7.2f |", ir + 1, Sdts, Sdtv, Ndts, Ndtv);
         printf("%8.4f%8.4f%8.4f%8.4f%8.4f", *t, com.kappa, com.omega, *dN, *dS);
         printf(" | %8.4f%8.4f\n", kappaS, kappaN);
      }
      if (fabs(*dS - dS0) < accu && fabs(*dN - dN0) < accu && fabs(com.omega - w0) < accu)
         break;
      dS0 = *dS;  dN0 = *dN;  w0 = com.omega;
   } /* for (ir) */
   if (ir == nround) status = -2;
   /* if(status) printf("\n\tstatus: %d\n", status); */
   return(status);
}



int Statistics(FILE* fout, double space[])
{
   /* This calculates base frequencies, using npatt & fpatt[]
   */
   int h, is, j, c[3], wname = 20;
   double f3x4tot[12], * fb3tot = com.pi, * fb3s = space;

   if (fout) {
      fprintf(fout, "\n\nns =%4d\tls =%4d", com.ns, com.ls);
      fprintf(fout, "\n\nCodon position x base (3x4) table for each sequence.");
   }
   zero(f3x4tot, 12);  zero(fb3s, 64 * com.ns);
   for (is = 0; is < com.ns; is++) 
      zero(com.f3x4s[is], 12);
   for (is = 0; is < com.ns; is++) {
      for (h = 0; h < com.npatt; h++) {
         j = FROM61[(int)com.z[is][h]];
         c[0] = j / 16; c[1] = (j % 16) / 4; c[2] = j % 4;
         fb3s[is * 64 + j] += com.fpatt[h];
         for (j = 0; j < 3; j++)
            com.f3x4s[is][j * 4 + c[j]] += com.fpatt[h] / com.ls;
      }
      for (j = 0; j < 12; j++) f3x4tot[j] += com.f3x4s[is][j] / com.ns;
      if (fout) {
         fprintf(fout, "\n\n%-*s", wname, com.spname[is]);
         for (j = 0; j < 3; j++) {
            fprintf(fout, "\nposition %2d:", j + 1);
            for (h = 0; h < 4; h++)
               fprintf(fout, "%5c:%7.5f", BASEs[h], com.f3x4s[is][j * 4 + h]);
         }
      }
   }
   if (fout) {
      fprintf(fout, "\n\nAverage");
      for (j = 0; j < 3; j++) {
         fprintf(fout, "\nposition %2d:", j + 1);
         for (h = 0; h < 4; h++)
            fprintf(fout, "%5c:%7.5f", BASEs[h], f3x4tot[j * 4 + h]);
      }
      for (is = 0, zero(fb3tot, 64); is < com.ns; is++)
         for (j = 0; j < 64; j++) fb3tot[j] += fb3s[is * 64 + j];
      fprintf(fout, "\n\nCodon usage for each species\n");
      printcums(fout, com.ns, fb3s, com.icode);
      fprintf(fout, "\nSums\n");
      printcums(fout, 1, fb3tot, com.icode);
   }

   return(0);
}

int GetFreqs(int is1, int is2, double f3x4[], double pi[])
{
   /* uses com.fcommon and com.f3x4s to calculate f3x4[] and pi[].
      Codon frequencies pi[] are calculated from the f3x4 table.
      The calculation is duplicated when com.fcommon=1.
   */
   int n = com.ncode, j, k, ic, b[3];

   if (com.fcommon)
      for (j = 0, zero(f3x4, 12); j < com.ns; j++)
         for (k = 0; k < 12; k++) f3x4[k] += com.f3x4s[j][k] / com.ns;
   else
      for (k = 0; k < 12; k++)
         f3x4[k] = (com.f3x4s[is1][k] + com.f3x4s[is2][k]) / 2;

   if (noisy >= 9)
      matout(F0, f3x4, 3, 4);
   for (j = 0; j < n; j++) {
      ic = FROM61[j]; b[0] = ic / 16; b[1] = (ic % 16) / 4; b[2] = ic % 4;
      pi[j] = f3x4[b[0]] * f3x4[4 + b[1]] * f3x4[8 + b[2]];
   }
   abyx(1 / sum(pi, n), pi, n);

   return (0);
}


int DistanceMatLWL85(FILE* fout)
{
   /* This implements 3 methods: LWL85 (Li, Wu & Luo 1985), LPB (Li 1993,
      Pamilo & Bianchi 1993), and LWL85m (equation 12 in book; check other refs).
      alpha is not used.
   */
   int i, j, k, h;
   char* codon1, * codon2;
   double L[3], sdiff[3], vdiff[3], Lt[3], sdifft[3], vdifft[3], A[3], B[3];
   double P[3], Q[3], a, b, dS, dN, pS2, S, N, Sd, Nd;

   for (i = 0; i < com.ns; i++) {
      for (j = 0; j < i; j++) {  /* pair i and j */
         for (k = 0; k < 3; k++) L[k] = sdiff[k] = vdiff[k] = 0;

         for (h = 0; h < com.npatt; h++) {
            codon1 = CODONs[(int)com.z[i][h]];
            codon2 = CODONs[(int)com.z[j][h]];
            difcodonLWL85(codon1, codon2, Lt, sdifft, vdifft, 0, com.icode);
            for (k = 0; k < 3; k++) {
               L[k] += Lt[k] * com.fpatt[h];
               sdiff[k] += sdifft[k] * com.fpatt[h];
               vdiff[k] += vdifft[k] * com.fpatt[h];
            }
         }

         for (k = 0; k < 3; k++) {
            P[k] = sdiff[k] / L[k];
            Q[k] = vdiff[k] / L[k];
            a = 1 - 2 * P[k] - Q[k];
            b = 1 - 2 * Q[k];
            A[k] = -log(a) / 2 + log(b) / 4;
            B[k] = -log(b) / 2;
         }
         if (fout) {
            fprintf(fout, "\n%d (%s) vs. %d (%s)\n\n", i + 1, com.spname[i], j + 1, com.spname[j]);
            fprintf(fout, "L(i):  %9.1f %9.1f %9.1f  sum=%9.1f\n", L[0], L[1], L[2], L[0] + L[1] + L[2]);
            fprintf(fout, "Ns(i): %9.4f %9.4f %9.4f  sum=%9.4f\n", sdiff[0], sdiff[1], sdiff[2], sdiff[0] + sdiff[1] + sdiff[2]);
            fprintf(fout, "Nv(i): %9.4f %9.4f %9.4f  sum=%9.4f\n", vdiff[0], vdiff[1], vdiff[2], vdiff[0] + vdiff[1] + vdiff[2]);
            fprintf(fout, "A(i):  %9.4f %9.4f %9.4f\n", A[0], A[1], A[2]);
            fprintf(fout, "B(i):  %9.4f %9.4f %9.4f\n", B[0], B[1], B[2]);

            Sd = L[1] * A[1] + L[2] * (A[2] + B[2]);
            Nd = L[1] * B[1] + L[0] * (A[0] + B[0]);
            pS2 = 1 / 3.;
            S = L[1] * pS2 + L[2];
            N = L[1] * (1 - pS2) + L[0];
            dS = Sd / S;
            dN = Nd / N;
            fprintf(fout, "LWL85:  dS = %7.4f dN = %7.4f w =%7.4f S =%7.1f N =%7.1f\n", dS, dN, dN / dS, S, N);
            pS2 = A[2] / (A[2] + B[2]);
            S = L[1] * pS2 + L[2];
            N = L[1] * (1 - pS2) + L[0];
            dS = Sd / S;
            dN = Nd / N;
            fprintf(fout, "LWL85m: dS = %7.4f dN = %7.4f w =%7.4f S =%7.1f N =%7.1f (rho = %.3f)\n", dS, dN, dN / dS, S, N, pS2);

            dS = (L[1] * A[1] + L[2] * A[2]) / (L[1] + L[2]) + B[2];
            dN = (L[0] * B[0] + L[1] * B[1]) / (L[0] + L[1]) + A[0];
            fprintf(fout, "LPB93:  dS = %7.4f dN = %7.4f w =%7.4f\n", dS, dN, dN / dS);
         }
      }
      if (noisy)  printf(" %3d", i + 1);
   }
   if (noisy)  printf("\n");
   if (fout)   fprintf(fout, "\n");
   return (0);
}



int GetKappa(void)
{
   /* This calculates mutational transition/transversion rate ratio kappa
      using 4-fold degenerate sites from pairwise comparisons
      under HKY85, weighting estimates by the numbers of sites
   */
   int is, js, j, k, h, i1, pos, c[2], aa[2], b[2][3], a, ndeg, by[3] = { 16,4,1 }, status = 0;
   double ka[2], F[2][16], S[2], wk[2], t, P, Q, pi[4];
   /* F&S&wk [0]: non-degenerate; [1]:4-fold;  S:sites */
   double kdefault = (com.kappa > 0 ? com.kappa : (com.icode == 1 ? 10 : 2));
   
   for (is = 0, com.kappa = 0; is < com.ns; is++) {
      for (js = 0; js < is; js++) {
         if (noisy >= 9) printf("\n%4d vs. %3d", is + 1, js + 1);
         for (k = 0; k < 2; k++) zero(F[k], 16);
         for (h = 0; h < com.npatt; h++) {
            c[0] = FROM61[(int)com.z[is][h]];
            c[1] = FROM61[(int)com.z[js][h]];
            for (k = 0; k < 2; k++) {
               b[k][0] = c[k] / 16;
               b[k][1] = (c[k] % 16) / 4;
               b[k][2] = c[k] % 4;
               aa[k] = GeneticCode[com.icode][c[k]];
            }

            /* find non-degenerate sites */
            for (pos = 0; pos < 3; pos++) {         /* check all positions */
               for (k = 0, ndeg = 0; k < 2; k++) {       /* two codons */
                  for (i1 = 0; i1 < 4; i1++) {
                     if (i1 == b[k][pos]) continue;
                     a = GeneticCode[com.icode][c[k] + (i1 - b[k][pos]) * by[pos]];
                     if (a == aa[k]) break;
                  }
                  if (i1 == 4) ndeg++;
               }
               if (ndeg == 2) {
                  F[0][b[0][pos] * 4 + b[1][pos]] += .5 * com.fpatt[h];
                  F[0][b[1][pos] * 4 + b[0][pos]] += .5 * com.fpatt[h];
               }

            }
            /* find 4-fold degenerate sites at 3rd positions */
            for (k = 0, ndeg = 0; k < 2; k++) {       /* two codons */
               for (j = 0, i1 = c[k] - b[k][2]; j < 4; j++)
                  if (j != b[k][2] && GeneticCode[com.icode][i1 + j] != aa[k]) break;
               if (aa[0] == aa[1] && j == 4) ndeg++;
            }
            if (ndeg < 2) continue;
            F[1][b[0][2] * 4 + b[1][2]] += .5 * com.fpatt[h];
            F[1][b[1][2] * 4 + b[0][2]] += .5 * com.fpatt[h];
         }  /* for (h) */
         for (k = 0; k < 2; k++) {  /* two kinds of sites */
            /*
            char *sitestr[2] = { "non-degenerate","4-fold" };
            if(noisy>3) printf("\n%s:\n",sitestr[k]);
            */
            S[k] = sum(F[k], 16);
            if (S[k] <= 0) { wk[k] = 0; continue; }
            for (j = 0; j < 16; j++) F[k][j] /= S[k];
            P = (F[k][0 * 4 + 1] + F[k][2 * 4 + 3]) * 2;
            Q = 1 - (F[k][0 * 4 + 0] + F[k][1 * 4 + 1] + F[k][2 * 4 + 2] + F[k][3 * 4 + 3]) - P;
            for (j = 0; j < 4; j++)
               pi[j] = sum(F[k] + j * 4, 4);
            DistanceF84(S[k], P, Q, pi, &ka[k], &t, NULL);
            wk[k] = (ka[k] > 0 ? S[k] : 0);

            /* matout(F0,F[k],4,4);  matout(F0,pi,1,4);  */
            /*
            if(noisy>3)
               printf("\nSPQkt:%9.4f%9.5f%9.5f%9.4f%9.4f\n",S[k],P,Q,ka[k],t);
            */
         }
         if (wk[0] + wk[1] == 0) {
            status = -1;
            ka[0] = kdefault;
            if (noisy > 3) printf("\ngot no kappa! fix it at %.4f\n", ka[0]);
         }
         else
            ka[0] = (ka[0] * wk[0] + ka[1] * wk[1]) / (wk[0] + wk[1]);
         com.kappa += ka[0] / (com.ns * (com.ns - 1.) / 2);
      }  /* for(js) */
   }     /* for(is) */

   return (status);
}


int CountSites(char z[], double pi[], double* Stot, double* Ntot, double fbS[], double fbN[])
{
   /* This calculates the total numbers of synonymous and nonsynonymous sites
      (Stot & Ntot) in the sequence z[] using com.kappa and pi[].
      It also count the base frequencies at the synonymous and nonsynonymous
      sites.  Total number of sites is scaled to be equal to sequence length
      even if some changes are to stop codons.  Since pi[] is scaled to sum
      to one, rates to stop codons are not considered.
      The counting goes through the sequence codon by codon, and so is different
      from the counting in codeml, which uses pi[] to count the sites.
   */
   int h, j, k, c[2], aa[2], b[3], by[3] = { 16,4,1 };
   double r, S, N, kappa = com.kappa;

   *Stot = *Ntot = 0;
   for (k = 0; k < 4; k++)
      fbS[k] = fbN[k] = 0;
   for (h = 0; h < com.npatt; h++) {
      c[0] = FROM61[(int)z[h]];
      b[0] = c[0] / 16; b[1] = (c[0] % 16) / 4; b[2] = c[0] % 4;
      aa[0] = GeneticCode[com.icode][c[0]];
      if (aa[0] == -1)
         zerror("stop codon");
      for (j = 0, S = N = 0; j < 3; j++) {
         for (k = 0; k < 4; k++) {    /* b[j] changes to k */
            if (k == b[j]) continue;
            c[1] = c[0] + (k - b[j]) * by[j];
            aa[1] = GeneticCode[com.icode][c[1]];
            if (aa[1] == -1) continue;
            r = pi[FROM64[c[1]]];
            if (k + b[j] == 1 || k + b[j] == 5)  r *= kappa; /* transition */
            if (aa[0] == aa[1]) { S += r; fbS[b[j]] += r * com.fpatt[h]; }
            else { N += r; fbN[b[j]] += r * com.fpatt[h]; }
         }
      }
      *Stot += com.fpatt[h] * S;
      *Ntot += com.fpatt[h] * N;
   }
   r = 3 * com.ls / (*Stot + *Ntot);  *Stot *= r;  *Ntot *= r;
   r = sum(fbS, 4);  for (k = 0; k < 4; k++) fbS[k] /= r;
   r = sum(fbN, 4);  for (k = 0; k < 4; k++) fbN[k] /= r;
   return (0);
}


int GetPMatCodon(double P[], double t, double kappa, double omega, double space[])
{
   /* Get PMat=exp(Q*t) for weighting pathways
   */
   int n = com.ncode, ic1, ic2, i, j, k, aa[2], ndiff, pos = 0, from[3], to[3];
   double* Q = P, * U = space + n * n, * V = U + n * n, * Root = V + n * n, mr, spacesqrt[NCODE];

   for (i = 0; i < n * n; i++) Q[i] = 0;
   for (i = 0; i < n; i++) {
      ic1 = FROM61[i]; from[0] = ic1 / 16; from[1] = (ic1 / 4) % 4; from[2] = ic1 % 4;
      for (j = 0; j < i; j++) {
         ic2 = FROM61[j];   to[0] = ic2 / 16;   to[1] = (ic2 / 4) % 4;   to[2] = ic2 % 4;
         aa[0] = GeneticCode[com.icode][ic1];
         aa[1] = GeneticCode[com.icode][ic2];
         if (aa[0] == -1 || aa[1] == -1)  continue;
         for (k = 0, ndiff = 0; k < 3; k++)
            if (from[k] != to[k]) { ndiff++; pos = k; }
         if (ndiff != 1)  continue;
         Q[i * n + j] = 1;
         if ((from[pos] + to[pos] - 1) * (from[pos] + to[pos] - 5) == 0)
            Q[i * n + j] *= kappa;
         if (aa[0] != aa[1])  Q[i * n + j] *= omega;
         Q[j * n + i] = Q[i * n + j];
      }
   }

   for (i = 0; i < n; i++) for (j = 0; j < n; j++)
      Q[i * n + j] *= com.pi[j];

   for (i = 0, mr = 0; i < n; i++) {
      Q[i * n + i] = -sum(Q + i * n, n);
      mr -= com.pi[i] * Q[i * n + i];
   }

   eigenQREV(Q, com.pi, n, Root, U, V, spacesqrt);
   for (i = 0; i < n; i++) Root[i] /= mr;
   PMatUVRoot(P, t, n, U, V, Root);
   return (0);
}



int CountDiffs(char z1[], char z2[], double* Sdts, double* Sdtv, double* Ndts, double* Ndtv, double PMat[])
{
   /* Count the numbers of synonymous and nonsynonymous differences between
      sequences z1 and z2, weighting pathways with PMat. No weighting if PMat=NULL
      Modified from difcodon()
      dmark[i] (=0,1,2) is the i_th different codon position (i=0,1,ndiff).
      step[j] (=0,1,2) is the codon position to be changed at step j (j=0,1,ndiff).
      b[i][j] (=0,1,2,3) is the nucleotide at position j (0,1,2) in codon i (0,1)
      sts,stv,nts,ntv are syn ts & tv and nonsyn ts & tv at a codon site.
      stspath[k] stvpath[k] ntspath[k] ntvpath[k] are syn ts & tv and
      nonsyn ts & tv differences on path k (k=2,6).
   */
   char str[4] = "   ";
   int n = com.ncode, h, i1, i2, i, k, transi, c[2], ct[2], aa[2], by[3] = { 16,4,1 };
   int dmark[3], step[3], b[2][3], bt1[3], bt2[3];
   int ndiff, npath, nstop, stspath[6], stvpath[6], ntspath[6], ntvpath[6];
   double sts, stv, nts, ntv; /* syn ts & tv, nonsyn ts & tv for 2 codons */
   double ppath[6], sump, p;

   *Sdts = *Sdtv = *Ndts = *Ndtv = 0;
   for (h = 0; h < com.npatt; h++) {
      c[0] = FROM61[(int)z1[h]];
      c[1] = FROM61[(int)z2[h]];
      if (c[0] == c[1]) continue;
      for (i = 0; i < 2; i++) {
         b[i][0] = c[i] / 16; b[i][1] = (c[i] % 16) / 4; b[i][2] = c[i] % 4;
         aa[i] = GeneticCode[com.icode][c[i]];
      }
      if (aa[0] == -1 || aa[1] == -1)
         zerror("stop codon in sequence.");
      ndiff = 0;  sts = stv = nts = ntv = 0;
      for (k = 0; k < 3; k++) dmark[k] = -1;
      for (k = 0; k < 3; k++) if (b[0][k] != b[1][k]) dmark[ndiff++] = k;
      npath = 1;
      if (ndiff > 1) npath = (ndiff == 2 ? 2 : 6);
      if (ndiff == 1) {
         transi = b[0][dmark[0]] + b[1][dmark[0]];
         transi = (transi == 1 || transi == 5);
         if (aa[0] == aa[1]) { if (transi) sts++; else stv++; }
         else { if (transi) nts++; else ntv++; }
      }
      else {   /* ndiff=2 or 3 */
         if (debug == DIFF) {
            printf("\n\nh=%d %s (%c) .. ", h + 1, getcodon(str, c[0]), AAs[aa[0]]);
            printf("%s (%c): ", getcodon(str, c[1]), AAs[aa[1]]);
         }
         nstop = 0;
         for (k = 0; k < npath; k++) {
            if (debug == DIFF) printf("\npath %d: ", k + 1);

            for (i1 = 0; i1 < 3; i1++)  step[i1] = -1;
            if (ndiff == 2) {
               step[0] = dmark[k];
               step[1] = dmark[1 - k];
            }
            else {
               step[0] = k / 2;
               step[1] = k % 2;
               if (step[0] <= step[1]) step[1]++;
               step[2] = 3 - step[0] - step[1];
            }
            for (i1 = 0; i1 < 3; i1++) bt1[i1] = bt2[i1] = b[0][i1];
            stspath[k] = stvpath[k] = ntspath[k] = ntvpath[k] = 0;
            /* mutations along each path */
            for (i1 = 0, ppath[k] = 1; i1 < ndiff; i1++) {
               bt2[step[i1]] = b[1][step[i1]];
               for (i2 = 0, ct[0] = ct[1] = 0; i2 < 3; i2++) {
                  ct[0] += bt1[i2] * by[i2];
                  ct[1] += bt2[i2] * by[i2];
               }
               ppath[k] *= PMat[FROM64[ct[0]] * n + FROM64[ct[1]]];
               for (i2 = 0; i2 < 2; i2++) aa[i2] = GeneticCode[com.icode][ct[i2]];

               if (debug == DIFF) printf("%s (%c) %.5f: ", getcodon(str, ct[1]), AAs[aa[1]], PMat[ct[0] * n + ct[1]]);

               if (aa[1] == -1) {
                  nstop++;  ppath[k] = 0; break;
               }
               transi = b[0][step[i1]] + b[1][step[i1]];
               transi = (transi == 1 || transi == 5);  /* transition? */

               if (aa[0] == aa[1]) { if (transi) stspath[k]++; else stvpath[k]++; }
               else { if (transi) ntspath[k]++; else ntvpath[k]++; }
               for (i2 = 0; i2 < 3; i2++) bt1[i2] = bt2[i2];
            }

            if (debug == DIFF) printf("  p =%.9f", ppath[k]);

         }  /* for(k,npath) */
         if (npath == nstop) {  /* all paths through stop codons */
            puts("all paths through stop codons..");
            if (ndiff == 2) { nts = .5; ntv = 1.5; }
            else { nts = .5; ntv = 2.5; }
         }
         else {
            sump = sum(ppath, npath);
            if (sump < 1e-20) {
               printf("\nsump=0, npath=%4d\nh=%2d ", npath, h + 1);
               printf("(%s ", getcodon(str, c[0]));
               printf("%s)", getcodon(str, c[1]));
               for (k = 0; k < npath; k++)
                  printf(" %9.6g", ppath[k]);
               printf("\n");
               matout(frub, PMat, n, n);
               exit(-1);
            }
            for (k = 0; k < npath; k++) {
               p = ppath[k] / sump;
               sts += stspath[k] * p;
               stv += stvpath[k] * p;
               nts += ntspath[k] * p;
               ntv += ntvpath[k] * p;
            }

            if (debug == DIFF) {
               for (k = 0; k < npath; k++)
                  printf("\n p =%.5f", ppath[k] / sump);
               printf("\n");
               printf(" syn ts & tv, nonsyn ts & tv:%9.5f%9.5f%9.5f%9.5f\n", sts, stv, nts, ntv);
            }
         }

         if (debug == DIFF) getchar();

      }     /* if (ndiff) */
      *Sdts += com.fpatt[h] * sts;
      *Sdtv += com.fpatt[h] * stv;
      *Ndts += com.fpatt[h] * nts;
      *Ndtv += com.fpatt[h] * ntv;
   }  /* for (h) */
   return (0);
}


int DistanceF84(double n, double P, double Q, double pi[],
   double* k_HKY, double* t, double* SEt)
{
   /* This calculates kappa and d from P (proportion of transitions) & Q
      (proportion of transversions) & pi under F84.
      When F84 fails, we try to use K80.  When K80 fails, we try
      to use JC69.  When JC69 fails, we set distance t to maxt.
      Variance formula under F84 is from Tateno et al. (1994), and briefly
      checked against simulated data sets.
   */
   int failF84 = 0, failK80 = 0, failJC69 = 0;
   double tc, ag, Y, R, a = 0, b = 0, A = -1, B = -1, C = -1, k_F84;
   double Qsmall = min2(1e-10, 0.1 / n), maxkappa = 999, maxt = 99;

   *k_HKY = -1;
   Y = pi[0] + pi[1];  R = pi[2] + pi[3];  tc = pi[0] * pi[1];  ag = pi[2] * pi[3];
   if (P + Q > 1) { *t = maxt; *k_HKY = 1; return(3); }
   if (P < -1e-10 || Q < -1e-10 || fabs(Y + R - 1)>1e-8) {
      printf("\nPQYR & pi[]: %9.5f%9.5f%9.5f%9.5f", P, Q, Y, R);
      matout(F0, pi, 1, 4);
      zerror("DistanceF84: input err.");
   }
   if (Q < Qsmall)  failF84 = failK80 = 1;
   else if (Y <= 0 || R <= 0 || (tc <= 0 && ag <= 0)) failF84 = 1;
   else {
      A = tc / Y + ag / R; B = tc + ag; C = Y * R;
      a = (2 * B + 2 * (tc * R / Y + ag * Y / R) * (1 - Q / (2 * C)) - P) / (2 * A);
      b = 1 - Q / (2 * C);
      if (a <= 0 || b <= 0) failF84 = 1;
   }
   if (!failF84) {
      a = -.5 * log(a); b = -.5 * log(b);
      if (b <= 0) failF84 = 1;
      else {
         k_F84 = a / b - 1;
         *t = 4 * b * (tc * (1 + k_F84 / Y) + ag * (1 + k_F84 / R) + C);
         *k_HKY = (B + (tc / Y + ag / R) * k_F84) / B; /* k_F84=>k_HKY85 */
         if (SEt) {
            a = A * C / (A * C - C * P / 2 - (A - B) * Q / 2);
            b = A * (A - B) / (A * C - C * P / 2 - (A - B) * Q / 2) - (A - B - C) / (C - Q / 2);
            *SEt = sqrt((a * a * P + b * b * Q - square(a * P + b * Q)) / n);
         }
      }
   }
   if (failF84 && !failK80) {  /* try K80 */
      if (noisy >= 9) printf("\na=%.5f  b=%.5f, use K80\n", a, b);
      a = 1 - 2 * P - Q;  b = 1 - 2 * Q;
      if (a <= 0 || b <= 0) failK80 = 1;
      else {
         a = -log(a); b = -log(b);
         if (b <= 0)  failK80 = 1;
         else {
            *k_HKY = (.5 * a - .25 * b) / (.25 * b);
            *t = .5 * a + .25 * b;
         }
         if (SEt) {
            a = 1 / (1 - 2 * P - Q); b = (a + 1 / (1 - 2 * Q)) / 2;
            *SEt = sqrt((a * a * P + b * b * Q - square(a * P + b * Q)) / n);
         }
      }
   }
   if (failK80) {
      if ((P += Q) >= .75) { failJC69 = 1; P = .75 * (n - 1.) / n; }
      *t = -.75 * log(1 - P * 4 / 3.);
      if (*t > maxt) *t = maxt;
      if (SEt) {
         *SEt = sqrt(9 * P * (1 - P) / n) / (3 - 4 * P);
      }
   }
   if (*k_HKY > maxkappa) *k_HKY = maxkappa;

   return(failF84 + failK80 + failJC69);
}



#if 0

double dsdnREV(int is, int js, double space[])
{
   /* This calculates ds and dn by recovering the Q*t matrix using the equation
         F(t) = PI * P(t) = PI * exp(Q*t)
      This is found not to work well and is not published.
      space[64*64*5]
      The code here is broken since I changed the coding.  Codons are now coded 0, 1, ..., 60.
   */
   int n = com.ncode, i, j, h;
   double* F = PMat, * Qt = F;
   double* Root = space + n * n, * pi = Root + n, * U = pi + n, * V = U + n * n;
   double* T1 = V + n * n, * T2 = T1 + n * n, t, small = 1e-6;

   fprintf(frst, "\npi in model\n");
   matout(frst, com.pi, 1, n);
   for(j=0; i<n*n; i++) 
      F[i] = 0;
   for(i=0; i<n; i++) {
      for (h = 0; h < com.npatt; h++) {
         F[com.z[is][h] * n + com.z[js][h]] += com.fpatt[h] / (2 * com.ls);
         F[com.z[js][h] * n + com.z[is][h]] += com.fpatt[h] / (2 * com.ls);
      }
   }
   if (fabs(1 - sum(F, n * n)) > 1e-6)
      zerror("Sum F != 1 in dsdnREV");

   for (i = 0; i < n; i++) {
      pi[i] = sum(F + i * n, n);
      /*
      if (F[i*n+i]<=small || F[i*n+i]<pi[i]/4)
      */
      if (F[i * n + i] <= small)  F[i * n + i] = 1 - pi[i] + F[i * n + i];
      else                  abyx(1 / pi[i], F + i * n, n);
   }
   if (eigen(1, F, n, Root, T1, U, V, T2)) zerror("eigen jgl");
   xtoy(U, V, n * n);
   matinv(V, n, n, T1);

   fprintf(frst, "\npi in data\n");
   matout(frst, pi, 1, n); 
   printf("\n");
   matout(frst, Root, 1, n);

   for(i=0; i<n; i++) {
      if (Root[i] <= 0)
         printf("  Root %d:%10.4f", i + 1, Root[i]);
      Root[i] = log(Root[i]);
   }
   for (i = 0; i < n; i++) 
      for (j = 0; j < n; j++)
         T1[i * n + j] = U[i * n + j] * Root[j];
   matby(T1, V, Qt, n, n, n);
   for (i = 0, t = 0; i < n; i++) t -= pi[i] * Qt[i * n + i];
   if (t <= 0) puts("err: dsdnREV");

   for (i = 0; i < n*n; i++)  /* remove negative numbers from rounding errors */
      Qt[i] += 1e-8; 

   matout(frst, Qt, n, n);
   printf("\nt = %.5f\n", t);

   return (0);
}


#endif
