/* paml.h 
*/

#if (!defined PAML_H)
#define PAML_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#define square(a) ((a)*(a))
#define FOR(i,n) for(i=0; i<n; i++)
#define FPN(file) fputc('\n', file)
#define F0 stdout
#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))
#define PI 3.141592653

#define beep putchar('\a')
#define spaceming2(n) ((n)*((n)*2+9+2)*(int)sizeof(double))

int ReadSeq (FILE *fout, FILE *fseq);
void EncodeSeqs (void);
void ReadPatternFreq (FILE* fout, char* fpattf);
int Initialize (FILE *fout);
int MoveCodonSeq (int ns, int ls, char *z[]);
int PatternWeight (FILE *fout);
int PatternJC69like (FILE *fout);
int PatternWeightSimple (int CollapsJC, double space[]);
int Site2Pattern (FILE *fout);

int f_and_x(double x[], double f[], int n, int fromx, int LastItem);
void SetSeed (unsigned int seed);
double rndu (void);
void randorder(int order[], int n, int space[]);
#define rnduab(a,b) ((a)+((b)-(a))*rndu())
double reflect(double x, double a, double b);
#define rndexp(mean) (-(mean)*log(rndu()))
double rndnormal(void);
double rndgamma(double s);
double rndbeta(double p, double q);
int rndpoisson(double m);
int rndNegBinomial(double shape, double mean);
int SampleCat (double P[], int n, double space[]);
int MultiNomial (int n, int nclass, double Prob[], int nobs[], double space[]);
int AutodGamma (double M[], double freqK[], double rK[], double *rho1,
    double alfa, double rho, int K);
int DiscreteGamma (double freqK[], double rK[], 
    double alfa, double beta, int K, int median);
double PointChi2 (double prob, double v);
#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))
double PDFGamma(double x, double alpha, double beta);
#define CDFGamma(x,alpha,beta) IncompleteGamma((beta)*(x),alpha,LnGamma(alpha))
double  PDF_IGamma(double x, double alpha, double beta);
#define CDF_IGamma(x,alpha,beta) (1-CDFGamma(1/(x),alpha,beta))
#define CDFChi2(x,v) CDFGamma(x,(v)/2.0,0.5)
int getab_gamma(double *a, double* b, double xmodetiles[3]);
double CDFBeta(double x, double p, double q, double lnbeta);
double InverseCDFBeta(double prob, double p, double q, double lnbeta);
double InverseCDF(double(*cdf)(double x,double par[]),
       double p,double x,double par[],double xb[2]);
double PointNormal (double prob);
double CDFNormal (double x);
double LnGamma (double alpha);
double DFGamma(double x, double alpha, double beta);
double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
double LBinormal(double h, double k, double r);
double CDFBinormal (double h, double k, double r);
double probBinomial (int n, int k, double p);
double probBetaBinomial (int n, int k, double p, double q);
long factorial(int n);
double Binomial(double n, int k, double *scale);

int GaussLegendreRule(double **x, double **w, int order);
int GaussLaguerreRule(double **x, double **w, int order);
double NIntegrateGaussLegendre(double(*fun)(double x), double a, double b, int order);

int ScatterPlot (int n, int nseries, int yLorR[], double x[], double y[],
    int nrow, int ncol, int ForE);
void rainbowRGB (double temperature, int *R, int *G, int *B);
void GetIndexTernary(int *ix, int *iy, double *x, double *y, int itriangle, int K);

int CodeChara (char b, int seqtype);
int dnamaker (char z[], int ls, double pi[]);
int picksite (char z[], int ls, int begin, int gap, char buffer[]);
int transform (char z[], int ls, int direction, int seqtype);
int RemoveIndel(void);
int f_mono_di (FILE *fout, char z[], int ls, int iring, 
    double fb1[], double fb2[], double CondP[]);
int PickExtreme (FILE *fout, char z[], int ls,int iring,int lfrag,int ffrag[]);

int print1seq (FILE*fout, char *z, int ls, int encoded, int pose[]);
void printSeqs(FILE *fout, int *pose, char keep[], int format);
void printSeqsMgenes (void);
int printsma(FILE*fout, char*spname[], char*z[],
    int ns, int l, int lline, int gap, int seqtype, 
    int transformed, int simple, int pose[]);
int printsmaCodon (FILE *fout,char * z[],int ns,int ls,int lline,int simple);
int zztox ( int n31, int l, char z1[], char z2[], double *x );
int testXMat (double x[]);
double SeqDivergence (double x[], int model, double alpha, double *kapa);
int symtest (double freq[], int n, int nobs, double space[], double *chisym, 
    double* chihom);
int dSdNNG1986(char *z1, char *z2, int lc, int icode, int transfed, 
    double *dS, double *dN, double *S, double *N);
int difcodonNG(char codon1[], char codon2[], double *SynSite,double *AsynSite, 
    double *SynDif, double *AsynDif, int transfed, int icode);
int testTransP (double P[], int n);
void pijJC69 (double pij[2], double t);
int PMatK80 (double P[], double t, double kapa);
int PMatT92 (double P[], double t, double kappa, double pGC);
int PMatTN93 (double P[], double a1t, double a2t, double bt, double pi[]);
int PMatUVRoot (double P[],double t,int n,double U[],double V[],double Root[]);
int PMatCijk (double PMat[], double t);
int PMatQRev(double P[], double pi[], double t, int n, double space[]);
int EvolveHKY85 (char source[], char target[], int ls, double t, 
    double rates[], double pi[], double kapa, int isHKY85);
double DistanceIJ (int is, int js, int model, double alpha, double *kappa);
int DistanceMatNuc (FILE *fout, FILE*f2base, int model, double alpha);
int EigenQREVbase (FILE* fout, double kappa[], double pi[], 
                   int *nR, double Root[], double Cijk[]);

int BootstrapSeq (char* seqfilename);
int rell(FILE*flnf, FILE*fout, int ntree);
int print1site (FILE*fout, int h);
int MultipleGenes (FILE* fout, FILE*fpair[], double space[]);
int lfunRates (FILE* fout, double x[], int np);
int AncestralSeqs (FILE *fout, double x[]);
void ListAncestSeq(FILE *fout, char *zanc);
int ChangesSites(FILE*fout, int coding, char *zanc);
int InitConditionalPNode (void);

int NucListall(char b, int *nb, int ib[4]);
char *getcodon(char codon[], int icodon);
char *getAAstr(char *AAstr, int iaa);
int Codon2AA(char codon[3], char aa[3], int icode, int *iaa);
int DNA2protein(char dna[], char protein[], int lc, int icode);
int printcu (FILE *f1, double fcodon[], int icode);
int printcums (FILE *fout, int ns, double fcodons[], int code);
int QtoPi (double Q[], double pi[], int n, double *space);
int PtoPi (double P[], double Pi[], int n, double *space);
int PtoX (double P1[], double P2[], double Pi[], double X[]);

void starttime(void);
char* printtime(char timestr[]);
void sleep2(int wait);
char *strc (int n, char c);
int printdouble(FILE*fout, double a);
void strcase (char *str, int direction);
void error2(char * message);
int  sort1(double x[], int n, int rank[], int descending, int space[]);
FILE *gfopen(char *filename, char *mode);
int  appendfile(FILE*fout, char*filename);

int zero (double x[], int n);
double sum (double x[], int n);
int fillxc (double x[], double c, int n);
int xtoy (double x[], double y[], int n);
int abyx (double a, double x[], int n);
int axtoy(double a, double x[], double y[], int n);
int axbytoz(double a, double x[], double b, double y[], double z[], int n);
int identity (double x[], int n);
double distance (double x[], double y[], int n);
double innerp(double x[], double y[], int n);
double norm(double x[], int n);

double Det3x3 (double x[3*3]);
int matby (double a[], double b[], double c[], int n,int m,int k);
int matIout (FILE *fout, int x[], int n, int m);
int matout (FILE *file, double x[], int n, int m);
int matout2 (FILE *fout, double x[], int n, int m, int wid, int deci);
int mattransp1 (double x[], int n);
int mattransp2 (double x[], double y[], int n, int m);
int matinv (double x[], int n, int m, double space[]);
int matexp (double Q[], double t, int n, int TimeSquare, double space[]);
int matsqrt (double A[], int n, double work[]);
int CholeskyDecomp (double A[], int n, double L[]);
int eigenQREV (double Q[], double pi[], int n, 
               double Root[], double U[], double V[], double spacesqrtpi[]);
int eigenRealSym(double A[], int n, double Root[], double work[]);

int MeanVar (double x[], int n, double *mean, double *var);
int variance (double x[], int n, int nx, double mx[], double vx[]);
int correl (double x[], double y[], int n, double *mx, double *my,
            double *v11, double *v12, double *v22, double *r);
int comparefloat  (const void *a, const void *b);
int comparedouble (const void *a, const void *b);
int DescriptiveStatistics(FILE *fout, char infile[], int nbin, int nrho);
int splitline (char line[], int fields[]);

double bound (int nx, double x0[], double p[], double x[],
    int (*testx) (double x[], int nx));
int gradient (int n, double x[], double f0, double g[], 
    double (* fun)(double x[],int n), double space[], int Central);
int Hessian (int nx, double x[], double f, double g[], double H[],
    double (*fun)(double x[], int n), double space[]);

int H_end (double x0[], double x1[], double f0, double f1,
    double e1, double e2, int n);
double LineSearch(double(*fun)(double x),double *f,double *x0,double xb[2],double step,double e);
double LineSearch2 (double(*fun)(double x[],int n), double *f, double x0[], 
    double p[], double h, double limit, double e, double space[], int n);

void xtoFreq(double x[], double freq[], int n);


int SetxBound (int np, double xb[][2]);
int ming1 (FILE *fout, double *f, double (* fun)(double x[], int n),
    int (*dfun) (double x[], double *f, double dx[], int n),
    int (*testx) (double x[], int n),
    double x0[], double space[], double e, int n);
int ming2 (FILE *fout, double *f, double (*fun)(double x[], int n),
    int (*dfun)(double x[], double *f, double dx[], int n),
    double x[], double xb[][2], double space[], double e, int n);
int minB (FILE*fout, double *lnL,double x[],double xb[][2],double e, double space[]);
int minB2 (FILE*fout, double *lnL,double x[],double xb[][2],double e, double space[]);


int Newton (FILE *fout, double *f, double (* fun)(double x[], int n),
    int (* ddfun) (double x[], double *fx, double dx[], double ddx[], int n),
    int (*testx) (double x[], int n),
    double x0[], double space[], double e, int n);

int nls2 (FILE *fout, double *sx, double * x0, int nx,
    int (* fun)(double x[], double y[], int nx, int ny),
    int (* jacobi)(double x[], double J[], int nx, int ny),
    int (*testx) (double x[], int nx),
    int ny, double e);

/* tree structure functions in treesub.c */
void NodeToBranch (void);
void BranchToNode (void);
void ClearNode (int inode);
int ReadaTreeN (FILE *ftree, int *haslength, int *haslabel, int copyname, int popline);
int ReadaTreeB (FILE *ftree, int popline);
int OutaTreeN (FILE *fout, int spnames, int printopt);
int OutaTreeB (FILE *fout);
int DeRoot (void);
int GetSubTreeN (int keep[], int space[]);
void printtree (int timebranches);
void PointconPnodes (void);
int SetBranch (double x[]);
int DistanceMat (FILE *fout, int ischeme, double alfa, double *kapa);
int LSDistance (double * ss, double x[], int (*testx)(double x[],int np));

int StepwiseAdditionMP (double space[]);
double MPScoreStepwiseAddition (int is, double space[], int save);
int AddSpecies (int species, int ib);
int StepwiseAddition (FILE* fout, double space[]);
int readx(double x[], int *fromfile);
double TreeScore(double x[], double space[]);

int PopEmptyLines (FILE* fseq, int lline, char line[]);
int hasbase (char *str);
int blankline (char *str);

void BranchLengthBD(int rooted, double birth, double death, double sample, 
     double mut);
int RandomLHistory (int rooted, double space[]);

void DescentGroup (int inode);
void BranchPartition (char partition[], int parti2B[]);
int NSameBranch (char partition1[],char partition2[], int nib1,int nib2,
    int IBsame[]);

int RootTN93 (int ischeme, double kapa1, double kapa2, double pi[], 
    double *f, double Root[]);
int EigenTN93 (int ischeme, double kapa1, double kapa2, double pi[],
    int *nR, double Root[], double Cijk[]);

int DownStatesOneNode (int ison, int father);
int DownStates (int inode);
int PathwayMP (FILE *fout, double space[]);
double MPScore (double space[]);
double RemoveMPNinfSites (double *nsiteNinf);

int MPInformSites (void);
int CollapsNode (int inode, double x[]);

int MakeTreeIb (int ns, int Ib[], int rooted);
int GetTreeI (int itree, int ns, int rooted);
int NumberofTrees (int ns, int rooted);
int ListTrees (FILE* fout, int ns, int rooted);
int GetIofTree (int rooted, int keeptree, double space[]);
void ReRootTree (int newroot);
int NeighborNNI (int i_tree);
int GetLHistoryI (int iLH);
int GetIofLHistory (void);
int CountLHistory(char LHistories[], double space[]);
int ReorderNodes (char LHistory[]);

int GetSubSeqs(int nsnew);

/* functions for evolving sequences */
int GenerateSeq (void);
int Rates4Sites (double rates[],double alfa,int ncatG,int ls, int cdf,
    double space[]);
void Evolve (int inode);
void EvolveJC (int inode);


/* dating using complex data */
int ReadTreeSeqs(FILE*fout);
int UseLocus (int locus, int copyconP, int setmodel, int setSeqName);
int GetGtree(int locus);
void printGtree(int printBlength);
void copySptree(void);
void printSptree(void);


#define  BASEseq     0
#define  CODONseq    1
#define  AAseq       2
#define  CODON2AAseq 3


#define FAST_RANDOM_NUMBER


#define RELEASE      0

#define VerStr "paml 3.14, January 2004"


/*
#define DEBUG 9
*/

#endif
