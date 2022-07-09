/* paml.h 
*/

#if (!defined PAML_H)
#define PAML_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <search.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <assert.h>

#define square(a) ((a)*(a))
#define F0 stdout
#if !defined(MAX)
#define MAX(a,b)                            (((a) > (b)) ? (a) : (b))
#endif
#if !defined(MIN)
#define MIN(a,b)                            (((a) < (b)) ? (a) : (b))
#endif

#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))
#define swap2(a,b,y) { y=a; a=b; b=y; }
#define Pi  3.1415926535897932384626433832795

int ReadSeq (FILE* fout, FILE* fseq, int cleandata, int locus, int read_seq_only);
int ScanFastaFile (FILE *f, int *ns, int *ls, int *aligned);
void EncodeSeqs (void);
void SetMapAmbiguity(int seqtype, int ModelAA2Codon);
void ReadPatternFreq (FILE* fout, char* fpattf);
int Initialize (FILE *fout);
int MoveCodonSeq (int ns, int ls, char *z[]);
int PatternWeight (void);
int PatternWeightJC69like (void);
int PatternWeightSimple (void);
int Site2Pattern (FILE *fout);

double getRoot(double (*f)(double), double (*df)(double), double initVal);
int f_and_x(double x[], double f[], int n, int fromx, int LastItem);
void bigexp(double lnx, double *a, double *b);
void SetSeed (int seed, int PrintSeed);
double rndu (void);
void rndu_vector (double r[], int n);
void randorder(int order[], int n, int space[]);
#define rnduab(a,b) ((a)+((b)-(a))*rndu())
double reflect(double x, double a, double b);
#define rndexp(mean) (-(mean)*log(rndu()))
double rnduM0V1 (void);
double rndNormal(void);
int rndNp(double x[], int n, int p, double mx[], double vx[], int isvroot);
int rndBinomial(int n, double p);
double rndBox(void);
double rndAirplane(void);
double rndStrawhat(void);
double rndBactrian(void);
double rndBactrianTriangle(void);
double rndBactrianLaplace(void);
double rndTriangle(void);
double rndLaplace (void);
double rndCauchy (void);
double rndt2 (void);
double rndt4 (void);
double rndlogt2 (double loc, double s);
double rndlogistic (void);
double rndloglogistic (double loc, double s);
double rndgamma(double alpha);
double rndbeta(double p, double q);
int rnddirichlet (double x[], double alpha[], int K);
int rndpoisson(double m);
int rndNegBinomial(double shape, double mean);
int SampleCat (double P[], int n, double space[]);
int MultiNomialAliasSetTable (int ncat, double prob[], double F[], int L[], double space[]);
int MultiNomialAlias (int n, int ncat, double F[], int L[], int nobs[]);
int MultiNomial2 (int n, int ncat, double prob[], int nobs[], double space[]);
int DiscreteBeta (double freq[], double x[], double p, double q, int K, int UseMedian);
int DiscreteGamma (double freqK[], double rK[], double alpha, double beta, int K, int UseMedian);
int AutodGamma (double Mmat[], double freqK[], double rK[], double *rho1, double alfa, double rho, int K);

double QuantileChi2 (double prob, double v);
#define QuantileGamma(prob,alpha,beta) QuantileChi2(prob,2.0*(alpha))/(2.0*(beta))
double  PDFGamma(double x, double alpha, double beta);
double logPDFGamma(double x, double alpha, double beta);
#define CDFGamma(x,alpha,beta) IncompleteGamma((beta)*(x),alpha,lgamma(alpha))
double logPriorRatioGamma(double xnew, double xold, double a, double b);
double  PDFinvGamma(double x, double alpha, double beta);
#define CDFinvGamma(x,alpha,beta) (1-CDFGamma(1/(x),alpha,beta))
#define CDFChi2(x,v) CDFGamma(x,(v)/2.0,0.5)
double logPriorRatioBeta(double xnew, double x, double p, double q);
double logPDFBeta(double x, double p, double q);
double PDFBeta(double x, double p, double q);
double CDFBeta(double x, double p, double q, double lnbeta);
double QuantileBeta(double prob, double p, double q, double lnbeta);
double Quantile(double(*cdf)(double x,double par[]), double p, double x, double par[], double xb[2]);
double QuantileNormal (double prob);
double PDFNormal (double x, double mu, double sigma2);
double logPDFNormal(double x, double mu, double sigma2);
double CDFNormal (double x);
double logCDFNormal (double x);
double PDFCauchy (double x, double m, double sigma);
double PDFloglogistic (double x, double loc, double s);
double PDFlogt2 (double x, double loc, double s);
double PDFt2 (double x, double m, double s);
double PDFt4 (double x, double m, double s);
double PDFt (double x, double loc, double scale, double df, double lnConst);
double CDFt (double x, double loc, double scale, double df, double lnbeta);
double PDFSkewT (double x, double loc, double scale, double shape, double df);
double PDFSkewN (double x, double loc, double scale, double shape);
double logPDFSkewN(double x, double loc, double scale, double shape);

int StirlingS2(int n, int k);
double lnStirlingS2(int n, int k);
double LnGamma(double alpha);
#define LnBeta(p,q) (lgamma(p) + lgamma(q) - lgamma(p+q))
double DFGamma(double x, double alpha, double beta);
double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
#define CDFBinormal(h,k,r)  LBinormal(-(h),-(k),r)   /* CDF for bivariate normal */
double LBinormal (double h, double k, double r);
double logLBinormal (double h, double k, double r);
double probBinomial (int n, int k, double p);
double probBetaBinomial (int n, int k, double p, double q);
double factorial (int n);
double Binomial(double n, int k, double *scale);
int BinomialK(double alpha, int n, double C[], double S[]);
int GaussLegendreRule(const double **x, const double **w, int order);
int GaussLaguerreRule(const double **x, const double **w, int order);
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
int f_mono_di (FILE *fout, char z[], int ls, int iring, double fb1[], double fb2[], double CondP[]);
int PickExtreme (FILE *fout, char z[], int ls, int iring, int lfrag, int ffrag[]);

int print1seq (FILE*fout, char *z, int ls, int pose[]);
void printSeqs(FILE *fout, char *z[], char *spnames[], int ns, int ls, int npatt, double fpatt[], int *pose, char keep[], int format);
int printPatterns(FILE *fout);
void printSeqsMgenes (void);
int printsma (FILE*fout, char*spname[], char*z[], int ns, int l, int lline, int gap, int seqtype, 
    int transformed, int simple, int pose[]);
int printsmaCodon (FILE *fout, char * z[], int ns, int ls, int lline, int simple);
int zztox ( int n31, int l, char z1[], char z2[], double *x );
int testXMat (double x[]);
double SeqDivergence (double x[], int model, double alpha, double *kapa);
int symtest (double freq[], int n, int nobs, double space[], double *chisym, 
    double* chihom);
int dSdNNG1986(char *z1, char *z2, int lc, int icode, int transfed, double *dS, double *dN, double *Ssites, double *Nsites);
int difcodonNG(char codon1[], char codon2[], double *SynSite,double *AsynSite, double *SynDif, double *AsynDif, int transfed, int icode);
int difcodonLWL85 (char codon1[], char codon2[], double sites[3], double sdiff[3], double vdiff[3], int transfed, int icode);
int testTransP (double P[], int n);
double testDetailedBalance (double P[], double pi[], int n);
void pijJC69 (double pij[2], double t);
int PMatK80 (double P[], double t, double kapa);
int PMatT92 (double P[], double t, double kappa, double pGC);
int PMatTN93 (double P[], double a1t, double a2t, double bt, double pi[]);
int PMatUVRoot (double P[],double t,int n,double U[],double V[],double Root[]);
int PMatCijk (double PMat[], double t);
int PMatQRev(double P[], double pi[], double t, int n, double space[]);
int EvolveHKY85 (char source[], char target[], int ls, double t, double rates[], double pi[], double kapa, int isHKY85);
double DistanceIJ (int is, int js, int model, double alpha, double *kappa);
int DistanceMatNuc (FILE *fout, FILE*f2base, int model, double alpha);
int EigenQGTRbase (FILE* fout, double kappa[], double pi[], int *nR, double Root[], double Cijk[]);
int  DistanceMatNG86 (FILE *fout, FILE*fds, FILE*fdn, FILE*dt, double alpha);
int  setmark_61_64 (void);

int BootstrapSeq (char* seqfilename);
int rell(FILE*flnf, FILE*fout, int ntree);
int print1site (FILE*fout, int h);
int MultipleGenes (FILE* fout, FILE* ftree, FILE*fpair[], double space[]);
int lfunRates (FILE* fout, double x[], int np);
int AncestralSeqs (FILE *fout, double x[]);
void ListAncestSeq(FILE *fout, char *zanc);
int ChangesSites(FILE*fout, int coding, char *zanc);

int NucListall(char b, int *nb, int ib[4]);
char *getcodon(char codon[], int icodon);
char *getAAstr(char *AAstr, int iaa);
int Codon2AA(char codon[3], char aa[3], int icode, int *iaa);
int DNA2protein(char dna[], char protein[], int lc, int icode);
int printcu (FILE *f1, double fcodon[], int icode);
int printcums (FILE *fout, int ns, double fcodons[], int code);
int QtoPi (double Q[], double pi[], int n, double *space);
int PtoPi (double P[], double pi[], int n, double *space);
int PtoX (double P1[], double P2[], double pi[], double X[]);

void starttimer(void);
char* printtime(char timestr[]);
void sleep2(int wait);
char *strc (int n, int c);
int printdouble(FILE*fout, double a);
void strcase (char *str, int direction);
void error2(char * message);
int  indexing(double x[], int n, int index[], int descending, int space[]);
int binarysearch (const void *key, const void *base, size_t n, size_t size, int(*compare)(const void *, const void *), int *found);
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
int matbytransposed (double a[], double b_transposed[], double c[], int n, int m, int k);
int matIout (FILE *fout, int x[], int n, int m);
int matout (FILE *file, double x[], int n, int m);
int matout2 (FILE *fout, double x[], int n, int m, int wid, int deci);
int mattransp1 (double x[], int n);
int mattransp2 (double x[], double y[], int n, int m);
int matinv (double x[], int n, int m, double space[]);
int matexp (double A[], int n, int nTaylorTerms, int nSquares, double space[]);
#ifdef USE_GSL
int matexpGSL (double A[], int n, double space[]);
#endif
int matsqrt (double A[], int n, double work[]);
int CholeskyDecomp (double A[], int n, double L[]);
int eigenQREV (double Q[], double pi[], int n, double Root[], double U[], double V[], double spacesqrtpi[]);
int eigenRealSym(double A[], int n, double Root[], double work[]);

int MeanVar (double x[], int n, double *mean, double *var);
int variance (double x[], int n, int nx, double mx[], double vx[]);
int correl (double x[], double y[], int n, double *mx, double *my, double *vxx, double *vxy, double *vyy, double *r);
int comparefloat  (const void *a, const void *b);
int comparedouble (const void *a, const void *b);
double Eff_IntegratedCorrelationTime (double x[], int n, double *mx, double *vx, double *rho1);
double Eff_IntegratedCorrelationTime2(double x[], int n, int nbatch, double *mx, double *vx);

int HPDinterval(double x[], int n, double HPD[2], double alpha);
int DescriptiveStatistics(FILE *fout, char infile[], int nbin, int propternary, int SkipColumns);
int DescriptiveStatisticsSimple (FILE *fout, char infile[], int SkipColumns);
int splitline(char line[], int nfields, int fields[]);
int scanfile (FILE*fin, int *nrecords, int *nx, int *HasHeader, char line[], int ifields[]);

double bound (int nx, double x0[], double p[], double x[], int (*testx) (double x[], int nx));
int gradient (int n, double x[], double f0, double g[], double (* fun)(double x[],int n), double space[], int Central);
int Hessian (int nx, double x[], double f, double g[], double H[], double (*fun)(double x[], int n), double space[]);
int HessianSKT2004 (double xmle[], double lnLm, double g[], double H[]);

int H_end (double x0[], double x1[], double f0, double f1, double e1, double e2, int n);
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

int ResetStepLengths(FILE *fout, double Pjump[], double finetune[], int nsteps);

/* tree structure functions in treesub.c */
void NodeToBranch (void);
void BranchToNode (void);
void ClearNode (int inode);
int ReadTreeN (FILE *ftree, int *haslength, int copyname, int popline);
int ReadTreeB (FILE *ftree, int popline);
int OutTreeN (FILE *fout, int spnames, int printopt);
int OutTreeB (FILE *fout);
int DeRoot (void);
int GetSubTreeN (int keep[], int space[]);
void PrintTree (int timebranches);
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

void BranchLengthBD(int rooted, double birth, double death, double sample, double mut);
int RandomLHistory (int rooted, double space[]);

void Tree2Partition (char partition[]);
int Partition2Tree (char splits[], int lsplit, int ns, int nsplit, double label[]);
void CladeSupport (FILE *fout, char treef[], int getSnames, char maintreef[], int pick1tree);
int GetNSfromTreeFile(FILE *ftree, int *ns, int *ntree);
int NSameBranch (char partition1[],char partition2[], int nib1,int nib2, int IBsame[]);

int QTN93 (int model, double Q[], double kappa1, double kappa2, double pi[]);
int RootTN93 (int ischeme, double kapa1, double kapa2, double pi[], double *scalefactor, double Root[]);
int EigenTN93 (int ischeme, double kapa1, double kapa2, double pi[], int *nR, double Root[], double Cijk[]);

int DownStatesOneNode (int ison, int father);
int DownStates (int inode);
int PathwayMP (FILE *fout, double space[]);
double MPScore (double space[]);
double RemoveMPNinfSites (double *nsiteNinf);
int MarkStopCodons(void);

int MPInformSites (void);
int CollapsNode (int inode, double x[]);

int MakeTreeIb (int ns, int Ib[], int rooted);
int GetTreeI (int itree, int ns, int rooted);
double CountTrees (int ns, int rooted);
double CountLHs (int ns);
int CountLHsTree (void);
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
int Rates4Sites (double rates[],double alfa,int ncatG,int ls, int cdf, double space[]);
void Evolve (int inode);
void EvolveJC (int inode);


int ReadTreeSeqs(FILE*fout);
int UseLocus (int locus, int copyconP, int setmodel, int setSeqName);
int GetGtree(int locus);
int printGtree(int printBlength);

void copySptree(void);
void printSptree(void);


extern int BASEseq, CODONseq, AAseq, CODON2AAseq, BINARYseq, BASE5seq;
extern int PrBranch, PrNodeNum, PrLabel, PrNodeStr, PrAge, PrOmega;


/* use mean (0; default) for discrete gamma instead of median (1) */
#define DGammaUseMedian 0

#define FAST_RANDOM_NUMBER

#define mBactrian  0.90
#define sBactrian  sqrt(1 - mBactrian*mBactrian)
#define aBox 0.5
#define bBox (sqrt(12 - 3*aBox*aBox) - aBox) / 2
#define aAirplane 1.0
#define aStrawhat 1.0

#define MAXNFIELDS 320000

#define PAML_RELEASE      0

#define pamlVerStr "paml version 4.10.5, March 2022"

#endif
