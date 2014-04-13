/*   tools.h */

#ifndef TOOLS_H
#define TOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#define square(a) ((a)*(a))
#define FOR(i,n) for(i=0; i<n; i++)
#define FPN(file) fputc('\n', file)
#define F0 stdout
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))
#define PI 3.141592653

int ReadSeq (FILE *fout, FILE *fseq, int seqtype);
int Initialize (FILE *fout, double space[], int seqtype);
int MoveCodonSeq (int ns, int ls, char *z[]);
int PatternWeight (FILE *fout, double space[]);
int PatternJC69like (FILE *fout);
int PatternWeightSimple (int CollapsJC, double space[], int* nconst);
int Site2Pattern (FILE *fout);

double rndu (void);
void SetSeed (int seed);
double rndnorm (void);
double rndgamma (double s);
int rndpoisson (double m);
int rndNegBinomial (double shape, double mean);
int SampleCat (double P[], int n, double space[]);
int MultiNomial (int n, int nclass, double Prob[], int nobs[], double space[]);
int AutodGamma (double M[], double freqK[], double rK[], double *rho1,
    double alfa, double rho, int K);
int DiscreteGamma (double freqK[], double rK[], 
    double alfa, double beta, int K, int median);
double PointChi2 (double prob, double v);
#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define CDFGamma(x,alpha,beta) IncompleteGamma((beta)*(x),alpha,LnGamma(alpha))
#define CDFChi2(x,v) CDFGamma(x,(v)/2.0,0.5)
double PointNormal (double prob);
double CDFNormal (double x);
double LnGamma (double alpha);
double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
double LBinormal(double h, double k, double r);
double CDFBinormal (double h, double k, double r);
double probBinomial (int n, int k, double p);
long factorial(int n);
double NchooseK(int n, int k);

char CodeChara (char b, int seqtype);
int dnamaker (char z[], int ls, double pi[]);
int picksite (char z[], int ls, int begin, int gap, char buffer[]);
int transform (char z[], int ls, int direction, int seqtype);
int RemoveIndel(char *z[], int *ns,int *ls,int posit[],int n31,int concensus);
int f_mono_di (FILE *fout, char z[], int ls, int iring, 
    double fb1[], double fb2[], double CondP[]);
int PickExtreme (FILE *fout, char z[], int ls,int iring,int lfrag,int ffrag[]);

int getseqma( char *seqfile, char  *z[], int *ns, int *l );
int printaSeq (FILE *fout, char z[], int ls, int lline, int gap);
int printsma (FILE *fout, char * z[], int ns, int l, int lline, int gap, 
    int transfed, int simple, int seqtype);
int printsmaPose (FILE *fout, char * z[], int ns, int l, int lline, int gap, 
    int transfed, int simple, int seqtype, int pose[]);
int zztox ( int n31, int l, char z1[], char z2[], double *x );
int testXMat (double x[]);
double SeqDivergence (double x[], int model, double alpha, double *kapa);
int symtest (double freq[], int n, int nobs, double space[], double *chisym, 
    double* chihom);
int ng1986_1(char z1[], char z2[], int lc, int code, int transfed, 
    double *Ks, double *Ka, double *S, double *N);
int difcodon (char codon1[], char codon2[], double *SynSite, double *AsynSite, 
    double *SynDif, double *AsynDif, int transfed, int icode);
int testTransP (double P[], int n);
int PMatUVRoot (double P[],double t,int n,double U[],double V[],double Root[]);
int PMatK80 (double P[], double t, double kapa);
int PMatTN93 (double P[], double a1t, double a2t, double bt, double pi[6]);
int EvolveHKY85 (char source[], char target[], int ls, double t, 
    double rates[], double pi[6], double kapa, int isHKY85);
int DistanceMatNuc (FILE *fout, int model, double alfa);
int EigenREV (FILE* fout, double kapa[], double pi[], 
              int *nR, double Root[], double Cijk[]);

int MultipleGenes (FILE* fout, double space[]);
int lfunAdG_rate (FILE* fout, double x[], int np);
int AncestralSeqs (FILE *fout, double x[], double space[]);

int pair(int i, int j);
char *getcode(int id);
char *getaa (int iaa, int id);
int getaa1 (char *dna, char *aa, int code, int iaa);
int printcu (FILE *f1, double *fb3, int code, double *space);
int printcums (FILE *fout, int ns, double *fb3, int code);
int PtoPi (double P[], double Pi[], int n, double *space);
int PtoX(double P1[], double P2[], double Pi[], double X[]);

char *strc (int n, char c);
void error(char * message);
int indexing (double x[], int n, int index[], int descending, double space[]);

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
int QRdecomp (double A[], int m, int n, double Q[]);
int CholeskyDecomp (double A[], int n, double L[]);

int eigen (int job, double A[], int n, double rr[], double ri[],
          double vr[], double vi[], double w[]);

int MeanVar (double x[], int n, double *mean, double *var);
int variance (double x[], int n, int nx, double mx[], double vx[]);
int correl (double x[], double y[], int n, double *mx, double *my,
            double *v11, double *v12, double *v22, double *r);

/* complex functions */
typedef struct { double re, im; } complex;
#define csize(a) (fabs(a.re)+fabs(a.im))

complex compl (double re,double im);
complex conj (complex a);
complex cplus (complex a, complex b);
complex cminus (complex a, complex b);
complex cby (complex a, complex b);
complex cdiv (complex a,complex b);
complex cexp (complex a);
complex cfactor (complex x, double a);
int cxtoy (complex x[], complex y[], int n);
int cmatby (complex a[], complex b[], complex c[], int n,int m,int k);
int cmatout (FILE * fout, complex x[], int n, int m);
int cmatinv( complex x[], int n, int m, double space[]);
int EigenUNREST (FILE *fout, double kapa[], double pi[], 
    int *nR, complex cRoot[], complex cU[], complex cV[]);
int cPMat (double P[],double t,int n,complex cU[],complex cV[],complex cRoot[]);

double bound (int nx, double x0[], double p[], double x[],
    int (*testx) (double x[], int nx));
int gradient (int n, double x[], double f0, double g[], 
    double (* fun)(double x[],int n), double space[], int Central);
int Hessian (int nx, double x[], double f, double g[], double H[],
    double (*fun)(double x[], int n), double space[]);

int H_end (double x0[], double x1[], double f0, double f1,
    double e1, double e2, int n);
double LineSearch2 (double(*fun)(double x[],int n), double *f, double x0[], 
    double p[], double h, double limit, double e, double space[], int n);

int SetxBound (int np, double xb[][2]);
int ming2 (FILE *fout, double *f, double (*fun)(double x[], int n),
    int (*dfun)(double x[], double *f, double dx[], int n),
    double x[], double xb[][2], double space[], double e, int n);

int ming1 (FILE *fout, double *f, double (* fun)(double x[], int n),
    int (* dfun) (double x[], double *f, double dx[], int n),
    int (*testx) (double x[], int n),
    double x0[], double space[], double e, int n);

int Newton (FILE *fout, double *f, double (* fun)(double x[], int n),
    int (* ddfun) (double x[], double *fx, double dx[], double ddx[], int n),
    int (*testx) (double x[], int n),
    double x0[], double space[], double e, int n);

int nls2 (FILE *fout, double *sx, double * x0, int nx,
    int (* fun)(double x[], double y[], int nx, int ny),
    int (* jacobi)(double x[], double J[], int nx, int ny),
    int (*testx) (double x[], int nx),
    int ny, double e1, double e2);

/* tree structure functions in treesub.c */
void NodeToBranch (void);
void BranchToNode (void);
void ClearNode (int inode);
int ReadaTreeN (FILE *ftree, int *hasbranch, int popline);
int ReadaTreeB (FILE *ftree, int popline);
int OutaTreeN (FILE *fout, int spnames, int branchlen);
int OutaTreeB (FILE *fout);
void PointLklnodes (void);
int SetBranch (double v[], int correct);
int DistanceMat (FILE *fout, int ischeme, double alfa, double *kapa);
int LSDistance (double * ss, double x[], int (*testx)(double x[],int np));

int StepwiseAdditionMP (double space[]);
double MPScoreStepwiseAddition (int is, double space[], int save);
int TreeSearch (FILE *fout, double space[]);
int AddSpecies (int species, int ib);
int StepwiseAddition (FILE* fout, double space[]);
double TreeScore(double x[], double space[]);

int PopEmptyLines (FILE* fseq, int lline, char line[]);
int hasbase (char *str);
int blankline (char *str);

void BranchLengthBD(int rooted, double birth, double death, double sample, 
     double mut);
int RandomLHistory (int rooted, double space[]);

void DescentGroup (int inode);
void BranchPartition (int partition[], int parti2B[]);
int NSameBranch (int partition1[],int partition2[], int nib1,int nib2);
void BranchPartitionLarge (char partition[], int parti2B[]);
int NSameBranchLarge (char partition1[],char partition2[], int nib1,int nib2);

int RootTN93 (int ischeme, double kapa1, double kapa2, double pi[], 
    double *f, double Root[]);
int EigenTN93 (int ischeme, double kapa1, double kapa2, double pi[6],
    int *nR, double Root[], double Cijk[]);

int DownStatesOneNode (int ison, int father);
int DownStates (int inode);
int PathwayMP (FILE *fout, double space[]);
double MPScore (double space[]);
double RemoveMPNinfSites (double *nsiteNinf);

int MPInformSites (void);
int CollapsNode (int inode, double x[]);
void OutSeqs (FILE *fseq, int seqtype, int pose[], int format);
int OutSeqsCodon (FILE* fseqtmp, int *pose);
int OutSeqsMgenes (void);

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

/* functions for evolving sequences */
int GenerateSeq (void);
int Rates4Sites (double rates[],double alfa,int ncatG,int ls, int cdf,
    double space[]);
void Evolve (int inode);
void EvolveJC (int inode);

extern char NUCs[], AAs[];
#define  NUCseq      0
#define  CODONseq    1
#define  AAseq       2
#define  CODON2AAseq 3

/*
#define FAST_RANDOM_NUMBER
*/

/*
#define DEBUG 9
*/

#endif
