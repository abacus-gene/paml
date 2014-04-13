/* yn00.c
   Pairwise estimation of dS and dN by the method of Yang & Nielsen 
   (2000 Mol. Biol. Evol. 17:32-43)

   Copyright, 1998, Ziheng Yang

                 cc -o yn00 -fast yn00.c tools.o -lm
                     yn00 <SequenceFileName>

  Codon sequences are encoded as 0,1,...,63, different from codeml.c.
  Routines in the two files should not be mixed.
*/
#include "tools.h"
#define NS            500
#define LSPNAME       30
#define NCODE         64
#define NGENE         1

int GetOptions (char *ctlf);
int ReadSeq1 (FILE *fout, FILE *fseq);
int EncodeSeq(void);
int Statistics(FILE *fout, double space[]);
int DistanceMatNG86 (FILE *fout, double alpha);
int DistanceYN00(int is, int js, double*S, double*N, double*dS,double*dN,
    double *t,double space[]);
int GetKappa (void);
int GetFreqs(int is1, int is2, double f3x4[], double pi[]);
int CountSites(char z[],double pi[],double*Stot,double*Ntot,
    double fbS[],double fbN[]);
int GetPMatCodon(double P[],double t, double kappa, double omega, double pi[],
    double e, double space[]);
int CountDiffs(char z1[],char z2[], 
               double*Sdts,double*Sdtv,double*Ndts,double*Ndtv,double PMat[]);
int DistanceF84(double n,double P,double Q,double pi[],double*k_HKY, double*t);
double dsdnREV (int is, int js, double space[]);

int ExpPattFreq(double t,double kappa,double omega,double pi[],double space[]);
int ConsistencyMC(void);
int InfiniteData(double t,double kappa,double omega,double f3x4_0[],
    double space[]);
void SimulateData2s64(FILE* fout, double f3x4_0[], double space[]);

struct common_info {
   char *z[NS], *spname[NS], seqf[96],outf[96];
   int ns,ls,npatt,codonf,icode,ncode,getSE,*pose,verbose, seqtype, cleandata;
   int fcommon,kcommon, iteration;
   double *fpatt, pi[NCODE], f3x4s[NS][12], kappa, omega;
   int ngene,posG[NGENE+1],lgene[NGENE],fix_rgene, model;
   double rgene[NGENE],piG[NGENE][NCODE], alpha;
}  com;


#define REALSEQUENCE
#include "treesub.c"


double PMat[NCODE*NCODE];
char *codonfreqs[]={"Fequal", "F1x4", "F3x4", "Fcodon"};
enum {Fequal, F1x4, F3x4, Fcodon} CodonFreqs;

FILE *frst, *frst1, *frub;
extern char BASEs[], AAs[], Nsensecodon[];
extern int GenetCode[][64];
extern int noisy;
enum {NODEBUG, KAPPA, SITES, DIFF} DebugFunctions;
int debug=0;

double omega_NG, dN_NG, dS_NG;


int main(int argc, char *argv[])
{
   char dsf[32]="2YN.dS",dnf[32]="2YN.dN",tf[32]="2YN.t";
   FILE *fout, *fseq, *fds,*fdn,*ft;
   char ctlf[96]="yn00.ctl";
   int n=com.ncode, is,js, j,k, wname=20, sspace;
   double t=0.4, dS=0.1,dN=0.1, S,N, f3x4[12], *space;

   /* ConsistencyMC(); */

   if (argc>1)  strcpy(ctlf, argv[1]); 
   com.seqtype=1;  com.cleandata=1;
   noisy=1; com.icode=0;  com.fcommon=0;  com.kcommon=1;
   GetOptions (ctlf);
   fout=fopen (com.outf,"w"); frst=fopen("rst","w"); frst1=fopen("rst1","w"); 
   frub=fopen ("rub","w");
   if (fout==NULL || frst==NULL) error2("outfile creation err.");
   fds=(FILE*)fopen(dsf,"w"); fdn=(FILE*)fopen(dnf,"w"); ft=(FILE*)fopen(tf,"w"); 
   if(fds==NULL || fdn==NULL || ft==NULL) error2("file open error");

   if((fseq=fopen (com.seqf,"r"))==NULL) {
      printf ("\n\nSequence file %s not found!\n", com.seqf);
      exit(-1);
   }
   ReadSeq ((com.verbose?fout:NULL), fseq);
   sspace=max2(200000,64*com.ns*sizeof(double));
   sspace=max2(sspace,64*64*5*sizeof(double));
   if ((space=(double*)malloc(sspace))==NULL) error2("oom space");

   com.kappa=4.6;  com.omega=1;

   k=max2(64*64,com.ls/3)*sizeof(double);
   if ((com.fpatt=(double*)malloc(k))==NULL)  error2("oom fpatt");
   EncodeSeq();   /* ls changed to ls/3 in this */
   com.npatt=com.ls;  FOR(k,com.ls) com.fpatt[k]=1;

   fprintf (fout,"YN00 %15s", com.seqf);
   Statistics(fout, space); 
   DistanceMatNG86(fout,0);  fflush(fout);

   if(noisy) puts("\nEstimation by the method of Yang & Nielsen (2000):");
   fputs("\nEstimation by the method of Yang & Nielsen (2000):\n",fout);
   if(!com.iteration) fputs("(equal weighting of pathways)\n",fout);

   if(com.fcommon) GetFreqs(-1,-1,f3x4,com.pi);
   matout(frst,f3x4,3,4);
   if(com.kcommon) {
      GetKappa();
      printf("kappa = %.2f\n\n",com.kappa);
      /* puts("kappa?"); scanf("%lf", &com.kappa); */
   }

   fputs("\nseq. seq.     S       N        t   kappa   omega     dN     dS\n\n",fout);
   fprintf(fds,"%6d\n", com.ns);  fprintf(fdn,"%6d\n", com.ns);
   fprintf(ft,"%6d\n", com.ns);
   FOR(is,com.ns) {
      if(noisy) printf (" %3d", is+1);
      fprintf(fds,"%-*s ", wname,com.spname[is]);
      fprintf(fdn,"%-*s ", wname,com.spname[is]);
      fprintf(ft,"%-*s ", wname,com.spname[is]);
      FOR (js,is) {
         if(noisy>3) printf(" vs. %3d", js+1);
         fprintf (fout, " %3d  %3d ", is+1,js+1);

         if(!com.fcommon) GetFreqs(is,js,f3x4,com.pi);
         if(!com.kcommon) GetKappa();
         j=DistanceYN00(is, js, &S, &N, &dS,&dN, &t,space);

         fprintf(fout,"%7.1f %7.1f %8.4f%8.4f%8.4f%7.4f%7.4f\n",S,N,t,com.kappa,com.omega,dN,dS);
	 /*
         fprintf(frst,"%8.4f%8.4f%8.4f%7.4f%7.4f\n",t,com.kappa,com.omega,dN,dS);
	 */
         fprintf(fds," %7.4f",dS); fprintf(fdn," %7.4f",dN); fprintf(ft," %7.4f",t);
      }    /* for (js) */
      FPN(fds); FPN(fdn); FPN(ft);   fflush(fds); fflush(fdn); fflush(ft);
   }       /* for (is) */
   fflush(frst);
   FPN(F0);
   return (0);
}



int GetOptions (char *ctlf)
{
   int i, nopt=7, lline=255; 
   char line[255], *pline, opt[20], comment='*';
   char *optstr[]={"seqfile","outfile", "verbose", "icode", 
        "weighting","commonkappa", "commonf3x4"};
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
      if ((pline=strstr(line, "="))==NULL) error2("option file.");

      for (i=0; i<nopt; i++) {
         if (strncmp(opt, optstr[i], 8)==0)  {
            if (noisy>2)
               printf ("\n%3d %15s | %-20s %6.2f", i+1,optstr[i],opt,t);
            switch (i) {
               case (0): sscanf(pline+2, "%s", com.seqf);    break;
               case (1): sscanf(pline+2, "%s", com.outf);    break;
               case (2): com.verbose=(int)t;     break;
               case (3): com.icode=(int)t;       break;
               case (4): com.iteration=(int)t;   break;
               case (5): com.kcommon=(int)t;     break;
               case (6): com.fcommon=(int)t;     break;
            }
            break;
         }
      }
      if (i==nopt)
         { printf ("\noption %s in %s\n", opt, ctlf);  exit (-1); }
   }
   com.ncode=Nsensecodon[com.icode];
   fclose (fctl);
   return (0);
}

int DistanceYN00(int is, int js, double*S, double*N, double*dS,double*dN,
    double *t,double space[])
{
/* calculates dS, dN, w, t by iteration.
   com.kappa & com.pi[] are calculated beforehand are not updated.
*/
   int j,k,ir,nround=10, status=0;
   double fbS[4],fbN[4],fbSt[4],fbNt[4], St,Nt, Sdts,Sdtv,Ndts,Ndtv,y;
   double w0=0,dS0=0,dN0=0, accu=5e-4,minomega=1e-5;

   if(*t==0) *t=1;  if(com.omega<=0) com.omega=1;
   FOR(k,4) fbS[k]=fbN[k]=0;
   if(debug) printf("\nCountSites\n");
   for(k=0,*S=*N=0; k<2; k++) {
      CountSites(com.z[k==0?is:js], com.pi, &St, &Nt, fbSt, fbNt);
      *S+=St/2; *N+=Nt/2; 
      FOR(j,4) { fbS[j]+=fbSt[j]/2; fbN[j]+=fbNt[j]/2; }
      if(noisy>3) printf("Seq. %d: S = %9.3f N=%9.3f\n",k+1,St,Nt);
   }
   if(noisy>3) {
      printf("Ave.  : S = %9.3f N=%9.3f\n\n",*S,*N);
      printf("Base freqs at syn & nonsyn sites");
      matout(F0,fbS,1,4);  matout(F0,fbN,1,4);
   }
   if(noisy>3) printf("\n      %-34s t k w dN dS\n","Diffs");
   /* initial values? */
   for (ir=0; ir<(com.iteration?nround:1); ir++) {   /* iteration */
      if(com.iteration) 
         GetPMatCodon(PMat,*t,com.kappa,com.omega,com.pi,1e-5,space);
      else
         FOR(j,NCODE*NCODE) PMat[j]=1;

      CountDiffs(com.z[is],com.z[js], &Sdts,&Sdtv,&Ndts,&Ndtv, PMat);

      if(DistanceF84(*S, Sdts/ *S, Sdtv/ *S, fbS, &y, dS)) status=-1;
      if(noisy>3) printf("   syn kappa=%8.4f\n",y);
      if(DistanceF84(*N, Ndts/ *N, Ndtv/ *N, fbN, &y, dN)) status=-1;
      if(noisy>3) printf("nonsyn kappa=%8.4f\n",y);

      if(*dS<1e-9) { status=-1; com.omega=1; }
      else         com.omega= max2(minomega, *dN/ *dS);
      *t = *dS * 3 * *S/(*S + *N) + *dN * 3 * *N/(*S + *N);
      if(noisy>3) {
         printf("%2d %7.2f%7.2f %7.2f%7.2f |",ir+1, Sdts,Sdtv,Ndts,Ndtv);
         printf("%8.4f%8.4f%8.4f%8.4f%8.4f\n",*t,com.kappa,com.omega,*dN,*dS);
      }
      if(fabs(*dS-dS0)<accu && fabs(*dN-dN0)<accu && fabs(com.omega-w0)<accu)
         break;
      dS0=*dS;  dN0=*dN;  w0=com.omega;
   } /* for (ir) */
   if(ir==nround) status=-2;
   /* if(status) printf("\n\tstatus: %d\n", status); */
   return(status);
}




int EncodeSeq(void)
{
/* encode sequences into 0, 1, ..., 63
   check for stop codons
*/
   char str[4]="";
   int is,h,j, it, indel=0, c[3];

   FOR (is, com.ns) {
      for (h=0; h<com.ls; h++) {
         for(j=0,it=0;j<3;j++) {
            c[j]=CodeChara(com.z[is][h*3+j],0);
            if(c[j]<0 || c[j]>3) { it=-1; indel=1; }
         }
         if (it)
            printf("Strange character at codon %d seq. %d",h+1,is+1);
         com.z[is][h]=(char)(it=c[0]*16+c[1]*4+c[2]);
         if (GenetCode[com.icode][it]==-1) {
           printf("\a\nStop codon %s at %d seq %d\n",getcodon(str,it),h+1,is+1);
           exit(-1);
	 }
      }
   }
   if (indel) error2("indel?");
   FOR(is,com.ns) com.z[is]=(char*)realloc(com.z[is], (com.ls+1)*sizeof(char));
   return(0);
}

int Statistics(FILE *fout, double space[])
{
/* This calculates base frequencies, using npatt & fpatt[]
*/
   int h, is,j, c[3], wname=20;
   double f3x4tot[12], *fb3tot=com.pi, *fb3s=space;

   if(fout) {
      fprintf(fout, "\n\nns =%4d\tls =%4d", com.ns, com.ls);
      fprintf(fout,"\n\nCodon position x base (3x4) table for each sequence.");
   }
   zero(f3x4tot,12);  zero(fb3s,64*com.ns);
   FOR(is,com.ns)  zero(com.f3x4s[is],3*4);
   for (is=0; is<com.ns; is++) {
      for (h=0; h<com.npatt; h++) {
         j=com.z[is][h]; c[0]=j/16; c[1]=(j%16)/4; c[2]=j%4;
         fb3s[is*64+j]+=com.fpatt[h];
         FOR(j,3) com.f3x4s[is][j*4+c[j]] += com.fpatt[h]/com.ls;
      }
      FOR(j,12) f3x4tot[j]+=com.f3x4s[is][j]/com.ns;
      if(fout) { 
         fprintf(fout,"\n\n%-*s", wname, com.spname[is]);
         FOR (j,3) {
            fprintf (fout, "\nposition %2d:", j+1);
            FOR(h,4) fprintf (fout,"%5c:%7.5f", BASEs[h], com.f3x4s[is][j*4+h]);
	 }
      }
   }
   if(fout) {
      fprintf (fout, "\n\nAverage");
      FOR (j,3) {
         fprintf (fout, "\nposition %2d:", j+1);
         FOR (h,4) fprintf (fout,"%5c:%7.5f", BASEs[h], f3x4tot[j*4+h]);
      }
      for(is=0,zero(fb3tot,64);is<com.ns;is++) 
         FOR(j,64) fb3tot[j]+=fb3s[is*64+j];
      fprintf (fout, "\n\nCodon usage for each species\n");
      printcums (fout, com.ns, fb3s, com.icode);
      fprintf (fout, "\nSums\n");
      printcums (fout, 1, fb3tot, com.icode);
   }

   return(0);
}

int GetFreqs(int is1, int is2, double f3x4[], double pi[])
{
/* uses com.fcommon and com.f3x4s to calculate f3x4[] and pi[].
   Codon frequencies pi[] are calculated from the f3x4 table.
   The calculation is duplicated when com.fcommon=1.
*/
   int j, k;
   double fstop=0;

   if (com.fcommon)
      for(j=0,zero(f3x4,12);j<com.ns;j++)
         FOR(k,12)  f3x4[k]+=com.f3x4s[j][k]/com.ns;
   else 
      FOR(k,12)  f3x4[k]=(com.f3x4s[is1][k]+com.f3x4s[is2][k])/2;

   if (noisy>=9)  matout(F0,f3x4,3,4);
   FOR(j,64) {
      pi[j] = f3x4[j/16] * f3x4[4+(j%16)/4] * f3x4[8+j%4];
      if (GenetCode[com.icode][j]==-1) { fstop+=pi[j]; pi[j]=0; }
   }
   FOR(j,64) pi[j]/=(1-fstop);
   if (fabs(1-sum(pi,64))>1e-6) error2("err GetFreqs()");

   return (0);
}


int DistanceMatNG86 (FILE *fout, double alpha)
{
/* Nei & Gojobori (1986)
   alpha used for nonsynonymous rates only.
   changed coding com.z[][] for this program.
*/
   int i,j, h, wname=15, status=0, ndiff;
   char codon1[4], codon2[4], str[4]="";
   double dS=-1,dN=-1,dN_dS=-1, ns, na, nst, nat, S,N, St,Nt, y;

   if(noisy) puts("\nEstimation by the method of Nei & Gojobori (1986):\n");
   if (fout) fprintf(fout,"\nNei & Gojobori 1986. dN/dS (dN, dS)");
   FOR (i,com.ns) {
      if (fout)  fprintf (fout, "\n%-*s", wname, com.spname[i]);
      FOR (j,i) {
         for (h=0,nst=nat=St=Nt=0; h<com.npatt; h++)  {
            strcpy(codon1, getcodon(str,com.z[i][h]));
            strcpy(codon2, getcodon(str,com.z[j][h]));
            ndiff=difcodonNG (codon1, codon2, &S, &N, &ns, &na, 0, com.icode);

            St+=S*com.fpatt[h];
            Nt+=N*com.fpatt[h];
            nst+=ns*com.fpatt[h];
            nat+=na*com.fpatt[h];
         }
         y=3*com.ls/(St+Nt); St*=y; Nt*=y;

         if (noisy>=9)
           printf("\n%3d%3d:Site%8.2f +%8.2f =%8.2f\tDiff%8.2f +%8.2f =%8.2f",
             i+1,j+1,St,Nt,St+Nt,nst,nat, nst+nat);

         dS = 1-4./3*nst/St;
         dN = 1-4./3*nat/Nt;
         if(dS<=0||dN<=0) { status=-1; if(noisy>=9) puts("\nNG86 distance."); }
         if(dS==1) dS=0;
         else      dS=(dS<0?-1:3./4*(-log(dS)));
         if(dN==1) dN=0;
         else dN=(dN<0?-1:-3./4*log(dN));

         dN_dS=(dS>0?dN/dS:-1);
         if(dS<=0||dN<=0||dN_dS<0) status=-1;
         if(fout) fprintf(fout, "%7.3f(%6.4f %6.4f)", dN_dS, dN, dS);

fprintf (frub, " NG %7.3f %6.4f %6.4f\n", dN_dS, dN, dS);
omega_NG=dN_dS;  dN_NG=dN; dS_NG=dS;
/* for infinite data */
fprintf (frst, " NG %7.3f %6.4f %6.4f\n", dN_dS, dN, dS);

      }
      if(noisy)  printf(" %3d",i+1);
   }
   if(noisy)  FPN(F0);
   if (fout) { FPN(fout);
   if (status) fputs("NOTE: -1 means the method is inapplicable.\n",fout); }
   return (status);
}



int GetKappa(void)
{
/* This calculates mutational transition/transversion rate ratio kappa 
   using 4-fold degenerate sites from pairwise comparisons 
   under HKY85, weighting estimates by the numbers of sites
*/
   int is,js,j,k,h, i1,pos,c[2],aa[2],b[2][3],a,ndeg,by[3]={16,4,1}, status=0;
   double ka[2], F[2][16],S[2],wk[2], t,P,Q,pi[4];
                 /* F&S&wk [0]: non-degenerate; [1]:4-fold;  S:sites */
   double kdefault=(com.kappa>0?com.kappa:(com.icode==1?10:2));
   char str1[4]="",str2[4]="", *metstr[2]={"non-degenerate","4-fold"};

   for(is=0,com.kappa=0;is<com.ns;is++) {
      FOR (js,is) {
         if(noisy>=9) printf ("\n%4d vs. %3d", is+1, js+1);
         FOR(k,2) zero(F[k],16);
         for(h=0; h<com.npatt; h++) {
            c[0]=com.z[is][h]; c[1]=com.z[js][h];
            FOR(k,2) {
               b[k][0]=c[k]/16; b[k][1]=(c[k]%16)/4; b[k][2]=c[k]%4;
               aa[k]=GenetCode[com.icode][c[k]];
            }
/*
printf("\nh=%2d: %s (%c) .. %s (%c) ", h+1,getcodon(str1,c[0]),AAs[aa[0]],
                                           getcodon(str2,c[1]),AAs[aa[1]]);
*/
            /* find non-degenerate sites */
            for(pos=0; pos<3; pos++) {         /* check all positions */
               for(k=0,ndeg=0;k<2;k++) {       /* two codons */
                  FOR(i1,4) {
                     if(i1==b[k][pos]) continue;
                     a=GenetCode[com.icode][c[k]+(i1-b[k][pos])*by[pos]];
                     if(a==aa[k]) break;
                  }
                  if(i1==4) ndeg++;
               }
               if(ndeg==2) {
                  F[0][b[0][pos]*4+b[1][pos]]+=.5*com.fpatt[h];
                  F[0][b[1][pos]*4+b[0][pos]]+=.5*com.fpatt[h];
               }

            }
            /* find 4-fold degenerate sites at 3rd positions */
            for(k=0,ndeg=0;k<2;k++) {       /* two codons */
               for(j=0,i1=c[k]-b[k][2]; j<4; j++) 
                  if(j!=b[k][2] && GenetCode[com.icode][i1+j]!=aa[k]) break;
               if(aa[0]==aa[1] && j==4) ndeg++;
            }
/*
fprintf(frub,"  %c ", (ndeg==2?'+':'-')); fflush(frub);
*/
            if (ndeg<2) continue;
            F[1][b[0][2]*4+b[1][2]]+=.5*com.fpatt[h]; 
            F[1][b[1][2]*4+b[0][2]]+=.5*com.fpatt[h];
         }  /* for (h) */
         FOR (k,2) {  /* two kinds of sites */
            if(noisy>3) printf("\n%s:\n",metstr[k]);
            S[k]=sum(F[k],16); 
            if(S[k]<=0) { wk[k]=0; continue; }
            FOR(j,16) F[k][j]/=S[k];
            P=(F[k][0*4+1]+F[k][2*4+3])*2;
            Q=1-(F[k][0*4+0]+F[k][1*4+1]+F[k][2*4+2]+F[k][3*4+3]) - P;
            FOR(j,4) pi[j]=sum(F[k]+j*4,4);
            DistanceF84(S[k], P,Q,pi, &ka[k], &t);
            wk[k]=(ka[k]>0?S[k]:0);

            /* matout(F0,F[k],4,4);  matout(F0,pi,1,4);  */
            if(noisy>3)
               printf("\nSPQkt:%9.4f%9.5f%9.5f%9.4f%9.4f\n",S[k],P,Q,ka[k],t);
         }
         if(wk[0]+wk[1]==0) {
            status=-1;  ka[0]=kdefault;
            if(noisy>3) printf("\ngot no kappa! fix it at %.4f\n",ka[0]);
         }
         else
             ka[0]=(ka[0]*wk[0]+ka[1]*wk[1])/(wk[0]+wk[1]);
         com.kappa+=ka[0]/(com.ns*(com.ns-1.)/2);
      }  /* for(js) */
   }     /* for(is) */

   return (status);
}


int CountSites(char z[],double pi[],double*Stot,double*Ntot,double fbS[],double fbN[])
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
   int h, j,k, c[2],aa[2], b[3], by[3]={16,4,1};
   double r, S,N, kappa=com.kappa;

   *Stot=*Ntot=0;  FOR(k,4) fbS[k]=fbN[k]=0;
   for (h=0; h<com.npatt; h++)  {
      c[0]=z[h]; b[0]=c[0]/16; b[1]=(c[0]%16)/4; b[2]=c[0]%4;
      aa[0]=GenetCode[com.icode][c[0]];
      if (aa[0]==-1) error2("stop codon");
      for (j=0,S=N=0; j<3; j++) {
         FOR (k,4) {    /* b[j] changes to k */
            if (k==b[j]) continue;
            c[1]=c[0]+(k-b[j])*by[j];
            aa[1]=GenetCode[com.icode][c[1]];
            if(aa[1]==-1) continue;
            r=pi[c[1]];
            if (k+b[j]==1 || k+b[j]==5) r*=kappa; /* transition */
            if (aa[0]==aa[1]) { S+=r; fbS[b[j]]+=r*com.fpatt[h]; }
            else              { N+=r; fbN[b[j]]+=r*com.fpatt[h]; }
         }
      }
      *Stot+=com.fpatt[h]*S;
      *Ntot+=com.fpatt[h]*N;
   }
   r=3*com.ls/(*Stot+*Ntot);  *Stot*=r;  *Ntot*=r;
   r=sum(fbS,4); FOR(k,4)fbS[k]/=r;
   r=sum(fbN,4); FOR(k,4)fbN[k]/=r;
   return (0);
}


int GetPMatCodon(double P[],double t, double kappa, double omega, double pi[],
    double e, double space[])
{
/* Get approximate PMat using polynomial approximation
   P(t) = I + Qt + (Qt)^2/2 + (Qt)^3/3!
*/
   int nterms=100, n=NCODE, i,j,k, c[2],ndiff,pos=0,from[3],to[3];
   double *Q=space, *Q1=Q+n*n, *Q2=Q1+n*n, mr, div;

   FOR (i,n*n) Q[i]=0;
   for (i=0; i<n; i++) FOR (j,i) {
      from[0]=i/16; from[1]=(i/4)%4; from[2]=i%4;
      to[0]=j/16;   to[1]=(j/4)%4;   to[2]=j%4;
      c[0]=GenetCode[com.icode][i];
      c[1]=GenetCode[com.icode][j];
      if (c[0]==-1 || c[1]==-1)  continue;
      for (k=0,ndiff=0; k<3; k++)  if (from[k]!=to[k]) { ndiff++; pos=k; }
      if (ndiff!=1)  continue;
      Q[i*n+j]=1;
      if ((from[pos]+to[pos]-1)*(from[pos]+to[pos]-5)==0)  Q[i*n+j]*=kappa;
      if(c[0]!=c[1])  Q[i*n+j]*=omega;
      Q[j*n+i]=Q[i*n+j];
   }
   FOR(i,n) FOR(j,n) Q[i*n+j]*=pi[j];
   for (i=0,mr=0;i<n;i++) { 
      Q[i*n+i]=-sum(Q+i*n,n); mr-=pi[i]*Q[i*n+i]; 
   }
   FOR(i,n*n) Q[i]*=t/mr;

   xtoy(Q,P,n*n);  FOR(i,n) P[i*n+i]++;   /* I + Qt */
   xtoy(Q,Q1,n*n);
   for (i=2,div=2;i<nterms;i++,div*=i) {  /* k is divisor */
      matby(Q, Q1, Q2, n, n, n);
      for(j=0,mr=0;j<n*n;j++) { P[j]+=Q2[j]/div; mr+=fabs(Q2[j]); }
      mr/=div;
      /* if(debug) printf("Pmat term %d: norm=%.6e\n", i, mr); */
      if (mr<e) break;
      xtoy(Q2,Q1,n*n);
   }

/*
testTransP(PMat, n);
fprintf(frub,"\nP(%.5f)", t);
matout (frub,PMat,NCODE,NCODE); fflush(frub);
*/
   FOR(i,n*n) if(P[i]<0) P[i]=0;

   return (0);
}



int CountDiffs(char z1[],char z2[], 
    double*Sdts,double*Sdtv,double*Ndts,double*Ndtv,double PMat[])
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
   char str[4]="";
   int h,i1,i2,i,k, transi, c[2],ct[2],aa[2], by[3]={16,4,1};
   int dmark[3], step[3], b[2][3], bt1[3], bt2[3];
   int ndiff, npath, nstop, stspath[6],stvpath[6],ntspath[6],ntvpath[6];
   double sts,stv,nts,ntv; /* syn ts & tv, nonsyn ts & tv for 2 codons */
   double ppath[6], sump,p;

   for (h=0,*Sdts=*Sdtv=*Ndts=*Ndtv=0; h<com.npatt; h++)  {
      c[0]=z1[h]; c[1]=z2[h];
      if (c[0]==c[1]) continue;
      FOR(i,2) {
         b[i][0]=c[i]/16; b[i][1]=(c[i]%16)/4; b[i][2]=c[i]%4;
         aa[i]=GenetCode[com.icode][c[i]];
      }
      if (aa[0]==-1||aa[1]==-1)
         error2("stop codon in sequence.");
      ndiff=0;  sts=stv=nts=ntv=0;
      FOR (k,3) dmark[k]=-1;
      FOR (k,3) if (b[0][k]!=b[1][k]) dmark[ndiff++]=k;
      npath=1;  if(ndiff>1) npath=(ndiff==2)?2:6;
      if (ndiff==1) {
         transi=b[0][dmark[0]]+b[1][dmark[0]];
         transi=(transi==1 || transi==5);
         if (aa[0]==aa[1])  { if (transi) sts++; else stv++; }
         else               { if (transi) nts++; else ntv++; }
      }
      else {   /* ndiff=2 or 3 */

if(debug==DIFF) {
   printf("\n\nh=%d %s (%c) .. ", h+1,getcodon(str,c[0]),AAs[aa[0]]);
   printf("%s (%c): ", getcodon(str,c[1]),AAs[aa[1]]);
}
         nstop=0;
         FOR(k,npath) {

if(debug==DIFF) printf("\npath %d: ", k+1);

            FOR (i1,3)  step[i1]=-1;
            if (ndiff==2) 
               { step[0]=dmark[k];   step[1]=dmark[1-k];  }
            else {
               step[0]=k/2;   step[1]=k%2;
               if (step[0]<=step[1]) step[1]++;
               step[2]=3-step[0]-step[1];
            }
            FOR(i1,3) bt1[i1]=bt2[i1]=b[0][i1];
            stspath[k]=stvpath[k]=ntspath[k]=ntvpath[k]=0;  
            /* mutations along each path */
            for (i1=0,ppath[k]=1; i1<ndiff; i1++) { 
               bt2[step[i1]] = b[1][step[i1]];
               for (i2=0,ct[0]=ct[1]=0; i2<3; i2++) {
                  ct[0]+=bt1[i2]*by[i2];
                  ct[1]+=bt2[i2]*by[i2];
               }
               ppath[k]*=PMat[ct[0]*64+ct[1]];
               FOR(i2,2) aa[i2]=GenetCode[com.icode][ct[i2]];

if(debug==DIFF) printf("%s (%c) %.5f: ", getcodon(str,ct[1]),AAs[aa[1]],PMat[ct[0]*64+ct[1]]);

               if (aa[1]==-1)       { nstop++;  ppath[k]=0; break; }
               transi=b[0][step[i1]]+b[1][step[i1]];
               transi=(transi==1 || transi==5);  /* transition? */

               if(aa[0]==aa[1]) { if(transi) stspath[k]++; else stvpath[k]++; }
               else             { if(transi) ntspath[k]++; else ntvpath[k]++; }
               FOR(i2,3) bt1[i2]=bt2[i2];
            }

if(debug==DIFF) printf("  p =%.9f", ppath[k]);

         }  /* for(k,npath) */
         if (npath==nstop) {  /* all paths through stop codons */
            puts ("all paths through stop codons..");
            if (ndiff==2) { nts=.5; ntv=1.5; }
            else          { nts=.5; ntv=2.5; }
         }
         else {
            sump=sum(ppath,npath);
            if(sump<1e-20) { 
               printf("\nsump=0, npath=%4d\nh=%d (%d,%d)\np:",
                    npath,h+1,c[0],c[1]);
               matout2(F0,ppath,1,npath,16,12);
               matout(frub,PMat,NCODE,NCODE); fflush(frub);

               getchar(); exit(-1);

               sump=1; FOR(k,npath) if(ppath[k]) ppath[k]=1./(npath-nstop);
            }
            for(k=0;k<npath;k++) { 
               p=ppath[k]/sump;
               sts+=stspath[k]*p; stv+=stvpath[k]*p;  
               nts+=ntspath[k]*p; ntv+=ntvpath[k]*p;
            }

if(debug==DIFF) {
   FOR(k,npath) printf("\n p =%.5f", ppath[k]/sump);  FPN(F0);
   printf(" syn ts & tv, nonsyn ts & tv:%9.5f%9.5f%9.5f%9.5f\n",sts,stv,nts,ntv);
}
         }

if(debug==DIFF) getchar();

      }     /* if (ndiff) */
      *Sdts+=com.fpatt[h]*sts;
      *Sdtv+=com.fpatt[h]*stv;
      *Ndts+=com.fpatt[h]*nts;
      *Ndtv+=com.fpatt[h]*ntv;
   }  /* for (h) */
   return (0);
}


int DistanceF84(double n,double P,double Q,double pi[],double*k_HKY, double*t)
{
/* This calculates kappa and d from P & Q & pi under F84.
   When F84 fails, we try to use K80.  When K80 fails, we try
   to use JC69.  When JC69 fails, we set distance t to maxt.

*/
   int failF84=0,failK80=0,failJC69=0;
   double tc,ag, Y,R, a=0,b=0, k_F84;
   double Qsmall=min2(1e-10,0.1/n), maxkappa=999,maxt=99;

   *k_HKY=-1;
   Y=pi[0]+pi[1];  R=pi[2]+pi[3];  tc=pi[0]*pi[1];  ag=pi[2]*pi[3];
   if (P+Q>1) { *t=maxt; *k_HKY=1; return(3); }
   if (P<-1e-10 || Q<-1e-10 || fabs(Y+R-1)>1e-8) {
      printf("\nPQYR & pi[]: %9.5f%9.5f%9.5f%9.5f",P,Q,Y,R);
      matout(F0,pi,1,4);
      error2("DistanceF84: input err.");
   }
   if(Q<Qsmall)          failF84=failK80=1;
   else if(Y<=0 || R<=0 || (tc<=0 && ag<=0)) failF84=1;
   else {
      a=(2*(tc+ag)+2*(tc*R/Y+ag*Y/R)*(1-Q/(2*Y*R)) - P) / (2*tc/Y+2*ag/R);
      b=1-Q/(2*Y*R);
      if (a<=0 || b<=0) failF84=1;
   }
   if (!failF84) {
      a=-.5*log(a); b=-.5*log(b);
      if(b<=0) failF84=1;
      else {
         k_F84=a/b-1;
         if(t) *t=4*b*(tc*(1+ k_F84/Y)+ag*(1+ k_F84/R)+Y*R);
         *k_HKY=(tc+ag + (tc/Y+ag/R)* k_F84)/(tc+ag); /* k_F84=>k_HKY85 */
      }
   }
   if(failF84 && !failK80) {  /* try K80 */
      if (noisy>=9) printf("\na=%.5f  b=%.5f, use K80\n", a,b);
      a=1-2*P-Q;  b=1-2*Q;
      if (a<=0 || b<=0) failK80=1;
      else {
         a=-log(a); b=-log(b);
         if(b<=0)  failK80=1;
         else {
            *k_HKY=(.5*a-.25*b)/(.25*b);
            if(*t) *t = .5*a+.25*b;
	 }
      }
   }
   if(failK80) {
      if((P+=Q)>=.75) { failJC69=1; P=.75*(n-1.)/n; }
      if(*t) { *t=-.75*log(1-P*4/3.); if(*t>maxt) *t=maxt; }
   }
   if(*k_HKY>maxkappa) *k_HKY=maxkappa;

   return(failF84 + failK80 + failJC69);
}



#if 0

double dsdnREV (int is, int js, double space[])
{
/* This calculates ds and dn by recovering the Q*t matrix using the equation
      F(t) = PI * P(t) = PI * exp(Q*t)
   space[64*64*5]
*/
   int n=NCODE, i,j, h;
   double *F=PMat, *Qt=F;
   double *Root=space+n*n,*pi=Root+n, *U=pi+n,*V=U+n*n;
   double *T1=V+n*n,*T2=T1+n*n, t, small=1e-6;
   
   fprintf(frst,"\npi in model\n");
   matout(frst,com.pi,1,n);
   FOR(i,n*n) F[i]=0;
   FOR (h,com.npatt) {
      F[com.z[is][h]*n+com.z[js][h]]+=com.fpatt[h]/(2*com.ls);
      F[com.z[js][h]*n+com.z[is][h]]+=com.fpatt[h]/(2*com.ls);
   }
   if(fabs(1-sum(F,n*n))>1e-6) error2("Sum F != 1 in dsdnREV");

   FOR (i,n) {
      pi[i]=sum(F+i*n, n);  
/*
      if (F[i*n+i]<=small || F[i*n+i]<pi[i]/4)
*/
      if (F[i*n+i]<=small)  F[i*n+i]=1-pi[i]+F[i*n+i];
      else                  abyx(1/pi[i], F+i*n, n); 
   }
   if (eigen (1, F, n, Root, T1, U, V, T2)) error2 ("eigen jgl");
   xtoy (U, V, n*n);
   matinv (V, n, n, T1);

fprintf(frst,"\npi in data\n");
matout (frst, pi, 1, n);   FPN(F0);
matout (frst, Root, 1, n);

   FOR (i,n) {
      if (Root[i]<=0) 
         printf ("  Root %d:%10.4f", i+1, Root[i]); 
      Root[i]=log(Root[i]);
   }
   FOR (i,n) FOR (j,n) T1[i*n+j]=U[i*n+j]*Root[j];
   matby (T1, V, Qt, n, n, n);
   for (i=0,t=0; i<n; i++) t-=pi[i]*Qt[i*n+i];
   if (t<=0) puts ("err: dsdnREV");

   FOR(i,n*n) Qt[i]+=1e-8;  /* remove negative numbers from rounding errors */

   matout(frst,Qt,n,n);
printf("\nt = %.5f\n", t);

   return (0);
}


int ExpPattFreq(double t,double kappa,double omega,double pi[],double space[])
{
/* This puts site pattern probabilities into com.fpatt[]
*/
   int i,j,h, n=NCODE;
   double y, sum=0;

   GetPMatCodon(PMat, t, kappa, omega, pi, 1e-8, space);
   for(i=0,h=0;i<n;i++) for(j=0;j<=i;j++) {
      if(GenetCode[com.icode][i]==-1 || GenetCode[com.icode][j]==-1) continue;
      com.z[0][h]=(char)i; com.z[1][h]=(char)j;
      y=pi[i]*PMat[i*n+j];
      if(i!=j) y+=pi[j]*PMat[j*n+i];
      com.fpatt[h++]=y;
      sum+=y;

      y=pi[i]*PMat[i*n+j]-pi[j]*PMat[j*n+i];
      if (fabs(y)>1e-5) { 
         printf("F(%.4f)[%d, %d] not sym %.6f?\n", t,i,j,y);
         getchar();
      }

   }
   if(h!=com.ncode*(com.ncode+1)/2) error2("ExpPattFreq: npatt");
   if(fabs(sum-1)>1e-5) printf("\nExpPattFreq: sum=1=%.6f?",sum);
   return(0);
}



int InfiniteData(double t,double kappa,double omega,double f3x4_0[],double space[])
{
/* This generates site pattern freqs, data when seq length is infinite
*/
   int j, n=com.ncode;
   double s, pi0[NCODE];

   for(j=0,s=0; j<64; j++) {
      pi0[j]=f3x4_0[j/16]*f3x4_0[4+(j%16)/4]*f3x4_0[8+j%4];
      if (GenetCode[com.icode][j]==-1) pi0[j]=0;
      s+=pi0[j];
   }
   FOR(j,64) pi0[j]=com.pi[j]=pi0[j]/s;

   printf("\nsum pi=1=%.6f\n",sum(pi0,NCODE));
   printf("\nt=%.3f\tkappa=%.3f\tomega=%.3f\n",t,kappa,omega);

   com.ns=2;
   com.ls=1; com.npatt=n*(n+1)/2;
   FOR(j,com.ns) sprintf(com.spname[j],"seq.%d",j+1);
   FOR(j,com.ns) if(com.z[j]) free(com.z[j]);
   FOR(j,com.ns) com.z[j]=(char*) malloc(com.npatt*sizeof(char));
   if (com.z[com.ns-1]==NULL) error2 ("oom");
   if (com.fpatt) free(com.fpatt);
   if ((com.fpatt=(double*)malloc(com.npatt*sizeof(double)))==NULL)  
      error2("oom fpatt");

   ExpPattFreq(t, kappa, omega, pi0, space);
   return(0);
}

int ConsistencyMC(void)
{
/* This does consistency analysis and computer simulation for 
   the YN00 method.  Needs testing since untested chagnes were made when 
   preparing a release in the paml package.
*/
   FILE *fout;
   char outf[64]="yn";
   int n=com.ncode, j, wname=30, sspace, realseq,fcodon;
   double t=0.4, dS=0.1,dN=0.1, S,N, f3x4[12], *space;
   double f3x4_data[][3*4]={
                            {0.25, 0.25, 0.25, 0.25,
                             0.25, 0.25, 0.25, 0.25,
                             0.25, 0.25, 0.25, 0.25},

                            {0.20517, 0.28293, 0.30784, 0.20406, /* mt pri */
                             0.40979, 0.27911, 0.18995, 0.12116,
                             0.15105, 0.43290, 0.37123, 0.04482},

                            {0.19020, 0.16201, 0.36655, 0.28124, /* hiv */
                             0.28889, 0.18805, 0.30179, 0.22127,
                             0.24875, 0.16894, 0.36822, 0.21410},

                            {0.14568, 0.24519, 0.33827, 0.27086,
                             0.35556, 0.18765, 0.24049, 0.21630, 
                             0.26444, 0.25728, 0.21012, 0.26815}, /* lysinNew*/
                           }; 

   noisy=1;        com.ns=2;
   com.icode=0;    com.ncode=Nsensecodon[com.icode];
   realseq=0;      fcodon=1;

   printf("\nrandom number seed & fcodon? ");
   scanf("%d%d", &j,&fcodon);
   SetSeed(j);
   printf("\nicode: %d fcommon: %d kcommon: %d & fcodon: %d\n",
      com.icode,com.fcommon,com.kcommon,fcodon);
   matout(F0,f3x4_data[fcodon],3,4);
   fout=fopen (outf,"w"); frst=fopen("rst","w"); frst1=fopen("rst1","w"); 
   frub=fopen ("rub","w");
   if (fout==NULL || frst==NULL) error2("outfile creation err.");

   sspace=max2(200000,64*com.ns*sizeof(double));
   sspace=max2(sspace,64*64*5*sizeof(double));
   if ((space=(double*)malloc(sspace))==NULL) error2("oom space");

   com.kappa=4.6;  com.omega=1;

/* SimulateData2s64(fout, f3x4_data[fcodon], space); */

   fputs("\nConsistency analysis: YN00\n",fout);
   if(!com.iteration) fputs("(equal weighting of pathways)\n",fout);

   for (; ; ) {
      printf("\nInfinite data\nt k w? ");
      scanf("%lf%lf%lf",&t,&com.kappa,&com.omega);
      if(t<=0 || com.kappa<=0 || com.omega<=0) break;
      InfiniteData(t, com.kappa, com.omega, f3x4_data[fcodon], space);
      fprintf(frst,"%7.3f%7.3f%7.3f%  ",t,com.kappa,com.omega);
      if(noisy) matout(F0,f3x4_data[fcodon],3,4);
      Statistics(fout, space); 
      GetKappa ();
      printf("\nEstimated kappa = %.2f\n\n", com.kappa);

      fputs("\n    S       N        t   kappa   omega     dN     dS\n\n",fout);
      GetFreqs(0,1,f3x4,com.pi);
      j=DistanceYN00(0,1, &S, &N, &dS,&dN, &t,space);
      fprintf(fout,"%7.1f %7.1f %8.4f%8.4f%8.4f%7.4f%7.4f\n",S,N,t,com.kappa,com.omega,dN,dS);
      fprintf(frst,"%8.4f%8.4f%8.4f%7.4f%7.4f\n",t,com.kappa,com.omega,dN,dS);
      fflush(frst);
   }
   return (0);
}

void SimulateData2s64(FILE* fout, double f3x4_0[], double space[])
{
/* This simulates two codon sequences and analyze using NG86 and YN00.
   It generates site pattern freqs and then samples from them to generate
   the seq data.  It is modified from a routine of a similar name from
   codeml.c, and differs from it in that codons are coded as 0,1,...,63.
   This should not be mixed with any thing in codeml.c
   f3x4 table is used to generate codon frequencies.
   x[0-2]: w,dN,dS for NG86; x[3-5]:w,dN,dS for YN00;
*/
   int il,ir,nr=2000, i,j,h, n=Nsensecodon[com.icode];
   int npatt0=n*(n+1)/2, nobs[NCODE*NCODE], nfail, nfail2;
   int ip,nil=1, ls[]={500, 100,500};
   double y, t=.5, x[6], S,N,dS=0,dN=0, f3x4[12];
   double t0,kappa0,omega0,pi0[NCODE],Efij[NCODE*NCODE], mx[6],vx[6],mse[2];
   char str[4]="";
   FILE*ftmp=NULL;

   puts("\nThis simulates data and analyses them using NG86 & YN00.");

   for(j=0,y=0; j<64; j++) {
      pi0[j]=f3x4_0[j/16]*f3x4_0[4+(j%16)/4]*f3x4_0[8+j%4];
      if (GenetCode[com.icode][j]==-1) pi0[j]=0;
      y+=pi0[j];
   }
   FOR(j,64) pi0[j]/=y;

   matout(frst,f3x4_0,3,4);
   com.icode=0; com.ns=2; com.ls=1;
   FOR(j,com.ns) com.z[j]=(char*) malloc(npatt0*sizeof(char));
   if(com.z[com.ns-1]==NULL) error2 ("oom");
   if((com.fpatt=(double*)malloc(npatt0*sizeof(double)))==NULL)
      error2("oom fpatt");

   printf("\nSite pattern probabilities.\nsum pi=1=%.6f\n",sum(pi0,NCODE));
   printf("%d categories in the multinomial.\n",npatt0);

   for(ip=0; ip<99; ip++) {
      t0=1; kappa0=1; omega0=1; FOR(i,6) x[i]=0;
/* com.kappa=-999; */ com.omega=-999; 
      printf("\nt0 & kappa0 & omega0? ");
      scanf("%lf%lf%lf",&t0,&kappa0,&omega0);
      printf("\nt & k & w: %8.2f%8.2f%8.2f\n", t0,kappa0,omega0);
      if(t0<=0 || kappa0<=0 || omega0<=0) exit(0);
      fprintf(frst, "\nt & k & w: %8.2f%8.2f%8.2f\n\n", t0,kappa0,omega0);
   
      fputs("Lc & t & k & w & dS & dN\n",frst);
      for(il=0; il<nil; il++,FPN(F0)) {
         ExpPattFreq(t0, kappa0, omega0, pi0, space);
         xtoy(com.fpatt,Efij,npatt0);
         for(h=1;h<npatt0; h++) Efij[h]+=Efij[h-1];
         if(fabs(1-Efij[npatt0-1])>1e-6) puts("Sum p_ij != 1.");
   
         com.ls=ls[il];
         printf("\nls = %d\n", com.ls);
         zero(mx,6); zero(vx,6); mse[0]=mse[1]=0;
         for(ir=0,nfail=nfail2=0; ir<nr; ir++){
            if((ir+1)%10==0) printf("\r%3d(%2d %2d)", ir+1,nfail,nfail2);
            if(noisy) printf("\n\nReplicate %3d\n", ir+1);
            MultiNomial (com.ls, npatt0, Efij, nobs, NULL);
   
            for(i=h=com.npatt=0; i<NCODE; i++) for(j=0;j<=i;j++) {
               if(GenetCode[com.icode][i]==-1 || GenetCode[com.icode][j]==-1) 
                  continue;
               if(nobs[h++]==0) continue;
               com.z[0][com.npatt]=(char)i; com.z[1][com.npatt]=(char)j;
               com.fpatt[com.npatt++]=nobs[h-1];
            }
            if(h!=npatt0) error2("h!=npatt0");
            j=DistanceMatNG86 (NULL, 0);
            if(j==-1 || omega_NG<1e-6 || dS_NG<1e-6 || dN_NG<1e-6) 
               { ir--; nfail++; continue; }
            x[0]=omega_NG; x[1]=dN_NG; x[2]=dS_NG; 

            if (GetKappa())  { ir--; nfail++; continue; }
            if(noisy) printf("\nEstimated kappa = %.2f\n\n", com.kappa);
            Statistics(NULL, space);
            GetFreqs(0,1,f3x4,com.pi);
            t=(.25*dS_NG+.75*dN_NG)*3;
            j=DistanceYN00(0, 1, &S, &N, &dS,&dN, &t,space);
            if(j==2) nfail2++;
            if(dS<1e-6 || dN<-1e-6) { ir--; continue; }
            x[3]=dN/dS; x[4]=dN; x[5]=dS;

if(noisy &&(x[0]<omega0/4 || x[0]>omega0*4 || x[3]<omega0/4 || x[3]>omega0*4))
  printf("w_NG=%.5f w_YN=%.5f\n",x[0],x[3]);

#if 0
            if((ftmp=fopen("codonseq.tmp","w"))==NULL) error2("seq file err");
            fprintf(ftmp," %6d %6d\n", com.ns,com.ls*3);
            for(i=0;i<2;i++,FPN(ftmp),fflush(ftmp)) {
               fprintf(ftmp,"seq.%-5d  ", i+1);
               FOR(h,com.npatt) FOR(j,(int)com.fpatt[h]) 
                  fprintf(ftmp,"%s", getcodon(str,com.z[i][h]));
            }
            fclose(ftmp);  exit(0);
#endif
            FOR(j,6) {
               vx[j]+=(x[j]-mx[j])*(x[j]-mx[j]);
               mx[j]=(mx[j]*ir+x[j])/(ir+1.);
            }
            mse[0]+=square(x[0]-omega0); mse[1]+=square(x[3]-omega0);
            if((ir+1)%10==0) {
               /* FOR(i,6) printf("%6.3f",x[i]);  printf(" "); */
               FOR(i,6) { printf("%6.3f",mx[i]); if(i==2) printf(" "); }
            }
         }  /* for(ir) */
         if(nr>1) {
            FOR(j,6) vx[j]=sqrt(vx[j]/(nr-1.)/nr); 
            FOR(j,2) mse[j]=sqrt(mse[j]/nr);
         }
         fprintf(frst,"%4d (%3d %2d)", com.ls,nfail,nfail2);
         FOR(i,6) fprintf(frst,"%7.3f +%6.3f", mx[i],vx[i]);
         FPN(frst); fflush(frst);
         fprintf(frst1,"%6d%6d %7.2f%7.2f%7.2f: %8.3f +%7.3f %8.3f +%7.3f\n",
             com.ls,nr, t0,kappa0,omega0, mx[0],vx[0],mx[3],vx[3]);
         fflush(frst1);
      }  /* for(il) */
   }     /* for(ip) */
   exit(0);
}

#endif
