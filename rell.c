/* RELL.c
   the RELL method of Kishino and Hasegawa (1989) and Kishino et al. (1990)
   for comparing models or trees

           gcc -o rell -O2 rell.c tools.o -lm 

           rell <lfhfile> <rstfile>


    The format of the lfhfile
       ntree   ls    npattern
       itree

       ipattern  weight   log{f(x)}   r   ***

*/


#include "tools.h"
#define NTREE 200

int main(int argc, char *argv[])
{
   char infile[32]="lfh", outfile[32]="r", line[250];
   int nboot=5000, h,i,j,k, it, ntree, ls;
   int npatt, weight, mltree;
   double P[NTREE], lfx[NTREE], *lfh, ml[NTREE], varl[NTREE*NTREE], max, y;
   FILE *flfh, *fout;


   if (argc>=2) strcpy (infile, argv[1]);
   if (argc>=3) strcpy (outfile, argv[2]);

   flfh=fopen (infile,"r");   fout=fopen (outfile,"w");
   if (flfh==NULL || fout==NULL) error ("file err");

   printf ("RELL.c %12s%12s\n", infile, outfile);
   fprintf (fout, "RELL.c %12s%12s\n", infile, outfile);

   fscanf (flfh, "%d%d%d", &ntree, &ls, & npatt);
   printf ("\nntree:%5d\tls:%5d\tnpatt:%5d\n", ntree,ls,npatt);
   lfh = (double *) malloc (ntree*ls*sizeof (double));
   if (lfh == NULL)  error ("oom");

   FOR (j,ntree) {
       printf ("\rReading tree # %d\n", j+1);
       fscanf (flfh, "%d", &k);
       for (h=0,npatt=0; h<ls; npatt++) {
           fscanf (flfh, "%d%d%lf", &k, &weight, &y);
           FOR (i, weight) lfh[j*ls+h++]=y;
           fgets (line, 240, flfh);
           printf ("\n%4d =%4d %6d", npatt+1, k, weight);
       }
       FPN (flfh);
   }

   printf ("\n\nBootstrapping..\n");
   zero (P, ntree);
   FOR (i,nboot) {
       for (j=0,zero(lfx,ntree); j<ls; j++) {
            it = (int) (rndu()*(double)ls);
            FOR (k, ntree) lfx[k] += lfh[k*ls+it];
       }
       for (k=1, mltree=0, y=lfx[0]; k<ntree; k++)
            if (y<lfx[k]) { y=lfx[k]; mltree=k; }
       P[mltree] += 1./(double)nboot;
       if (i%10==9) printf ("\rBootstrapping.. #%6d /%6d", i+1, nboot);    
   } 
   printf ("\n");

   variance (lfh,  ls, ntree, ml, varl);
   abyx ((double)ls, ml, ntree);    
   abyx ((double)ls, varl, ntree*ntree);    
   fprintf (fout,"\n\nMean and variance of likelihoods\n\n");
   FOR (i,ntree) fprintf (fout, "%12.4f", ml[i]);

   FPN (fout);   FPN (fout);
   FOR (i, ntree)     {
       for (j=0; j<=i; j++)  fprintf (fout, "%10.4f", varl[i*ntree+j]);
       FPN (fout);
   }

   FOR (i, ntree)  FOR (j,i)
       varl[i*ntree+j]=varl[j*ntree+i]=
        sqrt(varl[i*ntree+i]+varl[j*ntree+j]-2.0*varl[i*ntree+j]);
   FOR (i,ntree)  varl[i*ntree+i]=0;

   for (i=0,mltree=0,max=ml[0]; i<ntree; i++)
       if (ml[i]>max) { max=ml[i]; mltree=i; }
   fprintf (fout,"\n\n%12s%12s%12s%12s%12s\n\n", 
            "tree #", "li", "li-lML", "   SE", "P(RELL)");
   FOR (i, ntree)
       fprintf (fout,"%12d%12.4f%12.4f%12.4f%12.4f\n", i+1, ml[i], 
           ml[i]-max, varl[i*ntree+mltree], P[i]*100);
   printf ("\n");
   return(0);
}
