/* descriptive statistics

   cc -o ds -O3 ds.c tools.c -lm
   cl -Ox -W3 -D_CRT_SECURE_NO_WARNINGS ds.c tools.c
   ds in.txt
*/
#include "paml.h"

int DescriptiveStatistics(FILE* fout, char infile[], int nbin, int propternary, int SkipColumns);

int main(int argc, char *argv[])
{
   FILE *fp;
   char infile[96]="in.txt", outfile[96]="out.txt";
   int propternary=0;

   if(argc>1) strcpy(infile, argv[1]);
   if(argc>2) strcpy(outfile,argv[2]);
   fp = (FILE*)zfopen(outfile,"w");
   puts("results go into out.txt");
   starttimer();
   DescriptiveStatistics(fp, infile, 100, 0, 1);
   exit(0);
}
