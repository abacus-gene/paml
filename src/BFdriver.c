/* This drives the computation of the marginal likelihood (bayes factor) calculation
   using bpp and mcmctree.

   cc -o BFdriver -O3 BFdriver.c tools.c -lm
   cl -Ox BFdriver.c tools.c

   BFdriver <controlfilename> <npoints> <scriptname.sh>
   BFdriver mcmctree.ctl 16 tmp.sh
*/
#include "paml.h"

int main (int argc, char*argv[])
{
   int j, npoints=8, ixw, nline=1024;
   char resultsf[96]="betaweights.txt", ctlf[96]="mcmctree.ctl", scriptf[96]="tmp.sh";
   char ctlfi[96], tmpf[96], line[1024], *pline, *s;
   double beta, sign, weight;
   const double *xNI = NULL, *wNI = NULL;  /* Gauss-Legendre quadrature points */
   FILE  *fctl, *fctlb, *fcommand, *fresults;

   puts("Usage:\n\tBFdriver controlfilename npoints\n");
   puts("\tquadrature: log{M} = 0.5 * SUM w_b * E_b(log{f(X)})\n");
   if(argc<2) exit(-1);
   strcpy(ctlf, argv[1]);
   if(argc>2) npoints = atoi(argv[2]);
   if(argc>3) strcpy(scriptf, argv[3]);
   fctl = (FILE*)gfopen(ctlf, "r");
   if( s = strstr(ctlf, ".ctl") ) *s = '\0';

   fresults =  (FILE*)gfopen(resultsf, "w");
   fprintf(fresults, "%s\t%s\t%s\n", "beta", "weight", "ElnfX");
   GaussLegendreRule(&xNI, &wNI, npoints);
   for (j=0; j<npoints; j++) {
      if (j<npoints / 2) { ixw = npoints / 2 - 1 - j;  sign = -1; }
      else               { ixw = j - npoints / 2;    sign = 1; }
      beta = 0.5 + sign / 2 * xNI[ixw];
      weight = wNI[ixw];
      printf("b%02d: beta = %.4f  w = %8.6f\n", j+1, beta, weight );
      sprintf(ctlfi, "%s.b%02d.ctl\0", ctlf, j+1);
      fctlb = (FILE*)gfopen(ctlfi, "w");
      fprintf(fctlb, "BayesFactorBeta = %8.6f *  w=%8.6f.ctl\n", beta, weight);

      rewind(fctl);
      for ( ; ; ) {
         if (fgets(line, nline, fctl) == NULL) break;
         if (line[0] == '*') continue;
         if (strstr(line, "BayesFactorBeta")) continue;
         if (strstr(line, "outfile") || strstr(line, "mcmcfile")) {
            pline = strchr(line, '=');
            sscanf(pline + 1, "%s", tmpf);
            if( s = strstr(tmpf, ".txt") ) {
               *s = '\0';
            }
            sprintf(pline + 2, "%s.b%02d.txt\n\0", tmpf, j+1);
         }
         fputs(line, fctlb);
      }
      fclose(fctlb);
      fprintf(fresults, "%.6f\t%.6f\t\n", beta, weight);
   }
   fclose(fctl); 

   fcommand = (FILE*)gfopen("commands", "w");
   fprintf(fcommand, "#!/bin/bash\nfor I in {01..%02d}\n", npoints);
   fprintf(fcommand, "  do\n");
   fprintf(fcommand, "     echo \"#!/bin/bash\" > %s\n", scriptf);
   fprintf(fcommand, "     echo \"mcmctree %s.b$I.ctl > log.b$I.txt\" > %s\n", ctlf, scriptf);
   fprintf(fcommand, "     sleep 1\n");
   fprintf(fcommand, "     qsub -S /bin/bash -l h_vmem=4G -l tmem=4G -l h_rt=360:0:0 -cwd %s\n", scriptf);
   fprintf(fcommand, "  done\n");

   fputs("\n#grep BFbeta log.b*.txt\n", fcommand);
   fclose(fcommand);
   exit(0);
}
