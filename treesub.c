/* TREESUB.c
*/

#include "tools.h"

extern char NUCs[], AAs[], BINs[];
extern int noisy;

#ifdef  BASEML
#define EIGEN
#define REALSEQUENCE
#define NODESTRUCTURE
#define TREESEARCH
#define LSDISTANCE
#define LFUNCTIONS
#define RECONSTRUCTION
#endif 

#ifdef  CODEML
#define EIGEN
#define REALSEQUENCE
#define NODESTRUCTURE
#define TREESEARCH
#define LSDISTANCE
#define LFUNCTIONS
#define RECONSTRUCTION
#endif 

#ifdef  BASEMLG
#define REALSEQUENCE
#define NODESTRUCTURE
#define LSDISTANCE
#endif 

#ifdef  RECONSTRUCTION
#define PARSIMONY
#endif 

#define EqPartition(p1,p2,ns) (p1==p2||p1+p2+1==(1<<ns))

#ifdef REALSEQUENCE

int hasbase (char *str)
{
   char *p=str, *eqdel=".-";
   while (*p) if (*p==eqdel[0] || *p==eqdel[1] || isalpha(*p++)) return(1);
   return(0);
}

int blankline (char *str)
{
   char *p=str;
   while (*p) if (isalnum(*p++)) return(0);
   return(1);
}

int sreadi (int nums[], int n, char str[])
{
/* read n numbers into nums[] from str[]
*/
   int i;
   char *p=str;

   FOR (i,n) {
      while (*p && !isalnum(*p)) p++;
      sscanf (p, "%d", &nums[i]);
      while (*p && isalnum(*p)) p++;
   }
   return (0);
}

int ReadSeq (FILE *fout, FILE *fseq, int seqtype)
{
/* read in sequence, specify options: ngene, lgene[]
   seqtype: 0=nucleotides;  1=codons;  2:AAs, 3:BINs
   com.pose[] used to store gene marks

   char opt_c[]="AKGI";
      A:alpha given. K:kapa given
      G:many genes,  I:interlaved format
*/
   char *line, *p, *eqdel=".-?";
   int  i,j,k, ch, noptline, lspname=LSPNAME, indel=0, simple=1;
   int  lline=255,lt[NS], igroup, Sequential=1;

   fscanf (fseq, "%d%d", &com.ns, &com.ls);
   if (com.ns>NS) error ("too many sequences.. raise NS?");
   if (seqtype==CODONseq && com.ls%3!=0) {
      printf ("\n%d nucleotides", com.ls); error ("seq. len. err.");
   }
   if (noisy) printf ("\nns = %d  \tls = %d\n", com.ns, com.ls);

   FOR (j, com.ns)  if (com.z[j]) free(com.z[j]);
   if (com.pose) free(com.pose);   if (com.fpatt) free(com.fpatt);
   com.pose=(int*) malloc (com.ls*sizeof(int));
   com.fpatt=(int*) malloc (com.ls*sizeof(int));

   FOR (j, com.ns) com.z[j]=(char*) malloc ((com.ls+1)*sizeof(char));
   if (com.z[com.ns-1]==NULL) error ("oom");
   FOR (j,com.ns) com.z[j][com.ls]='\0';
   com.rgene[0]=1;   com.ngene=1;  
   com.pose=(int*) malloc (com.ls*sizeof(int));
   com.fpatt=(int*) malloc (com.ls*sizeof(int));
   lline=max(lline, com.ls+lspname+10);
   line=(char*) malloc (lline*sizeof(char));
   if (com.pose==NULL || com.fpatt==NULL || line==NULL) error ("oom");
   FOR (j, com.ls) com.pose[j]=0;      /* gene #1, default */

   /* first line */
   if (!fgets (line, lline, fseq)) error ("err ReadSeq: first line");
   for (j=0,noptline=0; j<lline && line[j] && line[j]!='\n'; j++) {
      if (!isalnum(line[j])) continue;
      line[j]=(char)toupper(line[j]);
      switch (line[j]) {
         case 'G': case 'A': case 'K':
                    noptline++; break;
         case 'S':  Sequential=1;  break;
         case 'I':  Sequential=0;  break;
         default :  printf (".Bad option %c\n", line[j]);  exit (-1);
      }
   }
   /* option lines */
   FOR (j, noptline) {
      line[0]=(char)fgetc(fseq);
      line[0]=(char)toupper(line[0]);
      switch (line[0]) {
      case ('G') :
         fscanf (fseq, "%d", &com.ngene);
         if (com.ngene > NGENE) error ("raise NGENE?");
         if (com.given_rgene) {    /* specified in the main program? */
            puts ("reading rates for genes from the seq. file, correct?");
            FOR (k, com.ngene) if (fscanf (fseq, "%lf", &com.rgene[k]) != 1)
               error ("err: given_rgene?");
         }
         fgets (line, lline, fseq);
         if (!blankline(line)) {
            sreadi (com.lgene, com.ngene, line);
            for (i=0,k=0; i<com.ngene; k+=com.lgene[i],i++)
               FOR (j, com.lgene[i]) com.pose[k+j]=i;
            for (i=0,k=0; i<com.ngene; i++) k+=com.lgene[i];
            if (k!=(seqtype==CODONseq?com.ls/3:com.ls)) error("ls for genes");
         }
         else {
            for(k=0; k<(seqtype==CODONseq||seqtype==CODON2AAseq?com.ls/3:com.ls);) {
               if (com.ngene>9)  fscanf(fseq,"%d", &ch);
               else {
                  do ch=fgetc(fseq); while (!isdigit(ch));
                  ch=ch-(int)'1'+1;
               }
               if (ch<1 || ch>com.ngene)
                  { printf("\ngene mark %d at %d?\n", ch, k+1);  exit (-1); }
               com.pose[k++]=ch-1;
            }
            if (!fgets (line, lline, fseq)) error ("err: option lines");
         }
         break;
      case ('A') :
         com.given_alfa=1;
         fscanf (fseq, "%lf", &com.alfa);
         if (!fgets (line, lline, fseq)) error ("err: option lines");
         break;
      case ('K') :
         com.given_kapa=1;
         fscanf (fseq, "%lf", &com.kapa);
         if (!fgets (line, lline, fseq)) error ("err: option lines");
         break;
      default :
         printf ("..Bad option %c\n", line[0]);  exit (-1);
      }
   }
   /* read sequence */
   if (Sequential)  {    /* sequential */
      if (noisy) printf ("\nSequential format..\n");
      FOR (j,com.ns) {
         FOR (i, 2*lspname) line[i]='\0';
         if (!fgets (line, lline, fseq)) error ("err: EOF?");
         if (blankline(line)) {
            if (!fgets (line, lline, fseq)) error ("err: EOF");
            if (!hasbase(line)) error("err reading seq.: delete blank lines");
         }
         p=line+(line[0]=='=' || line[0]=='>') ;
         strncpy (com.spname[j], p, lspname);
         p+=lspname;
         for (k=lspname; k>0; k--)
            if (!isgraph(com.spname[j][k]))   com.spname[j][k]=0;
            else    break;
         if (noisy) printf ("Reading seq #%2d: %s\n", j+1, com.spname[j]);
         for (k=0; k<com.ls; p++) {
            while (*p=='\n' || *p=='\0')  p=fgets(line, lline, fseq);
            if (isalpha(*p))  com.z[j][k++] = toupper(*p);
            else if (*p == eqdel[0]) {
               if (j==0) error ("err: . in 1st seq.?");
               com.z[j][k] = com.z[0][k];  k++;
            }
            else if (*p==eqdel[1] || *p==eqdel[2])
               { com.z[j][k++] = *p; indel=1; }
            else if (*p==EOF) error ("err: EOF?");
         }
      }                 /* for (j,com.ns) */
   }
   else   {              /* interlaved */
      if (noisy) printf ("\nInterlaved format..\n");
      FOR (j, com.ns) lt[j]=0;  /* temporary seq length */
      for (igroup=0; ; igroup++) {
         FOR (j, com.ns) if (lt[j]<com.ls) break;
         if (j==com.ns) break;
         FOR (j,com.ns) {
            if (!fgets (line, lline, fseq)) {
               printf ("\n%d, seq %d group %d", lt[j], j+1, igroup+1);
               error ("err: EOF?");
            }
            if (!hasbase(line)) {
               if (j)  {
                  printf ("\nspecies %d in group %d\n", j+1, igroup+1);
                  error ("err: empty line.");
	       }
               else {
                  if (!fgets (line, lline, fseq)) {
                     printf ("\n%d, seq %d group %d", lt[j], j+1, igroup);
                     error ("err: EOF?");
                  }
                  if (!hasbase(line)) {
                     printf ("\nspecies %d in group %d\n", j+1, igroup+1);
                     error ("err: empty line.");
		  }
               }
            }
            p=line;
            if (igroup==0) {
               strncpy (com.spname[j], p, lspname);
               k=strlen(com.spname[j]);
               p+=(k<lspname?k:lspname);

               for (k=lspname; k>0; k--)
                  if (!isgraph(com.spname[j][k]))   com.spname[j][k]=0;
                  else    break;
               if(noisy) printf ("Reading seq #%2d: %s\n",j+1,com.spname[j]);
            }
            for (; *p && *p!='\n'; p++) {
               if (lt[j] == com.ls) break;
               if (isalpha(*p))  com.z[j][lt[j]++] = toupper(*p);
               else if (*p == eqdel[0]) {
                  if (j==0) error ("err: . in 1st seq.");
                  com.z[j][lt[j]] = com.z[0][lt[j]];
                  lt[j]++;  continue;
               }
               else if (*p==eqdel[1] || *p==eqdel[2])
                  { com.z[j][lt[j]++] = toupper(*p); indel=1; }
            }         /* for (*p) */
         }            /* for (j,com.ns) */
      }               /* for (igroup) */
   }

   free (line);

   if (fout) {
      fprintf (fout,"\n\nlist of sequences, sequential\n");
      FOR (j,com.ns)  {
         fprintf (fout,"%s\n", com.spname[j]);
         printaSeq (fout, com.z[j], com.ls, 60, 10);
      }
   }
/****************/
   if (indel) {
      RemoveIndel (com.z,&com.ns,&com.ls, com.pose,(seqtype==CODONseq?3:1),0);
      if (fout) {
         fprintf (fout,"\nAfter deleting indels..  ls:%10d\n", com.ls);
         printsma (fout, com.z, com.ns, com.ls, 60,
            (seqtype==CODONseq?3:10), 0, simple, -1);
      }
   }
   if (fout) {
      fprintf (fout,"\n\nlist of sequences, sequential\n");
      FOR (j,com.ns)  {
         fprintf (fout,"%s\n", com.spname[j]);
         printaSeq (fout, com.z[j], com.ls, 60, 10);
      }
/*
      fprintf (fout,"\n\nlist of sequences, interlaved\n");
      printsma (fout, com.z, com.ns, com.ls, 60,
         (seqtype==CODONseq?3:10), 0, simple, -1);
      FPN (fout);
*/
   }
   if (noisy>2) printf ("\nSequences read..\n");
   return (0);
}

int PatternWeight (FILE *fout, double space[])
{
/* posC: positions for constant patterns
   space[com.ls*sizeof(int) + com.ns*com.ls*sizeof(char)]
*/
   int ncode=com.ncode, posC[NCODE], h, ht, j, same=0, ig, *pose1=(int*)space;
   char *z1[NS];

   FOR (j,com.ns) z1[j]=(char*)(pose1+com.ls)+j*com.ls;

   FOR (h,com.ls) com.fpatt[h]=0; 
   FOR (ig, com.ngene) com.lgene[ig]=0;
   for (ig=0,com.npatt=0; ig<com.ngene; ig++) {
      for (j=0,com.posG[ig]=com.npatt; j<ncode; j++) posC[j]=-1;
      for (h=0; h<com.ls; h++) {
         if (com.pose[h] != ig) continue;
         com.lgene[ig]++;

         /* look for pattern identical to that at h  */
         for (j=1; j<com.ns; j++)
            if (com.z[j][h]!=com.z[0][h]) break;
         if (j==com.ns) {                  /* constant pattern? */
            if (posC[com.z[0][h]-1]!=-1)
               { pose1[h]=posC[com.z[0][h]-1]; }
            else {
               FOR (j, com.ns) z1[j][com.npatt]=com.z[j][h];
               pose1[h]=posC[com.z[0][h]-1]=com.npatt++;
            }
         }
         else {
            for (ht=com.posG[ig],same=0; ht<com.npatt; ht++) {
               for (j=0,same=1; j<com.ns; j++)
                  if (z1[j][ht]!=com.z[j][h]) { same=0; break; }
               if (same) { pose1[h]=ht; break; }
            }
            if (!same) {
               FOR (j, com.ns) z1[j][com.npatt]=com.z[j][h];
               pose1[h]=com.npatt++;
            }
         }
         com.fpatt[pose1[h]]++;
      }     /* for (h)   */
   }        /* for (ig)  */
   com.posG[com.ngene]=com.npatt;
   for (j=1; j<com.ngene; j++) com.lgene[j]+=com.lgene[j-1];
   if (fout) {
      fprintf (fout, "\nns = %d  \tls = %d", com.ns,com.ls);
      if (com.ngene>1) {
         fprintf (fout,"\nngene =%3d, lengths =", com.ngene);
         FOR (j, com.ngene)
            fprintf (fout,"%7d", (j?com.lgene[j]-com.lgene[j-1]:com.lgene[j]));
         fprintf (fout,"\n  Starting position =");
         FOR (j, com.ngene)  fprintf (fout,"%7d", com.posG[j]+1);
      }
      if (noisy) printf ("\n# site patterns =%7d\n", com.npatt);
      fprintf (fout,"\n# of site patterns =%7d\n", com.npatt);
      if (com.seqtype!=CODONseq) {
         fprintf(fout,"\nSite pattern freqs (patterns shown below)\n");
         FOR (h,com.npatt) {
            fprintf (fout," %3d", com.fpatt[h]);
            if ((h+1)%15==0) FPN (fout);
	 }
      }
   }
   FOR (j,com.ns) {
      memcpy (com.z[j], z1[j], com.npatt);
      com.z[j][com.npatt]=0;
      com.z[j]=(char*)realloc (com.z[j], (com.npatt+1)*sizeof(char));
   }
   memcpy (com.pose, pose1, com.ls*sizeof(int));
   com.fpatt = (int*)realloc(com.fpatt, com.npatt*sizeof(int));

   return (0);
}

int Initialize (FILE *fout, double *space, int seqtype)
{
/* frequencies for different genes, site patterns and fpatts
   for nucleotide and amino-acid sequences
*/
   char *pch=(seqtype==0?NUCs:(seqtype==2?AAs:BINs));
   int wname=15, h, j,k, ig, lt, nconstp, nc=com.ncode;
   double fb[NS*NCODE], t;

   FOR (j, com.ns) 
      if (transform (com.z[j], com.ls,1,seqtype)) printf(" in seq #%2d\n",j+1);
   PatternWeight (fout, space);

   if (fout) {
      printsma (fout, com.z, com.ns, com.npatt, 60, 10, 1, 0, seqtype);
      fprintf (fout,"\nFrequencies..");
   }
   zero (fb, nc*com.ns);  zero (com.pi, nc);
   t=1./(double)com.lgene[0];
   for (h=0,ig=0,lt=0,nconstp=0; h<com.npatt; ) {
      for (j=1; j<com.ns; j++)  if(com.z[j][h]!=com.z[0][h]) break;
      if (j==com.ns) nconstp += com.fpatt[h];
      FOR (j,com.ns) fb[j*nc+com.z[j][h]-1]+=(double)com.fpatt[h]*t;

      lt+=com.fpatt[h++];
      if (lt==com.lgene[ig]) {
         if (fout && com.ngene>1)
            fprintf (fout,"\n\nGene # %2d (len %4d)",
                ig+1, com.lgene[ig]-(ig==0?0:com.lgene[ig-1]));
         if (fout) {
            fprintf (fout,"\n%*s", wname, "");
            FOR (k,nc) fprintf (fout,"%7c", pch[k]);
	 }
         FOR (j,com.ns) {
            if (fout) {
               fprintf (fout,"\n%-*s", wname, com.spname[j]);
               FOR (k,nc) fprintf (fout,"%7.4f", fb[j*nc+k]);
            }
            FOR (k,nc) com.pi[k]+=fb[j*nc+k]/t/((double)com.ls*com.ns);
         }
         for (j=0,zero(com.piG[ig],nc); j<com.ns; j++) 
            FOR (k,nc) com.piG[ig][k]+=fb[j*nc+k]/com.ns;
         if (com.seqtype==NUCseq) {
            com.piG[ig][4]=com.piG[ig][0]+com.piG[ig][1];
            com.piG[ig][5]=com.piG[ig][2]+com.piG[ig][3];
	 }
         if (fout && com.ngene>1) {
            fprintf (fout,"\n\n%-*s", wname, "Mean");
            FOR (k,nc) fprintf(fout, "%7.4f", com.piG[ig][k]);
         }
         ig++;
         zero (fb,nc*com.ns);
         if (ig<com.ngene) t=1./(com.lgene[ig]-com.lgene[ig-1]);
       } /* ig */
    }    /* for (h) */
   if (com.seqtype==NUCseq) {
      com.pi[4]=com.pi[0]+com.pi[1];  com.pi[5]=com.pi[2]+com.pi[3];
   }

   if (fout) {
      fprintf (fout, "\n\n%-*s", wname, "Average");
      FOR (j,nc) fprintf (fout,"%7.4f", com.pi[j]);
      fprintf (fout,"\n\n# constant sites: %6d (%6.2f%%)",
               nconstp, (double)nconstp*100./com.ls);
   }
   if (seqtype==NUCseq && com.model<=1) {
      fillxc (com.pi, 0.25, 4);
      com.pi[4]=com.pi[0]+com.pi[1]; com.pi[5]=com.pi[2]+com.pi[3]; 
      FOR (j, com.ngene) xtoy (com.pi, com.piG[j], 6);
   }
   else if (com.seqtype==AAseq) {
      for (j=0,lt=0; j<nc; j++) lt+=(com.pi[j]>0);
      if (lt<=4)  printf ("\n\a\t\tAre these a.a. sequences?");
   }
   if (seqtype!=NUCseq) {
      if (com.model==0)  FOR (j,nc) com.pi[j]=1./nc;
      if (com.model==0)  FOR (j,com.ngene) xtoy(com.pi, com.piG[j], nc);
   }
   for (h=0,com.lmax=-(double)com.ls*log((double)com.ls); h<com.npatt; h++)
      if (com.fpatt[h]>1)
          com.lmax+=(double)com.fpatt[h]*log((double)com.fpatt[h]);
   if (fout&&com.ngene==1)  fprintf(fout, "\nln{Lmax}:%16.6f\n", com.lmax);
   if (fout) fflush (fout);
   if (noisy>2) printf ("\n\nInitialized..\n");
   return(0);
}

int MPInformSites (void)
{
/* Outputs parsimony informative and noninformative sites into 
   two files named MPinf.seq and MPninf.seq
   Uses transformed sequences
*/
   char *imark, *pch=(com.seqtype==0?NUCs:(com.seqtype==2?AAs:BINs));
   int h, i, markb[NS], inf, lsinf;
   FILE *finf, *fninf;

   finf=fopen("MPinf.seq","w");
   fninf=fopen("MPninf.seq","w");

   puts ("\nSorting parsimony-informative sites: MPinf.seq & MPninf.seq");
   if ((imark=(char*)malloc(com.ls*sizeof(char)))==NULL) error("oom");
   for (h=0,lsinf=0; h<com.ls; h++) {
      for (i=0; i<com.ns; i++) markb[i]=0;
      for (i=0; i<com.ns; i++) markb[(int)com.z[i][com.pose[h]]-1]++;

      for (i=0,inf=0; i<com.ncode; i++)  if (markb[i]>=2)  inf++;
      if (inf>=2) { imark[h]=1; lsinf++; }
      else imark[h]=0;
   }
   fprintf (finf, "%6d%6d\n", com.ns, lsinf);
   fprintf (fninf, "%6d%6d\n", com.ns, com.ls-lsinf);
   for (i=0; i<com.ns; i++) {
      fprintf (finf, "\n%s\n", com.spname[i]);
      fprintf (fninf, "\n%s\n", com.spname[i]);
      for (h=0; h<com.ls; h++)
         fprintf ((imark[h]?finf:fninf), "%c", pch[com.z[i][com.pose[h]]-1]);
      FPN (finf); FPN(fninf);
   }
   free (imark);
   fclose(finf);  fclose(fninf);
   return (0);
}

int PatternJC69like (FILE *fout)
{
/* further collaps of site patterns for JC69-like models, called after
   PatternWeight().
*/
   char zh[NS],b;
   int npatt0=com.npatt, h, ht, j,k, same=0, ig;

   for (h=0,com.npatt=0,ig=-1; h<npatt0; h++) {
      if (ig<com.ngene-1 && h==com.posG[ig+1]) { com.posG[++ig]=com.npatt; }
      zh[0]=b=1;
      for (j=1; j<com.ns; j++) {
         for (k=0; k<j; k++) if (com.z[j][h]==com.z[k][h]) break;
         zh[j]=(k<j?zh[k]:++b);
      }
      for (ht=com.posG[ig],same=0; ht<com.npatt; ht++) {
         for (j=0,same=1; j<com.ns; j++)
            if (zh[j]!=com.z[j][ht]) { same=0; break; }
         if (same) break; 
      }
      if (same)  com.fpatt[ht]+=com.fpatt[h];
      else {
         FOR (j, com.ns) com.z[j][com.npatt]=zh[j];
         com.fpatt[com.npatt++]=com.fpatt[h];
      }
      FOR (k, com.ls) if (com.pose[k]==h) com.pose[k]=ht;
   }     /* for (h)   */
   com.posG[com.ngene]=com.npatt;
   FOR (j, com.ns) com.z[j][com.npatt]=0;

   if (noisy>=2) printf ("\nnew no. site patterns:%7d\n", com.npatt);
   if (fout) fprintf (fout, "\nno. site patterns:%7d\n", com.npatt);

   if (fout) {
      FOR (h,com.npatt) {
         fprintf (fout," %3d", com.fpatt[h]);
         if ((h+1)%15==0) FPN (fout);
      }
      printsma (fout, com.z, com.ns, com.npatt, 60, 10, 1, 0, com.seqtype);
   }
   return (0);
}

int Site2Pattern (FILE *fout)
{
   int h;
   fprintf(fout,"\n\nMapping site to pattern (i.e. site %d has pattern %d):\n",
      com.ls-1, com.pose[com.ls-2]+1);
   FOR (h, com.ls) {
      fprintf (fout, "%6d", com.pose[h]+1);
      if ((h+1)%10==0) FPN (fout);
   }
   FPN (fout);
   return (0);
}

#endif

#define gammap(x,alfa) (alfa*(1-pow(x,-1/alfa)))
/* DistanceREV () used to be here... 
*/

#if (NCODE==4) 

double SeqDivergence (double x[], int model, double alfa, double *kapa)
{
/* HKY85 model, following Tamura (1993, MBE, ..)
   alfa=0 if no gamma 
   return -1 if in error.
*/
   int i,j;
   double p[4], Y,R, a1,a2,b, P1,P2,Q,fd,tc,ag, a1t,a2t,bt, largek=999;

   if (testXMat (x)) error ("X err..");
   for (i=0,fd=1,zero(p,4); i<4; i++) {
      fd -= x[i*4+i];
      FOR (j,4) { p[i]+=x[i*4+j]/2;  p[j]+=x[i*4+j]/2; }
   }
   P1=x[0*4+1]+x[1*4+0];
   P2=x[2*4+3]+x[3*4+2];
   Q = fd-P1-P2;
   Y=p[0]+p[1];    R=p[2]+p[3];  tc=p[0]*p[1]; ag=p[2]*p[3];

   switch (model) {
   case (JC69):
      FOR (i,4) p[i]=.25;
   case (F81):
      for (i=0,b=0; i<4; i++)  b += p[i]*(1-p[i]);
      if (1-fd/b<=0) return (-1);
      if (alfa<=0) return (-b*log (1-fd/b));
      else return  (-b*gammap(1-fd/b,alfa));
   case (K80) :
      a1=1-2*(P1+P2)-Q;   b=1-2*Q;
      if (a1<=0 || b<=0) return (-1);
      if (alfa<=0)  { a1=-log(a1);  b=-log(b); }
      else          { a1=-gammap(a1,alfa); b=-gammap(b,alfa); }
      a1t=.5*a1-.25*b;  bt=.25*b;
      *kapa = min(a1t/bt, largek);
      return (a1t+2*bt);
   case (F84):
      a1=(2*(tc+ag)+2*(tc*R/Y+ag*Y/R)*(1-Q/(2*Y*R)) -P1-P2) / (2*tc/Y+2*ag/R);
      b = 1 - Q/(2*Y*R);
      if (a1<=0 || b<=0) return (-1);
      if (alfa<=0) { a1=-log(a1); b=-log(b); }
      else         { a1=-gammap(a1,alfa); b=-gammap(b,alfa); }
      a1t=.5*a1;  bt=.5*b;
      *kapa = a1t/bt-1;
      /* *kapa = min(*kapa, largek);   */    *kapa = max(*kapa, -.5);
      return  4*bt*(tc*(1+ *kapa/Y)+ag*(1+ *kapa/R)+Y*R);
   case (HKY85):         /* TN93  */
      a1=1-Y*P1/(2*tc)-Q/(2*Y);  a2=1-R*P2/(2*ag)-Q/(2*R);   b=1-Q/(2*Y*R);
      if (a1<=0 || a2<=0 || b<=0) return (-1);
      if (alfa<=0) { a1=-log(a1); a2=-log(a2); b=-log(b); }
      else   { a1=-gammap(a1,alfa); a2=-gammap(a2,alfa); b=-gammap(b,alfa);}
      a1t=.5/Y*(a1-R*b);  a2t=.5/R*(a2-Y*b);  bt=.5*b;
      *kapa = min((a1t+a2t)/2/bt, largek);
      return 4*p[0]*p[1]*a1t + 4*p[2]*p[3]*a2t + 4*Y*R*bt;
   }    /* switch  */
   return (-1);
}


#ifdef LSDISTANCE
#ifdef REALSEQUENCE
extern double *SeqDistance;

int DistanceMatNuc (FILE *fout, int model, double alfa, double *kapa)
{
   char *models[]={"JC69", "K80", "F81", "F84", "TN93"};
   int i,j, h, status=0;
   double x[16], kapat=0, t;

   if (model>=HKY85) model=4;
   if (fout) {
       fprintf(fout,"\nDistances:%5s ", models[model]);
       if (model!=JC69 && model!=F81) fprintf (fout, "(kapa) ");
       fprintf(fout," (alpha set at %.2f)\n", alfa);
   }
   for (i=0,*kapa=0; i<com.ns; i++) {
      if (fout) fprintf(fout,"\n%-12s", com.spname[i]);
      FOR (j, i) {
         for (h=0,zero(x,16); h<com.npatt; h++)
            x[(com.z[i][h]-1)*4 + com.z[j][h]-1] +=
               (double)com.fpatt[h]/(double)com.ls;
         if ((t=SeqDivergence (x, model, alfa, &kapat)) < 0)
            { puts ("\alarge distance..");   t = 5;  status=-1; }
         SeqDistance[i*(i-1)/2+j] = t;

         *kapa += kapat/(double)(com.ns*(com.ns-1)/2);
         if (fout) fprintf(fout,"%8.4f", SeqDistance[i*(i-1)/2+j]);
         if (fout && (model==K80 || model==F84 || model==HKY85))
            fprintf(fout,"(%7.4f)", kapat);
      }
   }
   if (fout && (model==K80 || model==F84 || model==HKY85)) 
      fprintf (fout, "\n\nMean kappa:%10.4f\n", *kapa);
   else 
      FPN (fout);
   return (status);
}

#endif
#endif

#ifdef BASEMLG
extern int CijkIs0[];
#endif

extern int nR;
extern double Cijk[], Root[];

int RootTN93 (int model, double kapa1, double kapa2, double pi[], 
    double *f, double Root[])
{
   double *p=pi;

   if (model==F84) { kapa2=1+kapa1/p[5]; kapa1=1+kapa1/p[4]; }
   else if (model==HKY85) kapa2=kapa1;

   *f=1/(2*p[0]*p[1]*kapa1+2*p[2]*p[3]*kapa2 + 2*p[4]*p[5]);

   Root[0] = 0;            /* Root[0]=0 for any model */
   Root[1] = - (*f);
   Root[2] = -(p[4]+p[5]*kapa2) * (*f);
   Root[3] = -(p[4]*kapa1+p[5]) * (*f);
   return (0);
}

int EigenTN93 (int model, double kapa1, double kapa2, double pi[6],
    int *nR, double Root[], double Cijk[])
{
/* initialize Cijk[] & Root[], which are the only part to be changed
   for a new substitution model
   for JC69, K80, F81, F84, HKY85, TN93
   Root: real Root divided by v, the number of nucleotide substitutions.
*/
   int i,j,k, nr;
   double *p=pi, f, U[16], V[16], t;

   if (model==JC69 || model==F81) kapa1=kapa2=com.kapa=1; 
   else if (com.model<=HKY85)     kapa2=kapa1;
   RootTN93 (model, kapa1, kapa2, p, &f, Root);

   *nR=nr = 2+(model==K80||model>=F84)+(model>=HKY85);
   U[0*4+0]=U[1*4+0]=U[2*4+0]=U[3*4+0]=1;
   U[0*4+1]=U[1*4+1]=1/p[4];   U[2*4+1]=U[3*4+1]=-1/p[5];
   U[0*4+2]=U[1*4+2]=0;  U[2*4+2]=p[3]/p[5];  U[3*4+2]=-p[2]/p[5];
   U[2*4+3]=U[3*4+3]=0;  U[0*4+3]=p[1]/p[4];  U[1*4+3]=-p[0]/p[4];

   xtoy (p, V, 4);
   V[1*4+0]=p[5]*p[0];   V[1*4+1]=p[5]*p[1];
   V[1*4+2]=-p[4]*p[2];  V[1*4+3]=-p[4]*p[3];
   V[2*4+0]=V[2*4+1]=0;  V[2*4+2]=1;   V[2*4+3]=-1;
   V[3*4+0]=1;  V[3*4+1]=-1;   V[3*4+2]=V[3*4+3]=0;

   FOR (i,4) FOR (j,4) {
      Cijk[i*4*nr+j*nr+0]=U[i*4+0]*V[0*4+j];
      switch (model) {
      case JC69:
      case F81:
         for (k=1,t=0; k<4; k++) t += U[i*4+k]*V[k*4+j];
         Cijk[i*4*nr+j*nr+1] = t;
         break;
      case K80:
      case F84:
         Cijk[i*4*nr+j*nr+1]=U[i*4+1]*V[1*4+j];
         for (k=2,t=0; k<4; k++) t += U[i*4+k]*V[k*4+j];
         Cijk[i*4*nr+j*nr+2]=t;
         break;
      case HKY85:   case TN93:
         for (k=1; k<4; k++)  Cijk[i*4*nr+j*nr+k] = U[i*4+k]*V[k*4+j];
         break;
      default:
         error ("err:model");
      }
   }
#ifdef BASEMLG
   FOR (i,64) CijkIs0[i] = (Cijk[i]==0);
#endif
   return(0);
}

#ifdef EIGEN

int EigenREV (FILE* fout, double kapa[], double pi[], 
              int *nR, double Root[], double Cijk[])
{
/* pi[] is constant
*/
   int i,j,k;
   double Q[16], U[16], V[16], T1[16], T2[16], mr;

   *nR=4;
   for (i=0,k=0; i<3; i++) for (j=i+1; j<4; j++)
      if (i*4+j!=11) Q[i*4+j]=Q[j*4+i]=kapa[k++];
   for (i=0,Q[3*4+2]=Q[2*4+3]=1; i<4; i++) FOR (j,4) Q[i*4+j] *= pi[j];
   for (i=0,mr=0; i<4; i++) 
      { Q[i*4+i]=0; Q[i*4+i]=-sum(Q+i*4, 4); mr-=pi[i]*Q[i*4+i]; }
   abyx (1/mr, Q, 16);
   if (fout) {
      mr=2*(pi[0]*Q[0*4+1]+pi[2]*Q[2*4+3]);
      fprintf (fout, "\nRate matrix Q, Average Ts/Tv =%9.4f", mr/(1-mr));
      matout (fout, Q, 4,4);
   }
   if ((k=eigen (1, Q, 4, Root, T1, U, V, T2))!=0) 
      error ("\ncomplex roots in EigenREV");
   xtoy (U, V, 16);
   matinv (V, 4, 4, T1);
   FOR (i,4) FOR(j,4) FOR(k,4) Cijk[i*4*4+j*4+k] = U[i*4+k]*V[k*4+j];
   return (0);
}
   
int EigenUNREST (FILE *fout, double kapa[], double pi[], 
    int *nR, complex cRoot[], complex cU[], complex cV[])
{
/* pi[] is changed
*/
   int i,j,k;
   double Q0[16], Q[16], U[16], V[16], T1[16], T2[16], mr;

   *nR=4;
   for (i=0,k=0; i<4; i++) FOR (j,4) 
      if (i!=j && i*4+j != 14)  Q[i*4+j]=kapa[k++];
   for (i=0,Q[14]=1; i<4; i++)  { Q[i*4+i]=0; Q[i*4+i]=-sum(Q+i*4, 4); }
   xtoy (Q, Q0, 16);

   i = eigen (1, Q, 4, Root, T1, U, V, T2) ; 

   FOR (i, 4)  { cRoot[i].re=Root[i], cRoot[i].im=T1[i]; }
   FOR (i, 16) { cU[i].re=cV[i].re=U[i]; cU[i].im=cV[i].im=V[i]; }
   cmatinv (cV, 4, 4, T1);

   FOR (i, 4) pi[i]=cV[0*4+i].re;
   abyx (1/sum(pi,4), pi, 4);

   for (i=0,mr=0; i<4; i++)  mr -= pi[i]*Q0[i*4+i];
   FOR (i,4) { cRoot[i].re /= mr; cRoot[i].im /= mr; }
   if (fout) {
      abyx (1/mr, Q0, 16);
      mr=pi[0]*Q0[0*4+1]+pi[1]*Q0[1*4+0]+pi[2]*Q0[2*4+3]+pi[3]*Q0[3*4+2];
      fprintf (fout,"\npi and rate matrix Q, Average Ts/Tv =%9.4f",mr/(1-mr));
      matout (fout, pi, 1, 4);
      matout (fout, Q0, 4, 4);
   }
   return (0);
}

int cPMat (double P[],double t,int n,complex cU[],complex cV[],complex cRoot[])
{
/* P(t) = cU * exp{cRoot*t} * cV
*/
   int i,j,k, status=0;
   complex cUd[NCODE*NCODE], cP, cY;
   double sum;

   if (t<1e-6) { identity (P, n); return(0); }
   FOR (i,n) cUd[i*n+0]=cU[i*n+0];
   for (j=1; j<n; j++) {
      cY.re=cRoot[j].re*t; cY.im=cRoot[j].im*t; cY=cexp(cY);
      for (i=0; i<n; i++)  cUd[i*n+j]=cby(cU[i*n+j],cY);
   }
   FOR (i,n)   {
      for (j=0,sum=0; j<n; j++) {
         for (k=0,cP=compl(0,0); k<n; k++) {
            cY = cby(cUd[i*n+k],cV[k*n+j]);
            cP.re+=cY.re;  cP.im+=cY.im;
         }
         P[i*n+j]=cP.re;
         sum+=P[i*n+j];
         if (P[i*n+j]<=0 || fabs(cP.im)>1e-4) status=-1;
      }
      if (fabs(sum-1)>1e-4) status=-1;
   }

   if (status)
      { printf ("\nerr cPMat.."); matout (F0, P, n, n); getchar(); }

   return (0);
}

#endif   /* ifdef BASEML  */

#endif   /* if (NCODE==4) */

#ifdef LSDISTANCE

static double *SeqDistance; 
static int *ancestor;

int fun_LS (double v[], double diff[], int np, int npair);

int fun_LS (double v[], double diff[], int np, int npair)
{
   int i,j, aa, it;
   double dexp;

   if (SetBranch(v, 0)) puts ("branch len err.");
   if (npair != com.ns*(com.ns-1)/2) error ("# seq pairs err.");
   FOR (i, com.ns) FOR (j, i) {
      it=i*(i-1)/2+j;
      for (aa=i,dexp=0; aa!=ancestor[it]; aa=nodes[aa].father)
         dexp += nodes[aa].branch;
      for (aa=j; aa!=ancestor[it]; aa=nodes[aa].father)
         dexp += nodes[aa].branch;
      diff[it]=SeqDistance[it]-dexp;
   }
   return(0);
}

int LSDistance (double *ss, double x[], int (*testx)(double x[],int np))
{
/* get Least Squares estimates of branch lengths for a given tree topology
*/
   int i,j, h, a1,a2;

   if ((*testx)(x, com.ntime)) {
      matout (F0, x, 1, com.ntime);
      error ("initial err in LSDistance()");
   }
   FOR (i, com.ns) FOR (j, i) {
      h=i*(i-1)/2+j;
      ancestor[h]=-1;
      for (a1=i; a1!=-1; a1=nodes[a1].father) {
         for (a2=j; a2!=-1; a2=nodes[a2].father)
            if (a1==a2) { ancestor[h]=a1; break; }
         if (ancestor[h] != -1) break;
      }
      if (ancestor[h] == -1) error ("no ancestor");
   }      
   i=nls2 (NULL, ss, x, com.ntime, fun_LS, NULL,
           testx, com.ns*(com.ns-1)/2, 1e-8, 1e-8);
   return (i);
}

#endif 

#ifdef NODESTRUCTURE

void BranchToNode (void)
{
   int i,j, from, to;
   
   tree.nnode=tree.nbranch+1;
   if (tree.origin<0 || tree.origin>com.ns*2-2) 
      { printf ("root at %d", tree.origin+1); error ("tree origin"); }
   FOR (j,tree.nnode) ClearNode (j);
   FOR (i,tree.nbranch) {
      from=tree.branches[i][0];
      to  =tree.branches[i][1];
      nodes[from].sons[nodes[from].nson++]=to;
      nodes[to].father=from;
      nodes[to].ibranch=i;
   }
/*
   printf("\nNode\n%7s%7s%7s%7s%7s\n","father","node","branch","nson:","sons");
   FOR (i, tree.nnode) {
      printf ("\n%7d%7d%7d:%7d",
         nodes[i].father+1, i+1, nodes[i].ibranch+1, nodes[i].nson);
      FOR (j, nodes[i].nson) printf ("%4d", nodes[i].sons[j]+1);
   }
*/
}

void NodeToBranchSub (int inode)
{
   int i, ison;

   FOR (i, nodes[inode].nson) {
      tree.branches[tree.nbranch][0]=inode;
      tree.branches[tree.nbranch][1]=ison=nodes[inode].sons[i];
      nodes[ison].ibranch=tree.nbranch++;
      if(nodes[ison].nson>0)  NodeToBranchSub(ison);
   }
}

void NodeToBranch (void)
{
   tree.nbranch=0;
   NodeToBranchSub (tree.origin);
   if (tree.nnode!=tree.nbranch+1)  puts ("nd!=nb+1");
}

void ClearNode (int inode)
{
   nodes[inode].father=nodes[inode].ibranch=-1;
   nodes[inode].nson=0;
/*   FOR (i, com.ns) nodes[inode].sons[i]=-1; */
}

int ReadaTreeB (FILE *ftree, int popline)
{
   char line[255];
   int i,j, state=0, YoungAncestor=0;

   fscanf (ftree, "%d", &tree.nbranch);
   FOR (j, tree.nbranch) {
      FOR (i,2) {
         if (fscanf (ftree, "%d", & tree.branches[j][i]) != 1) state=-1;
         tree.branches[j][i]--;
      }
      if (tree.branches[j][0]<com.ns) YoungAncestor=1;
/*
      printf ("\nBranch #%3d: %3d -> %3d",
         j+1, tree.branches[j][0]+1,tree.branches[j][1]+1);
*/
   }
   if (popline) fgets (line, 254, ftree);
   tree.origin=tree.branches[0][0];
   if (YoungAncestor && com.clock) error ("clock with young ancestors?");
   com.ntime = com.clock ? tree.nbranch+1-com.ns: tree.nbranch;

   BranchToNode ();
   return (state);
}

int OutaTreeB (FILE *fout)
{
   int j;
   char *fmt[]={" %3d..%-3d", " %2d..%-2d"};
   FOR (j, tree.nbranch)
      fprintf(fout, fmt[0], tree.branches[j][0]+1,tree.branches[j][1]+1);
   return (0);
}

int ReadaTreeN (FILE *ftree, int *hasbranch, int popline)
{
/* Read a tree from ftree, using the parenthesis node representation of trees.
   If branch lengths are provided, they are read in nodes[].branch, with
   hasbranch set to 1.  Both names and numbers for species are fine.  
   Species names are considered case-sensitive, with trailing blanks ignored;
   they have to match the names in the sequence data file.
*/
   int cnode, cfather=-1;  /* current node and father */
   int inodeb=0;  /* node number that will have the next branch length */
   int i,j, level=0, hasname, ch=' ';
   char check[NS], line[255], delimiters[]="(),:";

   while(isspace(ch)) ch=fgetc(ftree);
   ungetc(ch, ftree);
   if (isdigit(ch)) { ReadaTreeB (ftree, popline); return (0); }

   tree.nnode=com.ns;  tree.nbranch=0;  *hasbranch=0;
   FOR (i, 2*com.ns-1) ClearNode (i);
   FOR (i, com.ns) check[i]=0;
   for (;;) {
      ch = fgetc (ftree);
      if (ch==EOF) error ("err in user tree..");
      else if (!isgraph(ch)) continue;
      else if (ch=='(') {
         level++;
         cnode=tree.nnode++;
         if (cfather>=0) {
            nodes[cfather].sons[nodes[cfather].nson++] = cnode;
            nodes[cnode].father=cfather;
            tree.branches[tree.nbranch][0]=cfather;
            tree.branches[tree.nbranch][1]=cnode;
            nodes[cnode].ibranch=tree.nbranch++;
         }
         else
            tree.origin=cnode;
         cfather=cnode;
      }
      else if (ch==')')
         { level--;  inodeb=cfather; cfather=nodes[cfather].father; }
      else if (ch==':') {
         *hasbranch=1; fscanf(ftree, "%lf", &nodes[inodeb].branch);
      }
      else if (ch==',') ;
      else { /* read species name or number */
         line[0]=ch;  line[1]=fgetc(ftree);
         if (com.ns<10 && isdigit(line[0]) && isdigit(line[1])) 
           { ungetc(line[1], ftree); line[1]=0; }
         else 
            for (i=1; ; )  { 
               if (strchr(delimiters,line[i]))
                  { ungetc(line[i], ftree); line[i]=0; break; }
               line[++i]=fgetc(ftree);
	    }
         for (j=i-1; j>0; j--)    if (!isgraph(line[j])) line[j]=0;
         for (i=0,hasname=0; line[i]; i++)  if (!isdigit(line[i])) hasname=1;

         if (hasname) {   /* name */
            for (i=0; i<com.ns; i++) if (!strcmp(line,com.spname[i])) break;
            if ((cnode=i)==com.ns) printf("\nSpecies %s?\n", line);
	 }
         else {           /* number */
            sscanf(line, "%d", &cnode);   cnode--;
            if (cnode<0 || cnode>=com.ns) printf("\nsp number %d\n", cnode+1);
	 }
         nodes[cnode].father=cfather;
         nodes[cfather].sons[nodes[cfather].nson++]=cnode;
         tree.branches[tree.nbranch][0]=cfather;
         tree.branches[tree.nbranch][1]=cnode;
         nodes[cnode].ibranch=tree.nbranch++;
         check[inodeb=cnode]++;
      }
      if (level==0) break;
   }
   if (popline) fgets (line, 254, ftree);
   if (tree.nnode != tree.nbranch+1)
      printf ("\n??nnode%6d   nbranch%6d\n", tree.nnode, tree.nbranch);
   com.ntime = com.clock ? tree.nnode-com.ns+(tree.origin<com.ns)
                         : tree.nbranch;
   FOR (i, com.ns)  if (check[i]!=1) return (-1);

   return (0);
}

int OutSubTreeN (FILE *fout, int inode, int spnames, int branchlen);

int OutSubTreeN (FILE *fout, int inode, int spnames, int branchlen)
{
   int i;

   fputc ('(', fout);

   FOR (i, nodes[inode].nson) {
      if (nodes[nodes[inode].sons[i]].nson <= 0) {
         if (spnames) fprintf (fout, "%s", com.spname[nodes[inode].sons[i]]);
         else         fprintf (fout, "%d", nodes[inode].sons[i]+1);
      }
      else
         OutSubTreeN (fout, nodes[inode].sons[i], spnames, branchlen);
      if (branchlen)fprintf (fout,":%.5f", nodes[nodes[inode].sons[i]].branch);
      if (com.ns>9 || spnames || branchlen)
         if (i<nodes[inode].nson-1) fprintf (fout, (spnames ? ", " : ","));
   }
   fputc (')', fout);
   return (0);
}


int OutaTreeN (FILE *fout, int spnames, int branchlen)
{
   OutSubTreeN (fout, tree.origin, spnames, branchlen);
   return (0);
}

void PointLklnodes (void)
{ 
   int i, nintern=0, ncode=com.ncode;
   FOR (i, tree.nbranch+1)
      if (nodes[i].nson>0)
         nodes[i].lkl = com.chunk + ncode*com.npatt*nintern++;
/*
   if (nintern!=tree.nbranch-com.ns+1)
      printf ("\nerr in PointLklnodes, # intern nodes: %d", nintern);
*/
}

int SetBranch (double x[], int correct)
{
   static direction=1;
   int i, state=0, father;
   double fmid=4e-6;

   if (com.clock) {
       FOR (i,com.ntime) nodes[i+com.ns].divtime=x[i];
       FOR (i,tree.nnode) {
          if (i==tree.origin) continue;
          father=nodes[i].father;
          if (nodes[father].divtime<nodes[i].divtime && !correct) state=-1;
          else if (correct==1 && nodes[father].divtime<=nodes[i].divtime) {
             if (direction==1)
              nodes[i].divtime=x[i-com.ns]=nodes[father].divtime*(1-fmid);
             else
              nodes[father].divtime=x[father-com.ns]=nodes[i].divtime*(1+fmid);
             direction*=-1;
          }
          nodes[i].branch=nodes[father].divtime-nodes[i].divtime;
       }
   }
   else
      FOR (i,tree.nnode)
         if (i!=tree.origin && (nodes[i].branch=x[nodes[i].ibranch])<-1e-6)
            state=-1;
   return (state);
}


int CollapsNode (int inode, double x[]) 
{
/* Merge inode to its father. Update the first com.ntime elments of
   x[] only if (x!=NULL), by using either x[] if clock=1 or
   nodes[].branch if clock=0.  So when clock=0, the routine works
   properly only if SetBranch() is called before this routine, which
   is true if m.l. or l.s. has been used to estimate branch
   lengths.  
*/
   int i,j, ifather, ibranch, ison;

   if (inode==tree.origin || inode<com.ns) error("err CollapsNode");
   ibranch=nodes[inode].ibranch;   ifather=nodes[inode].father; 
   for (i=0; i<nodes[inode].nson; i++) {
      ison=nodes[inode].sons[i];
      tree.branches[nodes[ison].ibranch][0]=ifather;
   }
   for (i=ibranch+1; i<tree.nbranch; i++) 
      for (j=0; j<2; j++) tree.branches[i-1][j]=tree.branches[i][j];
   tree.nbranch--; com.ntime--;
   for (i=0; i<tree.nbranch; i++)  for (j=0; j<2; j++) 
        if (tree.branches[i][j]>inode)  tree.branches[i][j]--;
   BranchToNode();

   if (x) {
      if (com.clock) 
         for (i=inode+1; i<tree.nnode+1; i++) x[i-1-com.ns]=x[i-com.ns];
      else {
         for (i=ibranch+1; i<tree.nbranch+1; i++)  x[i-1]=x[i];
         SetBranch (x, 0);
      }
   }
/*
printf ("\ninode & ibranch: %d & %d\n", inode+1, ibranch+1);
SetBranch (x, 0);
OutaTreeN (F0, 0, 0);  FPN(F0);
OutaTreeB (F0);  FPN(F0);
FOR (i, com.ntime) printf ("%9.5f", x[i]);  FPN (F0);
OutaTreeN (F0, 1, 1);  FPN(F0);
*/
   return (0);
}


static int PARTITION=0;
void DescentGroup (int inode);
void DescentGroup (int inode)
{
   int i;
   for (i=0; i<nodes[inode].nson; i++) 
      if (nodes[inode].sons[i]<com.ns) PARTITION |= (1<<nodes[inode].sons[i]);
      else DescentGroup (nodes[inode].sons[i]);
}
void BranchPartition (int partition[], int parti2B[])
{
/* calculates branch partitions for tree.
   partition[nib] has bits on if the species are down that branch.
   parti2B[nib] maps nib to ibranch.
*/
   int i, nib=0;

   if (com.ns>32) error(">32 species");
   for (i=0,nib=0; i<tree.nbranch; i++) {
      if (tree.branches[i][1]>=com.ns){
         PARTITION=0;
         DescentGroup (tree.branches[i][1]);
         partition[nib]=PARTITION;
         parti2B[nib++]=i;
      }
   }  
}
int NSameBranch (int partition1[],int partition2[], int nib1,int nib2)
{
/* counts the number of correct bipartitions, for com.ns<32.
   nib1 and nib2 are the numbers of interior branches in the two trees
*/
   int i,j,count;

   for (i=0,count=0; i<nib1; i++)  for(j=0;j<nib2;j++)
       count+=EqPartition(partition1[i], partition2[j], com.ns);
   return (count);
}

void DescentGroupLarge (int inode);
void BranchPartitionLarge (char partition[], int parti2B[]);

static char *PARTITIONLARGE;

void DescentGroupLarge (int inode)
{
   int i;
   for (i=0; i<nodes[inode].nson; i++) 
      if (nodes[inode].sons[i]<com.ns) 
         PARTITIONLARGE[nodes[inode].sons[i]]=1;
      else 
         DescentGroupLarge (nodes[inode].sons[i]);
}

void BranchPartitionLarge (char partition[], int parti2B[])
{
/* calculates branch partitions for a large tree (ns>32).
   partition[nib*com.ns+j]=1 if species j is down that interior branch.
   parti2B[nib] maps nib to ibranch.
   partition[nib*com.ns].  nib: # of interior branches.
*/
   int i,j, nib;  /* number of internal branches */

   for (i=0,nib=0; i<tree.nbranch; i++) {
      if (tree.branches[i][1]>=com.ns){
         PARTITIONLARGE=partition+nib*com.ns;
         FOR (j,com.ns) PARTITIONLARGE[j]=0;
         DescentGroupLarge (tree.branches[i][1]);
         if (parti2B) parti2B[nib]=i;
         nib++;
      }
   }
   if (nib!=tree.nbranch-com.ns) error("err BranchPartitionLarge"); 
}


int NSameBranchLarge (char partition1[],char partition2[], int nib1,int nib2)
{
/* counts the number of correct bipartitions 
   nib1 and nib2 are the numbers of interior branches in the two trees
*/
   int i,j,k, count,nsame;

   for (i=0,count=0; i<nib1; i++)  for(j=0;j<nib2;j++) {
       for (k=0,nsame=0;k<com.ns;k++)
          nsame+=(partition1[i*com.ns+k]==partition2[j*com.ns+k]);
       if (nsame==0 || nsame==com.ns) count++;   /* room for improvement */
   }
   return (count);
}



#ifdef TREESEARCH

struct TREE {struct TREEB tree; struct TREEN nodes[2*NS-1]; double x[NP]; } 
       treebest, treestar;
int DecompTree (int inode, int ison1, int ison2);
#define hdID(i,j) (max(i,j)*(max(i,j)-1)/2+min(i,j))

int TreeSearch (FILE *fout, double space[]);

int TreeSearch (FILE *fout, double space[])
{
/* automatic tree search by star decomposition, nhomo<=1
   returns (0,1,2,3) for the 4s problem.
*/
   int status=0,stage=0, i,j, itree,ntree=0,ntreet,best=0,improve=1,collaps=0;
   int inode, nson=0, ison1,ison2, son1, son2;
   char treef[12]="trees. s";
   double x[NP], e=1e-7;
   FILE *ftree, *fsum=frst;
/*
fsum=NULL;
*/
   if (com.runmode==1) {   /* read the star-like tree from trees.*s */
      sprintf (treef, "trees.%ds", com.ns);
      if ((ftree=fopen (treef,"r"))==NULL) error ("no treefile");
      fscanf (ftree, "%d%d", &i, &ntree);
      if (ReadaTreeN(ftree, &i, 1)) error ("err tree file");
      fclose (ftree);
   }
   else {                  /* construct the star tree of ns species */
      tree.nnode = (tree.nbranch=tree.origin=com.ns)+1;
      for (i=0; i<tree.nbranch; i++)
         { tree.branches[i][0]=com.ns; tree.branches[i][1]=i; }
      com.ntime = com.clock?1:tree.nbranch;
      BranchToNode ();
   }
   if (noisy) { printf ("\n\nstage 0: ");  OutaTreeN (F0, 0, 0); }
   if (fsum) { fprintf (fsum,"\n\nstage 0: ");  OutaTreeN (fsum, 0, 0); }
   if (fout) { fprintf (fout,"\n\nstage 0: ");  OutaTreeN (fout, 0, 0); }
   PointLklnodes ();
   GetInitials (x, space);
   NFunCall=0;
   status=ming1(NULL,&tree.lnL,com.plfun,NULL,testx,x,space,e,com.np);

   if (noisy)  printf("\nlnL:%14.6f%6d", -tree.lnL, NFunCall);
   if (fsum) fprintf(fsum,"\nlnL:%14.6f%6d", -tree.lnL, NFunCall);
   if (fout) {
      fprintf(fout,"\nlnL(ntime:%3d  np:%3d):%14.6f\n",
         com.ntime, com.np, -tree.lnL);
      OutaTreeB (fout);  FPN(fout);
      FOR (i, com.np) fprintf (fout,"%9.5f", x[i]);  FPN (fout);
   }
   treebest.tree=tree;  memcpy (treebest.nodes, nodes, sizeof(nodes));
   FOR (i,com.np) treebest.x[i]=x[i];

   for (ntree=0,stage=1; ; stage++) {
      for (inode=treebest.tree.nnode-1; inode>=0; inode--) {
         nson=treebest.nodes[inode].nson;
         if (nson>3) break;
         if (com.clock) { if (nson>2) break; }
         else if (nson>2+(inode==treebest.tree.origin)) break;
      }
      if (inode==-1 || /*stage>com.ns-3+com.clock ||*/ !improve) { /* end */
         tree=treebest.tree;  memcpy (nodes, treebest.nodes, sizeof(nodes));
         if (noisy) {
            printf ("\n\nbest tree: ");  OutaTreeN (F0, 0, 0);
            printf ("   lnL:%14.6f\n", -tree.lnL);
         }
         if (fsum) {
            fprintf (fsum, "\n\nbest tree: ");  OutaTreeN (fsum, 0, 0);
            fprintf (fsum, "   lnL:%14.6f\n", -tree.lnL);
	 }
         if (fout) {
            fprintf (fout, "\n\nbest tree: ");  OutaTreeN (fout, 0, 0);
            fprintf (fout, "   lnL:%14.6f\n", -tree.lnL);
            OutaTreeN (fout, 1, 1);  FPN(fout);
         }
         break;
      }
      treestar=treebest;  memcpy (nodes, treestar.nodes, sizeof(nodes));

      if (collaps && stage) { 
         printf ("\ncollapsing nodes\n");
         OutaTreeN (F0, 1, 1);  FPN(F0);

         tree=treestar.tree;  memcpy (nodes, treestar.nodes, sizeof(nodes));
         for (i=com.ns,j=0; i<tree.nnode; i++)
            if (i!=tree.origin && nodes[i].branch<1e-7) 
               { CollapsNode (i, treestar.x);  j++; }
         treestar.tree=tree;  memcpy (treestar.nodes, nodes, sizeof(nodes));

         if (j)  { 
            fprintf (fout, "\n%d node(s) collapsed\n", j);
            OutaTreeN (fout, 1, 1);  FPN(fout);
	 }
         if (noisy) {
            printf ("\n%d node(s) collapsed\n", j);
            OutaTreeN (F0, 1, 1);  FPN(F0);
/*            if (j) getchar (); */
	 }
      }

      ntreet = nson*(nson-1)/2;
      if (!com.clock && inode==treestar.tree.origin && nson==4)  ntreet=3;
      com.ntime++;  com.np++;

      if (noisy) {
         printf ("\n\nstage %d:%6d trees, ntime:%3d  np:%3d\nstar tree: ",
            stage, ntreet, com.ntime, com.np);
         OutaTreeN (F0, 0, 0);
         printf ("  lnL:%10.3f\n", -treestar.tree.lnL);
      }
      if (fsum) {
       fprintf (fsum, "\n\nstage %d:%6d trees, ntime:%3d  np:%3d\nstar tree: ",
         stage, ntreet, com.ntime, com.np);
         OutaTreeN (fsum, 0, 0);
         fprintf (fsum, "  lnL:%10.6f\n", -treestar.tree.lnL);
      }
      if (fout) {
         fprintf (fout,"\n\nstage %d:%6d trees\nstar tree: ", stage, ntreet);
         OutaTreeN (fout, 0, 0);
         fprintf (fout, " lnL:%14.6f\n", -treestar.tree.lnL);
         OutaTreeN (fout, 1, 1);  FPN (fout);
      }

      for (ison1=0,itree=improve=0; ison1<nson; ison1++)
      for (ison2=ison1+1; ison2<nson&&itree<ntreet; ison2++,itree++,ntree++) {
         DecompTree (inode, ison1, ison2);
         son1=nodes[tree.nnode-1].sons[0];
         son2=nodes[tree.nnode-1].sons[1];

         for(i=com.np-1; i>0; i--)  x[i]=treestar.x[i-1];
         if (!com.clock)
            for (i=0; i<tree.nbranch; i++)
               x[i]=max(nodes[tree.branches[i][1]].branch*0.99, 0.0001);
         else
            for (i=com.ns; i<tree.nnode; i++)
               x[i-com.ns]=max(nodes[i].divtime, 0.0001);
         PointLklnodes ();
         NFunCall=0;

         if (noisy) {
            printf("\nS=%d:%3d/%d  T=%4d  ", stage,itree+1,ntreet,ntree+1);
            OutaTreeN (F0, 0, 0);
         }
         if (fsum) {
         fprintf(fsum, "\nS=%d:%3d/%d  T=%4d  ", stage,itree+1,ntreet,ntree+1);
            OutaTreeN (fsum, 0, 0);
         }
         if (fout) {
           fprintf(fout,"\nS=%d:%4d/%4d  T=%4d ",stage,itree+1,ntreet,ntree+1);
           OutaTreeN (fout, 0, 0);
         }
         status=ming1(NULL,&tree.lnL,com.plfun,NULL,testx,x,space,e,com.np);

         if (tree.lnL<treebest.tree.lnL) {
            treebest.tree=tree;  memcpy (treebest.nodes, nodes, sizeof(nodes));
            FOR (i, com.np) treebest.x[i]=x[i];
            best=itree+1;   improve=1;
         }
         if (noisy) printf("%6d%2c", NFunCall, (status?'?':'X'));
         if (fsum) {
            fprintf(fsum, "%6d%2c", NFunCall, (status?'?':'X'));
            for (i=com.ntime; i<com.np; i++)  fprintf(fsum, "%7.3f", x[i]);
            fflush(fsum);
	 }
         if (fout) {
            fprintf(fout,"\nlnL(ntime:%3d  np:%3d):%14.6f\n",
                         com.ntime, com.np, -tree.lnL);
            OutaTreeB (fout);   FPN(fout);
            FOR (i,com.np) fprintf(fout,"%9.5f", x[i]); 
            FPN(fout); fflush(fout);
         }

      }  /* for (itree) */
      son1=treebest.nodes[tree.nnode-1].sons[0];
      son2=treebest.nodes[tree.nnode-1].sons[1];
   }    /* for (stage) */

   if (com.ns<=4) return (best);
   else return (0);
}

int DecompTree (int inode, int ison1, int ison2)
{
/* decompose treestar at NODE inode into tree and nodes[]
*/
   int i, son1, son2;
   double bt, fmid=0.001, fclock=0.0001;

   tree=treestar.tree;  memcpy (nodes, treestar.nodes, sizeof(nodes));
   for (i=0,bt=0; i<tree.nnode; i++)
      if (i!=tree.origin) bt+=nodes[i].branch/tree.nbranch;

   nodes[tree.nnode].nson=2;
   nodes[tree.nnode].sons[0]=son1=nodes[inode].sons[ison1];
   nodes[tree.nnode].sons[1]=son2=nodes[inode].sons[ison2];
   nodes[tree.nnode].father=inode;
   nodes[son1].father=nodes[son2].father=tree.nnode;

   nodes[inode].sons[ison1]=tree.nnode;
   for (i=ison2; i<nodes[inode].nson; i++)
      nodes[inode].sons[i]=nodes[inode].sons[i+1];
   nodes[inode].nson--;

   tree.nnode++;
   NodeToBranch();
   if (!com.clock)
      nodes[tree.nnode-1].branch=bt*fmid;
   else
      nodes[tree.nnode-1].divtime=nodes[inode].divtime*(1-fclock);

   return(0);
}


#ifdef REALSEQUENCE


int MultipleGenes (FILE* fout, double space[])
{
/*   FILE *fseqtmp;
   char seqfs[20]="seq .tmp";  */
   int ig, j, ngene0, npatt0, lgene0[NGENE], posG0[NGENE+1];

   ngene0=com.ngene;  npatt0=com.npatt;
   FOR (ig, ngene0)   lgene0[ig]=com.lgene[ig];
   FOR (ig, ngene0+1) posG0[ig]=com.posG[ig];
   for (ig=0; ig<ngene0; ig++) {
      com.ngene=1; 
      com.ls=com.lgene[0]= ig==0?lgene0[0]:lgene0[ig]-lgene0[ig-1];
      com.npatt =  ig==ngene0-1 ? npatt0-posG0[ig] : posG0[ig+1]-posG0[ig];
      com.posG[0]=0;  com.posG[1]=com.npatt;
      FOR (j,com.ns) com.z[j]+=posG0[ig];   com.fpatt+=posG0[ig];
      xtoy (com.piG[ig], com.pi, com.ncode+(com.seqtype==NUCseq?2:0));

      printf ("\n\nGene # %d  ls: %d  npatt: %d\n",ig+1,com.ls,com.npatt);
      fprintf(fout,"\nGene # %d  ls: %d  npatt: %d\n",ig+1,com.ls,com.npatt);
      fprintf(frst,"\n\nGene # %d  ls: %d  npatt: %d\n",ig+1,com.ls,com.npatt);
      if (com.runmode==0)  Forestry (fout, space);
      else                 TreeSearch (fout, space);
      FOR (j,com.ns) com.z[j]-=posG0[ig];   com.fpatt-=posG0[ig];
#if 0
sprintf (seqfs, "seq%d.tmp", ig+1);
printf ("\nprinting sequence data for gene %d into %s\ndo this now!\n", ig+1,seqfs);
fseqtmp=fopen(seqfs, "w");
/* OutSeqData (fseqtmp, 0); */
fclose (fseqtmp);
#endif
   }
   return (0);
}


#endif   /* ifdef REALSEQUENCE */
#endif   /* ifdef TREESEARCH */
#endif   /* ifdef NODESTRUCTURE */


/***************************/
#ifdef PARSIMONY

void UpPassScoreOnly (int inode);
void UpPassScoreOnlyB (int inode);

static int *Nsteps, *chUB;   /* MM */
static char *Kspace, *chU, *NchU; 
/* Elements of chU are character states (there are NchU of them).  This 
   representation is used to speed up calculation for large trees.
   Bit operations on chUB are performed for binary trees
*/

void UpPassScoreOnly (int inode)
{
/* => VU, VL, & MM, theorem 2 */
   int ison, i, j;
   char *K=Kspace, maxK;  /* chMark (VV) not used in up pass */

   FOR (i,nodes[inode].nson)
      if (nodes[nodes[inode].sons[i]].nson>0)
          UpPassScoreOnly (nodes[inode].sons[i]);

   FOR (i,com.ncode) K[i]=0;
   FOR (i,nodes[inode].nson) 
      for (j=0,ison=nodes[inode].sons[i]; j<NchU[ison]; j++)
         K[chU[ison*com.ncode+j]]++;
   for (i=0,maxK=0; i<com.ncode; i++)  if (K[i]>maxK) maxK=K[i];
   for (i=0,NchU[inode]=0; i<com.ncode; i++)
      if (K[i]==maxK)  chU[inode*com.ncode+NchU[inode]++]=i;
   Nsteps[inode]=nodes[inode].nson-maxK;
   FOR (i, nodes[inode].nson)  Nsteps[inode]+=Nsteps[nodes[inode].sons[i]];
}

void UpPassScoreOnlyB (int inode)
{
/* uses bit operation, for binary trees only 
*/
   int ison1,ison2, i, change=0;

   FOR (i,nodes[inode].nson)
      if (nodes[nodes[inode].sons[i]].nson>0)
          UpPassScoreOnlyB (nodes[inode].sons[i]);

   ison1=nodes[inode].sons[0];  ison2=nodes[inode].sons[1];
   if ((chUB[inode]=(chUB[ison1] & chUB[ison2]))==0)
      { chUB[inode]=(chUB[ison1] | chUB[ison2]);  change=1; }
   Nsteps[inode]=change+Nsteps[ison1]+Nsteps[ison2];
}


int MPScore (double space[])
{
/* calculates MP score for a given tree using Hartigan's (1973) algorithm.
   sizeof(space) = nnode*sizeof(int)+(nnode+2)*ncode*sizeof(char).
   Uses Nsteps[nnode], chU[nnode*ncode], NchU[nnode].
*/
   int score, h, i, bintree,U[3],change;

   Nsteps=(int*)space;
   bintree=(tree.nnode==2*com.ns-1 - (nodes[tree.origin].nson==3));
   if (bintree)  chUB=Nsteps+tree.nnode;
   else {
      chU=(char*)(Nsteps+tree.nnode);
      NchU=chU+tree.nnode*com.ncode;  Kspace=NchU+tree.nnode;
   }
/*
printf (" |bintree=%d| ", bintree);
*/
   for (h=0,score=0; h<com.npatt; h++) {
      FOR (i,tree.nnode) Nsteps[i]=0;
      if (bintree) {           /* binary trees, use bit operation */
         FOR (i,com.ns)  chUB[i]=1<<(com.z[i][h]-1);
         UpPassScoreOnlyB (tree.origin);
         if (nodes[tree.origin].nson>2) {
            FOR (i,3) U[i]=chUB[nodes[tree.origin].sons[i]];
            change=2;
            if (U[0]&U[1]&U[2]) change=0;
            else if (U[0]&U[1] || U[1]&U[2] || U[0]&U[2]) change=1;
            for (i=0,Nsteps[tree.origin]=change; i<3; i++) 
               Nsteps[tree.origin]+=Nsteps[nodes[tree.origin].sons[i]];
	 }
      }
      else {                   /* polytomies, use characters */
         FOR (i,com.ns) {chU[i*com.ncode]=com.z[i][h]-1; NchU[i]=1; }
         for (i=com.ns; i<tree.nnode; i++)  NchU[i]=0;
         UpPassScoreOnly (tree.origin);
      }
      score+=Nsteps[tree.origin]*com.fpatt[h];
   }
   return (score);
}

int RemoveMPNinfSites (int *nsiteNinf)
{
/* Removes parsimony-noninformative sites and return the number of changes 
   at those sites.
   Changes .z[], .fpatt[], .npatt, etc.
*/   
   int  h,j, it, npatt0=com.npatt, markb[NCODE], nchange, identic;

   for (h=0,com.npatt=nchange=0,*nsiteNinf=0; h<npatt0; h++) {
      FOR (j, com.ncode) markb[j]=0;
      FOR (j, com.ns)  markb[com.z[j][h]-1]++;
      for (j=0,it=identic=0; j<com.ncode; j++) {
         if (markb[j]==com.ns) { identic=1; break; }
         if (markb[j]>=2) it++;
      }
      if (identic || it<2) { /* non-informative */
         FOR (j,com.ncode) if(markb[j]==1) nchange+=com.fpatt[h];
	 *nsiteNinf+=com.fpatt[h];
      }
      else {
         FOR (j, com.ns) com.z[j][com.npatt]=com.z[j][h];
         com.fpatt[com.npatt++]=com.fpatt[h];
      }
   }
   return (nchange);
}


#endif


#ifdef RECONSTRUCTION

static char *chMark, *chMarkU, *chMarkL; /* VV, VU, VL */
/* chMark, chMarkU, chMarkL (VV, VU, VL) have elements 0 or 1, marking
   whether the character state is present in the set */
static char *PATHWay, *NCharaCur, *ICharaCur, *CharaCur;
/* PATHWay, NCharaCur, ICharaCur, CharaCur are for the current 
   reconstruction.  
*/

int UpPass (int inode);
int DownPass (int inode);

int UpPass (int inode)
{
/* => VU, VL, & MM, theorem 2 */
   int n=com.ncode, i, j;
   char *K=chMark, maxK;   /* chMark (VV) not used in up pass */

   FOR (i,nodes[inode].nson)
      if (nodes[nodes[inode].sons[i]].nson>0) UpPass (nodes[inode].sons[i]);

   FOR (i, n) K[i]=0;
   FOR (i,nodes[inode].nson) 
      FOR (j, n)  if(chMarkU[nodes[inode].sons[i]*n+j]) K[j]++;
   for (i=0,maxK=0; i<n; i++)  if (K[i]>maxK) maxK=K[i];
   for (i=0; i<n; i++) {
      if (K[i]==maxK)         chMarkU[inode*n+i]=1; 
      else if (K[i]==maxK-1)  chMarkL[inode*n+i]=1;
   }
   Nsteps[inode]=nodes[inode].nson-maxK;
   FOR (i, nodes[inode].nson)  Nsteps[inode]+=Nsteps[nodes[inode].sons[i]];
   return (0);
}

int DownPass (int inode)
{
/* VU, VL => VV, theorem 3 */
   int n=com.ncode, i, j, ison;

   FOR (i,nodes[inode].nson) {
      ison=nodes[inode].sons[i];
      FOR (j,n) if (chMark[inode*n+j]>chMarkU[ison*n+j]) break;
      if (j==n) 
         FOR (j,n) chMark[ison*n+j]=chMark[inode*n+j];
      else 
         FOR (j,n)
            chMark[ison*n+j] = 
               chMarkU[ison*n+j]||(chMark[inode*n+j]&&chMarkL[ison*n+j]);
   }
   FOR (i,nodes[inode].nson)
      if (nodes[nodes[inode].sons[i]].nson>0) DownPass (nodes[inode].sons[i]);
   return (0);
}


int DownStates (int inode)
{
/* VU, VL => NCharaCur, CharaCur, theorem 4 */
   int i;

   FOR (i,nodes[inode].nson) 
      if (nodes[inode].sons[i]>=com.ns) 
         DownStatesOneNode (nodes[inode].sons[i], inode);
   return (0);
}

int DownStatesOneNode (int ison, int father)
{
/* States down inode, given father */
   char chi=PATHWay[father-com.ns];
   int n=com.ncode, j, in;

   if((in=ison-com.ns)<0) return (0);
   if (chMarkU[ison*n+chi]) {
      NCharaCur[in]=1;   CharaCur[in*n+0]=chi;
   }
   else if (chMarkL[ison*n+chi]) {
      for (j=0,NCharaCur[in]=0; j<n; j++) 
         if (chMarkU[ison*n+j] || j==chi) CharaCur[in*n+NCharaCur[in]++]=j;
   }
   else {
      for (j=0,NCharaCur[in]=0; j<n; j++) 
         if (chMarkU[ison*n+j]) CharaCur[in*n+NCharaCur[in]++]=j;
   }
   PATHWay[in]=CharaCur[in*n+(ICharaCur[in]=0)];
   FOR (j, nodes[ison].nson)  if (nodes[ison].sons[j]>=com.ns) break;
   if (j<nodes[ison].nson) DownStates (ison);

   return (0);
}

int InteriorStatesMP (int job, int h, int *nchange, char NChara[NS-1], 
    char Chara[(NS-1)*NCODE], double space[]);

int InteriorStatesMP (int job, int h, int *nchange, char NChara[NS-1], 
    char Chara[(NS-1)*NCODE], double space[])
{
/* sizeof(space) = nnode*sizeof(int)+3*nnode*ncode*sizeof(char)
   job: 0=# of changes; 1:equivocal states
*/
   int n=com.ncode, i,j;

   Nsteps=(int*)space;            chMark=(char*)(Nsteps+tree.nnode);
   chMarkU=chMark+tree.nnode*n;   chMarkL=chMarkU+tree.nnode*n;
   FOR (i,tree.nnode) Nsteps[i]=0;
   memset (chMark, 0, 3*n*tree.nnode*sizeof(char));
   FOR (i,com.ns)  chMark[i*n+com.z[i][h]-1]=chMarkU[i*n+com.z[i][h]-1]=1;
   UpPass (tree.origin);
   *nchange=Nsteps[tree.origin];
   if (job==0) return (0);
   FOR (i,n) chMark[tree.origin*n+i]=chMarkU[tree.origin*n+i];
   DownPass (tree.origin);
   FOR (i,tree.nnode-com.ns) 
      for (j=0,NChara[i]=0; j<n; j++) 
         if (chMark[(i+com.ns)*n+j])  Chara[i*n+NChara[i]++] = j;
   return (0);     
}

int PathwayMP (FILE *fout, double space[])
{
/* Hartigan, JA.  1973.  Minimum mutation fits to a given tree. 
   Biometrics, 29:53-65.
*/
   char *pch=(com.seqtype==0?NUCs:AAs), visit[NS-1];
   int n=com.ncode, nid=tree.nbranch-com.ns+1, it, i,j,k, h, npath;
   int nchange, nchange0;
   char nodeb[NNODE], Equivoc[NS-1];

   PATHWay=(char*)malloc(nid*(n+3)*sizeof(char));
   NCharaCur=PATHWay+nid;  ICharaCur=NCharaCur+nid;  CharaCur=ICharaCur+nid;

   for (j=0,visit[i=0]=tree.origin-com.ns; j<tree.nbranch; j++) 
     if (tree.branches[j][1]>=com.ns) 
        visit[++i]=tree.branches[j][1]-com.ns;
/*
   printf ("\nOrder in nodes: ");
   FOR (j, nid) printf ("%4d", visit[j]+1+com.ns); FPN(F0);
*/
   for (h=0; h<com.npatt; h++) {
      fprintf (fout, "\n%4d%6d  ", h+1, com.fpatt[h]);
      FOR (j, com.ns) fprintf (fout, "%c", pch[com.z[j][h]-1]);
      fprintf (fout, ":  ");

      FOR (j,com.ns) nodeb[j]=com.z[j][h]-1;

      InteriorStatesMP (1, h, &nchange, NCharaCur, CharaCur, space); 
      ICharaCur[j=tree.origin-com.ns]=0;  PATHWay[j]=CharaCur[j*n+0];
      FOR (j,nid) Equivoc[j]=(NCharaCur[j]>1);
      DownStates (tree.origin);
/*
printf ("\n%4d%6d  ", h+1, com.fpatt[h]);
FOR (j, com.ns) printf ("%c", pch[com.z[j][h]-1]);
printf ("  ");
FOR (j,nid) {
   printf ("\nnode %2d: %4d  ", j+com.ns+1, NCharaCur[j]);
   FOR (i,n) if (chMark[(j+com.ns)*n+i]) printf ("%c", pch[i]);
   printf ("  ");
   FOR (i,n) if (chMarkU[(j+com.ns)*n+i]) printf ("%c", pch[i]);
   printf ("  \t\t  ");
   FOR (i,NCharaCur[j]) printf ("%c", pch[CharaCur[j*n+i]]);
}
*/
      for (npath=0; ;) {
         for (j=0,k=visit[nid-1]; j<NCharaCur[k]; j++) {
            PATHWay[k]=CharaCur[k*n+j]; npath++; 
            FOR (i, nid) fprintf (fout, "%c", pch[PATHWay[i]]);
            fprintf (fout, "  ");

            FOR (i,nid) nodeb[i+com.ns]=PATHWay[i];
            for (i=0,nchange0=0; i<tree.nbranch; i++) 
            nchange0+=(nodeb[tree.branches[i][0]]!=nodeb[tree.branches[i][1]]);
            if (nchange0!=nchange) 
               { puts("\a\nerr:PathwayMP"); fprintf(fout,".%d. ", nchange0); }

         }
         for (j=nid-2; j>=0; j--) {
            if(Equivoc[k=visit[j]] == 0) continue;
            if (ICharaCur[k]+1<NCharaCur[k]) {
               PATHWay[k] = CharaCur[k*n + (++ICharaCur[k])];
               DownStates (k+com.ns);
/*
printf ("\nDownstates for node %4d", k+com.ns+1);
FOR (j,nid) {
   printf ("\nnode %2d: %4d  ", j+com.ns+1, NCharaCur[j]);
   FOR (i,n) if (chMark[(j+com.ns)*n+i]) printf ("%c", pch[i]);
   printf ("  \t\t ");
   FOR (i,NCharaCur[j]) printf ("%c", pch[CharaCur[j*n+i]]);
}
getchar();
*/
               break;
            }
            else { /* if (next equivocal node is not ancestor) update node k */
               for (i=j-1; i>=0; i--) if (Equivoc[(int)visit[i]]) break;
               if (i>=0) { 
                  for (it=k+com.ns,i=visit[i]+com.ns; ; it=nodes[it].father)
                     if (it==tree.origin || nodes[it].father==i) break;
                  if (it==tree.origin) {
/*
                     printf ("\n\a\t\t\t\t*** node %6d%6d", k+1+com.ns, i+1);
*/
                     DownStatesOneNode (k+com.ns, nodes[k+com.ns].father);
		  }
	       }
	    }
         }
         if (j<0) break;
       }
       fprintf (fout, " |%4d (%d)", npath, nchange);
   }   /* for (h) */
   free (PATHWay);
   return (0);
}

#endif

#ifdef LFUNCTIONS
#ifdef RECONSTRUCTION



int AncestralSeqs (FILE *fout, double x[], double space[])
{
/* fhs[]: fh[] for each site pattern
   pnode[nid*npatt]: MAP prob at each node for all patterns
   pnsite[nid][n]: prob for each chara for all nodes at one site
*/
   char *pch=(com.seqtype==0?NUCs:(com.seqtype==2?AAs:BINs));
   char NChara[NS-1], Chara[(NS-1)*NCODE], *zz[NNODE];
   int n=com.ncode, i,j,h,ig, npathMP,nid=tree.nbranch-com.ns+1,it,ipath,npath;
   int nodeb[NNODE], markb[NCODE], icat, nchange,nchange0;
   double lnL=0, *fhs, y, *Ptb, pmin=0.01, ppath, ppathMP,pMAP[3],pMP[3],fs[3];
   double *pnode, pnsite[(NS-1)*NCODE], pMAPnode[NS-1], pMAPnode2[NS-1];
#if 0   
   double smallp=.2;
#endif
   if (com.alfa) error ("gamma rates for AncestralSeqs");
   fhs=(double*)malloc(com.npatt*sizeof(double));
   Ptb=(double*)malloc(tree.nbranch*n*n*sizeof(double));
   pnode=(double*)malloc((nid*com.npatt+1)*(sizeof(double)+sizeof(char)));
   FOR (j,nid) zz[com.ns+j]=(char*)(pnode+nid*com.npatt)+j*com.npatt;
   FOR (j,com.ns) zz[j]=com.z[j];
   if (fhs==NULL || Ptb==NULL || pnode==NULL) error ("oom");

   if (noisy) printf ("\nReconstructing ancestral states\n");
   fprintf (fout, "\nReconstructing ancestral states for the tree..\n");
   OutaTreeN (fout, 0, 0);  
   fprintf (fout, "\nThis can also be represented as \n");
   OutaTreeB (fout);  
   fprintf (fout, "\nwhich you use to deduce ancestral nodes %d to %d. ", 
      com.ns+1, tree.nnode);
   fprintf (fout, "\nSome results are shown for site patterns only, i.e.,\n");
   fprintf (fout, "sites with the same data, and a table maps sites to patterns.\n");
   fprintf (fout, "\n(1) # & freq. of site pattern, reconstruction (prob.)\n");
   fprintf (fout,
      "parsimony reconstructions marked with *, |# paths(# changes)\n");

   if (SetParameters (x)) puts ("par err.");
   /* calculate fh[] for sites */
   for (ig=0; ig<com.ngene; ig++) {
      if (com.Mgene>1) SetPGene(ig, 1, 1, 0, x);
      PartialLikelihood (tree.origin, ig);
      for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
         for (i=0,fhs[h]=0; i<com.ncode; i++) 
            fhs[h] += com.pi[i]*nodes[tree.origin].lkl[h*n+i];
         if (fhs[h]<=0) {
            matout (F0, x, 1, com.np);
            printf ("\a\nAncestralSeqs: h=%4d  fh=%9.4f \n", h, fhs[h]);
         }
         lnL -= log(fhs[h])*(double)com.fpatt[h];
      }
   }

   /* reconstruct sequences */
   zero(pMAPnode,nid);  fillxc (pMAPnode2, 1., nid);
   for (h=0,ig=-1,zero(pMAP,3),zero(pMP,3),zero(fs,3); h<com.npatt; h++) {
      if (ig<com.ngene-1 && h==com.posG[ig+1]) { /* get all transP's */
         ig++;
         if (com.Mgene>1) SetPGene(ig, 1, 1, 0, x);
#ifdef CODEML
         for (j=0,y=com.rgene[ig]; j<tree.nbranch; j++) 
          PMatUVRoot(Ptb+j*n*n,nodes[tree.branches[j][1]].branch*y,n,U,V,Root);
#elif BASEML
         for (j=0,y=com.rgene[ig]; j<tree.nbranch; j++) 
           PMatCijk(Ptb+j*n*n, nodes[tree.branches[j][1]].branch*y);
#elif BINML
         PMatBranches (x, n, Ptb); 

FOR (i, tree.nbranch) {
   fprintf (fout, "\nbranch %d", i+1);
   matout (fout, Ptb+i*n*n, n, n);
}
#endif
      }
      fprintf (fout,"\n%4d%6d  ", h+1, com.fpatt[h]);
      FOR (j, com.ns) fprintf (fout, "%c", pch[com.z[j][h]-1]);
      fprintf (fout, ":  ");

      FOR (j,nid*n) pnsite[j]=0;
      FOR (j, n) markb[j]=0;  /* icat: variable, informative sites etc. */
      FOR (j, com.ns) { nodeb[j]=com.z[j][h]-1; markb[nodeb[j]]++; }
      for (j=0,icat=3,it=0; j<n; j++) {
         if (markb[j]==com.ns) { icat=1; break; }
         else if (markb[j]>=2)   it++;
      }
      if (icat==3 && it<2) icat=2;

      if (0 && ((com.seqtype!=AAseq && nid<15) || nid<7)) {
         InteriorStatesMP (0, h, &nchange, NULL, NULL, space); 
         FOR (j,nid) NChara[j]=0;
         FOR (j,n) if(markb[j])  FOR(it,nid) Chara[it*n+NChara[it]++]=j;
      }
      else {
         InteriorStatesMP (1, h, &nchange, NChara, Chara, space);
#ifdef BINML
         for (i=0; i<nid; i++) {
#if 0
            if (NChara[i]==n) continue;
            if (com.ns+i==tree.origin) {
               if (com.pi[(int)Chara[i*n]]>smallp) continue;
            }
            else {
               j=nodes[com.ns+i].ibranch;
               if ((Chara[i*n]==0 && Ptb[j*n*n+0*n+1]<1-smallp) ||
                   (Chara[i*n]==1 && Ptb[j*n*n+1*n+0]<1-smallp) )  continue;
	    }
#endif
            if (NChara[i]==n || Chara[i*n]==0) continue;
            NChara[i]=n; Chara[i*n]=0; Chara[i*n+1]=1; printf("\nmodi.");
	 }
#endif
      }
      for (j=0,npath=1; j<nid; j++)  npath*=NChara[j]; 

      for (ipath=0,ppath=ppathMP=0,npathMP=0; ipath<npath; ipath++) {
         for (j=nid-1,it=ipath; j>=0; it/=NChara[j--])
            nodeb[com.ns+j]=Chara[j*n+it%NChara[j]];
         for (j=0,nchange0=0; j<tree.nbranch; j++) 
            nchange0+=(nodeb[tree.branches[j][0]]!=nodeb[tree.branches[j][1]]);
         if (ipath && nchange0>nchange+4) continue;  

         for (j=0,y=com.pi[nodeb[tree.origin]]; j<tree.nbranch; j++) 
         y*=Ptb[j*n*n+nodeb[tree.branches[j][0]]*n+nodeb[tree.branches[j][1]]];

         FOR (j,nid) pnsite[j*n+nodeb[j+com.ns]]+=y;
         if (nchange0==nchange)  { npathMP++; ppathMP+=y; }
         if (y>ppath) ppath=y;
         if ((y/fhs[h]>pmin /* && y/fhs[h]<1-pmin */ )|| nchange0==nchange) {
            FOR (j, nid)  fprintf (fout, "%c", pch[nodeb[com.ns+j]]);
            fprintf (fout, nchange0==nchange?"* ":"  ");
            if (y/fhs[h]>pmin) fprintf (fout, "(%.3f) ", y/fhs[h]);
         }
/*
         if (noisy && (ipath==0 || ipath==npath-1 || (ipath+1)%1000==0))
           printf("pattern%4d/%4d\tpath%9d/%9d\n",h+1,com.npatt,ipath+1,npath);
*/
      }  /* for (ipath) */
      fprintf (fout, "|%3d (%d)", npathMP, nchange);
      for (j=0,y=(double)com.fpatt[h]/com.ls/fhs[h]; j<icat; j++) {
         pMAP[j]+=y*ppath;   pMP[j]+=y*ppathMP/npathMP;
         fs[j]+=(double)com.fpatt[h]/com.ls;
      }
      FOR (j, nid) {
         for (i=1,it=0; i<n; i++) if (pnsite[j*n+i]>pnsite[j*n+it]) it=i;
         zz[com.ns+j][h]=it+1;   pnode[j*com.npatt+h]=pnsite[j*n+it]/fhs[h];
         pMAPnode[j]+=com.fpatt[h]*pnode[j*com.npatt+h]/com.ls;
         pMAPnode2[j]*=pow(pnode[j*com.npatt+h], (double)com.fpatt[h]);

      }
   }  /* for (h) */
   if (noisy)  printf("\n");
   fprintf (fout,"\n\nlnL =%12.6f", -lnL);
   FOR (j, 2) {
      fprintf(fout,"\naccuracy at all, variable and informative sites (%s):\n",
         j==0?"likelihood":"parsimony");
      FOR (it, 3)  fprintf (fout,"%9.5f", (j==0?pMAP[it]:pMP[it])/fs[it]);
   }
   fprintf (fout,"\n\n(2) Accuracy at each node for each site pattern\n");
   FOR (h, com.npatt) {
      fprintf (fout,"\n%4d%6d  ", h+1, com.fpatt[h]);
      FOR (j, com.ns) fprintf (fout, "%c", pch[com.z[j][h]-1]);
      fprintf (fout, ":  ");
      FOR (j, nid) fprintf (fout, "%c", pch[zz[com.ns+j][h]-1]); 
      fprintf (fout, ":  ");
      FOR (j, nid) 
       fprintf(fout,"%c (%5.3f) ",pch[zz[com.ns+j][h]-1],pnode[j*com.npatt+h]);
   }
   Site2Pattern (fout);
   fprintf (fout,"\n\n(3) list of extant and reconstructed sequences");
   FOR (j,com.ns) fprintf (fout,"\n%s", com.spname[j]);
   FOR (j,tree.nnode-com.ns) fprintf (fout,"\nnode #%d", j+com.ns+1);
   printsmaPose (fout, zz, tree.nnode, com.ls, 60,10,1,0,com.seqtype,com.pose);
   fprintf(fout,"\n\nOverall accuracy of the %d ancestral sequences:", nid);
   matout2 (fout, pMAPnode, 1, nid, 9, 5);  fputs("for a site.\n",fout);
   matout2 (fout, pMAPnode2, 1, nid, 9, 5); fputs("for the sequence.\n", fout);

   free (fhs);  free (Ptb);  free (pnode);

   return (0);
}

#endif

int lfunAdG_rate (FILE* fout, double x[], int np)
{
/* for dG, AdG or similar non-parametric models
   fhK[<npatt] stores rates for conditional mean (re), and 
   fhK[<2*npatt] stores the most probable rate category number
*/
   char *pch=(com.seqtype==NUCseq?NUCs:AAs);
   int  ig, h, ir, il, it=0, i,j, direction=-1;
   double lnL, fh, t=0, re, mre, vre, b1[NCATG], b2[NCATG];

   if (SetParameters (x)) puts ("par err. lfunAdG_rate");
   fprintf (fout,"\nFrequencies and rates for categories (K=%d)", com.ncatG);
   matout (fout, com.freqK, 1, com.ncatG);
   matout (fout, com.rK, 1, com.ncatG);
   if (com.rho) {
      fprintf (fout, "\nTransition prob matrix over sites");
      matout2 (fout, com.MK, com.ncatG, com.ncatG, 8, 4);
   }
   zero (com.fhK, com.npatt*com.ncatG);
   for (ig=0; ig<com.ngene; ig++) {
     if (com.Mgene>1) SetPGene(ig, 1, 1, com.nalfa>1, x);
      for (ir=0; ir<com.ncatG; ir++) {
         FOR (i, tree.nnode)
            nodes[i].branch *= (ir==0?com.rK[ir]:com.rK[ir]/com.rK[ir-1]);
         PartialLikelihood (tree.origin, ig);
         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            for (i=0,fh=0; i<com.ncode; i++)
               fh += com.pi[i]*nodes[tree.origin].lkl[h*com.ncode+i];
            com.fhK[ir*com.npatt+h] += fh;
         }
      }
      FOR (i, tree.nnode)  nodes[i].branch/=com.rK[com.ncatG-1];
   }
   if (com.rho==0) {     /* dG model */
      for (h=0,lnL=0,mre=vre=0; h<com.npatt; h++) {
         for (ir=0,fh=0,it=0,t=re=0; ir<com.ncatG; ir++) {
            if (com.fhK[ir*com.npatt+h]>t)
               { t=com.fhK[ir*com.npatt+h]; it=ir; }
            fh+=com.freqK[ir]*com.fhK[ir*com.npatt+h];
            re+=com.freqK[ir]*com.fhK[ir*com.npatt+h]*com.rK[ir];
         }
         lnL -= com.fpatt[h]*log(fh);

         re/=fh;
         mre+=com.fpatt[h]*re/com.ls;     vre+=com.fpatt[h]*re*re/com.ls;
         com.fhK[h]=re;                   com.fhK[com.npatt+h]=(double)it;

         fprintf (fout,"\n%4d%6d%14.8f%9.2f%9.4f%6d  ",
                 h+1, com.fpatt[h], log(fh), com.ls*fh, re, it+1);
         if (com.seqtype!=CODONseq) 
            FOR (j,com.ns) fprintf (fout, "%c", pch[com.z[j][h]-1]);
      }
      vre-=mre*mre;
      fprintf (fout, "\n\nrates along the sequence\n");
      for (h=0; h<com.ls; h++) {
         it = com.pose[h];
        fprintf(fout,"\n%4d%9.4f%6.0f ",h+1,com.fhK[it],com.fhK[com.npatt+it]);
         if (com.seqtype!=CODONseq)  FOR (j,com.ns)
               fprintf (fout,"%c", pch[com.z[j][it]-1]);
      }
   }
   else {      /* Auto-dGamma model */
      fprintf (fout, "\n\nrates along the sequence\n");
      h = (direction==1?com.ls-1:0);
      for (il=0,lnL=0,mre=vre=0; il<com.ls; h-=direction,il++) {
         if (il==0)
            FOR(ir,com.ncatG) b1[ir]=com.fhK[ir*com.npatt+com.pose[h]];
         else {
            for (ir=0; ir<com.ncatG; ir++) {
               for (j=0,fh=0; j<com.ncatG; j++)
                  fh+=com.MK[ir*com.ncatG+j]*b1[j];
               b2[ir] = fh*com.fhK[ir*com.npatt+com.pose[h]];
            }
            xtoy (b2, b1, com.ncatG);
         }
         if (il%5==4)
            { t=sum(b1,com.ncatG); abyx(1/t,b1,com.ncatG); lnL-=log(t); }
         for (ir=0,it=0,re=fh=t=0; ir<com.ncatG; ir++) {
            re+=com.freqK[ir]*b1[ir]*com.rK[ir];
            fh+=com.freqK[ir]*b1[ir];
            if (b1[ir]>t) {it=ir; t=b1[ir]; }
         }
         re/=fh;  mre+=re/com.ls;   vre+=re*re/com.ls;

         fprintf (fout,"\n%4d%9.4f%6d ", h+1, re, it+1);
         if (com.seqtype!=CODONseq)   FOR (j,com.ns)
            fprintf (fout, "%c", pch[com.z[j][com.pose[h]]-1]);
      }
      vre-=mre*mre;
      for (ir=0,fh=0; ir<com.ncatG; ir++)  fh += com.freqK[ir]*b1[ir];
      lnL -= log(fh);
   }
   fprintf (fout,"\n\nlnL=%14.6f", -lnL);
   fprintf (fout,"\nmean(r^)=%9.4f  var(r^)=%9.4f", mre, vre);
   fprintf (fout,"\nAccuracy of rate prediction: corr(r^,r) =%9.4f", 
      sqrt(com.alfa*vre));

   FPN(fout);
   return (0);
}

double lfunAdG (double x[], int np)
{
/* Auto-Discrete-Gamma rates for sites
*/
   int  nscale=1, h, ir,ig, i,j;
   int  direction=-1;  /* 1: n->1;  -1: 1->n */
   double lnL, b1[NCATG], b2[NCATG], fh;

   NFunCall++;
   if (SetParameters (x)) puts ("par err. AdG");
   FOR (i, com.npatt*com.ncatG) com.fhK[i]=0;
   for (ig=0; ig<com.ngene; ig++) {
      SetPGene(ig, com.Mgene>1, com.Mgene>1, com.nalfa>1, x);
      for (ir=0; ir<com.ncatG; ir++) {
         FOR (i, tree.nnode)
            nodes[i].branch *= (ir==0?com.rK[ir]:com.rK[ir]/com.rK[ir-1]);
         PartialLikelihood (tree.origin, ig);
         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            for (i=0,fh=0; i<com.ncode; i++)
               fh += com.pi[i]*nodes[tree.origin].lkl[h*com.ncode+i];
            com.fhK[ir*com.npatt+h] += fh;
         }
      }
      FOR (i, tree.nnode)  nodes[i].branch/=com.rK[com.ncatG-1];
   }
   h = (direction==1?com.ls-1:0);
   for (ir=0,lnL=0; ir<com.ls; h-=direction,ir++) {
      if (ir==0)
         FOR(i,com.ncatG) b1[i]=com.fhK[i*com.npatt+com.pose[h]];
      else {
         for (i=0; i<com.ncatG; i++) {
            for (j=0,fh=0; j<com.ncatG; j++)
               fh+=com.MK[i*com.ncatG+j]*b1[j];
            b2[i] = fh*com.fhK[i*com.npatt+com.pose[h]];
         }
         xtoy (b2, b1, com.ncatG);
      }
      if((ir+1)%nscale==0) {
         fh=sum(b1,com.ncatG);
         if (fh<1e-40) printf ("h,fh%6d%12.4e\n", h, fh);
         abyx(1/fh,b1,com.ncatG); lnL-=log(fh);
      }
   }
   for (i=0,fh=0; i<com.ncatG; i++)  fh += com.freqK[i]*b1[i];
   lnL-=log(fh);
   return (lnL);
}

double lfundG (double x[], int np)
{
/* discrete gamma rates for sites
*/
   char *pch=(com.seqtype==0?NUCs: (com.seqtype==2?AAs:NULL) );
   int  h, ir, i, ig;
   double lnL, fh=0;

   NFunCall++;
   if (SetParameters (x)) puts ("\npar err..");
   zero (com.fhK, com.npatt);
   for (ig=0; ig<com.ngene; ig++) {
      SetPGene(ig, com.Mgene>1, com.Mgene>1, com.nalfa>1, x);
      for (ir=0; ir<com.ncatG; ir++) {
         FOR (i, tree.nnode)
            nodes[i].branch *= (ir==0?com.rK[ir]:com.rK[ir]/com.rK[ir-1]);
         PartialLikelihood (tree.origin, ig);
         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            if (com.fpatt[h]==0) continue;
            for (i=0,fh=0; i<com.ncode; i++)
               fh += com.pi[i]*nodes[tree.origin].lkl[h*com.ncode+i];
            com.fhK[h] += com.freqK[ir]*fh;
         }
      }
      FOR (i, tree.nnode)  nodes[i].branch/=com.rK[com.ncatG-1];
   }
   for (h=0,lnL=0; h<com.npatt; h++) {
      if (com.fpatt[h]==0) continue;
      if (com.fhK[h]<=0) {
         matout (F0, x, 1, np);
         printf ("\nlfundG: h=%4d  fhK=%9.4f", h, com.fhK[h]);
      }
      lnL-=log(com.fhK[h])*(double)com.fpatt[h];

      if (com.print<0){
         fprintf (flfh, "\n%6d%6d%16.10f%9.2f  ", 
            h+1,com.fpatt[h],log(com.fhK[h]),com.ls*com.fhK[h]);
            if (pch) FOR (i,com.ns) fprintf (flfh, "%c", pch[com.z[i][h]-1]);
      }
   }
   return (lnL);
}

double lfun (double x[], int np)
{
   char *pch=(com.seqtype==0?NUCs: (com.seqtype==2?AAs:NULL) );
   int  h, i, ig;
   double lnL=0, fh;

   NFunCall++;
   if (SetParameters (x)) puts ("\npar err..");
   for (ig=0; ig<com.ngene; ig++) {
      if (com.Mgene>1) SetPGene(ig, 1, 1, 0, x);
      PartialLikelihood (tree.origin, ig);
      for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
         if (com.fpatt[h]==0) continue;
         for (i=0,fh=0; i<com.ncode; i++) 
            fh += com.pi[i]*nodes[tree.origin].lkl[h*com.ncode+i];
         if (fh<=0) {
            matout (F0, x, 1, np);
            printf ("\a\nlfun: h=%4d  fh=%9.4f \n", h, fh);
         }
         lnL-=log(fh)*(double)com.fpatt[h];

         if (com.print<0){
            fprintf (flfh, "\n%6d%6d%16.10f%9.2f  ", 
               h+1, com.fpatt[h], log(fh), com.ls*fh);
            if (pch) FOR (i,com.ns) fprintf (flfh, "%c", pch[com.z[i][h]-1]);
	 }
      }
   }
   return (lnL);
}

#endif         /* #ifdef LFUNCTIONS */

void BranchLengthBD (double birth, double death, double sample, double mut)
{
/* Generate random branch lengths (nodes[].branch) using the birth and
   death process with species sampling
   Note: older interior node have larger node numbers, so root is 
   node com.ns*2-2
*/
   int i,j, imin,fixt0=1;
   double la=birth, mu=death, rho=sample, tmin, r, t[NS-1];
   double phi, eml;

   if (fixt0) t[com.ns-2]=1;

   if (fabs(la-mu)>1e-6) {
      eml=exp(mu-la);  phi=(rho*la*(eml-1)+(mu-la)*eml)/(eml-1);
      for (i=0; i<com.ns-1-(fixt0); i++) 
         { r=rndu(); t[i]=log((phi-r*rho*la)/(phi-r*rho*la+r*(la-mu)))/(mu-la); }
   }
   else  
      for (i=0; i<com.ns-1-(fixt0); i++) 
         { r=rndu();  t[i]=r/(1+la*rho*(1-r)); }

   /* bubble sort */
   for (i=0; i<com.ns-1-1; i++) {
      for (j=i+1,tmin=t[i],imin=i; j<com.ns-1; j++) 
         if (tmin>t[j]) { tmin=t[j]; imin=j; }
      t[imin]=t[i];  t[i]=tmin;
   }
   FOR (i,com.ns) nodes[i].divtime=0;
   for (i=com.ns; i<tree.nnode; i++) nodes[i].divtime=t[i-com.ns];
   for (i=0; i<tree.nnode; i++) 
      if (i!=tree.origin) 
         nodes[i].branch=(nodes[nodes[i].father].divtime-nodes[i].divtime)*mut;
}


#ifdef EVOLVE

int RandomLHistory (int rooted, double space[])
{
/* random coalescence tree, with each labeled history having equal probability.
   interior nodes are numbered ns, ns+1, ..., 2*ns-1-!rooted
*/
   int ns=com.ns, i, j, it=0, *nodea=(int*)space;
   double t;

   for (i=0; i<2*ns-1-!rooted; i++) ClearNode(i);
   for (i=0; i<ns; i++) nodea[i]=i;
   for (i=ns,t=0; i>(1+!rooted); i--) {
      nodes[it=2*ns-i].nson=2;
      j=i*rndu();     nodes[nodea[j]].father=it;  nodes[it].sons[0]=nodea[j];
      nodea[j]=nodea[i-1];
      j=(i-1)*rndu(); nodes[nodea[j]].father=it;  nodes[it].sons[1]=nodea[j];
      nodea[j]=it;
      if (!rooted && i==3) {
         nodes[it].nson++; 
         nodes[nodea[1-j]].father=it; nodes[it].sons[2]=nodea[1-j];
      }
   }
   tree.origin=it;  tree.nnode=ns*2-1-!rooted;
   NodeToBranch();
   return (0);
}


int Rates4Sites (double rates[],double alfa,int ncatG,int ls, int cdf,
    double space[])
{
/* Rates for sites from the gamma (ncatG=0) or discrete-gamma (ncatG>1)
   Rates are converted into the c.d.f. if cdf=1, which is useful for
   simulation under JC69-like models. 
*/
   int h, ir,j, *counts=(int*)(space+2*ncatG);
   double *rK=space, *freqK=space+ncatG;

   if (alfa==0) 
      for (h=0; h<ls; h++) rates[h]=1;
   else {
      if (ncatG>1) {
         DiscreteGamma (freqK, rK, alfa, alfa, ncatG, 0);
         MultiNomial (ls, ncatG, freqK, counts, space+3*ncatG);
         for (ir=0,h=0; ir<ncatG; ir++) 
            for (j=0; j<counts[ir]; j++)  rates[h++]=rK[ir];
      }
      else 
         for (h=0; h<ls; h++) rates[h]=rndgamma(alfa)/alfa;
      if (cdf) {
         for (h=1; h<ls; h++) rates[h]+=rates[h-1];
         abyx (1/rates[ls-1], rates, ls);
      }
   }
   return (0);
}


static double PMat[16], Qfactor, kapa1, kapa2;

int GenerateSeq (void)
{
/* makes sequence com.z[tree.origin] with com.pi[] and evolves it
   along the tree, using nodes[].branch, com.ns, com.ls, com.ncode, com.model,
   com.kapa, com.alfa, com.rates
   z[ns*2-1]
*/
   int i;
   double *p=com.pi;

   if (com.ncode!=4) error ("ncode");
   if (com.ns<2 || com.ls<1) error ("err GenerateSeq()");

   if (com.model==F84) { kapa1=1+com.kapa/p[4]; kapa2=1+com.kapa/p[5]; }
   else                  kapa1=kapa2=com.kapa;
   Qfactor=1/(2*p[0]*p[1]*kapa1+2*p[2]*p[3]*kapa2 + 2*p[4]*p[5]);
   if (com.model>K80) 
      dnamaker (com.z[tree.origin], com.ls, p);
   else {
      FOR (i,com.ls) com.z[tree.origin][i]=1+rndu()*4;
      com.z[tree.origin][i]=0;
   }
   if (com.model==JC69) EvolveJC (tree.origin);
   else                 Evolve   (tree.origin);
   com.npatt=0;
   return (0);
}

void Evolve (int inode)
{
/* evolve sequence com.z[tree.origin] along the tree to generate com.z[], 
   using nodes[].branch and  com.model
   Needs com.z[nnode], while com.z[0] -- com.z[ns-1] constitute the data.
   Nucleotides only.
*/
   int is, h,i,j, ison, ib, n=com.ncode;
   double t, r;
   
   for (is=0; is<nodes[inode].nson; is++) {
      ison=nodes[inode].sons[is];
      strcpy (com.z[ison], com.z[inode]);

      if (com.alfa==0) {
         if (com.model<=K80) PMatK80 (PMat, nodes[ison].branch, com.kapa);
         else {
            t=nodes[ison].branch*Qfactor;
            PMatTN93 (PMat, t*kapa1, t*kapa2, t, com.pi);
         }
         FOR (i, n) for(j=1; j<n; j++)  PMat[i*n+j]+=PMat[i*n+j-1];
         FOR (h, com.ls) {
            for (j=0,ib=com.z[ison][h]-1,r=rndu(); j<n; j++)
               if (r<PMat[ib*n+j]) break;
            if (j-ib) com.z[ison][h]=j+1;
         }
      }
      else {
         FOR (h, com.ls) {
            if (h==0 || com.rates[h]!=com.rates[h-1]) {
               if (com.model<=K80) 
                  PMatK80 (PMat, nodes[ison].branch*com.rates[h], com.kapa);
               else {
                  t=nodes[ison].branch*Qfactor*com.rates[h];
                  PMatTN93 (PMat, t*kapa1, t*kapa2, t, com.pi);
               }
               FOR (i, n) for(j=1; j<n; j++)  PMat[i*n+j]+=PMat[i*n+j-1];
            }
            for (j=0,ib=com.z[ison][h]-1,r=rndu(); j<n; j++)
               if (r<PMat[ib*n+j]) break;
            if (j-ib) com.z[ison][h]=j+1;
         }
      }
      if (nodes[ison].nson) Evolve (ison); 
   }  /* for (is) */
}

void EvolveJC (int inode)
{
/* Special version of Evolve that works with JC69-like (poisson) model only.
   This can be used to generate amino acid sequences also.
   For each branch in the tree, this determines the number of mutations 
   and then assigns mutations to sites at random.
   When alfa>0, com.rates[] are the accumulative probabilities.
*/
   int is,j, h, nmut, imut, ison;
   double r;

   if (com.alfa && fabs(com.rates[com.ls-1]-1)>1e-4) 
      printf ("rates c.d.f.: 1 = %.6f?\n", com.rates[com.ls-1]);
   for (is=0; is<nodes[inode].nson; is++) {
      ison=nodes[inode].sons[is];
      strcpy (com.z[ison], com.z[inode]);
      nmut = rndpoisson (nodes[ison].branch*com.ls);
      for (imut=0; imut<nmut; imut++) {
         if (com.alfa==0) h=(int)(rndu()*com.ls);
         else 
            for (h=0,r=rndu(); h<com.ls; h++) if (r<com.rates[h]) break;
         j=(int)(rndu()*(com.ncode-1))+1;
         if (j>=com.z[ison][h]) j++;
         com.z[ison][h]=j;
      }
      if (nodes[ison].nson) EvolveJC (ison); 
   }
}

int PatternWeightSimple (int CollapsJC, double space[])
{
/* Counts site patterns, changes .z[], .fpatt[], .npatt, etc.
   This is modified from PatternWeight(), with options for multiple genes 
   deleted.
*/   
   int  h, ht, b,j,k, same=0;
   char *z1[NS], zh[NS];

   FOR (j,com.ns) z1[j]=(char*)space+j*com.ls;
   FOR (h,com.ls) com.fpatt[h]=0; 

   for (h=0,com.npatt=0; h<com.ls; h++) {
      if (CollapsJC) {
         zh[0]=b=1;
         for (j=1; j<com.ns; j++) {
            for (k=0; k<j; k++) if (com.z[j][h]==com.z[k][h]) break;
            zh[j]=(k<j?zh[k]:++b);
         }
      }
      else 
         for (j=0; j<com.ns; j++) zh[j]=com.z[j][h];
      for (ht=0,same=0; ht<com.npatt; ht++) {
         for (j=0,same=1; j<com.ns; j++)
	    if (z1[j][ht]!=zh[j]) { same=0; break; }
         if (same)  break; 
      }
      if (same)   com.fpatt[ht]++;
      else {
	 FOR (j, com.ns) z1[j][com.npatt]=zh[j];
         com.fpatt[com.npatt++]=1;
      }
   }
   FOR (j,com.ns) {
      memcpy (com.z[j], z1[j], com.npatt);
      com.z[j][com.npatt]=0;
   }
   return (0);
}


void OutSeqData (FILE *fseq, int format)
{
/* format=0:PAML &PHYLIP format
          1:NEXUS, with help from John Huelsenbeck.  thank john.
*/
   int h, i,j, sites=(com.npatt==0);

   if (format==0) {  /* PAML & PHYLIP format */
      fprintf (fseq, "%8d%8d\n", com.ns, com.ls);
      for (j=0; j<com.ns; FPN(fseq),j++) {
         fprintf (fseq, "seq.%d\n", j+1);
         for (h=0; h<(sites?com.ls:com.npatt); h++) 
            FOR (i,(sites?1:com.fpatt[h])) 
               fprintf (fseq, "%c", NUCs[com.z[j][h]-1]);
      }
   }
   else if (format==1) {  /* NEXUS */
      fprintf(fseq, "\nbegin data;\n   dimensions ntax=%d nchar=%d;\n",
               com.ns,com.ls);
      fprintf (fseq, "   format datatype=dna;\n   matrix\n");

      for (j=0; j<com.ns; FPN(fseq),j++) {
         fprintf (fseq, "%6d  ", j+1);
         for (h=0; h<(sites?com.ls:com.npatt); h++) 
            FOR (i,(sites?1:com.fpatt[h])) 
               fprintf (fseq, "%c", NUCs[com.z[j][h]-1]);
      }
      fprintf (fseq, "   ;\nend;\n");

fputs("\nbegin paup;\n   set autoclose=yes warnreset=no notifybeep=no;", fseq);
fputs("\n   set crit=parsimony;\n   hsearch nocollapse swap=none;\n", fseq);
fputs("   savetrees fmt=altnex append;\nend;\n\n",fseq);
   }
   fflush (fseq);
}

#endif


int gradient2 (int n, double x[], double f0, double g[], 
    double (*fun)(double x[],int n), double space[], int xmark[]);

extern int noisy, Iround, NFunCall;
extern double SIZEp;

int gradient2 (int n, double x[], double f0, double g[], 
    double (*fun)(double x[],int n), double space[], int xmark[])
{
/* f0=fun(x) is provided.
   xmark=0: central; 1: upper; -1: down
*/
   int i,j;
   double *x0=space, *x1=space+n, eh0=1e-7, eh;

   FOR (i,n)  {
      eh=eh0*(fabs(x[i])+1);
      if (xmark[i]==0 && SIZEp<1) {   /* central */
         FOR (j, n)  x0[j]=x1[j]=x[j];
	 eh=pow(eh,.67);  x0[i]-=eh; x1[i]+=eh;
	 g[i]=((*fun)(x1,n)-(*fun)(x0,n))/(eh*2.0);
      }
      else  {                         /* forward */
	 FOR (j, n)  x1[j]=x[j];
	 if (xmark[i]) eh*=-xmark[i];
	 x1[i] += eh;
	 g[i]=((*fun)(x1,n)-f0)/eh;
      }
   }
   return(0);
}

int ming2 (FILE *fout, double *f, double (*fun)(double x[], int n),
    int (*dfun)(double x[], double *f, double dx[], int n),
    double x[], double xb[][2], double space[], double e, int n)
{
/* n-D minimization with bounds using the BFGS algorithm

   g0[n] g[n] p[n] x[n] y[n] s[n] z[n] H[n*n] tv[2*n]
   using bound()
*/
   int i,j, i1,i2,it, maxround=1000, fail=0, *xmark, *ix, nfr=n;
   double small=1.e-20;     /* small value for checking |w|=0   */
   double f0, *g0, *g, *p, *x0, *y, *s, *z, *H, *C, *tv;
   double w,v, alfa, am, h;

   f0 = *f = (*fun)(x,n);
   if (noisy>2) {
      printf ("\n\nIterating by ming2\nInitial: fx= %12.6f\nx=", f0);
      FOR (i,n) printf ("%8.4f", x[i]);   FPN (F0);
   }
   g0=space;   g=g0+n;  p=g+n;   x0=p+n;
   y=x0+n;     s=y+n;   z=s+n;   H=z+n;  C=H+n*n, tv=C+n*n;
   xmark=(int*)(tv+2*n);  ix=xmark+n;
   FOR (i,n)  { xmark[i]=0; ix[i]=i; }
   xtoy (x, x0, n);

   if (dfun)  (*dfun) (x0, &f0, g0, n);
   else       gradient2 (n, x0, f0, g0, fun, tv, xmark);

   SIZEp=0; identity (H,nfr);
   FOR (Iround, maxround) {
      for (i=0,zero(p,n); i<nfr; i++)  FOR (j,nfr)
         p[ix[i]] -= H[i*nfr+j]*g0[ix[j]];
      SIZEp=norm(p,n);
/*
matout (F0, g0, 1, n);
matout (F0, p, 1, n);
*/
      for (i=0,am=20; i<n; i++) {  /* max step length */
         if (p[i]>0 && (xb[i][1]-x0[i])/p[i]<am) am=(xb[i][1]-x0[i])/p[i];
         else if (p[i]<0 && (xb[i][0]-x0[i])/p[i]<am) am=(xb[i][0]-x0[i])/p[i];
      }

      if (Iround==0)    h=fabs(2*f0*.01/innerp(g,p,n));
      else              h=norm(s,n)/SIZEp;
      h=max(h,1e-5);    h=min(h,am/8);
      *f=f0;
      alfa = LineSearch2 (fun, f, x0, p, h, am, .001, tv, n);

      if (alfa<=0) {
	 if (fail) {
	    if (SIZEp>.1) printf("\nSIZEp:%9.4f  Iround:%5d", SIZEp, Iround+1);
	    Iround=maxround;  break;
	 }
	 else   { /* printf ("\a .. "); */ identity(H, nfr); fail=1; }
      }
      else  {
	 fail=0;
	 FOR(i,n)  x[i]=x0[i]+alfa*p[i];

	 if (fout) {
	    fprintf (fout, "\n%3d %7.4f%14.6f  x", Iround+1, SIZEp, *f);
	    FOR (i,n) fprintf (fout, "%8.5f  ", x[i]);
	    fflush (fout);
	 }
	 if (SIZEp<0.0001 && H_end (x0,x,f0,*f,e,e,n))
	    { xtoy(x,x0,n); break; }
      }

      if (dfun)  (*dfun) (x, f, g, n);
      else       gradient2 (n, x, *f, g, fun, tv, xmark);

      /* modify the working set */
      FOR (i, n) {         /* add constraints, reduce H */
         if (xmark[i]) continue;
         if (fabs(x[i]-xb[i][0])<1e-6 && -g[i]<0)  xmark[i]=-1;
         else if (fabs(x[i]-xb[i][1])<1e-6 && -g[i]>0)  xmark[i]=1;
         if (xmark[i]==0) continue;
         xtoy (H, C, nfr*nfr);
         FOR (it, nfr) if (ix[it]==i) break;
         for (i1=it; i1<nfr-1; i1++) ix[i1]=ix[i1+1];
         for (i1=0,nfr--; i1<nfr; i1++) FOR (i2,nfr)
            H[i1*nfr+i2]=C[(i1+(i1>=it))*(nfr+1) + i2+(i2>=it)];
      }
      for (i=0,it=0,w=0; i<n; i++) {  /* delete a constraint, enlarge H */
         if (xmark[i]==-1 && -g[i]>w)     { it=i; w=-g[i]; }
         else if (xmark[i]==1 && -g[i]<-w) { it=i; w=g[i]; }
      }
      if (w>10*SIZEp/nfr) {          /* *** */
         xtoy (H, C, nfr*nfr);
         FOR (i1,nfr) FOR (i2,nfr) H[i1*(nfr+1)+i2]=C[i1*nfr+i2];
         FOR (i1,nfr+1) H[i1*(nfr+1)+nfr]=H[nfr*(nfr+1)+i1]=0;
         H[(nfr+1)*(nfr+1)-1]=1;
         xmark[it]=0;   ix[nfr++]=it;
      }

      if (noisy>2) {
         printf (" |%5d:", nfr);
         FOR (i,n)  if (xmark[i]) printf ("%4d", i+1);
      }
      for (i=0,f0=*f; i<nfr; i++)
        {  y[i]=g[ix[i]]-g0[ix[i]];  s[i]=x[ix[i]]-x0[ix[i]]; }
      FOR (i,n) { g0[i]=g[i]; x0[i]=x[i]; }

      for (i=0,w=v=0.; i<nfr; i++) {
	 for (j=0,z[i]=0.; j<nfr; j++) z[i]+=H[i*nfr+j]*y[j];
	 w+=y[i]*z[i];    v+=y[i]*s[i];
      }
      if (fabs(v)<small)   { identity(H,nfr); fail=1; continue; }
      FOR (i,nfr)  FOR (j,nfr)
	 H[i*nfr+j] += ((1+w/v)*s[i]*s[j]-z[i]*s[j]-s[i]*z[j])/v;
   }    /* for (Iround,maxround)  */

   if (Iround==maxround) {
      if (fout) fprintf (fout,"\ncheck convergence!\n");
      return(-1);
   }
   return(0);
}
