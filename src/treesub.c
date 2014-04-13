/* TREESUB.c
   subroutines that operates on trees, inserted into other programs 
   such as baseml, basemlg, codeml, and pamp.
*/

#include "tools.h"

extern char BASEs[], nBASEs[], *EquateNUC[], AAs[], BINs[];
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

#ifdef  EVOLVE
#define BIRTHDEATH
#endif 

#define EqPartition(p1,p2,ns) (p1==p2||p1+p2+1==(1<<ns))

#ifdef REALSEQUENCE

int PopEmptyLines (FILE* fseq, int lline, char line[])
{
/* pop out empty lines in the sequence data file.
   returns -1 if EOF.
*/
   char *eqdel=".-?", *p;
   int i;

   for (i=0; ;i++) {
      p=fgets (line, lline, fseq);
      if (p==NULL) return(-1);
      while (*p) 
         if (*p==eqdel[0] || *p==eqdel[1] || *p==eqdel[2] || isalpha(*p)) 
            return(0);
         else p++;
   }
}

int hasbase (char *str)
{
   char *p=str, *eqdel=".-?";
   while (*p) 
      if (*p==eqdel[0] || *p==eqdel[1] || *p==eqdel[2] || isalpha(*p++)) 
         return(1);
   return(0);
}

int blankline (char *str)
{
   char *p=str;
   while (*p) if (isalnum(*p++)) return(0);
   return(1);
}


int GetSeqFileType(FILE *fseq, int *paupseq);
int IdenticalSeqs(void);

int GetSeqFileType(FILE *fseq, int *paupseq)
{
/* paupstart="begin data" and paupend="matrix" identify paup seq files.
   Modify if necessary.
*/
   int lline=1000;
   char line[1000], *paupstart="begin data",*paupend="matrix", *p;
   char *ntax="ntax",*nchar="nchar";

   if(fscanf(fseq,"%d%d",&com.ns,&com.ls)==2) { *paupseq=0; return(0); }
   *paupseq=1;
   puts("\nseq file is not paml format.  Trying to read it as a paup file.");

   for ( ; ; ) {
      if(fgets(line,lline,fseq)==NULL) error2("seq err1: EOF");
      strcase(line,0);
      if(strstr(line,paupstart)) break;
   }
   for ( ; ; ) {
      if(fgets(line,lline,fseq)==NULL) error2("seq err2: EOF");
      strcase(line,0);
      if((p=strstr(line,ntax))!=NULL) {
         while (*p != '=') { if(*p==0) error2("seq err"); p++; }
         sscanf(p+1,"%d", &com.ns);
         if((p=strstr(line,nchar))==NULL) error2("err: expect nchar");
         while (*p != '=') { if(*p==0) error2("err"); p++; }
         sscanf(p+1,"%d", &com.ls);
         break;
      } 
   }
   printf("\nns: %d\tls: %d\n", com.ns, com.ls); 
   for ( ; ; ) {
      if(fgets(line,lline,fseq)==NULL) error2("seq err1: EOF");
      strcase(line,0);
      if (strstr(line,paupend)) break;
   }
   return(0);
}

int ReadSeq (FILE *fout, FILE *fseq)
{
/* read in sequence, translate into protein (CODON2AAseq), and 
   This counts ngene but does not initialize lgene[].
   com.seqtype: 0=nucleotides; 1=codons; 2:AAs; 3:CODON2AAs; 4:BINs
   com.pose[] used to store gene marks.  ls/3 gene marks for codon sequences
   char opt_c[]="AKGI";
      A:alpha given. K:kappa given
      G:many genes,  I:interlaved format
   lcode: # of characters in the final sequence 
*/
   char *p,*p1, eq='.', *line;
   int i,j,k, ch, noptline=0, lspname=LSPNAME, miss=0, paupseq;
   int lline=10000,lt[NS], igroup, lcode, Sequential=1,basecoding=0;
   int n31=(com.seqtype==CODONseq||com.seqtype==CODON2AAseq?3:1);
   int gap=(n31==3?3:10), nchar=(com.seqtype==AAseq?20:4);
   char *pch=((com.seqtype<=1||com.seqtype==CODON2AAseq)?BASEs:(com.seqtype==2?AAs:BINs));

   if (com.seqtype==4) error2("seqtype==BINs, check with author");
   if (noisy>=9 && (com.seqtype<=CODONseq||com.seqtype==CODON2AAseq)) {
      puts("\n\nAmbiguity character definition table:\n");
      FOR (i,(int)strlen(BASEs)) {
         printf("%c (%d): ", BASEs[i],nBASEs[i]);
         FOR(j,nBASEs[i])  printf("%c ", EquateNUC[i][j]);
         FPN(F0);
      }
   }
   GetSeqFileType(fseq, &paupseq);

   if (com.ns>NS) error2 ("too many sequences.. raise NS?");
   if (com.ls%n31!=0) {
      printf ("\n%d nucleotides", com.ls); error2 ("seq. len. err.");
   }
   lcode=com.ls/n31;
   if (noisy) printf ("\nns = %d  \tls = %d\n", com.ns, com.ls);

   FOR(j,com.ns) {
      if (com.spname[j]) free(com.spname[j]);
      if (com.z[j]) free(com.z[j]);
      com.spname[j]=(char*)malloc((LSPNAME+1)*sizeof(char));
      com.z[j]=(char*)malloc(com.ls*sizeof(char));
      FOR(i,lspname) com.spname[j][i]=0;
   }
   if(com.pose) free(com.pose);   
   if(com.fpatt) free(com.fpatt);
   com.pose=(int*) malloc(com.ls*sizeof(int));
   com.fpatt=(double*) malloc(com.ls*sizeof(double));
   if(com.z[com.ns-1]==NULL || com.pose==NULL||com.fpatt==NULL) error2("oom zpf");
   com.rgene[0]=1;   com.ngene=1;  

   lline=max2(lline, com.ls/n31*(n31+1)+lspname+50);
   if((line=(char*)malloc(lline*sizeof(char)))==NULL) error2("oom line");
   FOR (j,com.ls) com.pose[j]=0;      /* gene #1, default */
   if(paupseq) goto readseq;

   /* first line */
   if (!fgets(line,lline,fseq)) error2 ("err ReadSeq: first line");
   for (j=0; j<lline && line[j] && line[j]!='\n'; j++) {
      if (!isalnum(line[j])) continue;
      line[j]=(char)toupper(line[j]);
      switch (line[j]) {
         case 'G': noptline++;   break;
         case 'C': basecoding=1; break;
         case 'S': Sequential=1; break;
         case 'I': Sequential=0; break;
         default : printf ("\nBad option %c\n", line[j]);  exit (-1);
      }
   }
   if (strchr(line,'C')) {   /* protein-coding DNA sequences */
      if(com.seqtype==2) error2("option C?");
      else if (com.seqtype==0) {
         if (com.ls%3!=0 || noptline<1)  error2("option C?");
         com.ngene=3; FOR(i,3) com.lgene[i]=com.ls/3;
         for (i=0;i<com.ls;i++) com.pose[i]=i%3;
#if BASEML
         com.coding=1;
#endif
      }
      noptline--;
   }

   /* option lines */
   FOR (j, noptline) {
      line[0]=0;
      while (!isalnum(line[0]))  line[0]=(char)fgetc(fseq); 
      line[0]=(char)toupper(line[0]);
      switch (line[0]) {
      case ('G') :
         if(basecoding) error2("incorrect option format, use GC?\n");
         if (fscanf(fseq,"%d",&com.ngene)!=1) error2 ("expecting #gene here..");
         if (com.ngene>NGENE) error2 ("raise NGENE?");
         if (com.fix_rgene) {    /* specified in the main program? */
            puts ("reading rates for genes from the seq. file, correct?");
            FOR (k, com.ngene) if (fscanf (fseq, "%lf", &com.rgene[k]) != 1)
               error2 ("err: fix_rgene?");
         }
         fgets (line, lline, fseq);
         if (!blankline(line)) {      /* #sites in each gene on the 2nd line */
            for (i=0,p=line; i<com.ngene; i++) {
               while (*p && !isalnum(*p)) p++;
               if (sscanf(p,"%d",&com.lgene[i])!=1) break;
               while (*p && isalnum(*p)) p++;
            }
            for (; i<com.ngene; i++)
               if (fscanf(fseq,"%d", &com.lgene[i])!=1) error2("EOF at lgene");

            for (i=0,k=0; i<com.ngene; k+=com.lgene[i],i++)
               FOR (j, com.lgene[i]) com.pose[k+j]=i;
            for (i=0,k=0; i<com.ngene; i++) k+=com.lgene[i];
            if (k!=lcode) {
               matIout(F0, com.lgene, 1, com.ngene);
               printf("\n%6d != %d", lcode, k);  error2("ls for genes");
            }
         }
         else {                   /* site marks on later line(s)  */
            for(k=0; k<lcode;) {
               if (com.ngene>9)  fscanf(fseq,"%d", &ch);
               else {
                  do ch=fgetc(fseq); while (!isdigit(ch));
                  ch=ch-(int)'1'+1;  /* assumes 1,2,...,9 are consecutive */
               }
               if (ch<1 || ch>com.ngene)
                  { printf("\ngene mark %d at %d?\n", ch, k+1);  exit (-1); }
               com.pose[k++]=ch-1;
            }
            if(!fgets(line,lline,fseq)) error2("err: option lines");
         }
         break;
      default :
         printf ("..Bad option %c\n", line[0]);  exit (-1);
      }
   }

   readseq:
   /* read sequence */
   if (Sequential)  {    /* sequential */
      if (noisy) printf ("\nSequential format..\n");
      FOR (j,com.ns) {
         lspname=LSPNAME;
         FOR (i, 2*lspname) line[i]='\0';
         if (!fgets (line, lline, fseq)) error2("err: EOF?");
         if (blankline(line)) {
            if (!fgets (line, lline, fseq)) error2("err: EOF");
            if (!hasbase(line))
               { printf("err: empty line (seq %d)\n",j+1); exit(-1); }
         }
         p=line+(line[0]=='=' || line[0]=='>') ;
         while(isspace(*p)) p++;
         if ((ch=strstr(p,"  ")-p)<lspname && ch>0) lspname=ch;
         strncpy (com.spname[j], p, lspname);
         k=strlen(com.spname[j]);
         p+=(k<lspname?k:lspname);

         for (; k>0; k--) /* trim spaces */
            if (!isgraph(com.spname[j][k]))   com.spname[j][k]=0;
            else    break;

         if (noisy>=2) printf ("Reading seq #%2d: %s\n", j+1, com.spname[j]);
         for (k=0; k<com.ls; p++) {
            while (*p=='\n' || *p=='\0') {
               p=fgets(line, lline, fseq);
               if(p==NULL)
                  { printf("\nEOF at site %d, seq %d\n", k+1,j+1); exit(-1); }
            }
            *p=(char)toupper(*p);  
            p1=strchr(pch,*p);
            if (p1 && p1-pch>=nchar)  
               miss=1;
            if (*p==eq) {
               if (j==0) error2 ("err: . in 1st seq.?");
               com.z[j][k]=com.z[0][k];  k++;
            }
            else if (p1) 
               com.z[j][k++]=*p;
            else if (isalpha(*p)) 
               { printf("\nerr: %c at %d seq %d.",*p,k+1,j+1); exit(0); }
            else if (*p==(char)EOF) error2 ("err: EOF?");
         }           /* for(k) */
      }              /* for (j,com.ns) */
   }
   else {            /* interlaved */
      if (noisy) printf ("\nInterlaved format..\n");
      FOR (j, com.ns) lt[j]=0;  /* temporary seq length */
      for (igroup=0; ; igroup++) {
         /*
         printf ("\nreading block %d ", igroup+1);  matIout(F0,lt,1,com.ns);*/

         FOR (j, com.ns) if (lt[j]<com.ls) break;
         if (j==com.ns) break;
         FOR (j,com.ns) {
            if (!fgets(line,lline,fseq)) {
               printf("\nerr reading site %d, seq %d group %d",
                  lt[j]+1,j+1,igroup+1);
               error2 ("err: EOF?");
            }
            if (!hasbase(line)) {
               if (j) {
                  printf ("\n%d, seq %d group %d", lt[j]+1, j+1, igroup+1);
                  error2 ("err: empty line.");
               }
               else 
                  if (PopEmptyLines(fseq,lline,line)==-1) {
                     printf ("\n%d, seq %d group %d", lt[j]+1, j+1, igroup+1);
                     error2 ("err: EOF?");
                  }
            }
            p=line;
            if (igroup==0) {
               lspname=LSPNAME;
               while(isspace(*p)) p++;
               if ((ch=strstr(p,"  ")-p)<lspname && ch>0) lspname=ch;
               strncpy (com.spname[j], p, lspname);
               k=strlen(com.spname[j]);
               p+=(k<lspname?k:lspname);

               for (; k>0; k--)   /* trim spaces */
                  if (!isgraph(com.spname[j][k]))  com.spname[j][k]=0;
                  else   break;
               if(noisy>=2) printf("Reading seq #%2d: %s\n",j+1,com.spname[j]);
            }
            for (; *p && *p!='\n'; p++) {
               if (lt[j]==com.ls) break;
               *p=(char)toupper(*p);
               p1=strchr(pch,*p);
               if (p1 && p1-pch>=nchar) 
                  miss=1;
               if (*p==eq) {
                  if (j==0) {
                     printf("err: . in 1st seq, group %d.\n",igroup);
                     exit (-1);
                  }
                  com.z[j][lt[j]] = com.z[0][lt[j]];
                  lt[j]++;
               }
               else if (p1)
                  com.z[j][lt[j]++]=*p;
               else if (isalpha(*p)) {
                  printf("\nerr:%c at %d seq %d block %d.",
                          *p,lt[j]+1,j+1,igroup+1);
                  exit(-1);
               }
               else if (*p==(char)EOF) error2 ("err: EOF");
            }         /* for (*p) */
         }            /* for (j,com.ns) */
      }               /* for (igroup) */
   }
   free (line);
   if(!miss) com.cleandata=1;
   else if (com.cleandata) {  /* remove ambiguity characters */
      if(noisy)  puts("\nSites with gaps or missing data will be removed.");
      if(fout) {
         fprintf(fout,"\nBefore deleting alignment gaps. %d sites\n",com.ls);
         printsma(fout,com.spname,com.z,com.ns,com.ls,com.ls,gap,0,NULL);
      }
      RemoveIndel ();
      lcode=com.ls/n31;
      if(fout) fprintf(fout,"\nAfter deleting gaps. %d sites\n",com.ls);
   }

   if(fout) printsma(fout,com.spname,com.z,com.ns,com.ls,com.ls,gap,0,NULL);
/*
   if(fout) printsma(fout,com.spname,com.z,com.ns,com.ls,60,0,0,NULL);
*/
   if (noisy>=2) printf ("\nSequences read..\n");
   if(n31==3) com.ls=lcode;

   /*  IdenticalSeqs(); */

#ifdef CODEML
   if(com.seqtype==CODON2AAseq) {
      if (noisy>2) puts("\nTranslating into AA sequences\n");
      FOR(i,com.ns) {
         if (noisy>2) printf("Translating sequence %d\n",i+1);
         DNA2protein(com.z[i], com.z[i], com.ls,com.icode);
      }
      com.seqtype=AAseq;
      if(fout) {
         fputs("\nTranslated AA Sequences\n",fout);
         fprintf(fout,"%4d  %6d",com.ns,com.ls);
         printsma(fout,com.spname,com.z,com.ns,com.ls,com.ls,10,0,NULL);
      }
   }
#endif

   return (0);
}

int IdenticalSeqs(void)
{
/* This checks for identical sequences and create a data set of unique 
   sequences.  The file name is <SeqDataFile>.tmp.  This is casually 
   written and need more testing.
   The routine is called right after the sequence data are read.
   For codon sequences, com.ls has the number of codons, which are NOT
   coded.
*/
   char tmpf[96], keep[NS];
   FILE *ftmp;
   int is,js,h, same,nkept=com.ns;
   int ls1=com.ls*(com.seqtype==CODONseq||com.seqtype==CODON2AAseq?3:1);

   puts("\nIdenticalSeqs: not tested\a");
   FOR(is,com.ns) keep[is]=1;
   FOR(is,com.ns)  { 
      if(!keep[is]) continue;
      FOR(js,is) {
         for(h=0,same=1; h<ls1; h++)
            if(com.z[is][h]!=com.z[js][h]) break;
         if(h==ls1) {
            printf("Seqs. %3d & %3d (%s & %s) are identical!\n",
               js+1,is+1,com.spname[js],com.spname[is]);
            keep[is]=0;
         }
      }
   }
   FOR(is,com.ns) if(!keep[is]) nkept--;
   if(nkept<com.ns) {
      strcpy(tmpf,com.seqf);  strcat(tmpf,".tmp");
      if((ftmp=fopen(tmpf,"w"))==NULL) error2("IdenticalSeqs: file error");
      printSeqs(ftmp,NULL,keep,1);
      fclose(ftmp);
      printf("\nUnique sequences collected in %s.\n",tmpf);
   }
   return(0);
}


int PatternWeight (FILE *fout, double space[])
{
/* This collaps sites into patterns, by looking at nc characters at a time 
   so that nb=1 for nucleotide and aa seqs and =3 for codon seqs.
   With codon seqs, com.ls is the # of codons, while com.z[] is 3*com.ls
   nucleotides long.
   This should work whether or not the data are encoded.
   poset & zt[] (space) is temporary and is copied back to com.z[]
   com.pose[i] takes value from 0 to npatt and maps sites to patterns
           (i=0,...,ls)
   space[com.ls*sizeof(int) + com.ns*com.ls*nb*sizeof(char)]

   This is rather inefficient for codon sequences when cleandata=1.  Think 
   about transforming data first before this routine.  Trouble is how 
   to deal with InitializeCodon, which now takes non-coded sequences

*/
   int h,ht,j,k=-1,newpatt,ig, *poset=(int*)space;
   int nb=(com.seqtype==CODONseq?3:1), gap=(nb==3?3:10);
   char *zt[NS], b[NS][3], miss[]="-?"; /* b[][] holds data at site h */

   FOR (j,com.ns) zt[j]=(char*)(poset+com.ls)+j*com.ls*nb;
   FOR (h,com.ls) com.fpatt[h]=0;
   FOR (ig, com.ngene) com.lgene[ig]=0;
   for(ig=0,com.npatt=0; ig<com.ngene; ig++) {
      com.posG[ig]=com.npatt;
      for(h=0; h<com.ls; h++) {
         if(com.pose[h] != ig) continue;
         FOR(j,com.ns) FOR(k,nb) b[j][k]=com.z[j][h*nb+k];

         FOR(j,com.ns) FOR(k,nb) 
            if(b[j][k]!=miss[0]&&b[j][k]!=miss[1]) goto NotAllMiss;
         NotAllMiss:
         if(noisy>=9 && j==com.ns && k==nb) puts("Some sites have no data (no harm).");

         com.lgene[ig]++;

         for(ht=com.posG[ig],newpatt=1; ht<com.npatt; ht++) {
            FOR(j,com.ns) FOR(k,nb)
               if(b[j][k]!=zt[j][ht*nb+k]) goto CheckNextSite;
            poset[h]=ht; newpatt=0; break;    /* h has same data as ht */
            CheckNextSite: ;
         }
         if(newpatt) {
            FOR(j,com.ns) FOR(k,nb) zt[j][com.npatt*nb+k]=b[j][k];
            poset[h]=com.npatt++;
            if(noisy && (h+1)%2000==0) 
               printf("\r%d patterns at %d sites", com.npatt,h+1);
         }
         com.fpatt[poset[h]]++;
      }     /* for (h)  */
   }        /* for (ig) */

   if(com.seqtype==CODONseq && com.ngene==3 &&com.lgene[0]==com.ls/3) {
      puts("\nCheck the G option in data file? (Enter)\n");
      getchar();
   }

   com.posG[com.ngene]=com.npatt;
   for(j=1; j<com.ngene; j++) com.lgene[j]+=com.lgene[j-1];

   if(com.lgene[com.ngene-1]!=com.ls) { puts("\algene vs. ls, missing data"); }

   FOR(j,com.ns) {
      memcpy (com.z[j], zt[j], com.npatt*nb);
      com.z[j]=(char*)realloc(com.z[j],com.npatt*nb*sizeof(char));
   }
   memcpy(com.pose, poset, com.ls*sizeof(int));
   com.fpatt=(double*)realloc(com.fpatt, com.npatt*sizeof(double));

   if (fout) {
      fprintf (fout, "\nns = %d  \tls = %d", com.ns,com.ls);
      if (com.ngene>1) {
         fprintf (fout,"\nngene =%3d, lengths =", com.ngene);
         FOR(j,com.ngene)
            fprintf (fout,"%7d", (j?com.lgene[j]-com.lgene[j-1]:com.lgene[j]));
         fprintf(fout,"\n  Starting at pattern");
         FOR(j,com.ngene) fprintf(fout,"%7d", com.posG[j]+1);
      }
      if(noisy) printf("\n# site patterns = %d\n", com.npatt);
      fprintf(fout,"\n# site patterns = %d\n", com.npatt);
      FOR (h,com.npatt) {
         fprintf(fout," %4.0f", com.fpatt[h]);
         if((h+1)%15==0) FPN(fout);
      }
      FPN(fout);
      newpatt=com.npatt*nb;
      printsma(fout,com.spname,com.z,com.ns,newpatt,newpatt,gap,1,NULL);
   }
   return (0);
}


int Initialize (FILE *fout, double *space, int seqtype)
{
/* Count site patterns (com.fpatt) and calculate base or amino acid frequencies 
   in genes and species.
   Ambiguity characters in sequences are resolved by iteration. 
   For frequencies in each species, they are resolved within that sequence.
   For average base frequencies among species, they are resolved over all 
   species.

   This routine is called by baseml and aaml.  codonml uses another
   routine InitializeCode()
*/
   char *pch=(seqtype==0?BASEs:(seqtype==2?AAs:BINs)), indel[]="-?";
   int wname=30, h,j,k, ig, nconstp, nc=com.ncode;
   int irf,nrf=20, miss=0, b,nb,ib[4];
   double pi0[20],t,lmax=0;

   PatternWeight(fout,space);

   if(fout)  fprintf(fout,"\nFrequencies..");
   for(h=0,nconstp=0; h<com.npatt; h++) {
      for (j=1;j<com.ns;j++)  if(com.z[j][h]!=com.z[0][h]) break;
      if (j==com.ns && com.z[0][h]!=indel[0] && com.z[0][h]!=indel[1])  
         nconstp+=(int)com.fpatt[h];
   }
   for (ig=0,zero(com.pi,nc); ig<com.ngene; ig++) {
      if (com.ngene>1)  fprintf (fout,"\n\nGene # %2d (len %4d)",
                             ig+1, com.lgene[ig]-(ig==0?0:com.lgene[ig-1]));
      fprintf(fout,"\n%*s",wname, ""); FOR(k,nc) fprintf(fout,"%7c",pch[k]);

      /* freqs in species first, com.pi and pi0 are new and current freqs 
         Those are just output, not used later anymore.
      */
      for(j=0,zero(com.piG[ig],nc); j<com.ns; j++) {
         FOR(k,nc) { pi0[k]=1./nc; com.pi[k]=0; }
         for(irf=0; irf<nrf; irf++) {
            for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
               b=strchr(pch,com.z[j][h])-pch;
               if(b<0) error2("Initialize: strange char");
               else if (b<nc)
                  com.pi[b]+=com.fpatt[h];
               else if(com.seqtype==BASEseq) {
                  NucListall(com.z[j][h], &nb, ib);
                  for(k=0,t=0; k<nb; k++) t+=pi0[ib[k]];
                  FOR(k,nb) com.pi[ib[k]]+=pi0[ib[k]]/t*com.fpatt[h];
               }
               else if(com.seqtype==AAseq)
                  FOR(k,20) com.pi[k]+=pi0[k]*com.fpatt[h];
            }
            abyx(1/sum(com.pi,nc), com.pi, nc);
            if(distance(com.pi,pi0,nc)<1e-6) break;
            xtoy(com.pi,pi0,nc);
         }   /* for(irf) */
         fprintf(fout,"\n%-*s",wname,com.spname[j]);
         FOR(k,nc) fprintf(fout,"%7.4f",com.pi[k]);
      }    /* for(j) */

      /* freqs in genes com.piG[], averaged over species.  These are used in 
         ML calculations 
      */
      FOR(k,nc) { pi0[k]=1./nc; com.piG[ig][k]=0; }
      for(irf=0; irf<nrf; irf++) {
         FOR(j,com.ns) {
            for(h=com.posG[ig]; h<com.posG[ig+1]; h++) {
               b=strchr(pch,com.z[j][h])-pch;
               if(b<nc)
                  com.piG[ig][b]+=com.fpatt[h];
               else if(com.seqtype==BASEseq) {
                  NucListall(com.z[j][h], &nb, ib);
                  for(k=0,t=0,miss=1; k<nb; k++) t+=pi0[ib[k]];
                  FOR(k,nb) com.piG[ig][ib[k]]+=pi0[ib[k]]/t*com.fpatt[h];
               }
               else if(com.seqtype==AAseq)
                  for(k=0,miss=1;k<20;k++) com.piG[ig][k]+=pi0[k]*com.fpatt[h];
            }   /* for(h) */
         }      /* for(j) */
         abyx(1/sum(com.piG[ig],nc), com.piG[ig], nc);
         if(distance(com.piG[ig],pi0,nc)<1e-6) break;
         xtoy(com.piG[ig],pi0,nc);
      }         /* for(irf) */
      if(com.ngene>1) {
         fprintf(fout,"\n\n%-*s",wname,"Mean");
         FOR(k,nc) fprintf(fout,"%7.4f",com.piG[ig][k]);
      }
   }            /* for(ig) */
   for (ig=0,zero(com.pi,nc); ig<com.ngene; ig++) {
      t=(ig==0?com.lgene[0]:com.lgene[ig]-com.lgene[ig-1])/(double)com.ls;
      FOR(k,nc)  com.pi[k]+=com.piG[ig][k]*t;
   }    /* for(ig) */
   fprintf (fout, "\n\n%-*s", wname, "Average");
   FOR (j,nc) fprintf(fout,"%7.4f", com.pi[j]);
   if(miss) fputs("\n(Ambiguity characters are used to calculate freqs.)\n",fout);

   fprintf (fout,"\n\n# constant sites: %6d (%.2f%%)",
            nconstp, (double)nconstp*100./com.ls);
   if (com.model==0 || (seqtype==BASEseq && com.model==1)) {
      fillxc (com.pi, 1./nc, nc);
      FOR (j,com.ngene) xtoy (com.pi, com.piG[j], nc);
   }

   if (com.seqtype==AAseq) {
      for (j=0,t=0; j<nc; j++) t+=(com.pi[j]>0);
      if (t<=4)  printf ("\n\a\t\tAre these a.a. sequences?");
   }
   if (!miss && com.ngene==1) {
      for (h=0,lmax=-(double)com.ls*log((double)com.ls); h<com.npatt; h++)
         if (com.fpatt[h]>1) lmax+=com.fpatt[h]*log((double)com.fpatt[h]);
   }
   if (fout) {
     fprintf(fout, "\nln Lmax (unconstrained) = ");
     if (lmax)  fprintf(fout, " %.6f\n", lmax);
     else     fputs(" unavailable due to missing or partitioned data.\n",fout);
     fflush (fout);
   }

   if(com.cleandata)  /* for aa and base sequences */
     FOR(j,com.ns) if(transform(com.z[j],com.npatt,1,com.seqtype)) error2("uh?");

   if (noisy>2) printf ("\n\nInitialized..\n");
   return(0);
}

int RemoveIndel(void)
{
/* Remove ambiguity characters and indels in the untranformed sequences, 
   Changing com.ls and com.pose[] (site marks for multiple genes).
   For codonml, com.ls is still 3*#codons
   Called at the end of ReadSeq, when com.pose[] are still site marks.
   All characters in com.z[][] not found in the character string pch are
   considered ambiguity characters and are removed.
*/
   int  h,k, j,js,lnew,nindel, n31,nchar;
   char b, *pch, *miss;  /* miss[h]=1 if site (codon) h is missing, 0 otherwise */

   if(com.seqtype==CODONseq||com.seqtype==CODON2AAseq)
      { n31=3; nchar=4; pch=BASEs; }
   else {
      n31=1;
      if(com.seqtype==AAseq)        { nchar=20; pch=AAs; }
      else if(com.seqtype==BASEseq) { nchar=4; pch=BASEs; }
      else                          { nchar=2; pch=BINs; }
    }

   if (com.ls%n31) error2 ("ls in RemoveIndel.");
   if((miss=(char*)malloc(com.ls/n31 *sizeof(char)))==NULL)
      error2("oom miss");
   FOR (h,com.ls/n31) miss[h]=0;
   for (js=0; js<com.ns; js++) {
      for (h=0,nindel=0; h<com.ls/n31; h++) {
         for (k=0; k<n31; k++) {
            b=(char)toupper(com.z[js][h*n31+k]);
            FOR(j,nchar) if(b==pch[j]) break;
            if(j==nchar) { miss[h]=1; nindel++; }
         }
      }
      if (nindel) printf ("\n%6d ambiguity characters in seq. %d", nindel,js+1);
   }
   if(noisy) puts("\nThe following sites are removed: ");
   for (h=0; h<com.ls/n31; h++)  if(miss[h]) printf(" %2d", h+1);

   for (h=0,lnew=0; h<com.ls/n31; h++)  {
      if(miss[h]) continue;
      for (js=0; js<com.ns; js++) {
         for (k=0; k<n31; k++)
            com.z[js][lnew*n31+k]=com.z[js][h*n31+k];
      }
      com.pose[lnew++]=com.pose[h];
   }
   com.ls=lnew*n31;
   free(miss);
   return (0);
}



int MPInformSites (void)
{
/* Outputs parsimony informative and noninformative sites into 
   two files named MPinf.seq and MPninf.seq
   Uses transformed sequences
*/
   char *imark, *pch=(com.seqtype==0?BASEs:(com.seqtype==2?AAs:BINs));
   int h, i, markb[NS], inf, lsinf;
   FILE *finf, *fninf;

puts("\nMPInformSites: missing data not dealt with yet?\n");

   finf=fopen("MPinf.seq","w");
   fninf=fopen("MPninf.seq","w");
   if (finf==NULL || fninf==NULL) error2("MPInformSites: file creation error");

   puts ("\nSorting parsimony-informative sites: MPinf.seq & MPninf.seq");
   if ((imark=(char*)malloc(com.ls*sizeof(char)))==NULL) error2("oom imark");
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

if(!com.cleandata) {
  puts("PatternJC69like: does not work with missing data\n");
  return(0);
}

   for (h=0,com.npatt=0,ig=-1; h<npatt0; h++) {
      if (ig<com.ngene-1 && h==com.posG[ig+1]) { com.posG[++ig]=com.npatt; }
      zh[0]=b=0;
      for (j=1; j<com.ns; j++) {
         FOR(k,j) if (com.z[j][h]==com.z[k][h]) break;
         zh[j]= (char)(k<j?zh[k]:++b);
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
   if (noisy>=2) printf ("\nnew no. site patterns:%7d\n", com.npatt);

   if (fout) {
      fprintf(fout,"\nnew number of site patterns = %d\n", com.npatt);
      FOR (h,com.npatt) {
         fprintf(fout," %4.0f", com.fpatt[h]);
         if((h+1)%15==0) FPN(fout);
      }
/*    for baseml only, to print out the patterns.  comment out.

      fprintf (fout, "\nno. site patterns:%7d\n", com.npatt);
      FOR(j,com.ns) FOR(h,com.npatt) com.z[j][h]=BASEs[com.z[j][h]];
      printsma(fout,com.spname,com.z,com.ns,com.npatt,com.npatt,10,0,NULL);
      exit(0);
*/
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


int print1seq (FILE*fout, char *z, int ls, int encoded, int pose[])
{
/* This prints out one sequence, which may be coded.  Codon seqs are coded
   as 0,1,2,...60 if called from codeml.c or 0,1,2,...,63 otherwise.
   This may be risking.  Check when use.
   z[] contains patterns if (pose!=NULL)
   This uses com.seqtype.
*/
   int i, h,hp, gap=10;
   char *pch=(com.seqtype==0?BASEs:(com.seqtype==2?AAs:BINs)), str[4]="";
   int nb=(com.seqtype==CODONseq?3:1);


   if(!encoded) {  /* raw data, not coded */
      FOR(h,ls) {
         hp=(pose?pose[h]:h);
         FOR(i,nb) 
            fputc(z[hp*nb+i],fout);
/*
         FOR(i,nb) fputc(z[(pose?pose[h]:h) *nb+i],fout);
*/
         if(com.seqtype==CODONseq || (h+1)%gap==0) fputc(' ',fout);
      }
   }
   else {    /* cleandata, coded */
      FOR(h,ls) {
         hp=(pose?pose[h]:h);
         if(com.seqtype!=CODONseq) fputc(pch[z[hp]],fout);
         else {
#ifdef CODEML
            fprintf(fout,"%s",getcodon(str,FROM61[z[hp]])); /* 0,1,...,60 */
#else
            fprintf(fout,"%s",getcodon(str,z[hp]));         /* 0,1,...,63 */
#endif
         }
         if(com.seqtype==CODONseq || (h+1)%gap==0) fputc(' ',fout);
      }
   }
   return(0);
}

void printSeqs(FILE *fout, int *pose, char keep[], int format)
{
/* Print sequences into fout, using paml (format=0) or paup (format=1) formats.
   Use pose=NULL if called before site patterns are collapsed.  
   Use NULL for keep if all sequences are to be printed.
   Sequences may (com.cleandata==1) and may not (com.cleandata=0) be coded.
   com.z[] has site patterns if pose!=NULL.
   This uses com.seqtype, and com.ls is the number of codons for codon seqs.
   See notes in print1seq()
   
   format=0:PAML & PHYLIP format
          1:NEXUS, with help from John Huelsenbeck.  Thanks to John.
*/
   int j,ls1=(com.seqtype==CODONseq?3*com.ls:com.ls), nskept=com.ns, wname=30;
   char *dt=(com.seqtype==AAseq?"protein":"dna");

   if(keep) FOR(j,com.ns) nskept -= !keep[j];
   if (format==1) {  /* NEXUS format */
      fprintf(fout,"\nbegin data;\n");
      fprintf(fout,"   dimensions ntax=%d nchar=%d;\n", nskept,ls1);
      fprintf(fout,"   format missing=? gap=- datatype=%s;\n   matrix\n",dt);
   }
   else    
      fprintf(fout,"%8d %8d\n",nskept,ls1);

   for(j=0; j<com.ns; j++,FPN(fout)) {
      if(keep && !keep[j]) continue;
      fprintf(fout,"%s%-*s  ", (format?"      ":""),wname,com.spname[j]);
      print1seq(fout,com.z[j],com.ls,com.cleandata, pose);
   }
   if (format==1) fprintf(fout, "   ;\nend;");
   FPN(fout);
   fflush(fout);
}

#define gammap(x,alpha) (alpha*(1-pow(x,-1/alpha)))
/* DistanceREV () used to be here, moved to pamp. 
*/

#if (NCODE==4) 

double SeqDivergence (double x[], int model, double alpha, double *kappa)
{
/* HKY85 model, following Tamura (1993, MBE, ..)
   alpha=0 if no gamma 
   return -1 if in error.
*/
   int i,j;
   double p[4], Y,R, a1,a2,b, P1,P2,Q,fd,tc,ag, a1t,a2t,bt;
   double small=1e-6,largek=999;

   if (testXMat(x)) error2 ("X err..");
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
      if (1-fd/b<=0) return (fd);

      if (alpha<=0) return (-b*log (1-fd/b));
      else return  (-b*gammap(1-fd/b,alpha));
   case (K80) :
      a1=1-2*(P1+P2)-Q;   b=1-2*Q;
/*      if (a1<=0 || b<=0) return (-1); */
      if (a1<=0 || b<=0) return (fd);
      if (alpha<=0)  { a1=-log(a1);  b=-log(b); }
      else          { a1=-gammap(a1,alpha); b=-gammap(b,alpha); }
      a1t=.5*a1-.25*b;  bt=.25*b;
      if(bt>small) *kappa = a1t/bt; else *kappa=largek;
/*
fprintf (frst, "%9.3e %9.3e %9.4f", P1+P2, Q, (P1+P2)/Q);
*/
      return (a1t+2*bt);
   case (F84):
      a1=(2*(tc+ag)+2*(tc*R/Y+ag*Y/R)*(1-Q/(2*Y*R)) -P1-P2) / (2*tc/Y+2*ag/R);
      b = 1 - Q/(2*Y*R);
/*      if (a1<=0 || b<=0) return (-1); */
      if (a1<=0 || b<=0) return (fd);
      if (alpha<=0) { a1=-log(a1); b=-log(b); }
      else         { a1=-gammap(a1,alpha); b=-gammap(b,alpha); }
      a1t=.5*a1;  bt=.5*b;
      *kappa = a1t/bt-1;
      /* *kappa = min2(*kappa, largek);   */    *kappa = max2(*kappa, -.5);
      return  4*bt*(tc*(1+ *kappa/Y)+ag*(1+ *kappa/R)+Y*R);
   case (HKY85):         /* TN93  */
      a1=1-Y*P1/(2*tc)-Q/(2*Y);  a2=1-R*P2/(2*ag)-Q/(2*R);   b=1-Q/(2*Y*R);
/*      if (a1<=0 || a2<=0 || b<=0) return (-1); */
      if (a1<=0 || a2<=0 || b<=0) return (fd);
      if (alpha<=0) { a1=-log(a1); a2=-log(a2); b=-log(b); }
      else   { a1=-gammap(a1,alpha); a2=-gammap(a2,alpha); b=-gammap(b,alpha);}
      a1t=.5/Y*(a1-R*b);  a2t=.5/R*(a2-Y*b);  bt=.5*b;
      *kappa = largek;
      if (bt>0) *kappa = min2((a1t+a2t)/2/bt, largek);
      return 4*p[0]*p[1]*a1t + 4*p[2]*p[3]*a2t + 4*Y*R*bt;
   }    /* switch  */
   return (-1);
}


#ifdef LSDISTANCE
#ifdef REALSEQUENCE
extern double *SeqDistance;

int DistanceMatNuc (FILE *fout, int model, double alpha)
{
/* This calculates pairwise distances.  The data may be clean and coded 
   (com.cleandata==1) or not.  
   In the latter case, ambiguity sites are not used (pairwise deletion).
   Site patterns are used.
*/
   char *models[]={"JC69", "K80", "F81", "F84", "TN93"}, distf[32]="2base.t";
   int is,js,k1,k2, h, miss=0, status=0, nc=4;
   double x[16], kappat=0, t,bigD=5;
   FILE *fdist=(FILE*)fopen(distf,"w");
   
   if(fdist==NULL) error2("DistnaceMatNuc: file creation error");
   fprintf(fdist,"%6d\n", com.ns);
   if (model>=HKY85) model=4; /* TN93 here */
   if (fout) {
      fprintf(fout,"\nDistances:%5s", models[model]);
      if (model!=JC69 && model!=F81) fprintf (fout, "(kappa) ");
      fprintf(fout," (alpha set at %.2f)\n", alpha);
      fprintf(fout,"\nThis matrix is not used in later m.l. analysis.\n");
      if(!com.cleandata) fputs("\n(Pairwise deletion.)",fout);
   }
   for (is=0; is<com.ns; is++) {
      if (fout) fprintf(fout,"\n%-12s", com.spname[is]);
      if(fdist) fprintf(fdist,"%-10s ", com.spname[is]);
      FOR (js, is) {
         if(com.cleandata)
            for (h=0,zero(x,16); h<com.npatt; h++)
               x[com.z[is][h]*nc+com.z[js][h]] += com.fpatt[h];
         else 
            for (h=0,zero(x,16); h<com.npatt; h++) {
               for(k1=0;k1<nc;k1++) if(BASEs[k1]==com.z[is][h]) break;
               for(k2=0;k2<nc;k2++) if(BASEs[k2]==com.z[js][h]) break;
               if(k1<nc && k2<nc)   x[k1*nc+k2] += com.fpatt[h];
               else                 miss=1;
            }
         abyx(1./sum(x,16),x,16);
         if ((t=SeqDivergence(x, model, alpha, &kappat)) < 0)
            { t=-1; status=-1; }
         SeqDistance[is*(is-1)/2+js] = (t>=0?t:bigD);
         if(fdist) fprintf(fdist," %7.4f", t);
         if (fout) fprintf(fout,"%8.4f", t);
         if (fout && (model==K80 || model==F84 || model==HKY85))
            fprintf(fout,"(%7.4f)", kappat);
/*
fprintf(frst,"%5d%5d %8.4f %8.4f\n", is+1,js+1, t, kappat);
*/
       }
       if(fdist) FPN(fdist);
   }
/*
fflush(frst);
*/
   FPN(fout);
   fclose(fdist);
   if(status) puts("\ndistance formula sometimes inapplicable.."); 
   return(status);
}

#endif
#endif

#ifdef BASEMLG
extern int CijkIs0[];
#endif

extern int nR;
extern double Cijk[], Root[];

int RootTN93 (int model, double kappa1, double kappa2, double pi[], 
    double *f, double Root[])
{
   double T=pi[0],C=pi[1],A=pi[2],G=pi[3],Y=T+C,R=A+G;

   if (model==F84) { kappa2=1+kappa1/R; kappa1=1+kappa1/Y; }
   else if (model==HKY85) kappa2=kappa1;

   *f=1/(2*T*C*kappa1+2*A*G*kappa2 + 2*Y*R);

   Root[0] = 0;
   Root[1] = - (*f);
   Root[2] = -(Y+R*kappa2) * (*f);
   Root[3] = -(Y*kappa1+R) * (*f);
   return (0);
}


int EigenTN93 (int model, double kappa1, double kappa2, double pi[],
    int *nR, double Root[], double Cijk[])
{
/* initialize Cijk[] & Root[], which are the only part to be changed
   for a new substitution model
   for JC69, K80, F81, F84, HKY85, TN93
   Root: real Root divided by v, the number of nucleotide substitutions.
*/
   int i,j,k, nr;
   double f, U[16],V[16], t;
   double T=pi[0],C=pi[1],A=pi[2],G=pi[3],Y=T+C,R=A+G;

   if (model==JC69 || model==F81) kappa1=kappa2=com.kappa=1; 
   else if (com.model<=HKY85)     kappa2=kappa1;
   RootTN93 (model, kappa1, kappa2, pi, &f, Root);

   *nR=nr = 2+(model==K80||model>=F84)+(model>=HKY85);
   U[0*4+0]=U[1*4+0]=U[2*4+0]=U[3*4+0]=1;
   U[0*4+1]=U[1*4+1]=1/Y;   U[2*4+1]=U[3*4+1]=-1/R;
   U[0*4+2]=U[1*4+2]=0;  U[2*4+2]=G/R;  U[3*4+2]=-A/R;
   U[2*4+3]=U[3*4+3]=0;  U[0*4+3]=C/Y;  U[1*4+3]=-T/Y;

   xtoy (pi, V, 4);
   V[1*4+0]=R*T;   V[1*4+1]=R*C;
   V[1*4+2]=-Y*A;  V[1*4+3]=-Y*G;
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
         error2 ("err:model");
      }
   }
#ifdef BASEMLG
   FOR (i,64) CijkIs0[i] = (Cijk[i]==0);
#endif
   return(0);
}

#ifdef EIGEN

int EigenREV (FILE* fout, double rates[], double pi[], 
              int *nR, double Root[], double Cijk[])
{
/* pi[] is constant
*/
   int i,j,k;
   double Q[16], U[16], V[16], T1[16], T2[16], mr;

   *nR=4;
   zero (Q, 16);
   for (i=0,k=0; i<3; i++) for (j=i+1; j<4; j++)
      if (i*4+j!=11) Q[i*4+j]=Q[j*4+i]=rates[k++];
   for (i=0,Q[3*4+2]=Q[2*4+3]=1; i<4; i++) FOR (j,4) Q[i*4+j] *= pi[j];
   for (i=0,mr=0; i<4; i++) 
      { Q[i*4+i]=0; Q[i*4+i]=-sum(Q+i*4, 4); mr-=pi[i]*Q[i*4+i]; }
   abyx (1/mr, Q, 16);
   if (fout) {
      mr=2*(pi[0]*Q[0*4+1]+pi[2]*Q[2*4+3]);
      fprintf (fout, "\nRate matrix Q, Average Ts/Tv =%9.4f", mr/(1-mr));
      matout (fout, Q, 4,4);
   }
   if ((k=eigen (1,Q,4,Root,T1,U,V,T2))!=0) 
      error2 ("\ncomplex roots in EigenREV?");
   xtoy (U, V, 16);
   matinv (V, 4, 4, T1);
   FOR (i,4) FOR(j,4) FOR(k,4) Cijk[i*4*4+j*4+k] = U[i*4+k]*V[k*4+j];

   return (0);
}



int EigenUNREST (FILE *fout, double kappa[], double pi[], 
    int *nR, complex cRoot[], complex cU[], complex cV[])
{
/* pi[] is changed
*/
   int i,j,k;
   double Q0[16], Q[16], U[16], V[16], T1[16], T2[16], mr;

   *nR=4;
   for (i=0,k=0; i<4; i++) FOR (j,4) 
      if (i!=j && i*4+j != 14)  Q[i*4+j]=kappa[k++];
   for (i=0,Q[14]=1; i<4; i++)  { Q[i*4+i]=0; Q[i*4+i]=-sum(Q+i*4, 4); }
   xtoy(Q,Q0,16);

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

double *SeqDistance=NULL; 
int *ancestor=NULL;

int fun_LS (double x[], double diff[], int np, int npair);

int fun_LS (double x[], double diff[], int np, int npair)
{
   int i,j, aa, it=-np;
   double dexp;

   if (SetBranch(x) && noisy>2) puts ("branch len.");
   if (npair != com.ns*(com.ns-1)/2) error2 ("# seq pairs err.");
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

int LSDistance (double *ss,double x[],int (*testx)(double x[],int np))
{
/* get Least Squares estimates of branch lengths for a given tree topology
*/
   int i,j, h, a1,a2;

   if ((*testx)(x, com.ntime)) {
      matout (F0, x, 1, com.ntime);
      puts ("initial err in LSDistance()");
   }
   FOR (i, com.ns) FOR (j, i) {
      h=i*(i-1)/2+j;
      ancestor[h]=-1;
      for (a1=i; a1!=-1; a1=nodes[a1].father) {
         for (a2=j; a2!=-1; a2=nodes[a2].father)
            if (a1==a2) { ancestor[h]=a1; break; }
         if (ancestor[h] != -1) break;
      }
      if (ancestor[h] == -1) error2 ("no ancestor");
   }      
   i=nls2((com.ntime>20&&noisy>=3?F0:NULL),
      ss,x,com.ntime,fun_LS,NULL,testx,com.ns*(com.ns-1)/2,1e-6);
   return (i);
}

#endif 

#ifdef NODESTRUCTURE

void BranchToNode (void)
{
/* tree.root need to be specified before calling this
*/
   int i, from, to;
   
   tree.nnode=tree.nbranch+1;

/*
   if (tree.root<0 || tree.root>com.ns*2-2) 
      { printf ("root at %d", tree.root+1); error2 ("tree root"); }
*/
   FOR (i,tree.nnode)
      { nodes[i].father=nodes[i].ibranch=-1;  nodes[i].nson=0; }
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

void NodeToBranchSub (int inode);

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
   NodeToBranchSub (tree.root);
   if (tree.nnode!=tree.nbranch+1)  puts ("nd!=nb+1");
}


void ClearNode (int inode)
{
/* a source of confusion. Try not to use this routine.
*/
   nodes[inode].father=nodes[inode].ibranch=-1;
   nodes[inode].nson=0;
   nodes[inode].branch=nodes[inode].divtime=0;
   /* nodes[inode].label=0; */
   /* nodes[inode].branch=0; clear node structure only, not branch lengths */
   /* FOR (i, com.ns) nodes[inode].sons[i]=-1; */
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
         if(tree.branches[j][i]<0 || tree.branches[j][i]>com.ns*2-1) 
            error2("ReadaTreeB: node numbers out of range");
      }
      if (tree.branches[j][0]<com.ns) YoungAncestor=1;
/*
      printf ("\nBranch #%3d: %3d -> %3d",
         j+1, tree.branches[j][0]+1,tree.branches[j][1]+1);
*/
   }
   if(popline) fgets(line, 254, ftree);
   tree.root=tree.branches[0][0];

   if(YoungAncestor)  puts("Ancestors in the data?  Check lkl carefully.");

/*
   com.ntime = com.clock ? (tree.nbranch+1)-com.ns+(tree.root<com.ns)
                         : tree.nbranch;
*/

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

int GetTreeFileType(FILE *ftree, int *ntree, int *pauptree, int shortform);

int GetTreeFileType(FILE *ftree, int *ntree, int *pauptree, int shortform)
{
/* paupstart="begin trees" and paupend="translate" identify paup tree files.
   paupch=";" will be the last character before the list of trees.
   Modify if necessary.
*/
   int i,k, lline=10000, ch=0, paupch=';';
   char line[10000];
   char *paupstart="begin tree", *paupend="translate";

   *pauptree=0;
   k=fscanf(ftree,"%d%d",&i,ntree);
   if(k==2 && i==com.ns) return(0);                 /* old paml style */
   else if(k==1) { *ntree=i; return(0); }           /* phylip & molphy style */
   while(ch!='(' && !isalnum(ch)) ch=fgetc(ftree);  /* treeview style */
   if(ch=='(') { *ntree=-1; ungetc(ch,ftree); return(0); }

   puts("\n# seqs in tree file does not match.  Read as the nexus format.");
   for ( ; ; ) {
      if(fgets(line,lline,ftree)==NULL) error2("tree err1: EOF");
      strcase(line,0);
      if (strstr(line,paupstart)) { *pauptree=1; *ntree=-1; break; }
   }
   if(shortform) return(0);
   for ( ; ; ) {
      if(fgets(line,lline,ftree)==NULL) error2("tree err2: EOF");
      strcase(line,0);
      if (strstr(line,paupend)) break;
   }
   for ( ; ; ) {
      if((ch=fgetc(ftree))==EOF) error2("tree err3: EOF");
      if (ch==paupch) break;
   }
   if(fgets(line,lline,ftree)==NULL) error2("tree err4: EOF");

   return(0);
}

int PaupTreeRubbish(FILE *ftree);
int PaupTreeRubbish(FILE *ftree)
{
/* This pops out the rubbish, typically "tree PAUP_1 = [&U]" with
   "[&U]" optional, before each tree representation in a paup tree file
*/
   int ch;
   char line[10000], *paup1="tree", treename[64], c;

   fscanf(ftree,"%s",line);
   strcase(line,0);
   if (strstr(line,paup1)==NULL) return(-1);
   fscanf(ftree, "%s%s%c", treename, line, &c);

   for (; ;) {
      if((ch=fgetc(ftree))==EOF) { 
         puts("err or end of treefile: expect ]");
         return(-1);
      }
      else if (ch==']') break;
      else if (ch=='(') error2("err treefile: strange");
   }
   return(0);
}



int ReadaTreeN (FILE *ftree, int *length_label, int popline)
{
/* Read a tree from ftree, using the parenthesis node representation of trees.
   *length_label = 1 if the tree has lengths and 
                 = 2 if the tree has labels
                 = 3 if the tree has lengths and labels
   Branch lengths are read in nodes[].branch, and branch (node) labels (integers) are 
   preceeded by # and read in nodes[].label.
   Both names and numbers for species are accepted.  
   Species names are considered case-sensitive, with trailing blanks ignored;
   they have to match the names in the sequence data file.
*/
   int cnode, cfather=-1;  /* current node and father */
   int inodeb=0;  /* node number that will have the next branch length */
   int i,j, level=0, hasname, haslength=0, haslabel=0, ch=' ', lline=255;
   char check[NS], line[255], delimiters[]="(),:#$;";

   *length_label=0;
   tree.nnode=com.ns;  tree.nbranch=0;
   FOR(i,2*com.ns-1) {
      nodes[i].father=nodes[i].ibranch=-1;
      nodes[i].nson=nodes[i].label=0;
      nodes[i].branch=0;
      /* nodes[i].divtime=0; */
      /* at least for tips, note clock=3 */
   }
   while(isspace(ch)) ch=fgetc(ftree);  /* skip spaces */
   ungetc(ch,ftree);
   if (isdigit(ch)) { ReadaTreeB(ftree,popline); return(0); }

   FOR(i,com.ns) check[i]=0;
   for (;;) {
      ch = fgetc (ftree);
      if (ch==EOF) return(-1);
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
            tree.root=cnode;
         cfather=cnode;
      }
      else if (ch==')')
         { level--;  inodeb=cfather; cfather=nodes[cfather].father; }
      else if (ch==':') {
         haslength=1; fscanf(ftree,"%lf",&nodes[inodeb].branch);
      }
      else if (ch=='#'||ch=='$') {
         haslabel=1; fscanf(ftree,"%d",&nodes[inodeb].label);
      }
      else if (ch==',') ;
      else if (ch==';' && level!=0) error2("; in treefile");
      else { /* read species name or number */
         line[0]=(char)ch;  line[1]=(char)fgetc(ftree);
/*         if(line[1]==(char)EOF) error2("eof in tree file"); */
         if (com.ns<10 && isdigit(line[0]) && isdigit(line[1])) 
           { ungetc(line[1], ftree); line[1]=0; }
         else 
            for (i=1; i<lline; )  { /* read species name until a delimiter */
               if (strchr(delimiters,line[i]) || line[i]==(char)EOF || line[i]=='\n')
                  { ungetc(line[i],ftree); line[i]=0; break; }
            line[++i]=(char)fgetc(ftree);
         }
         for(j=i-1;j>0;j--) if(!isgraph(line[j])) line[j]=0; /* trim spaces*/
         for(i=0,hasname=0; line[i]; i++)  /* name or number? */
            if(!isdigit(line[i])) { hasname=1; break; }  

         if (hasname) {   /* name */
            for (i=0; i<com.ns; i++) if (!strcmp(line,com.spname[i])) break;
            if ((cnode=i)==com.ns) printf("\nSpecies %s?\n", line);
         }
         else {           /* number */
            sscanf(line, "%d", &cnode);   cnode--;
            if (cnode<0 || cnode>=com.ns) {
               printf("\nReadaTreeN: species %d not in data.\n",cnode+1); 
               exit(-1); 
            }
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
      printf ("\nnnode%6d !=  nbranch%6d + 1\n", tree.nnode, tree.nbranch);

/* check that it is o.k. to comment out this line
   com.ntime = com.clock ? (tree.nbranch+1)-com.ns+(tree.root<com.ns)
                         : tree.nbranch;
*/
   FOR(i,com.ns)  if(check[i]!=1) return(-1);
   if(tree.nbranch>2*com.ns-2) { 
      printf("nbranch %d", tree.nbranch); error2("too many branches in tree?");
   }
   *length_label=haslength+2*haslabel;
   return (0);
}


int IncludeNodeLabel=0;

int OutSubTreeN (FILE *fout, int inode, int spnames, int branchlen);

int OutSubTreeN (FILE *fout, int inode, int spnames, int branchlen)
{
   int i,ison;

   fputc ('(', fout);
   FOR (i, nodes[inode].nson) {
      ison=nodes[inode].sons[i];
      if(nodes[ison].nson==0) {
         if(spnames) {
            if(IncludeNodeLabel) fprintf(fout, "%d_",ison+1);
            fprintf(fout,"%s",com.spname[ison]);
         }
         else
            fprintf(fout,"%d",ison+1);
      }
      else
         OutSubTreeN(fout,ison,spnames,branchlen);
      if(IncludeNodeLabel && nodes[ison].nson) fprintf(fout," %d ",ison+1);
      if(branchlen) { 

/*  Add branch labels to be read by Rod Page's TreeView. */
#ifdef CODEML
         if(com.verbose>1 && com.seqtype==1 && com.model && !com.NSsites) {
         fprintf(fout," #%.2f ", nodes[ison].omega);
}
#endif
         fprintf(fout,":%.6f", nodes[ison].branch);
      }

      if(com.ns>9 || spnames || branchlen)
         if(i<nodes[inode].nson-1) fprintf(fout,", ");
   }
   fputc (')', fout);

   /* output node times for treeview */
/*
   if(noisy>=9 && com.ns>9 && com.clock && branchlen)
      fprintf(fout," #%.1f ", nodes[inode].divtime);
*/
   if(noisy>=9 && com.ns>9 && branchlen)
      fprintf(fout," #%d ", inode+1);

   return (0);
}


int OutaTreeN (FILE *fout, int spnames, int branchlen)
{
/* IncludeNodeLabel=1 will label the nodes
*/
   OutSubTreeN(fout,tree.root,spnames,branchlen);  
   if(IncludeNodeLabel) fprintf(fout," %d ", tree.root+1);
   fputc(';',fout);
   return(0);
}

void PointLklnodes (void)
{
/* This points the nodes[com.ns+inode].lkl to the right space in com.lkl.
   The space is different depending on com.cleandata (0 or 1)
   This routine updates internal nodes com.lkl only.  
   End nodes (com.lkl0) are updated in InitPartialLikelihood().
*/
   int nintern=0, i;

   FOR(i,tree.nbranch+1)
      if(nodes[i].nson>0)  /* more thinking */
         nodes[i].lkl = com.lkl + com.ncode*com.npatt*nintern++;
}


double ScaleTimes_clock3=1, YoungDate=0, divtime_oldSon[2*NS-1];
/* for clock=3: mutation=mut/ScaleTimes_clock3; 
                divtime=divtime*ScaleTimes_clock3 */
void SetDivTime(int inode, double x[]);
void SetDivTime_OldSon(int inode);
static int innode_time=0;  
/* number of internal node times, usd to deal with young ancestors.  
   This might be broken.
*/


int GetSeqTimes(void)
{
/* Known divergence times are rescaled by ScaleTimes_clock3
   Andrew Rambaut's TipDate models
*/
   int i,k;
   double young=-1,old=-1,date;

   if(com.clock!=3) return(-1);
   FOR(i,com.ns) {
      k=strlen(com.spname[i])-1;
      for(; k>0 &&(com.spname[i][k]=='.'||isdigit(com.spname[i][k])); k--);
      sscanf(com.spname[i]+k+1, "%lf", &date);
      nodes[i].divtime=date;
      if(date<0) error2("date<0");
      if(i==0) young=old=date;
      else { old=min2(old,date); young=max2(young,date); }
      printf("%9.2f for species %d (%s)\n", date,i+1,com.spname[i]);
   }
   YoungDate=young;
   /* ScaleTimes_clock3=(YoungDate-old)*5; */
   ScaleTimes_clock3=100;
   FOR(i,com.ns) {
      nodes[i].divtime=(YoungDate-nodes[i].divtime)/ScaleTimes_clock3;
      divtime_oldSon[i]=nodes[i].divtime;
   }
   if(noisy) printf("\nDates range: (%.2f, %.2f), (0, %.2f) after scale\n",
                 young,old,(young-old)/ScaleTimes_clock3);

   return(0);
}


void SetDivTime_OldSon(int inode)
{
/* This is for (clock==3), and finds the oldest son for each internal node, 
   using nodes[].divtime for tips (which are initialized in GetSeqTimes).
   Call this in GetInitial as this does not change for the same tree:
   SetDivTime_OldSon(tree.root).
   ntime include the node times as well mutation rate so that ntime==com.ns.
*/
   int i,ison;

   if(com.clock!=3) error2("SetDivTime_OldSon");
   divtime_oldSon[inode]=0;
   for(i=0; i<nodes[inode].nson; i++) {
      ison=nodes[inode].sons[i];
      if(nodes[ison].nson) {
         if(ison<com.ns) puts("ison<ns for clock=3. not implemented.");
         SetDivTime_OldSon(ison);
      }
      divtime_oldSon[inode]=max2(divtime_oldSon[inode],divtime_oldSon[ison]);
   }
   if(inode<com.ns) divtime_oldSon[inode]=nodes[inode].divtime;
}


void SetDivTime(int inode, double x[])
{
/* This is called from SetBranch(), to set up divtime for nodes under clock 
   models
   When clock=3, this routine sets up nodes[].divtime and then SetBranch() sets
   up branch lengths:
      divtime_node i = divtime_oldson+(divtime_father-divtime_oldson)*x[i]
*/
   int i,ison;

   FOR (i,nodes[inode].nson) {
      ison=nodes[inode].sons[i];
      if(nodes[ison].nson) {
         if(com.clock<=2) 
            nodes[ison].divtime=nodes[inode].divtime*x[innode_time++];
         else if(com.clock==3) 
            nodes[ison].divtime=divtime_oldSon[ison]+
               (nodes[inode].divtime-divtime_oldSon[ison]) * x[innode_time++];
         SetDivTime(ison,x);
      }
   }

}


int SetBranch (double x[])
{
   int i, status=0;
   double small=-1.0e-6, mut=(com.clock==3? x[com.ntime-1] : 1);

   if(com.clock==0)
      FOR(i,tree.nnode) {
         if(i!=tree.root) 
            if((nodes[i].branch=x[nodes[i].ibranch])<small)  status=-1;
      }
   else {
      if(LASTROUND==0) {
         innode_time=0;
         nodes[tree.root].divtime=x[0]; innode_time=1;
         SetDivTime(tree.root, x);
      }
      else
         for(i=0; i<tree.nnode-com.ns; i++) 
            nodes[i+com.ns].divtime=x[i];

      FOR (i,tree.nnode) {  /* nodes[].divtime to nodes[].branch */
         /* if(i==tree.root) { innode_time++; continue; } */
         if(i==tree.root) 
            { if(tree.root<com.ns) innode_time++;  continue; }
         nodes[i].branch=(nodes[nodes[i].father].divtime-nodes[i].divtime)*mut;
         if(nodes[i].branch<small)  status=-1;
      }
   }

#if (BASEML || CODEML) /* multiply by rates in local clock models */
   if(com.clock==2)
      FOR(i,tree.nbranch) if(rateBranch[i]) 
         nodes[tree.branches[i][1]].branch*=x[innode_time+rateBranch[i]-1];
#endif

   return(status);
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

   if (inode==tree.root || inode<com.ns) error2("err CollapsNode");
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
         SetBranch (x);
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


void DescentGroup (int inode);
void BranchPartition (char partition[], int parti2B[]);

static char *PARTITION;

void DescentGroup (int inode)
{
   int i;
   for (i=0; i<nodes[inode].nson; i++) 
      if (nodes[inode].sons[i]<com.ns) 
         PARTITION[nodes[inode].sons[i]]=1;
      else 
         DescentGroup (nodes[inode].sons[i]);
}

void BranchPartition (char partition[], int parti2B[])
{
/* calculates branch partitions.
   partition[0,...,ns-1] marks the species bi-partition by the first interior 
   branch.  It uses 0 and 1 to indicate which side of the branch each species 
   is.  
   partition[ns,...,2*ns-1] marks the second interior branch.
   parti2B[0] maps the partition (internal branch) to the branch in tree.
   Use NULL for parti2B if this information is not needed.
   partition[nib*com.ns].  nib: # of interior branches.
*/
   int i,j, nib;  /* number of internal branches */

   for (i=0,nib=0; i<tree.nbranch; i++) {
      if (tree.branches[i][1]>=com.ns){
         PARTITION=partition+nib*com.ns;
         FOR (j,com.ns) PARTITION[j]=0;
         DescentGroup (tree.branches[i][1]);
         if (parti2B) parti2B[nib]=i;
         nib++;
      }
   }
   if (nib!=tree.nbranch-com.ns) error2("err BranchPartition"); 
}


int NSameBranch (char partition1[],char partition2[], int nib1,int nib2,
    int IBsame[])
{
/* counts the number of correct (identical) bipartitions.
   nib1 and nib2 are the numbers of interior branches in the two trees
   correctIB[0,...,(correctbranch-1)] lists the correct interior branches, 
   that is, interior branches in tree 1 that is also in tree 2.
   IBsame[i]=1 if interior branch i is correct.
*/
   int i,j,k, nsamebranch,nsamespecies;

   for (i=0,nsamebranch=0; i<nib1; i++)  for(j=0,IBsame[i]=0; j<nib2; j++) {
      for (k=0,nsamespecies=0;k<com.ns;k++)
         nsamespecies+=(partition1[i*com.ns+k]==partition2[j*com.ns+k]);
      if (nsamespecies==0 || nsamespecies==com.ns)
         { nsamebranch++;  IBsame[i]=1;  break; } 
   }
   return (nsamebranch);
}



int AddSpecies (int is, int ib)
{
/* Add species (is) to tree at branch ib.  The tree currently has 
   is+1-1 species.  Interior node numbers are increased by 2 to make 
   room for the new nodes.
   if(com.clock && ib==tree.nbranch), the new species is added as an
   outgroup to the rooted tree.
*/
   int i,j, it;

   if(ib>tree.nbranch+1 || (ib==tree.nbranch && !com.clock)) return(-1);

   if(ib==tree.nbranch && com.clock) { 
      FOR(i,tree.nbranch) FOR(j,2)
         if (tree.branches[i][j]>=is) tree.branches[i][j]+=2;
      it=tree.root;  if(tree.root>=is) it+=2;
      FOR(i,2) tree.branches[tree.nbranch+i][0]=tree.root=is+1;
      tree.branches[tree.nbranch++][1]=it;
      tree.branches[tree.nbranch++][1]=is;
   }
   else {
      FOR(i,tree.nbranch) FOR(j,2)
         if (tree.branches[i][j]>=is) tree.branches[i][j]+=2;
      it=tree.branches[ib][1];
      tree.branches[ib][1]=is+1;
      tree.branches[tree.nbranch][0]=is+1;
      tree.branches[tree.nbranch++][1]=it;
      tree.branches[tree.nbranch][0]=is+1;
      tree.branches[tree.nbranch++][1]=is;
      if (tree.root>=is) tree.root+=2;
   }
   BranchToNode ();
   return (0);
}


#ifdef TREESEARCH

static struct TREE
  {struct TREEB tree; struct TREEN nodes[2*NS-1]; double x[NP]; } 
  treebest, treestar;
/*
static struct TREE 
  {struct TREEB tree; struct TREEN nodes[2*NS-1];} treestar;
*/

int Perturbation(FILE* fout, int initialMP, double space[]);

int Perturbation(FILE* fout, int initialMP, double space[])
{
/* heuristic tree search by the NNI tree perturbation algorithm.  
   Some trees are evaluated multiple times as no trees are kept.
   This needs more work.
*/
   int step=0, ntree=0, nmove=0, improve=0, ineighb, i;
   int sizetree=(2*com.ns-1)*sizeof(struct TREEN);
   double *x=treestar.x;
   FILE *ftree;

if(com.clock) error2("\n\aerr: pertubation does not work with a clock yet.\n");
if(initialMP&&!com.cleandata)
   error2("\ncannot get initial parsimony tree for gapped data yet.");

   fprintf(fout, "\n\nHeuristic tree search by NNI perturbation\n");
   if (initialMP) {
      if (noisy) printf("\nInitial tree from stepwise addition with MP:\n");
      fprintf(fout, "\nInitial tree from stepwise addition with MP:\n");
      StepwiseAdditionMP (space);
   }
   else {
      if (noisy) printf ("\nInitial tree read from file %s:\n", com.treef);
      fprintf(fout, "\nInitial tree read from file.\n");
      if ((ftree=fopen (com.treef,"r"))==NULL) error2 ("treefile not exist?");
      fscanf (ftree, "%d%d", &i, &ntree);
      if (i!=com.ns) error2 ("ns in the tree file");
      if(ReadaTreeN (ftree, &i, 1)) error2 ("err tree..");
      fclose(ftree);
   }
   if (noisy) { FPN (F0);  OutaTreeN (F0,0,0);  FPN(F0); }
   tree.lnL=TreeScore(x, space);
   if (noisy) { OutaTreeN (F0,0,1);  printf("\n lnL = %.4f\n",-tree.lnL); }
   OutaTreeN (fout,1,1);  fprintf(fout, "\n lnL = %.4f\n",-tree.lnL);
   if (com.np>com.ntime) {
      fprintf(fout, "\tparameters:"); 
      for(i=com.ntime; i<com.np; i++) fprintf(fout, "%9.5f", x[i]);
      FPN(fout);
   }
   fflush(fout);
   treebest.tree=tree;  memcpy(treebest.nodes, nodes, sizetree);

   for (step=0; ; step++) {
      for (ineighb=0,improve=0; ineighb<(tree.nbranch-com.ns)*2; ineighb++) {
         tree=treebest.tree; memcpy (nodes, treebest.nodes, sizetree);
         NeighborNNI (ineighb);
         if(noisy) {
            printf("\nTrying tree # %d (%d move[s]) \n", ++ntree,nmove);
            OutaTreeN (F0,0,0);  FPN(F0);
         }
         tree.lnL=TreeScore(x, space);
         if (noisy) { OutaTreeN(F0,1,1); printf("\n lnL = %.4f\n",-tree.lnL);}
         if (noisy && com.np>com.ntime) {
            printf("\tparameters:"); 
            for(i=com.ntime; i<com.np; i++) printf("%9.5f", x[i]);
            FPN(F0);
         }
         if (tree.lnL<=treebest.tree.lnL) {
            treebest.tree=tree;  memcpy (treebest.nodes, nodes, sizetree);
            improve=1; nmove++;
            if (noisy) printf(" moving to this tree\n");
            if (fout) {
               fprintf(fout, "\nA better tree:\n");
               OutaTreeN(fout,0,0); FPN(fout); OutaTreeN(fout,1,1); FPN(fout); 
               fprintf(fout, "\nlnL = %.4f\n", tree.lnL);
               if (com.np>com.ntime) {
                  fprintf(fout,"\tparameters:"); 
                  for(i=com.ntime; i<com.np; i++) fprintf(fout,"%9.5f", x[i]);
                  FPN(fout);
               }
               fflush(fout);
          }
         }
      }
      if (!improve) break;
   }
   tree=treebest.tree;  memcpy (nodes, treebest.nodes, sizetree);
   if (noisy) {
      printf("\n\nBest tree found:\n");
      OutaTreeN(F0,0,0);  FPN(F0);  OutaTreeN(F0,1,1);  FPN(F0); 
      printf("\nlnL = %.4f\n", tree.lnL);
   }
   if (fout) {
      fprintf(fout, "\n\nBest tree found:\n");
      OutaTreeN(fout,0,0);  FPN(fout);  OutaTreeN(fout,1,1);  FPN(fout); 
      fprintf(fout, "\nlnL = %.4f\n", tree.lnL);
   }
   return (0);
}


static int *_U0, *_step0, _mnnode;
/* up pass characters and changes for the star tree: each of size npatt*nnode*/

int StepwiseAdditionMP (double space[])
{
/* tree search by species addition.
*/
   char *z0[NS];
   int  ns0=com.ns, is, i,j,h, tiestep=0,tie,bestbranch=0;
   int sizetree=(2*com.ns-1)*sizeof(struct TREEN);
   double bestscore=0,score;

   _mnnode=com.ns*2-1;
   _U0=(int*)malloc(com.npatt*_mnnode*sizeof(int));
   _step0=(int*)malloc(com.npatt*_mnnode*sizeof(int));
   if (noisy>2) 
     printf("\n%9ld bytes for MP (U0 & N0)\n", 2*com.npatt*_mnnode*sizeof(int));
   if (_U0==NULL || _step0==NULL) error2 ("oom U0&step0");

   FOR (i,ns0)  z0[i]=com.z[i];
   tree.nbranch=tree.root=com.ns=3;
   FOR (i, tree.nbranch) { tree.branches[i][0]=com.ns; tree.branches[i][1]=i; }
   BranchToNode ();
   FOR (h, com.npatt)
      FOR (i,com.ns)
        { _U0[h*_mnnode+i]=1<<(com.z[i][h]-1); _step0[h*_mnnode+i]=0; }
   for (is=com.ns,tie=0; is<ns0; is++) {
      treestar.tree=tree;  memcpy (treestar.nodes, nodes, sizetree);

      for (j=0; j<treestar.tree.nbranch; j++,com.ns--) {
         tree=treestar.tree;  memcpy (nodes, treestar.nodes, sizetree);
         com.ns++;
         AddSpecies (is, j);
         score=MPScoreStepwiseAddition(is, space, 0);
/*
OutaTreeN(F0, 0, 0); 
printf(" Add sp %d (ns=%d) at branch %d, score %.0f\n", is+1,com.ns,j+1,score);
*/
         if (j && score==bestscore) tiestep=1;
         if (j==0 || score<bestscore || (score==bestscore&&rndu()<.1)) {
            tiestep=0;
            bestscore=score; bestbranch=j;
         }
      }
      tie+=tiestep;
      tree=treestar.tree;  memcpy (nodes, treestar.nodes, sizetree);
      com.ns=is+1;
      AddSpecies (is, bestbranch);
      score=MPScoreStepwiseAddition(is, space, 1);

      if (noisy)
       { printf("\r  Added %d [%5.0f steps]",is+1,-bestscore); fflush(F0);}
/*
OutaTreeN(F0,0,0);
FOR (i,tree.nnode) {
  printf ("\nnode %d: ", i+1);
  FOR (j,com.npatt) printf ("%3d", _step0[j*_mnnode+i]);
}
getchar ();
*/
   }  
   if (noisy>2) printf("  %d stages with ties, ", tie);
   tree.lnL=bestscore;
   free(_U0); free(_step0);
   return (0);
}

double MPScoreStepwiseAddition (int is, double space[], int save)
{
/* this changes only the part of the tree affected by the newly added 
   species is.
   save=1 for the best tree, so that _U0 & _step0 are updated
*/
   int *U,*N,U3[3], h,ist, i,father,son2,*pU0=_U0,*pN0=_step0;
   double score;

   U=(int*)space;  N=U+_mnnode;
   for (h=0,score=0; h<com.npatt; h++,pU0+=_mnnode,pN0+=_mnnode) {
      FOR (i, tree.nnode) { U[i]=pU0[i-2*(i>=is)]; N[i]=pN0[i-2*(i>=is)]; }
      U[is]=1<<(com.z[is][h]-1);  N[is]=0;
      for (ist=is; (father=nodes[ist].father)!=tree.root; ist=father) {
         if ((son2=nodes[father].sons[0])==ist)  son2=nodes[father].sons[1];
         N[father]=N[ist]+N[son2];
         if ((U[father]=U[ist]&U[son2])==0)
            { U[father]=U[ist]|U[son2];  N[father]++; }
      }
      FOR (i,3) U3[i]=U[nodes[tree.root].sons[i]];
      N[tree.root]=2;
      if (U3[0]&U3[1]&U3[2]) N[tree.root]=0;
      else if (U3[0]&U3[1] || U3[1]&U3[2] || U3[0]&U3[2]) N[tree.root]=1;
      FOR(i,3) N[tree.root]+=N[nodes[tree.root].sons[i]];

      if (save) {
         memcpy (pU0, U, tree.nnode*sizeof(int));
         memcpy (pN0, N, tree.nnode*sizeof(int));
      }
      score+=N[tree.root]*com.fpatt[h];
   }
   return (score);
}


int readx(double x[], int *fromfile)
{
/* this reads parameters from file, used as initial values
   if(runmode>0), this reads common substitution parameters only, used for 
   heuristic tree search.
   fromfile=0: if nothing read from file, 1: read from file, -1:fix parameters
*/
   static int times=0;
   int i, npin;
   double *xin;

   times++;  *fromfile=0;
   if(fin==NULL || (com.runmode>0 && times>1)) return(0);
   if(com.runmode<=0) { npin=com.np; xin=x; }
   else               { npin=com.np-com.ntime; xin=xcom; }

   if(npin<=0) return(0);
   if(com.runmode>0&&com.seqtype==1&&com.model) error2("option or in.codeml");
   printf("\nReading initials/paras from file (np=%d). Stop if wrong.\n",npin);
   fscanf(fin,"%lf",&xin[i=0]);
   *fromfile=1;
   if(xin[0]==-1) { *fromfile=-1; LASTROUND=1; }
   else           i++;
   for( ;i<npin;i++) if(fscanf(fin,"%lf",&xin[i])!=1) break;
   if(i<npin)
      { printf("err at #%d. Edit or remove it.\n",i+1); exit(-1); }
   if(com.runmode>0) {
      for(i=0;i<npin;i++) x[com.ntime+i]=xcom[i];
      fclose(fin);  fin=NULL;
      matout(F0,xcom,1,npin);
      puts("Those are fixed for tree search.  Stop if wrong.");
   }
   return(0);
}


double TreeScore(double x[], double space[])
{
   static int fromfile=0;
   int i;
   double xb[NP][2], e=1e-6, lnL=0;

   if(com.clock==2) error2("local clock in TreeScore");
   com.ntime = com.clock ? tree.nnode-com.ns : tree.nbranch;

   GetInitials(x, &i);  /* this shoulbe be improved??? */
   if(i) fromfile=1;
   PointLklnodes();

   if(com.method==0 || !fromfile) SetxBound(com.np, xb);
#ifndef SIMULATE
   if(!com.cleandata) InitPartialLikelihood ();
#endif

   if(fromfile) {
      for(i=com.ntime;i<com.np;i++) x[i]=xcom[i-com.ntime];
      lnL=com.plfun(x,com.np);
      com.np=com.ntime;
   }
   NFunCall=0;
   if(com.method==0 || com.ntime==0)
      ming2(NULL,&lnL,com.plfun,NULL,x,xb, space,e,com.np);
   else
      minB(NULL, &lnL,x,xb, space);

   return(lnL);
}


int StepwiseAddition (FILE* fout, double space[])
{
/* heuristic tree search by species addition.  Species are added in the order 
   of occurrence in the data.
   Try to get good initial values.
*/
   char *z0[NS], *spname0[NS];
   int ns0=com.ns, is, i,j, bestbranch=0, randadd=0, order[NS];
   int sizetree=(2*com.ns-1)*sizeof(struct TREEN);
   double bestscore=0,score, *x=treestar.x;

   if (noisy) printf("\n\nHeuristic tree search by stepwise addition\n");
   if (fout) fprintf(fout, "\n\nHeuristic tree search by stepwise addition\n");
   FOR (i,ns0)  { z0[i]=com.z[i]; spname0[i]=com.spname[i]; }
   tree.nbranch=tree.root=com.ns=(com.clock?2:3);  

   FOR(i,ns0) order[i]=i;
   if(randadd) {
      FOR(i,ns0)
         { j=(int)(ns0*rndu()); is=order[i]; order[i]=order[j]; order[j]=is; }
      if(noisy) FOR(i,ns0) printf(" %d", order[i]+1);
      if(fout) { 
         fputs("\nOrder of species addition:\n",fout); 
         FOR(i,ns0)fprintf(fout,"%3d  %-s\n", order[i]+1,com.spname[order[i]]);
      }
      FOR(i,ns0) { com.z[i]=z0[order[i]]; com.spname[i]=spname0[order[i]]; }
   }

   FOR (i,tree.nbranch) { tree.branches[i][0]=com.ns; tree.branches[i][1]=i; }
   BranchToNode ();
   for (is=com.ns; is<ns0; is++) {                  /* add the is_th species */
      treestar.tree=tree;  memcpy (treestar.nodes, nodes, sizetree);

      for (j=0; j<treestar.tree.nbranch+(com.clock>0); j++,com.ns--) { 
         tree=treestar.tree;  memcpy(nodes, treestar.nodes, sizetree);
         com.ns++;
         AddSpecies(is,j);
         /* recalculate com.pi ? */
         score=TreeScore(x, space);
         if (noisy>1)
            { printf("\n "); OutaTreeN(F0, 0, 0); printf("%12.3f",-score); }

         if (j==0 || score<bestscore || (score==bestscore&&rndu()<.2)) {
            treebest.tree=tree;  memcpy(treebest.nodes, nodes, sizetree);
            xtoy (x, treebest.x, com.np);
            bestscore=score; bestbranch=j;
         }
      }
      tree=treebest.tree;  memcpy(nodes,treebest.nodes, sizetree);
      xtoy (treebest.x, x, com.np);
      com.ns=is+1;

      if (noisy) {
         printf("\n\nAdded sp. %d, %s [%.3f]\n",is+1,com.spname[is],-bestscore);
         OutaTreeN(F0,0,0);  FPN(F0);  OutaTreeN(F0,1,0);  FPN(F0);
         if (com.np>com.ntime) {
            printf("\tparameters:"); 
            for(i=com.ntime; i<com.np; i++) printf("%9.5f", x[i]);
            FPN(F0);
         }
      }
      if (fout) {
         fprintf(fout,"\n\nAdded sp. %d, %s [%.3f]\n",
                 is+1, com.spname[is], -bestscore);
         OutaTreeN(fout,0,0); FPN(fout);
         OutaTreeN(fout,1,1); FPN(fout);
         if (com.np>com.ntime) {
            fprintf(fout, "\tparameters:"); 
            for(i=com.ntime; i<com.np; i++) fprintf(fout, "%9.5f", x[i]);
            FPN(fout);
         }
         fflush(fout);
      }
   }  
   tree.lnL=bestscore;
   /*
fprintf(frst1," %12.5f  ",-bestscore); 
OutaTreeN(frst1,1,1); FPN(frst1); fflush(frst1);
*/
   return (0);
}


int DecompTree (int inode, int ison1, int ison2);
#define hdID(i,j) (max2(i,j)*(max2(i,j)-1)/2+min2(i,j))

int StarDecomposition (FILE *fout, double space[])
{
/* automatic tree search by star decomposition, nhomo<=1
   returns (0,1,2,3) for the 4s problem.
*/
   int status=0,stage=0, i,j, itree,ntree=0,ntreet,best=0,improve=1,collaps=0;
   int inode, nson=0, ison1,ison2, son1, son2;
   int sizetree=(2*com.ns-1)*sizeof(struct TREEN);
   double x[NP];
   FILE *ftree, *fsum=frst;
/*
fsum=NULL;
*/
   if (com.runmode==1) {   /* read the star-like tree from tree file */
      if ((ftree=fopen (com.treef,"r"))==NULL) error2 ("no treefile");
      fscanf (ftree, "%d%d", &i, &ntree);
      if (ReadaTreeN(ftree, &i, 1)) error2 ("err tree file");
      fclose (ftree);
   }
   else {                  /* construct the star tree of ns species */
      tree.nnode = (tree.nbranch=tree.root=com.ns)+1;
      for (i=0; i<tree.nbranch; i++)
         { tree.branches[i][0]=com.ns; tree.branches[i][1]=i; }
      com.ntime = com.clock?1:tree.nbranch;
      BranchToNode ();
   }
   if (noisy) { printf("\n\nstage 0: ");  OutaTreeN(F0,0,0); }
   if (fsum) { fprintf(fsum,"\n\nstage 0: ");  OutaTreeN(fsum,0,0); }
   if (fout) { fprintf(fout,"\n\nstage 0: ");  OutaTreeN(fout,0,0); }

   tree.lnL=TreeScore(x,space);

   if (noisy)  printf("\nlnL:%14.6f%6d", -tree.lnL, NFunCall);
   if (fsum) fprintf(fsum,"\nlnL:%14.6f%6d", -tree.lnL, NFunCall);
   if (fout) {
      fprintf(fout,"\nlnL(ntime:%3d  np:%3d):%14.6f\n",
         com.ntime, com.np, -tree.lnL);
      OutaTreeB (fout);  FPN(fout);
      FOR (i, com.np) fprintf (fout,"%9.5f", x[i]);  FPN (fout);
   }
   treebest.tree=tree;  memcpy(treebest.nodes,nodes,sizetree);
   FOR (i,com.np) treebest.x[i]=x[i];
   for (ntree=0,stage=1; ; stage++) {
      for (inode=treebest.tree.nnode-1; inode>=0; inode--) {
         nson=treebest.nodes[inode].nson;
         if (nson>3) break;
         if (com.clock) { if (nson>2) break; }
         else if (nson>2+(inode==treebest.tree.root)) break;
      }
      if (inode==-1 || /*stage>com.ns-3+com.clock ||*/ !improve) { /* end */
         tree=treebest.tree;  memcpy (nodes, treebest.nodes, sizetree);

         if (noisy) {
            printf("\n\nbest tree: ");  OutaTreeN(F0,0,0);
            printf("   lnL:%14.6f\n", -tree.lnL);
         }
         if (fsum) {
            fprintf(fsum, "\n\nbest tree: ");  OutaTreeN(fsum,0,0);
            fprintf(fsum, "   lnL:%14.6f\n", -tree.lnL);
         }
         if (fout) {
            fprintf(fout, "\n\nbest tree: ");  OutaTreeN(fout,0,0);
            fprintf(fout, "   lnL:%14.6f\n", -tree.lnL);
            OutaTreeN(fout,1,1);  FPN(fout);
         }
         break;
      }
      treestar=treebest;  memcpy(nodes,treestar.nodes,sizetree);

      if (collaps && stage) { 
         printf ("\ncollapsing nodes\n");
         OutaTreeN (F0, 1, 1);  FPN(F0);

         tree=treestar.tree;  memcpy(nodes, treestar.nodes, sizetree);
         for (i=com.ns,j=0; i<tree.nnode; i++)
            if (i!=tree.root && nodes[i].branch<1e-7) 
               { CollapsNode (i, treestar.x);  j++; }
         treestar.tree=tree;  memcpy(treestar.nodes, nodes, sizetree);

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
      if (!com.clock && inode==treestar.tree.root && nson==4)  ntreet=3;
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
               x[i]=max2(nodes[tree.branches[i][1]].branch*0.99, 0.0001);
         else
            for (i=1,x[0]=max2(x[0],.01); i<com.ntime; i++)  x[i]=.5;

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
         tree.lnL=TreeScore(x, space);

         if (tree.lnL<treebest.tree.lnL) {
            treebest.tree=tree;  memcpy (treebest.nodes, nodes, sizetree);
            FOR(i,com.np) treebest.x[i]=x[i];
            best=itree+1;   improve=1;
         }
         if (noisy) printf("%6d%2c %+8.2f", 
                       NFunCall,(status?'?':'X'),treestar.tree.lnL-tree.lnL);
         if (fsum) {
            fprintf(fsum, "%6d%2c", NFunCall, (status?'?':'X'));
            for (i=com.ntime; i<com.np; i++)  fprintf(fsum, "%7.3f", x[i]);
            fprintf(fsum, " %+8.2f", treestar.tree.lnL-tree.lnL);
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
   int sizetree=(2*com.ns-1)*sizeof(struct TREEN);
   double bt, fmid=0.001, fclock=0.0001;

   tree=treestar.tree;  memcpy (nodes, treestar.nodes, sizetree);
   for (i=0,bt=0; i<tree.nnode; i++)
      if (i!=tree.root) bt+=nodes[i].branch/tree.nbranch;

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
/* This does the separate analysis of multiple-gene data.
   Note that com.pose[] is not correct and so RateAncestor = 0 should be set
   in baseml and codeml.
*/
   int ig=0, j, ngene0, npatt0, lgene0[NGENE], posG0[NGENE+1];
   int nb = ((com.seqtype==1 && !com.cleandata) ? 3 : 1);

   ngene0=com.ngene;  npatt0=com.npatt;
   FOR (ig, ngene0)   lgene0[ig]=com.lgene[ig];
   FOR (ig, ngene0+1) posG0[ig]=com.posG[ig];

   ig=0;
/*
   printf("\nStart from gene (1-%d)? ", com.ngene);
   scanf("%d", &ig); 
   ig--;
*/

   for ( ; ig<ngene0; ig++) {
      com.ngene=1;
      com.ls=com.lgene[0]= ig==0?lgene0[0]:lgene0[ig]-lgene0[ig-1];
      com.npatt =  ig==ngene0-1 ? npatt0-posG0[ig] : posG0[ig+1]-posG0[ig];
      com.posG[0]=0;  com.posG[1]=com.npatt;
      FOR (j,com.ns) com.z[j]+=posG0[ig]*nb;   com.fpatt+=posG0[ig];
      xtoy (com.piG[ig], com.pi, com.ncode);

      printf ("\n\nGene #%4d  ls:%4d  npatt:%4d\n",ig+1,com.ls,com.npatt);
      fprintf(fout,"\nGene #%4d  ls:%4d  npatt:%4d\n",ig+1,com.ls,com.npatt);
      fprintf(frst,"\nGene #%4d  ls:%4d  npatt:%4d\n",ig+1,com.ls,com.npatt);

#ifdef CODEML
      if(com.seqtype==CODONseq) DistanceMatNG86(fout,0);
#endif

      if (com.runmode==0)  Forestry(fout);
#ifdef CODEML
      else if (com.runmode==-2) {
         if(com.seqtype==CODONseq) PairwiseCodon(fout,space);
         else                      PairwiseAA(fout);
      }
#endif
      else                         StepwiseAddition(fout, space);

      FOR (j,com.ns) com.z[j]-=posG0[ig]*nb;
      com.fpatt-=posG0[ig];
   }
   com.ngene=ngene0;  com.npatt=npatt0;  com.ls=lgene0[ngene0-1];
   FOR(ig,ngene0)   com.lgene[ig]=lgene0[ig];
   FOR(ig,ngene0+1) com.posG[ig]=posG0[ig];
   return (0);
}


int printSeqsMgenes (void)
{
/* separate sites from different partitions (genes) into different files.
   called before sequences are coded.
   In case of codons, com.ls is the number of codons.
*/
   FILE *fseq;
   char seqf[20];
   int ig, lg, i,j,h;
   int n31=(com.seqtype==CODONseq||com.seqtype==CODON2AAseq?3:1);

   puts("Separating sites in genes into different files.\n");
   for (ig=0, FPN(F0); ig<com.ngene; ig++) {
      for (h=0,lg=0; h<com.ls; h++)  if (com.pose[h]==ig) lg++;
      sprintf (seqf, "Gene%d.seq\0", ig+1);
      if ((fseq=fopen(seqf,"w"))==NULL) error2("file creation err.");
      printf ("%d sites in gene %d go to file %s\n", lg, ig+1,seqf);

      fprintf (fseq, "%8d%8d\n", com.ns, lg*n31);
      for (j=0; j<com.ns; FPN(fseq),j++) {
         fprintf(fseq,"%-20s  ", com.spname[j]);
         if (n31==1)  {       /* nucleotide or aa sequences */
            FOR (h,com.ls)   
               if (com.pose[h]==ig) fprintf (fseq, "%c", com.z[j][h]);
         }
         else {               /* codon sequences */
            FOR (h,com.ls)
               if (com.pose[h]==ig) {
                  FOR (i,3) fprintf(fseq,"%c", com.z[j][h*3+i]);
                  fputc(' ',fseq);
               }
         }
      }
      fclose (fseq);
   }
   return (0);
}

int printSeqsMgenes2 (void)
{
/* This print sites from certain genes into one file.
   called before sequences are coded.
   In case of codons, com.ls is the number of codons.
*/
   FILE *fseq;
   char seqf[20]="newseqs";
   int ig, lg, i,j,h;
   int n31=(com.seqtype==CODONseq||com.seqtype==CODON2AAseq?3:1);
   
   int ngenekept=0;
   char *genenames[44]={"atpa", "atpb", "atpe", "atpf", "atph", "petb", "petg", "psaa",
"psab", "psac", "psaj", "psba", "psbb", "psbc", "psbd", "psbe",
"psbf", "psbh", "psbi", "psbj", "psbk", "psbl", "psbn", "psbt",
"rl14", "rl16", "rl2", "rl20", "rl36", "rpob", "rpoc", "rpod", "rs11",
"rs12", "rs14", "rs18", "rs19", "rs2", "rs3", "rs4", "rs7", "rs8",
"ycf4", "ycf9"};
   int wantgene[44]={0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                     0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                     0, 0, 0, 0};
/*
for(ig=0,lg=0; ig<com.ngene; ig++) wantgene[ig]=!wantgene[ig];
*/

   if(com.ngene!=44) error2("ngene!=44");
   FOR(h,com.ls) { 
      printf("%3d",com.pose[h]); 
      if((h+1)%20==0) FPN(F0); if((h+1)%500==0) getchar();
   }
   matIout(F0,com.lgene,1,com.ngene);
   matIout(F0,wantgene,1,com.ngene);

   for(ig=0,lg=0; ig<com.ngene; ig++) 
      if(wantgene[ig]) { ngenekept++; lg+=com.lgene[ig]; }

   if((fseq=fopen(seqf,"w"))==NULL) error2("file creation err.");
   fprintf(fseq,"%4d %4d  G\nG  %d  ", com.ns, lg*n31, ngenekept);
   FOR(ig,com.ngene) if(wantgene[ig]) fprintf(fseq," %3d", com.lgene[ig]);
   FPN(fseq);

   for (j=0; j<com.ns; FPN(fseq),j++) {
      fprintf(fseq,"%-20s  ", com.spname[j]);
      if (n31==1)  {       /* nucleotide or aa sequences */
         FOR (h,com.ls)   
            if(wantgene[ig=com.pose[h]]) fprintf(fseq,"%c",com.z[j][h]);
      }
      else {               /* codon sequences */
         FOR (h,com.ls)
            if (wantgene[ig=com.pose[h]]) {
               FOR (i,3) fprintf(fseq,"%c", com.z[j][h*3+i]);
               fputc(' ', fseq);
            }
      }
   }
   FPN(fseq); 
   FOR(ig,com.ngene) if(wantgene[ig]) fprintf(fseq," %s", genenames[ig]);
   FPN(fseq);
   fclose (fseq);

exit(0);

   return (0);
}

#endif   /* ifdef REALSEQUENCE */
#endif   /* ifdef TREESEARCH */
#endif   /* ifdef NODESTRUCTURE */



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
      if (K[i]==maxK)  chU[inode*com.ncode+NchU[inode]++]=(char)i;
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


double MPScore (double space[])
{
/* calculates MP score for a given tree using Hartigan's (1973) algorithm.
   sizeof(space) = nnode*sizeof(int)+(nnode+2)*ncode*sizeof(char).
   Uses Nsteps[nnode], chU[nnode*ncode], NchU[nnode].
   if(BitOperation), bit operations are used on binary trees.
*/
   int h,i, BitOperation,U[3],change;
   double score;

   Nsteps=(int*)space;
   BitOperation=(tree.nnode==2*com.ns-1 - (nodes[tree.root].nson==3));
   BitOperation=(BitOperation&&com.ncode<32);
   if (BitOperation)  chUB=Nsteps+tree.nnode;
   else {
      chU=(char*)(Nsteps+tree.nnode);
      NchU=chU+tree.nnode*com.ncode;  Kspace=NchU+tree.nnode;
   }
   for (h=0,score=0; h<com.npatt; h++) {
      FOR (i,tree.nnode) Nsteps[i]=0;
      if (BitOperation) { 
         FOR (i,com.ns)  chUB[i]=1<<(com.z[i][h]);
         UpPassScoreOnlyB (tree.root);
         if (nodes[tree.root].nson>2) {
            FOR (i,3) U[i]=chUB[nodes[tree.root].sons[i]];
            change=2;
            if (U[0]&U[1]&U[2]) change=0;
            else if (U[0]&U[1] || U[1]&U[2] || U[0]&U[2]) change=1;
            for (i=0,Nsteps[tree.root]=change; i<3; i++) 
               Nsteps[tree.root]+=Nsteps[nodes[tree.root].sons[i]];
       }
      }
      else {                   /* polytomies, use characters */
         FOR(i,com.ns)
            {chU[i*com.ncode]=(char)(com.z[i][h]); NchU[i]=(char)1; }
         for (i=com.ns; i<tree.nnode; i++)  NchU[i]=0;
         UpPassScoreOnly (tree.root);
      }
      score+=Nsteps[tree.root]*com.fpatt[h];
/*
printf("\nh %3d:    ", h+1);
FOR(i,com.ns) printf("%2d  ", com.z[i][h]);
printf(" %6d ", Nsteps[tree.root]);
if((h+1)%10==0) exit(1);
*/
   }

   return (score);
}

double RemoveMPNinfSites (double *nsiteNinf)
{
/* Removes parsimony-noninformative sites and return the number of changes 
   at those sites.
   Changes .z[], .fpatt[], .npatt, etc.
*/
   int  h,j, it, npatt0=com.npatt, markb[NCODE], gt2;
   double MPScoreNinf;

   for (h=0,com.npatt=0,MPScoreNinf=0,*nsiteNinf=0; h<npatt0; h++) {
      FOR (j, com.ncode) markb[j]=0;
      FOR (j, com.ns)  markb[com.z[j][h]]++;
      for (j=0,it=gt2=0; j<com.ncode; j++)
         if (markb[j]>=2) { it++; gt2=1; }
      if (it<2) {                         /* non-informative */
       *nsiteNinf+=com.fpatt[h];
         FOR (j,com.ncode) if(markb[j]==1) MPScoreNinf+=com.fpatt[h];
         if (!gt2) MPScoreNinf-=com.fpatt[h];
      }
      else {
         FOR (j, com.ns) com.z[j][com.npatt]=com.z[j][h];
         com.fpatt[com.npatt++]=com.fpatt[h];
      }
   }
   return (MPScoreNinf);
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
             (char)(chMarkU[ison*n+j]||(chMark[inode*n+j]&&chMarkL[ison*n+j]));
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
         if (chMarkU[ison*n+j] || j==chi) CharaCur[in*n+NCharaCur[in]++]=(char)j;
   }
   else {
      for (j=0,NCharaCur[in]=0; j<n; j++) 
         if (chMarkU[ison*n+j]) CharaCur[in*n+NCharaCur[in]++]=(char)j;
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
   FOR (i,com.ns)  chMark[i*n+com.z[i][h]]=chMarkU[i*n+com.z[i][h]]=1;
   UpPass (tree.root);
   *nchange=Nsteps[tree.root];
   if (job==0) return (0);
   FOR (i,n) chMark[tree.root*n+i]=chMarkU[tree.root*n+i];
   DownPass (tree.root);
   FOR (i,tree.nnode-com.ns) 
      for (j=0,NChara[i]=0; j<n; j++) 
         if (chMark[(i+com.ns)*n+j])  Chara[i*n+NChara[i]++]=(char)j;
   return (0);     
}

int PathwayMP (FILE *fout, double space[])
{
/* Hartigan, JA.  1973.  Minimum mutation fits to a given tree. 
   Biometrics, 29:53-65.
*/
   char *pch=(com.seqtype==0?BASEs:AAs), visit[NS-1];
   int n=com.ncode, nid=tree.nbranch-com.ns+1, it, i,j,k, h, npath;
   int nchange, nchange0;
   char nodeb[NNODE], Equivoc[NS-1];

   PATHWay=(char*)malloc(nid*(n+3)*sizeof(char));
   NCharaCur=PATHWay+nid;  ICharaCur=NCharaCur+nid;  CharaCur=ICharaCur+nid;

   for (j=0,visit[i=0]=(char)(tree.root-com.ns); j<tree.nbranch; j++) 
     if (tree.branches[j][1]>=com.ns) 
        visit[++i]=(char)(tree.branches[j][1]-com.ns);
/*
   printf ("\nOrder in nodes: ");
   FOR (j, nid) printf ("%4d", visit[j]+1+com.ns); FPN(F0);
*/
   for (h=0; h<com.npatt; h++) {
      fprintf (fout, "\n%4d%6.0f  ", h+1, com.fpatt[h]);
      FOR (j, com.ns) fprintf (fout, "%c", pch[com.z[j][h]]);
      fprintf (fout, ":  ");

      FOR (j,com.ns) nodeb[j]=(char)(com.z[j][h]);

      InteriorStatesMP (1, h, &nchange, NCharaCur, CharaCur, space); 
      ICharaCur[j=tree.root-com.ns]=0;  PATHWay[j]=CharaCur[j*n+0];
      FOR (j,nid) Equivoc[j]=(char)(NCharaCur[j]>1);
      DownStates (tree.root);

      for (npath=0; ;) {
         for (j=0,k=visit[nid-1]; j<NCharaCur[k]; j++) {
            PATHWay[k]=CharaCur[k*n+j]; npath++; 
            FOR (i, nid) fprintf (fout, "%c", pch[PATHWay[i]]);
            fprintf (fout, "  ");

            FOR (i,nid) nodeb[i+com.ns]=PATHWay[i];
            for (i=0,nchange0=0; i<tree.nbranch; i++) 
            nchange0+=(nodeb[tree.branches[i][0]]!=nodeb[tree.branches[i][1]]);
            if (nchange0!=nchange) 
               { puts("\a\nerr:PathwayMP"); fprintf(fout,".%d. ", nchange0);}

         }
         for (j=nid-2; j>=0; j--) {
            if(Equivoc[k=visit[j]] == 0) continue;
            if (ICharaCur[k]+1<NCharaCur[k]) {
               PATHWay[k] = CharaCur[k*n + (++ICharaCur[k])];
               DownStates (k+com.ns);
               break;
            }
            else { /* if (next equivocal node is not ancestor) update node k */
               for (i=j-1; i>=0; i--) if (Equivoc[(int)visit[i]]) break;
               if (i>=0) { 
                  for (it=k+com.ns,i=visit[i]+com.ns; ; it=nodes[it].father)
                     if (it==tree.root || nodes[it].father==i) break;
                  if (it==tree.root)
                     DownStatesOneNode(k+com.ns, nodes[k+com.ns].father);
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



#if(BASEML || CODEML)


int BootstrapSeq (char* seqf)
{
/* Generates bootstrap samples.
   This should be used right after ReadSeq() and before the data are coded
   jackknife if(lsb<com.ls && com.ngene==1)
   gmark[] lists sites gene by gene
*/
   int ib,nboot=500, h,is,ig,lg,j, start, lsb=com.ls, n31=1,gap=10;
   int *sites=(int*)malloc(com.ls*sizeof(int)), *gmark=NULL;
   int gcounts[NGENE];
   FILE *fseq=(FILE*)fopen(seqf,"w");

   printf("\nGenerating bootstrap samples in file %s\n", seqf);
   if(com.seqtype==CODONseq||com.seqtype==CODON2AAseq) { n31=3; gap=1; }
   if(sites==NULL) error2("oom in BootstrapSeq");
   if(fseq==NULL) error2("BootstrapSeq: file open error2");
   if(com.ngene>1) {
      if(lsb<com.ls) error2("jackknife when #gene>1");
      if((gmark=(int*)malloc(com.ls*sizeof(int)))==NULL) error2("oom in boot");
      if(noisy) puts("Stratefied sampling");

      FOR(ig,com.ngene) com.lgene[ig]=gcounts[ig]=0;
      FOR(h,com.ls) com.lgene[com.pose[h]]++;
      for(j=1; j<com.ngene; j++) com.lgene[j]+=com.lgene[j-1];

      FOR(h,com.ls) {     /* create gmark[] */
         ig=com.pose[h];
         start=(ig==0?0:com.lgene[ig-1]);
         gmark[start + gcounts[ig]++] = h;
      }
   }

   for (ib=0; ib<nboot; ib++,FPN(fseq)) {
      if(com.ngene<=1)
         FOR(h,lsb) sites[h]=(int)(rndu()*com.ls);
      else {
         FOR(ig,com.ngene) {
            start=(ig==0?0:com.lgene[ig-1]);
            lg=(ig==0?com.lgene[ig]:com.lgene[ig]-com.lgene[ig-1]);
            FOR(h,lg)  sites[start+h]=gmark[start+(int)(rndu()*lg)];
         }
      }

      fprintf(fseq,"%6d %6d\n", com.ns, lsb*n31);
      for(is=0;is<com.ns;is++,FPN(fseq)) {
         fprintf(fseq,"%-20s  ", com.spname[is]);
         FOR(h,lsb) {
            FOR(j,n31) fprintf(fseq,"%c", com.z[is][sites[h]*n31+j]);
            if((h+1)%gap==0) fprintf(fseq," ");
         }
      }
      if(noisy && (ib+1)%10==0) printf("\rdid sample #%d", ib+1);
   }
   free(sites);  if(com.ngene>1) free(gmark);
   return(0);
}


int rell(FILE*flnf, FILE*fout, int ntree)
{
/* This implements three methods for tree topology comparison.  The first 
   tests the log likelihood difference using a normal approximation 
   (Kishino and Hasegawa 1989).  The second does approximate bootstrap sampling
   (the RELL method, Kishino and Hasegawa 1989, 1993).  The third is a 
   modification of the K-H test with a correction for multiple comparison 
   (Shimodaira and Hasegawa 1999) .
   The routine reads input from the file lnf.
   pattg[lgene] lists site patterns for gene ig for stratefied sampling
  
   com.space[ntree*(npatt+nr+4)]: 
    lnf[ntree*npatt] lnL0[ntree] lnL[ntree*nr] pRELL[ntree] pSH[ntree] vdl[ntree]
*/
   char line[1024];
   int lline=1023, ntree0,ls0,npatt0;
   int nr=10000, itree, h,ir,j,k, ig,lgene,npattg, btree, *pattg, status=0;
   double *lnf, *lnL0, *lnL, *pRELL, *pSH, *vdl, y, mdl;

   puts("\nTree comparisons (Kishino & Hasegawa 1989; Shimodaira & Hasegawa 1999)");
   fputs("\nTree comparisons (Kishino & Hasegawa 1989; Shimodaira & Hasegawa 1999)\n",fout);
   fprintf(fout,"Number of replicates: %d\n", nr);

   fscanf(flnf,"%d%d%d", &ntree0, &ls0, & npatt0);
   if (ntree0!=ntree)  error2("rell: input data file strange.  Check.");
   if (ls0!=com.ls || npatt0!=com.npatt)
      error2("rell: input data file incorrect.");
   k = ntree*(com.npatt+nr+4)*sizeof(double) + com.ls*sizeof(int);
   if(com.sspace<k) {
      if((com.space=(double*)realloc(com.space,com.sspace=k))==NULL)
         error2("oom space");
      com.sspace=k;
      if(noisy) printf("space reset to %d bytes in rell.\n",k);
   }
   lnf=com.space; lnL0=lnf+ntree*com.npatt; lnL=lnL0+ntree; pRELL=lnL+ntree*nr;
   pSH=pRELL+ntree; vdl=pSH+ntree; pattg=(int*)(vdl+ntree);

   /* read lnf from file flnf, calculates lnL0[] & find ML tree */
   for(itree=0,btree=0,FPN(F0);itree<ntree; itree++) {
      printf("\rReading tree # %d", itree+1);
      fscanf(flnf,"%d",&k);
      if(k!=itree+1) { printf("\nerr: lnf, reading tree %d.",itree+1); return(-1); }
      for(h=0,lnL0[itree]=0; h<com.npatt; h++) {
         fscanf (flnf, "%d%d%lf", &j, &k, &y);
         if(j!=h+1) { printf("\nlnf, patt %d.",h+1); return(-1); }
         fgets(line,lline,flnf);
         lnL0[itree]+=com.fpatt[h]*(lnf[itree*com.npatt+h]=y);
      }
      if(itree && lnL0[itree]>lnL0[btree]) btree=itree;
   }
   FPN(F0);
   /* calculates SEs (vdl) by sitewise comparison */
   FOR(itree,ntree) {
      if(itree==btree) { vdl[itree]=0; continue; }
      mdl=(lnL0[itree]-lnL0[btree])/com.ls;
      for(h=0,vdl[itree]=0; h<com.npatt; h++) {
         y=lnf[itree*com.npatt+h]-lnf[btree*com.npatt+h];
         vdl[itree]+=com.fpatt[h]*(y-mdl)*(y-mdl);
      }
      vdl[itree]=sqrt(vdl[itree]);
   }
   /* bootstrap resampling */
   if(com.ngene==1)
      for(h=0,k=0;h<com.npatt;h++) FOR(j,com.fpatt[h]) pattg[k++]=h;
   zero(pRELL,ntree); zero(pSH,ntree); zero(lnL,ntree*nr);
   for(ir=0; ir<nr; ir++) {
      for(ig=0; ig<com.ngene; ig++) {
         lgene=(ig?com.lgene[ig]-com.lgene[ig-1]:com.lgene[ig]);
         npattg=com.posG[ig+1]-com.posG[ig];
         if(com.ngene>1) {
            for(h=com.posG[ig],k=0;h<com.posG[ig+1];h++)
               FOR(j,com.fpatt[h]) pattg[k++]=h;
            if(k!=lgene) error2("rell: # sites for gene");
         }
         FOR(j,lgene) {
            h=pattg[(int)(lgene*rndu())];
            FOR(itree,ntree) lnL[itree*nr+ir]+=lnf[itree*com.npatt+h];
         }
      }
      for(j=1,itree=0; j<ntree; j++) 
         if(lnL[j*nr+ir]>lnL[itree*nr+ir]) itree=j;
         else if(fabs(lnL[j*nr+ir]-lnL[itree*nr+ir])<1e-5 && rndu()>.5) itree=j;
      pRELL[itree]+=1./nr;
      if((ir+1)%500==0) printf("\rBootstrapping.. #%6d /%6d",ir+1,nr);
   }

   /* Shimodaira & Hasegawa correction (1999), working on lnL[ntree*nr] */
   for(j=0; j<ntree; j++)  /* step 3: centering */
      for(ir=0,y=sum(lnL+j*nr,nr)/nr; ir<nr; ir++) lnL[j*nr+ir]-=y;
   for(itree=0; itree<ntree; itree++) {  /* steps 4 & 5 */
      for(ir=0; ir<nr; ir++) {
         /* y is the best lnL in replicate ir.  The calculation is repeated */
         for(j=1,y=lnL[ir]; j<ntree; j++) if(lnL[j*nr+ir]>y) y=lnL[j*nr+ir];
         if(y-lnL[itree*nr+ir]>lnL0[btree]-lnL0[itree]) pSH[itree]+=1./nr;
      }
   }
   
   fprintf(fout,"\n%6s %12s %9s %9s%8s%10s%9s\n\n",
      "tree","li","Dli"," +- SE","pKH","pSH","pRELL");
   FOR(j,ntree) {
      mdl=lnL0[j]-lnL0[btree]; 
      if(j==btree || fabs(vdl[j])<1e-6) { y=-1; pSH[j]=-1; status=-1; }
      else y=1-CDFNormal(-mdl/vdl[j]);
      fprintf(fout,"%6d%c%12.3f %9.3f %9.3f%8.3f%10.3f%9.3f\n",
           j+1,(j==btree?'*':' '),lnL0[j],mdl,vdl[j],y,pSH[j],pRELL[j]);
   }
   fputs("\npKH: P value for KH normal test (Kishino & Hasegawa 1989)\n",fout);
   fputs("pRELL: RELL bootstrap proportions (Kishino & Hasegawa 1989) (watch out for ties)\n",fout);
   fputs("pSH: P value with multiple-comparison correction (MC in table 1 of Shimodaira & Hasegawa 1999)\n",fout);
   if(status) fputs("(-1 for P values means N/A)\n",fout);

   FPN(F0);
   return(0);
}

#endif


#ifdef LFUNCTIONS
#ifdef RECONSTRUCTION


int ProbSitePattern (double x[], double *lnL, double fhs[]);
int AncestralMarginal (FILE *fout, double x[], double fhs[]);
int AncestralJointPPSG2000 (FILE *fout, double x[], double fhs[]);

int ProbSitePattern (double x[], double *lnL, double fhs[])
{
/* This calculates probabilities for observed site patterns.  
   It works with constant or gamma rates.
   This does not deal with underflows.
*/
   int ig, i,h, ir;
   double fh;

   if (SetParameters(x)) puts ("par err.");
   zero (fhs, com.npatt);
   if (com.alpha==0) {
      for (ig=0,*lnL=0; ig<com.ngene; ig++) {
         if(com.Mgene>1) SetPGene(ig, 1, 1, 0, x+com.ntime);
         PartialLikelihood (tree.root, ig);
         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            for (i=0; i<com.ncode; i++) 
               fhs[h] += com.pi[i]*nodes[tree.root].lkl[h*com.ncode+i];
            *lnL -= log(fhs[h])*com.fpatt[h];
         }
      }
   }
   else {
      for (ig=0; ig<com.ngene; ig++) {
         if(com.Mgene>1)
            SetPGene(ig, com.Mgene>1, com.Mgene>1, com.nalpha>1, x+com.ntime);
         for (ir=0; ir<com.ncatG; ir++) {
            SetPSiteClass(ir, x+com.ntime);
            PartialLikelihood (tree.root, ig);
            for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
               for (i=0,fh=0; i<com.ncode; i++)
                  fh += com.pi[i]*nodes[tree.root].lkl[h*com.ncode+i];
               fhs[h] += com.freqK[ir]*fh;
            }
         }
      }
      for (h=0,*lnL=0; h<com.npatt; h++)
         *lnL-=log(fhs[h])*com.fpatt[h];
   }
   return (0);
}


void ListAncestSeq(FILE *fout, char *zanc)
{
/* zanc[nintern*com.npatt] holds ancestral sequences.
   Extant sequences are coded if cleandata.
*/
   int wname=15, j;

   fputs("\n\n\nList of extant and reconstructed sequences\n\n",fout);
   for(j=0;j<com.ns;j++,FPN(fout)) {
      fprintf(fout,"%-*s ", wname,com.spname[j]);
      print1seq(fout,com.z[j],com.ls,com.cleandata, com.pose);
   }
   for(j=0;j<tree.nnode-com.ns;j++,FPN(fout)) {
      fprintf(fout,"node #%-*d", wname-5,com.ns+j+1);
      print1seq(fout, zanc+j*com.npatt,com.ls,1, com.pose);
   }
}


int AncestralMarginal (FILE *fout, double x[], double fhs[])
{
/* Ancestral reconstruction for each interior node.  This works under both 
   the one rate and gamma rates models.
   pnode[npatt*nid] stores the prob for the best chara at a node and site.
   The best character is kept in za[], coded as 0,...,n-1.
   The data may be coded (com.cleandata==1) or not (com.cleandata=0).
   Call ProbSitePatt() before running this routine.
   pMAPnode[NS-1], pMAPnodeA[] stores the MAP probabilities (accuracy)
   for a site and for the entire sequence, respectively.
   pChar1node[npatt*ncode] stores prob for each char at each pattern at 
   the current node (in).
*/
   char *pch=(com.seqtype==0?BASEs:(com.seqtype==2?AAs:BINs)), *za;
   char str[4]="",codon[4]="", aa[4]="";
   int n=com.ncode, inode,ir, ic=0,b[3],i,j,k1=-1,k2,k3;
   int h,hp,ig, best, oldroot=tree.root;
   int nid=tree.nnode-com.ns, nchange;
   double lnL=0, fh, y, pbest, *pChar1node, *pnode, p1,p2=-1;
   double pMAPnode[NS-1], pMAPnodeA[NS-1], smallp=0.0001;

   char coding=0, *bestAA=NULL;
   double pAA[21], *pbestAA=NULL;
    /* bestAA[nid*npatt], pbestAA[nid*npatt]: 
       To reconstruct aa seqs using codon or nucleotide seqs, universal code */

   if(_nnodeScale) { puts("Ancestral: not implemented with scaling."); return(-1); }
   if(noisy) puts("\nMarginal reconstruction of ancestral sequences\n");

   fprintf (fout,"\n(1) Marginal reconstruction of ancestral sequences ");
   fprintf (fout,"(eq. 4 in Yang et al. 1995).\n    This algorithm is OK.\n");
   pChar1node=(double*)malloc(com.npatt*n*sizeof(double));
   pnode=(double*)malloc((nid*com.npatt+1)*(sizeof(double)+sizeof(char)));
   if (pnode==NULL||pChar1node==NULL) error2 ("oom pnode");
   za=(char*)(pnode+nid*com.npatt);

#ifdef BASEML
   if(com.seqtype==0 && com.ls%3==0 && com.coding) { coding=1; k1=com.ls/3; }
#endif
   if(com.seqtype==1) { coding=1; k1=com.npatt; }
   if(coding==1) {
      if((pbestAA=(double*)malloc(nid*k1*2*sizeof(double)))==NULL) 
         error2("oom pbestAA");
      bestAA=(char*)(pbestAA+nid*k1);
   }

   if(SetParameters(x)) puts("par err.");
   if(com.verbose>1) 
      fprintf(fout,"\nProb distribs at nodes, those with p<%.4f not listed\n",
         smallp);

   /* This loop reroots the tree at inode & reconstructs sequence at inode */
   for (inode=com.ns; inode<tree.nnode; inode++) {
      zero (pChar1node,com.npatt*n);
      if(noisy) printf ("Node %3d: ", inode+1);
      ReRootTree (inode);
      if (com.alpha==0) {   /* constant rate for sites (ngene>1 OK) */
         for (ig=0,lnL=0; ig<com.ngene; ig++) {
            if(com.Mgene>1) SetPGene(ig, 1, 1, 0, x+com.ntime);
            PartialLikelihood (tree.root, ig);
            for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
               for (i=0,fh=0,pbest=0,best=-1; i<n; i++) {
                  y=com.pi[i]*nodes[tree.root].lkl[h*n+i];
                  fh +=  y;
                  if ((y=y/fhs[h])>pbest) { pbest=y; best=i; }
                  pChar1node[h*n+i]=y;
               }
               if (fabs(fh-fhs[h])>min2(1e-6,fh*.001)) error2("err: fh 1?");
               lnL -= log(fh)*(double)com.fpatt[h];
               za[(inode-com.ns)*com.npatt+h]=(char)best;
               pnode[(inode-com.ns)*com.npatt+h]=pbest;
            }
         }
      }
      else {                /* gamma rates among sites */
         /* sum over the rate distribution (ir) to calculate pChar1node[] */
         zero(com.fhK,com.npatt); /* calculates lnL for error2 checking */
         for (ig=0; ig<com.ngene; ig++) {
            if(com.Mgene>1)
               SetPGene(ig,com.Mgene>1,com.Mgene>1,com.nalpha>1,x+com.ntime);
            for (ir=0; ir<com.ncatG; ir++) {
               SetPSiteClass(ir, x+com.ntime);
               PartialLikelihood (tree.root, ig);
               for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
                  for (i=0,fh=0; i<n; i++) {
                     fh += (y=com.pi[i]*nodes[tree.root].lkl[h*n+i]);
                     pChar1node[h*n+i] += com.freqK[ir]*y;
                  }
                  com.fhK[h] += com.freqK[ir]*fh;
               }
            }
         }
         for (h=0,lnL=0; h<com.npatt; h++) {
            FOR (i,n) pChar1node[h*n+i]/=fhs[h];
            if (fabs(1-sum(pChar1node+h*n,n))>1e-5) error2("Err: sum!=1");

            for (i=0,best=-1,pbest=-1; i<n; i++)
               if (pChar1node[h*n+i]>pbest) {best=i; pbest=pChar1node[h*n+i];}
            za[(inode-com.ns)*com.npatt+h]=(char)best;
            pnode[(inode-com.ns)*com.npatt+h]=pbest;
            if (fabs(com.fhK[h]-fhs[h])>1e-6) error2("err: fh 2?");
            lnL-=log(com.fhK[h])*com.fpatt[h];
         }
      }     /* if (alpha) */
      if(noisy) printf ("lnL = %12.6f\n", -lnL);

      /* print Prob distribution at inode if com.verbose>1 */
      if (com.verbose>1) {
         fprintf(fout,"\nProb distribution at node %d, by site\n",inode+1);
         if (com.ngene>1) fputs("\nSite (g) Freq  Data: \n\n",fout);
         else             fputs("\nSite   Freq   Data: \n\n",fout);
         for (h=0;h<com.ls;h++,FPN(fout)) {
            hp=com.pose[h];
            fprintf (fout,"%4d%7.0f   ", h+1,com.fpatt[hp]);
            print1site(fout,hp);   fputs(": ",fout);
            FOR(j,n)
               if(pChar1node[hp*n+j]>smallp) {
                  if (com.seqtype!=CODONseq) { str[0]=pch[j]; str[1]=0; }
                  else {
#ifdef CODEML
if(com.seqtype==1 &&(FROM61[j]<0||FROM61[j]>63)) { puts("strange 1"); getchar();}
                     getcodon(str,FROM61[j]);
#endif
                  }
                  fprintf(fout,"%s(%5.3f) ",str,pChar1node[hp*n+j]);
               }
         }
      }     /* if (verbose) */


      /* find the best amino acid from codon & nuc seqs */
#ifdef CODEML
      if(com.seqtype==CODONseq)
         FOR(h,com.npatt) {
            FOR(j,20) pAA[j]=0; 
            FOR(j,n) {
               i=GenetCode[com.icode][FROM61[j]];
               pAA[i]+=pChar1node[h*n+j];
            }
            /* matout(F0,pAA,1,20); */
            for(j=0,best=0,pbest=0; j<20; j++) 
               if(pAA[j]>pbest) { pbest=pAA[j]; best=j; }
            bestAA[(inode-com.ns)*com.npatt+h]=(char)best;
            pbestAA[(inode-com.ns)*com.npatt+h]=pbest;
         }
#endif

      if(com.seqtype==0 && coding) { /* coding seqs analyzed by baseml */
         FOR(h,com.ls/3) {  /* h-th codon */
            /* sums up probs for the 20 AAs for each node. Stop codons are 
               ignored, and so those probs are approxiamte. */
            for(j=0,y=0;j<20;j++) pAA[j]=0;
            FOR(k1,4) FOR(k2,4) FOR(k3,4) {
               ic=k1*16+k2*4+k3;
               b[0]=com.pose[h*3+0]*n+k1; 
               b[1]=com.pose[h*3+1]*n+k2; 
               b[2]=com.pose[h*3+2]*n+k3;
               fh=pChar1node[b[0]]*pChar1node[b[1]]*pChar1node[b[2]];
               if((ic=GenetCode[com.icode][ic])==-1) y+=fh;
               else                                  pAA[ic]+=fh;
            }
            if(fabs(1-y-sum(pAA,20))>1e-6) error2("AncestralMarginal strange?");

            for(j=0,best=0,pbest=0; j<20; j++) 
               if(pAA[j]>pbest) { pbest=pAA[j]; best=j; }

            bestAA[(inode-com.ns)*com.ls/3+h]=(char)best;
            pbestAA[(inode-com.ns)*com.ls/3+h]=pbest/(1-y);
         }
      }
   }        /* for (inode), This closes the big loop */
   ReRootTree(oldroot);

   if(com.seqtype==0 && coding) { /* coding seqs analyzed by baseml */
      fputs("\nBest amino acids reconstructed from nucleotide model.\n",fout);
      fputs("Prob at each node listed by amino acid (codon) site\n",fout);
      fputs("(Please ignore if not relevant)\n\n",fout);
      for(h=0;h<com.ls/3;h++,FPN(fout)) {
         fprintf(fout,"%4d ", h+1);
         FOR(j,com.ns) {
            if(com.cleandata)
               FOR(i,3) codon[i]=BASEs[com.z[j][com.pose[h*3+i]]];
            else 
               FOR(i,3) codon[i]=com.z[j][com.pose[h*3+i]];
            Codon2AA(codon, aa, com.icode, &i);
            fprintf(fout," %s(%c)",codon,AAs[i]);
         }
         fprintf(fout,": ");
         for (j=0; j<tree.nnode-com.ns; j++) {
            fprintf(fout," %1c (%5.3f)",
               AAs[bestAA[j*com.ls/3+h]],pbestAA[j*com.ls/3+h]);
         }
      }
   }

   /* calculate accuracy measures */
   zero(pMAPnode,nid);  fillxc(pMAPnodeA, 1., nid);
   for (inode=0; inode<tree.nnode-com.ns; inode++) {
      FOR (h,com.npatt) {
         pMAPnode[inode]+=com.fpatt[h]*pnode[inode*com.npatt+h]/com.ls;
         pMAPnodeA[inode]*=pow(pnode[inode*com.npatt+h], com.fpatt[h]);
      }
   }

   fputs("\n\nProb of best character at each node, listed by site ",fout);
   /* if (com.verbose==0) fputs("(choose verbose to list constant sites.)",fout); */
   if (com.ngene>1) fputs("\n\nSite (g) Freq  Data: \n",fout);
   else             fputs("\n\nSite   Freq   Data: \n",fout);

   FOR (h, com.ls) {
      hp=com.pose[h];
      fprintf(fout,"\n%4d ",h+1);

#ifdef CODEML
      /* fprintf (frst1,"\n%4d ", h+1); */
#endif

      if (com.ngene>1) {  /* which gene the site is from */
         for(ig=1;ig<com.ngene;ig++) if(hp<com.posG[ig]) break;
         fprintf(fout,"(%d)",ig);
      }
      fprintf(fout," %5.0f   ", com.fpatt[hp]);
      print1site(fout,hp);   fputs(": ",fout);
      FOR(j,nid) {
         if (com.seqtype!=CODONseq)
            fprintf(fout,"%c(%5.3f) ",pch[za[j*com.npatt+hp]],pnode[j*com.npatt+hp]);
         else {

#ifdef CODEML
            i=GenetCode[com.icode][ic=FROM61[za[j*com.npatt+hp]]];
            fprintf(fout," %s %1c %5.3f (%1c %5.3f)",
               getcodon(str,ic),AAs[i],pnode[j*com.npatt+hp], 
               AAs[bestAA[j*com.npatt+hp]],pbestAA[j*com.npatt+hp]);

/*
fprintf(frst1," %1c %5.3f",
   AAs[bestAA[j*com.npatt+hp]],pbestAA[j*com.npatt+hp]);
*/
#endif
         }
      }
   }


   /* Map changes onto branches 
      k1 & k2 are the two characters; p1 and p2 are the two probs. */
   fputs("\n\nSummary of changes along branches.\n",fout);
   fputs("Check root for directions of change.\n",fout);
   fputs("Ignore branches/sites with alignment gaps?\n",fout);
   for(j=0;j<tree.nbranch;j++,FPN(fout)) {
      inode=tree.branches[j][1];  nchange=0;
      fprintf(fout,"\nBranch %d:%5d..%-2d",j+1,tree.branches[j][0]+1,inode+1);
      if(inode<com.ns) fprintf(fout," (%s)\n\n",com.spname[inode]);
      else fputs("\n\n",fout);
      for(h=0;h<com.ls;h++) {
         hp=com.pose[h];
         if (com.seqtype!=CODONseq) {
            if(inode<com.ns)
               k2=(com.cleandata?pch[com.z[inode][hp]]:com.z[inode][hp]);
            else {
               k2=pch[za[(inode-com.ns)*com.npatt+hp]]; 
               p2=pnode[(inode-com.ns)*com.npatt+hp];
            }
            k1=pch[za[(tree.branches[j][0]-com.ns)*com.npatt+hp]];
            p1=pnode[(tree.branches[j][0]-com.ns)*com.npatt+hp];
         }
         else {
#ifdef CODEML
            if(inode<com.ns) {
               if(!com.cleandata)
                  { Codon2AA(com.z[inode]+hp*3,aa,com.icode,&k2); k2=AAs[k2];}
               else
                  k2=AAs[GenetCode[com.icode][FROM61[com.z[inode][hp]]]];
            }
            else {
               k2=AAs[bestAA[(inode-com.ns)*com.npatt+hp]];
               p2=pbestAA[(inode-com.ns)*com.npatt+hp];
            }
            k1=AAs[bestAA[(tree.branches[j][0]-com.ns)*com.npatt+hp]];
            p1=pbestAA[(tree.branches[j][0]-com.ns)*com.npatt+hp];
#endif
         }
         if(k1==k2) continue;
         fprintf(fout,"\t%4d ",h+1);

#ifdef SITELABELS
         if(sitelabels) fprintf(fout," %5s   ",sitelabels[h]);
#endif
         if(inode<com.ns) fprintf(fout,"%c %.3f -> %1c\n",k1,p1,k2);
         else             fprintf(fout,"%c %.3f -> %1c %.3f\n",k1,p1,k2,p2);
         ++nchange;
      }
   }

   ListAncestSeq(fout, za);
   fprintf(fout,"\n\nOverall accuracy of the %d ancestral sequences:", nid);
   matout2(fout,pMAPnode, 1, nid, 9,5);  fputs("for a site.\n",fout);
   matout2(fout,pMAPnodeA, 1, nid, 9,5); fputs("for the sequence.\n", fout);

   /* best amino acid sequences from codonml */
#ifdef CODEML
   if(com.seqtype==1) {
      fputs("\n\nAmino acid sequences inferred by codonml.\n",fout);
      if(!com.cleandata) 
         fputs("Results unreliable for sites with alignment gaps.\n",fout);
      for(inode=0; inode<nid; inode++) {
         fprintf(fout,"\nNode #%-10d  ",com.ns+inode+1);
         for(h=0;h<com.ls;h++) {
            fprintf(fout,"%c",AAs[bestAA[inode*com.npatt+com.pose[h]]]);
            if((h+1)%10==0) fputc(' ',fout);
         }
      }
      FPN(fout);
   }
#endif
   ChangesSites(fout,za);

   free(pnode); free(pChar1node);
   if(coding) free(pbestAA);

   return (0);
}


int ChangesSites(FILE*fout, char *zanc)
{
/* this lists and counts changes at sites from reconstructed ancestral sequences
   com.z[] has the data, and zanc[] has the ancestors
   For codon sequences (codonml or baseml with com.coding), synonymous and 
   nonsynonymous changes are counted separately.
   Added in Nov 2000.
*/
   char *pch=(com.seqtype==0?BASEs:(com.seqtype==2?AAs:BINs));
   char codon[2][4]={"   ","   "};
   int  i,h,hp,inode,k1,k2,d, coding=0, ls1=com.ls;
   double St,Nt,ns,na, nst,nat;

   if(!com.cleandata) {
      puts("\a\nthink about counting changes for unclean data\n");
      return(0);
   }
#if BASEML
   if(com.coding) { coding=1; ls1=com.ls/3; }
#endif
   if(com.seqtype==CODONseq) coding=1;

   if(coding) {
      fputs("\n\nChanges at sites (syn nonsyn).\n\n",fout);
      fputs("syn and nonsyn changes at each site\n",frst1);

      for(h=0; h<ls1; h++) {
         hp=com.pose[h];
         fprintf(fout,"%4d ",h+1);
         for(inode=0,nst=nat=0; inode<tree.nnode; inode++) {

FOR(i,3) codon[0][i]=codon[1][i]=-1;  /*  to force crashes */

            if(inode==tree.root) continue;
            if(com.seqtype==CODONseq) {
#ifdef CODEML
               k1=zanc[(nodes[inode].father-com.ns)*com.npatt+hp];
               getcodon(codon[0],FROM61[k1]);
               if(inode>=com.ns)
                  getcodon(codon[1],FROM61[zanc[(inode-com.ns)*com.npatt+hp]]);
               else {
                  if(com.cleandata)
                     getcodon(codon[1],FROM61[com.z[inode][hp]]);
                  else
                     FOR(i,3) codon[1][i]=com.z[inode][hp*3+i];
               }
#endif
            }
            else { /* baseml coding reconstruction */
               FOR(i,3) {
                  k1=zanc[(nodes[inode].father-com.ns)*com.npatt+com.pose[h*3+i]];
                  codon[0][i]=BASEs[k1];
               }
               if(inode>=com.ns)
                  FOR(i,3)
                     codon[1][i]=BASEs[zanc[(inode-com.ns)*com.npatt+com.pose[h*3+i]]];
               else {
                  if(com.cleandata)
                     FOR(i,3) codon[1][i]=BASEs[com.z[inode][com.pose[h*3+i]]];
                  else  
                     FOR(i,3) codon[1][i]=com.z[inode][com.pose[h*3+i]];
               }
            }

            FOR(i,3)
               if(codon[0][i]!=codon[1][i]) break;
            if(i<3) {
               difcodonNG(codon[0],codon[1],&St,&Nt,&ns,&na,0,com.icode);
               nst+=ns;  nat+=na;
               fprintf(fout," %3s.%3s ",codon[0],codon[1]);
            }
         }
         fprintf(fout," (%.1f  %.1f)\n", nst,nat);
         fprintf(frst1," %4d  %.1f  %.1f\n", h+1,nst,nat);
      }
   }
   else {  /* noncoding nucleotide or aa sequences */
      fputs("\n\nCounts of changes at sites\n\n",fout);
      for(h=0;h<com.ls;h++) {
         hp=com.pose[h];
         fprintf(fout,"%4d ",h+1);
         for(inode=0,d=0;inode<tree.nnode;inode++) {
            if(inode==tree.root) continue;
            k1=pch[ zanc[(nodes[inode].father-com.ns)*com.npatt+hp] ];
            if(inode<com.ns) {
               if(com.cleandata) k2=pch[com.z[inode][hp]];
               else              k2=com.z[inode][hp];
            }
            else  
               k2=pch[ zanc[(inode-com.ns)*com.npatt+hp] ];
            if(k1!=k2) {
               d++;
               fprintf(fout," %c%c", k1,k2);
            }
         }
         fprintf(fout," (%d)\n", d);
      }
   }
   return(0);
}



#if BASEML
extern complex cU[], cV[], cRoot[];
#endif

int GetPMatBranch (double Pt[], double t, int inode)
{
/* P(t) for branch leading to inode, used in baseml and codeml

   For branch&site models, this uses IClass to set U&V&Root to the right place 
   in the 5 sets of U&V&Root vectors (w0b,w0f,w1b,w1f,w2f):
       iclass=0: w0b (0), w0f (1)
       iclass=1: w1b (2), w1f (3)
       iclass=2: w0b (0), w2f (4)
       iclass=3: w1b (2), w2f (4)
*/
   int updateUVRoot=0;
   int iUVR=0, iUVR0[4][2]={{0,1}, {2,3}, {0,4}, {2,4}};


#if CODEML
   if(com.seqtype==1 && com.NSsites && com.model) { /* branch&site models */
      iUVR=iUVR0[IClass][nodes[inode].label];
      U=_UU[iUVR]; V=_VV[iUVR]; Root=_Root[iUVR]; 
   }
   else if (com.seqtype==CODONseq && com.model==2 && N_OMEGA<=5) { /* branch model */
      iUVR=nodes[inode].label;  U=_UU[iUVR]; V=_VV[iUVR]; Root=_Root[iUVR];
   }
   else if (com.seqtype==CODONseq && com.model) { /* branch model */
      updateUVRoot=(com.omega!=nodes[inode].omega);
      if(com.NSsites && Qfactor_NS!=Qfactor_NS_branch[nodes[inode].label])
         { Qfactor_NS=Qfactor_NS_branch[nodes[inode].label]; updateUVRoot=1; }
      if(updateUVRoot)
         EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, com.kappa, nodes[inode].omega, Pt);
   }
   if (com.seqtype==AAseq && com.model==Poisson)
      PMatJC69like(Pt,t,com.ncode);
   else
      PMatUVRoot(Pt,t,com.ncode,U,V,Root);

#elif BASEML

   if (com.nhomo==2)
      EigenTN93(com.model,nodes[inode].kappa,1,com.pi,&nR,Root,Cijk);
   else if (com.nhomo>2)
      EigenTN93(com.model,nodes[inode].kappa,1,nodes[inode].pi,&nR,Root,Cijk);
   if (com.model<=K80)
      PMatK80(Pt, t, (com.nhomo==2?nodes[inode].kappa:com.kappa));
   else if (com.model<=REV)  PMatCijk(Pt, t);
   else                      cPMat (Pt, t, com.ncode, cU, cV, cRoot);

#endif

   return(0);
}


double *PMatTips;
char *nodeChar, *ancSeq;  /* C_z(i) in PPSG2000; lkl stores L_z(i) */
/* nodeChar[(z-com.ns)*com.npatt*n+h*n+i] is C_z(i), probability 
   of the best reconstruction upto node z, given that z's father
   has character i.
   nodeChar[nintern*n*com.npatt]
*/

void UpPassPPSG2000 (int inode, int igene)
{
/* Steps 2 & 4 in PPSG2000, modified from PartialLikelihood().  
   It calculates the L's in PPSG for node inode (z), storing them in 
   nodes[inode].lkl.  
   nodes[inode].lkl[h*n+i] and nodeChar[(inode-com.ns)*com.npatt*n+h*n+i] are 
   prob and chara at inode, given that father has ichar.
   inode has daughter nodes x (that is, x and y in PPSG).
   ichar is for father, and jchar is for inode.
   nodes[root].lkl[h] is the prob for best reconstruction at the site

   Need more thoughts about ambiguity data (error2 in running stewart.aa).
*/
   int n=com.ncode, ichar,jchar, jbest=-1, x,ix,h;
   double t,pj,pjbest=-1;

   FOR(ix,nodes[inode].nson) 
      if(nodes[x=nodes[inode].sons[ix]].nson>0) /* if son x is not a tip */
         UpPassPPSG2000(x,igene);

   if(inode!=tree.root)
      GetPMatBranch(PMat,nodes[inode].branch*com.rgene[igene]*_rateSite,inode);
   for(h=com.posG[igene]; h<com.posG[igene+1]; h++) { /* to save space */
      for(ichar=0; ichar<(inode!=tree.root?n:1); ichar++) { /* ichar for father */
         for(jchar=0; jchar<n; jchar++) { /* jchar for inode, for finding best*/
            pj=(inode!=tree.root?PMat[ichar*n+jchar]:com.pi[jchar]);
            for(ix=0; ix<nodes[inode].nson; ix++) {
               x=nodes[inode].sons[ix];  /* this loops through inode's sons */
               if(com.cleandata && nodes[x].nson<1)  /* x is tip */
                  pj*=PMatTips[x*n*n+jchar*n+com.z[x][h]];
               else                            /* x is internal node or tip */
                  pj*=nodes[x].lkl[h*n+jchar];
            }
            if(jchar==0 || pjbest<pj)
               { pjbest=pj; jbest=jchar; } 
         }
         nodes[inode].lkl[h*n+ichar]=pjbest; 
         nodeChar[(inode-com.ns)*com.npatt*n+h*n+ichar]=(char)jbest;
      }

      if(_nnodeScale && _nodeScale[inode]) {  /* scaling to avoid underflow */
         for(ichar=0,t=0;ichar<n;ichar++)
            if(nodes[inode].lkl[h*n+ichar]>t) t=nodes[inode].lkl[h*n+ichar];
         for(ichar=0;ichar<n;ichar++)  nodes[inode].lkl[h*n+ichar]/=t;
         _nodeScaleF[h]+=log(t);
      }
   }
}

void DownPassPPSG2000(int inode, int igene)
{
/* this reads out the best chara for inode from nodeChar[] into ancSeq[].
   c0 is best chara for father (if inode is not root).
*/
   int i,ison, h;
   char c0=0;

   for(h=com.posG[igene]; h<com.posG[igene+1]; h++) {
      if(inode!=tree.root) 
         c0=ancSeq[(nodes[inode].father-com.ns)*com.npatt+h];
      ancSeq[(inode-com.ns)*com.npatt+h]
         =nodeChar[(inode-com.ns)*com.npatt*com.ncode+h*com.ncode+c0];
   }
   FOR(i,nodes[inode].nson)
      if(nodes[ison=nodes[inode].sons[i]].nson>1) 
         DownPassPPSG2000(ison,igene);
}


int AncestralJointPPSG2000 (FILE *fout, double x[], double fhs[])
{
/* Z. Yang, 8 June 2000.  
   Joint ancestral reconstruction, taking character states for all nodes at a 
   site as one entity, based on the algorithm of Pupko et al. (2000) 
   Molecular Biology and Evolution 17:890-896.

   fhs[]: fh[] for each site pattern
   nodes[].lkl[] are used to store the probabilities (L's in PSSG), so that 
   partial likelihood values are destroyed.
   ancSeq[nintern*com.npatt] stores the ancestral seqs, from UpPassPPSG2000.
*/
   char *pch=(com.seqtype==0?BASEs:(com.seqtype==2?AAs:BINs)),codon[4]="";
   int n=com.ncode,nintern=tree.nnode-com.ns, i,h,hp,ig;
   double *psites;

   nodeChar=(char*)malloc(nintern*com.npatt*(n+1)*sizeof(char));
   ancSeq=nodeChar+nintern*com.npatt*n;
   if(nodeChar==NULL) error2("oom nodeChar in AncestralJointPPSG2000");
   if(com.cleandata && com.sspace<com.ns*n*n*(int)sizeof(double)) {
      com.sspace=com.ns*n*n*sizeof(double);
      if((com.space=(double*)realloc(com.space,com.sspace))==NULL)
         error2("oom PPSG");
      if(noisy) printf("\nspace reset to %d bytes in Ancestral.\n",com.sspace);
   }
   if(_nnodeScale)  zero(_nodeScaleF,com.npatt);

   PMatTips=psites=com.space;
   for(ig=0; ig<com.ngene; ig++) {
      if(com.Mgene>1) SetPGene(ig,1,1,0,x+com.ntime);
      FOR(i,com.ns)
         GetPMatBranch(PMatTips+i*n*n,nodes[i].branch*com.rgene[ig]*_rateSite,i);
      UpPassPPSG2000(tree.root, ig);
      for (h=com.posG[ig]; h<com.posG[ig+1]; h++)
         psites[h]=nodes[tree.root].lkl[h*n+0]/fhs[h];
      DownPassPPSG2000(tree.root,ig); /* enumerate best characters */
   }
   if(_nnodeScale)  FOR(h,com.npatt) fhs[h]-=_nodeScaleF[h];

   fputs("\n\n(2) Joint reconstruction of ancestral sequences ", fout);
   fputs("(eqn. 2 in Yang et al. 1995),\nusing the algorithm of ",fout);
   fputs("Pupko et al. (2000 Mol Biol Evol 17:890-896).\n",fout);

   fputs("\nListed by site, reconstruction (prob.)\n",fout);
   if (com.ngene>1) fputs("\n\nSite (g) Freq  Data: \n",fout);
   else             fputs("\n\nSite   Freq   Data: \n",fout);

   /* Output ancestral sequences */
   for (h=0; h<com.ls; h++) {
      hp=com.pose[h];
      for(ig=1;ig<com.ngene;ig++) if (hp<com.posG[ig]) break;
      ig--;

      fprintf (fout,"\n%4d ", h+1);
      if (com.ngene>1) fprintf(fout,"(%d) ", ig+1);
      fprintf (fout," %6.0f  ",com.fpatt[hp]);
      print1site(fout,hp); fputs(": ",fout);
      FOR(i,nintern) {
         if(com.seqtype==1) {
#ifdef CODEML
            fprintf(fout,"%s ",getcodon(codon,FROM61[ancSeq[i*com.npatt+hp]]));
#endif
         }
         else
            fprintf(fout,"%c",pch[ancSeq[i*com.npatt+hp]]);
      }
      fprintf(fout,"  (%.4f)",psites[hp]);
   }  /* for (h) */

   ListAncestSeq(fout, ancSeq);
   free(nodeChar);
   if(noisy)  FPN(F0);
   return (0);
}


int AncestralSeqs (FILE *fout, double x[])
{
/* Ancestral sequence reconstruction using likelihood (Yang et al. 1995).
   Marginal works with constant rate and variable rates among sites.
   Joint works only with constant rate among sites.
*/
   double lnL=0, *fhs;

   if(com.Mgene==1) {
      puts("When Mgene=1, use RateAncestor = 0.");
      return(0);
   }

   if(_nnodeScale) 
      { puts("\nAncestralSeq with scaling, not done yet."); return(-1); }
   if (tree.nnode==com.ns) 
      { puts("\nNo ancestral nodes to reconstruct..\n");  return(0); }
   if (noisy) printf ("\nReconstructed ancestral states go into file rst.\n");
   fprintf(fout, "\nAncestral reconstruction by %sML.\n",
          (com.seqtype==0?"BASE":(com.seqtype==1?"CODON":"AA")));
   FPN(fout);  OutaTreeN(fout,1,1);  FPN(fout);  FPN(fout);
   OutaTreeN(fout,0,0);  FPN(fout);  FPN(fout);
   OutaTreeB (fout);     FPN(fout);

   fputs("\ntree with node labels for Rod Page's TreeView\n",fout);
   IncludeNodeLabel=1;
   OutaTreeN(fout,1,1);  FPN(fout);
   IncludeNodeLabel=0;

   fprintf (fout, "\nNodes %d to %d are ancestral\n", com.ns+1,tree.nnode);
   if ((fhs=(double*)malloc(com.npatt*sizeof(double)))==NULL) error2("oom fhs");
   
   ProbSitePattern (x, &lnL, fhs);
   if (com.alpha) {
      puts("\nRates are variable among sites.");
      puts("Only marginal reconstructions are implemented.");
   }
   if(!com.verbose) fputs("Constant sites not listed for verbose=0\n",fout);
   if(!com.cleandata) fputs("Unreliable at sites with alignment gaps\n",fout);

   AncestralMarginal(fout,x,fhs);  fflush(fout);
   if(com.alpha==0 && tree.nnode>com.ns+1 && com.cleandata)
      AncestralJointPPSG2000 (fout, x, fhs);
   FPN(fout);
   free (fhs);
   return (0);
}

#endif

int SetnScale(int inode);

void InitializeNodeScale(void)
{
/* This allocates memory to hold scale factors for nodes and also decide on the 
   internal nodes for scaling by calling SetnScale().  The memory required is 
   _nnodeScale*com.npatt*sizeof(double).
*/
   int k;

   _nnodeScale=0; FOR(k,tree.nnode) _nodeScale[k]=0;
   SetnScale(tree.root);
   k=_nnodeScale*com.npatt*sizeof(double);
   if(lklSiteClass) k*=com.ncatG;
   if(_nnodeScale) {
      if((_nodeScaleF=(double*)realloc(_nodeScaleF,k))==NULL) 
         error2("oom nscale");
      memset(_nodeScaleF,k,0);

      if(noisy) {
         printf("\n%d nodes used for scaling: ",_nnodeScale);
         FOR(k,tree.nnode) if(_nodeScale[k]) printf(" %2d",k+1);
         /* sleep(2000); FPN(F0); */
      }
   }

}

int SetnScale(int inode)
{
/* This marks nodes for applying scaling factors when calculating f[h].
*/
   int i,d=0, every;

   if(com.seqtype==0)       every=201;  /* baseml */
   else if(com.seqtype==1)  every=21;   /* codonml */
   else                     every=101;  /* aaml */

   FOR(i,nodes[inode].nson)
      d+=(nodes[nodes[inode].sons[i]].nson?SetnScale(nodes[inode].sons[i]):1);
   if(inode!=tree.root && d>every)
      { _nodeScale[inode]=1;  d=1; _nnodeScale++; }
   return(d);
}

int fx_r(double x[], int np);

int lfunRates (FILE* fout, double x[], int np)
{
/* for dG, AdG or similar non-parametric models
   This distroys com.fhK[], and in return,
   fhK[<npatt] stores rates for conditional mean (re), and 
   fhK[<2*npatt] stores the most probable rate category number.
   fhs[npatt] stores fh=log(fh).
*/
   int  ir,il,it=0, h,hp,j, nscale=1, direction=-1;
   double lnL=0,fh,fh1, t, re,mre,vre, b1[NCATG],b2[NCATG],*fhs;

   if (noisy) printf("\nEstimated rates fo3r sites go into file %s\n",ratef);
   if (SetParameters(x)) puts ("par err. lfunAdG_rate");

   fprintf(fout, "\nEstimated rates for sites from %sML.\n",
          (com.seqtype==0?"BASE":(com.seqtype==1?"CODON":"AA")));
   OutaTreeN(fout,1,1); FPN(fout);
   fprintf (fout,"\nFrequencies and rates for categories (K=%d)", com.ncatG);
   matout (fout, com.freqK, 1, com.ncatG);
   matout (fout, com.rK, 1, com.ncatG);
   if (com.rho) {
      fprintf(fout,"\nTransition prob matrix over sites");
      matout2(fout,com.MK,com.ncatG,com.ncatG,8,4);
   }

   if((fhs=(double*)malloc(com.npatt*sizeof(double)))==NULL) error2("oom fhs");
   fx_r(x,np);
   if(_nnodeScale)
      FOR(h,com.npatt) {
         t=com.fhK[0*com.npatt+h];  lnL-=t*com.fpatt[h];
         for(ir=1,fh=com.freqK[ir]*(com.fhK[h]=1); ir<com.ncatG; ir++)
            fh+=com.freqK[ir]*
               (com.fhK[ir*com.npatt+h]=exp(com.fhK[ir*com.npatt+h]-t));
         fhs[h]=fh*exp(t);
      }
   if (com.rho==0) {     /* dG model */
      fputs("\nSite Freq  Data    ln(f)    Nexp.   Rates\n\n",fout);
      for (h=0,mre=vre=0; h<com.npatt; h++) {
         for (ir=0,fh=0,it=0,t=re=0; ir<com.ncatG; ir++) {
            if((fh1=com.freqK[ir]*com.fhK[ir*com.npatt+h])>t) 
               { t=fh1; it=ir; }
            fh+=fh1;  re+=fh1*com.rK[ir];
         }

         if(!_nnodeScale) fhs[h]=fh;
         lnL -= com.fpatt[h]*log(fh);

         re/=fh;
         mre+=com.fpatt[h]*re/com.ls;     vre+=com.fpatt[h]*re*re/com.ls;
         com.fhK[h]=re;                   com.fhK[com.npatt+h]=it+1.;
      }
      vre-=mre*mre;
      FOR(h,com.ls) {
         hp=com.pose[h];
         fprintf(fout,"%4d %5.0f  ",h+1,com.fpatt[hp]);  print1site(fout,hp);
         fprintf(fout," %8.3f%8.2f%8.3f%6.0f\n",
           log(fhs[hp]),com.ls*fhs[hp],com.fhK[hp],com.fhK[com.npatt+hp]);
      }
      if(!com.cleandata || com.ngene>1) 
         fputs("\n(With ambiguity or multigene data, ignore Nexp.)\n",fout);
   }
   else {      /* Auto-dGamma model */
      fputs("\nSite Freq  Data  Rates\n\n",fout);
      h = (direction==1?com.ls-1:0);
      for (il=0,mre=vre=0; il<com.ls; h-=direction,il++) {
         hp=com.pose[h];
         if (il==0)
            FOR(ir,com.ncatG) b1[ir]=com.fhK[ir*com.npatt+hp];
         else {
            for (ir=0; ir<com.ncatG; ir++) {
               for (j=0,fh=0; j<com.ncatG; j++)
                  fh+=com.MK[ir*com.ncatG+j]*b1[j];
               b2[ir] = fh*com.fhK[ir*com.npatt+hp];
            }
            xtoy (b2, b1, com.ncatG);
         }
         if ((il+1)%nscale==0)
            { fh=sum(b1,com.ncatG); abyx(1/fh,b1,com.ncatG); lnL-=log(fh); }

         for (ir=0,it=-1,re=fh1=t=0; ir<com.ncatG; ir++) {
            re+=com.freqK[ir]*b1[ir]*com.rK[ir];
            fh1+=com.freqK[ir]*b1[ir];
            if (b1[ir]>t) {it=ir; t=b1[ir]; }
         }
         re/=fh1;  mre+=re/com.ls;   vre+=re*re/com.ls;

         fprintf(fout,"%4d %5.0f  ",h+1,com.fpatt[hp]);  print1site(fout,hp);
         fprintf(fout," %8.3f%6.0f\n", re,it+1.);
      }  /* for(il) */
      vre-=mre*mre;
      for (ir=0,fh=0; ir<com.ncatG; ir++)  fh += com.freqK[ir]*b1[ir];
      lnL -= log(fh);
   }
   if (noisy) printf ("lnL=%14.6f\n", -lnL);
   fprintf (fout,"\nlnL=%14.6f\n", -lnL);
   if(com.ngene==1) {
      fprintf (fout,"\nmean(r^)=%9.4f  var(r^)=%9.4f", mre, vre);
      fprintf (fout,"\nAccuracy of rate prediction: corr(r^,r) =%9.4f\n", 
               sqrt(com.alpha*vre));
   }
   free(fhs);
   return (0);
}


double lfunAdG (double x[], int np)
{
/* Auto-Discrete-Gamma rates for sites
   See notes in lfundG().
*/
   int  nscale=1, h,il, ir, j, FPE=0;
   int  direction=-1;  /* 1: n->1;  -1: 1->n */
   double lnL=0, b1[NCATG], b2[NCATG], fh;

   NFunCall++;
   fx_r(x,np);
   if(_nnodeScale)
      FOR(h,com.npatt) {
         fh=com.fhK[0*com.npatt+h];
         lnL-=fh*com.fpatt[h];
         for(ir=1,com.fhK[h]=1; ir<com.ncatG; ir++) 
            com.fhK[ir*com.npatt+h]=exp(com.fhK[ir*com.npatt+h]-fh);
      }
   h = (direction==1?com.ls-1:0);
   for (il=0; il<com.ls; h-=direction,il++) {
      if (il==0)
         FOR(ir,com.ncatG) b1[ir]=com.fhK[ir*com.npatt+com.pose[h]];
      else {
         for (ir=0; ir<com.ncatG; ir++) {
            for (j=0,fh=0; j<com.ncatG; j++)
               fh+=com.MK[ir*com.ncatG+j]*b1[j];
            b2[ir]=fh*com.fhK[ir*com.npatt+com.pose[h]];
         }
         xtoy(b2,b1,com.ncatG);
      }
      if((il+1)%nscale==0) {
         fh=sum(b1,com.ncatG);
         if(fh<1e-40) {
            if(!FPE) {
               FPE=1; printf ("h,fh%6d %12.4e\n", h+1,fh);
               print1site(F0,h);  FPN(F0);
            }
            fh=1e-70;
         }
         abyx(1/fh,b1,com.ncatG); lnL-=log(fh);
      }
   }
   for (ir=0,fh=0; ir<com.ncatG; ir++)  fh+=com.freqK[ir]*b1[ir];
   lnL-=log(fh);
   return (lnL);
}


double lfundG (double x[], int np)
{
/* discrete rates for sites.  
   This deals with scaling for nodes to avoid underflow if(_nnodeScale).
   The routine calls fx_r() to calculate com.fhK[], which holds log{f(x|r)} 
   when scaling or f(x|r) when not.  Scaling factors are set and used for each 
   site class (ir) to calculate log(f(x|r).  When scaling is used, the routine 
   converts com.fhK[] into f(x|r), by collecting scaling factors into lnL.  
   The rest of the calculation then becomes the same and relies on f(x|r).  
   Check notes in fx_r.
   This is also used for NSsites models in codonml.  
   Note that scaling is used between fx_r() and PartialLikelihood()
*/
   int h,ir, it, FPE=0;
   double lnL=0, fh=0;

   NFunCall++;
   fx_r(x,np);
   if(_nnodeScale)  /* com.fhK[] has log{f(x|r}.  Note the scaling for nodes */
      FOR(h,com.npatt) {
         for(ir=1,it=0; ir<com.ncatG; ir++) /* select term for scaling */
            if(com.fhK[ir*com.npatt+h]>com.fhK[it*com.npatt+h]) it=ir;

         fh=com.fhK[it*com.npatt+h];
         lnL-=fh*com.fpatt[h];
         for(ir=0,com.fhK[it*com.npatt+h]=1; ir<com.ncatG; ir++) 
            if(ir!=it) com.fhK[ir*com.npatt+h]=exp(com.fhK[ir*com.npatt+h]-fh);
      }

   FOR(h,com.npatt) {
      if (com.fpatt[h]==0) continue;
      for(ir=0,fh=0; ir<com.ncatG;ir++) 
         fh+=com.freqK[ir]*com.fhK[ir*com.npatt+h];

      if(fh<=0) {
         if(!FPE) {
            FPE=1;  matout(F0,x,1,np);
            printf("\nlfundG: h=%4d  fhK=%9.6e\ndata: ", h, fh);
            print1site(F0,h);  FPN(F0);
         }
         fh=1e-70;
      }
      lnL-=log(fh)*com.fpatt[h];

      if(com.print<0){
         fprintf(flnf,"\n%6d %6.0f %12.6f %9.2f  ",
            h+1,com.fpatt[h],log(fh),com.ls*fh);
         print1site(flnf,h);
      }
   }
   return(lnL);
}


int SetPSiteClass(int iclass, double xcom[])
{
/* This sets parameters for the iclass-th site class
   This is used by PartialLikelihood() and also updatelkl in both algorithms
   For method=0 and 1.
*/
   int k=com.nrgene+!com.fix_kappa;

   _rateSite=com.rK[iclass];
#if CODEML
   if(com.seqtype==CODONseq && com.NSsites) {
      _rateSite=1;
      if (!com.model) {
         if(com.aaDist) {
            if(com.aaDist<10)       POMEGA=xcom+k+com.ncatG-1+2*iclass;
            else if(com.aaDist==11) POMEGA=xcom+k+com.ncatG-1+4*iclass;
            else if(com.aaDist==12) POMEGA=xcom+k+com.ncatG-1+5*iclass;
         }
         EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, com.kappa,com.rK[iclass],PMat);
      }
   }
#endif
   return (0);
}


int fx_r (double x[], int np)
{
/* This calculates f(x|r) if(_nnodeScale) or log{f(x|r)} if otherwise, 
   that is, the (log) probability of observing data x at a site, given the 
   rate r for the site.
   The results are stored in com.fhK[com.ncatG*com.npatt].
   This deals with underflows with large trees and uses global variables 
   _nodeScale and _nodeScaleF[_nnodeScale*com.npatt].
*/
   int  h, ir, i,k, ig, FPE=0;
   double fh;

   if(SetParameters(x)) puts("\npar err..");
   for(ig=0; ig<com.ngene; ig++) { /* alpha may differ over ig */
      if(com.Mgene>1)
         SetPGene(ig, com.Mgene>1, com.Mgene>1, com.nalpha>1, x+com.ntime);
      for(ir=0; ir<com.ncatG; ir++) {

#ifdef CODEML
         IClass=ir;
#endif
         if(ir && lklSiteClass) {  /* shift com.lkl */
            if(_nnodeScale) _nodeScaleF+=com.npatt*_nnodeScale;
            for(i=com.ns;i<tree.nnode;i++)
               nodes[i].lkl+=(tree.nnode-com.ns)*com.ncode*com.npatt;
         }
         SetPSiteClass(ir,x+com.ntime);
         PartialLikelihood(tree.root,ig);

         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            if (com.fpatt[h]==0) continue;
            for (i=0,fh=0; i<com.ncode; i++)
               fh += com.pi[i]*nodes[tree.root].lkl[h*com.ncode+i];

            if (fh<0) {
               if(fh<-1e-10 && !FPE) { /* note that 0 is o.k. here */
                  FPE=1; matout(F0,x,1,np);
                  printf ("\nfx_r: h=%3d  fhK=%9.6e\nData: ", h+1,fh);
                  if(com.seqtype==0||com.seqtype==2)
                     { print1site(F0,h);  FPN(F0); }
               }
               fh=1e-70;
            }
            if(!_nnodeScale)
               com.fhK[ir*com.npatt+h]=fh;
            else
               for(k=0,com.fhK[ir*com.npatt+h]=log(fh); k<_nnodeScale; k++)
                  com.fhK[ir*com.npatt+h]+=_nodeScaleF[k*com.npatt+h];
         }  /* for (h) */
      }   /* for (ir) */
      if(lklSiteClass) {
         if(_nnodeScale) _nodeScaleF-=(com.ncatG-1)*_nnodeScale*com.npatt;
         for(i=com.ns;i<tree.nnode;i++)
            nodes[i].lkl-=(com.ncatG-1)*(tree.nnode-com.ns)*com.ncode*com.npatt;
      }
   }  /* for(ig) */
   return(0);
}


double lfun (double x[], int np)
{
   int  h,i,k, ig, FPE=0;
   double lnL=0, fh;

   NFunCall++;
   if(SetParameters(x)) puts ("\npar err..");
   for(ig=0; ig<com.ngene; ig++) {
      if(com.Mgene>1) SetPGene(ig,1,1,0,x+com.ntime);
      PartialLikelihood (tree.root, ig);
      for(h=com.posG[ig]; h<com.posG[ig+1]; h++) {
         if(com.fpatt[h]==0) continue;
         for(i=0,fh=0; i<com.ncode; i++) 
            fh += com.pi[i]*nodes[tree.root].lkl[h*com.ncode+i];
         if(fh<=0) {

            if(!FPE) {
               FPE=1;  matout(F0,x,1,np);
               printf("\alfun: h=%4d  fh=%9.6e\nData: ", h+1,fh);
               print1site(F0,h);  FPN(F0);
            }

            fh=1e-70;
         }
         fh=log(fh);
         FOR(k,_nnodeScale) fh+=_nodeScaleF[k*com.npatt+h];

         lnL-=fh*com.fpatt[h];
         if (com.print<0) {
            fprintf(flnf,"\n%6d %6.0f %12.6f %9.2f ",
               h+1,com.fpatt[h],fh,com.ls*exp(fh));
            print1site(flnf,h);
         }
      }  /* for(h) */
   }     /* for(ig) */
   return (lnL);
}



int print1site (FILE*fout,int h)
{
/* This print out one site in the data com.z[].  It may be the h-th 
   site in the original data file or the h-th pattern.
   The data are coded if com.cleandata==1 or not if otherwise.
*/
   char str[4]="";
   char *pch=(com.seqtype==0?BASEs:(com.seqtype==2?AAs:BINs));
   int i, ib,nb=(com.seqtype==CODONseq?3:1), iaa=0,ic=0;

   if(!com.cleandata) 
      FOR(i,com.ns) {
         FOR(ib,nb) fprintf(fout,"%c", com.z[i][h*nb+ib]);
#ifdef CODEML
         if(nb==3)  {
            Codon2AA(com.z[i]+h*nb, str, com.icode, &iaa);
            fprintf(fout," (%c) ",AAs[iaa]);
	         /*fprintf(frst1,"%c", AAs[iaa]);*/
         }
#endif
      }
   else {
      FOR(i,com.ns) {
         if(nb==1) fprintf(fout,"%c", pch[com.z[i][h]]);
#ifdef CODEML
         else {
            ic=FROM61[com.z[i][h]];  iaa=GenetCode[com.icode][ic];
            fprintf(fout,"%s (%c) ", getcodon(str,ic), AAs[iaa]);

            /*fprintf(frst1,"%c", AAs[iaa]);*/
         }
#endif
      }
   }
   return(0);
}
   

int InitPartialLikelihood (void)
{
/* set partial likelihood at tips of the tree, considering missing data.
   Need testing if sequences in the data are ancestors.
*/
   static int times=0;
   int n=com.ncode, is,j,k,h;
   int i=0, nb[3]={-1,-1,-1},ib[3][4]={-1,-1,-1}, i0=0,i1=0,i2=0,ic=0,nsense=0;
   char *pch=(com.seqtype==0?BASEs:(com.seqtype==2?AAs:BINs));

   /*
   if(times++) { puts("InitPartialLikelihood: times>0"); sleep(3000); }
   */
   FOR(is,com.ns) nodes[is].lkl=com.lkl0 + is* n*com.npatt;
   for (is=0;is<com.ns;is++) {
      zero(nodes[is].lkl, com.npatt*n);
      for (h=0;h<com.npatt;h++) {
         if(com.seqtype==CODONseq) {
#ifdef CODEML
            for (h=0;h<com.npatt;h++) {
               FOR(i,3)  NucListall(com.z[is][h*3+i],&nb[i],ib[i]);
               for(i0=0,nsense=0; i0<nb[0]; i0++) FOR(i1,nb[1]) FOR(i2,nb[2]){
                  ic=ib[0][i0]*16+ib[1][i1]*4+ib[2][i2];         
                  ic=FROM64[ic];
                  if(ic>-1) { nsense++;  nodes[is].lkl[h*n+ic]=1; }
               }
            }
            if(nsense==0) {
               printf("\nseq. %d pattern %d",is+1,h+1);
               error2("Initlkl: stop codon");
            }
#endif
         }
         else {
            k=strchr(pch,com.z[is][h])-pch; 
            if(k<0) { printf("Character %d\n",com.z[is][h]); exit(-1); }
            if(k<n)                       nodes[is].lkl[h*n+k]=1;
            else if (com.seqtype==AAseq)  FOR(k,n) nodes[is].lkl[h*n+k]=1;
            else if (com.seqtype==BASEseq)
               FOR(j,nBASEs[k]) 
                  nodes[is].lkl[h*n+ (strchr(BASEs,EquateNUC[k][j])-BASEs)]=1;
         }
      }  /* for (h) */
   }     /* for (is) */
   return(0);
}


/* November, 1999, Minimization branch by branch ***/
int minbranches(double *lnL,double x[],double tb[2],double e,double space[]);
int lfunt(double t, int a,int b,double x[],double *l, double space[]);
int lfuntdd(double t, int a,int b,double x[], double *l,double*dl,double*ddl,
    double space[]);
int lfunt_SiteClass(double t, int a,int b,double x[],double *l,double space[]);
int lfuntdd_SiteClass(double t, int a,int b,double x[],
    double *l,double*dl,double*ddl,double space[]);

int minB (FILE*fout, double *lnL,double x[],double xb[][2],double space[])
{
/* optimizing one branch at one time
   Z. Yang, November 1999
   This calls minbranches to estimate branch lengths and ming2 to 
   estimate other paramters.
   xb[] are destroyed.
   Note that minbranches() does not update nodes[].lkl.
*/
   int i, npcom=com.np-com.ntime, ir,maxr=(npcom?500:1);
   double *xcom=x+com.ntime, lnL0= *lnL,dl, e=.1, tb[2];

   if(com.ntime==0) error2("minB: should not come here");

   /* tb[0]=xb[0][0]; tb[1]=xb[0][1]; */
   tb[0]=1e-8, tb[1]=50;
   if(npcom)
      for(i=com.ntime,xcom=x+com.ntime; i<com.np; i++)
         { xb[i-com.ntime][0]=xb[i][0]; xb[i-com.ntime][1]=xb[i][1]; }
   else
      e=.5e-3;
   com.ntime=0;  com.fix_branch=2;
   if(*lnL<=0) *lnL=com.plfun(x,com.np);
   FOR(ir,maxr) {
      if(npcom) {
         if(noisy>2) 
            printf("\n\nRound %d phase I: Paras (%d) (e=%.2e)",ir+1,npcom,e);
         ming2(NULL,lnL,com.plfun,NULL,xcom,xb, space,e,npcom);
         if(noisy>2) { FPN(F0); matout(F0,xcom,1,npcom); }
      }

      if(noisy>2)
         printf("\nRound %d phase II: Blengths (%d, e=%.3e)\n",ir+1,tree.nbranch,e);

#if(CODEML)
      if(com.NSsites==0) POMEGA=xcom+com.nrgene+!com.fix_kappa;
#endif
      minbranches(lnL,x,tb, e, space);

      if((dl= lnL0- *lnL)<.4e-6) break;
      if(dl<0.5) e=.05;

      lnL0= *lnL;
      if(fout) {
         fprintf(fout,"%4d %12.5f x ", ir+1,*lnL);
         FOR(i,com.np) fprintf(fout,"%9.5f",x[i]);
         FPN(fout);  fflush(fout);
      }
   }
   FOR(i,tree.nnode) 
      if(i!=tree.root) x[nodes[i].ibranch]=nodes[i].branch;
   com.ntime=tree.nbranch;  com.fix_branch=0;
   /* *lnL=-com.plfun(x,com.np); */
   if(noisy>2) printf("\nlnL  = %12.6f\n",- *lnL);

   return(0);
}

int updatelkl(double x[], int inode)
{
/* update lkl for inode.  
   Note that _nodeScaleF and nodes[].lkl are shifted if(lklSiteClass).
   called from minbranches()
*/
   int ig,i,ir;

   if(lklSiteClass==0)
      FOR(ig,com.ngene) {
         if(com.Mgene>1)
            SetPGene(ig,com.Mgene>1,com.Mgene>1,com.nalpha>1,x+tree.nbranch);
         PartialLikelihood(inode,ig);
      }
   else {
      FOR(ir,com.ncatG) {
#ifdef CODEML
         IClass=ir;
#endif
         if(ir) {
            _nodeScaleF+=_nnodeScale*com.npatt;
            for(i=com.ns;i<tree.nnode;i++)
               nodes[i].lkl+=(tree.nnode-com.ns)*com.ncode*com.npatt;
         }
         SetPSiteClass(ir, x+tree.nbranch);
         FOR(ig,com.ngene) {
            if(com.Mgene>1)
               SetPGene(ig,com.Mgene>1,com.Mgene>1,com.nalpha>1,x+tree.nbranch);
            if(com.nalpha>1) SetPSiteClass(ir, x+tree.nbranch);
            PartialLikelihood(inode,ig);
         }
      }  /* for(ir) */

      _nodeScaleF-=(com.ncatG-1)*_nnodeScale*com.npatt;
      for(i=com.ns;i<tree.nnode;i++)
         nodes[i].lkl-=(com.ncatG-1)*(tree.nnode-com.ns)*com.ncode*com.npatt;
   }

   return(0);
}


int minbranches (double *lnL,double x[],double tb[2],double e,double space[])
{
/* Z. Yang, November 1999.
   optimizing one branch at a time
   
   for each branch a..b, reroot the tree at b, and 
   then calculate partial likelihood for node a.
   For each branch, this routine determines the Newton search direction 
   p=-dl/dll.  It then halves the steplength to make sure -lnL is decreased.
   When the Newton solution is correct, this strategy will waste one 
   extra call to lfunt.  It does not seem possible to remove calculation of l 
   in lfuntddl().
   lfun or lfundG and thus SetParameters is called once beforehand to set up 
   global points like POMEGA.
   This works with NSsites and NSbranch models.
   At the end of the routine, nodes[].lkl are not updated.
*/
   int ib,oldroot=tree.root, a,b;
   int icycle, maxcycle=1000, icycleb,ncycleb=10, i;
   double lnL0=0, l0,l,dl,ddl, t,t0,t00, diffb, p,step=1, small=1e-15,y;

   lnL0=l0=l=*lnL;
   if(noisy>2) printf("lnL0 =    %14.6f\n",-lnL0);
   FOR(icycle,maxcycle) {
      if(noisy>2) printf("Cycle %2d: ",icycle+1);
      for(ib=0,diffb=0;ib<tree.nbranch;ib++) {
         t=t0=t00=nodes[tree.branches[ib][1]].branch; 
         l0=l;
         a=tree.branches[ib][0]; b=tree.branches[ib][1];
         FOR(i,tree.nnode) _oldlkl[i]=1;
         ReRootTree(b);
         updatelkl(x,a);
         for(icycleb=0; icycleb<ncycleb; icycleb++) {  /* iterating a branch */
            if(!lklSiteClass)
               lfuntdd(t,a,b,x, &y,&dl,&ddl,space);
            else
               lfuntdd_SiteClass(t,a,b,x, &y,&dl,&ddl,space);

            if((fabs(y-l)>1e-4 && noisy>2) || (fabs(y-l)>1e-5 && noisy>=9)) {
               printf("\nwarning rounding error? b=%d cycle=%d lnL=%12.7f != %12.7f\n",ib,icycleb,l,y);
               /* exit(-1); */
            }

            p=-dl/fabs(ddl);
            /* p=-dl/ddl; newton direction */
            if (fabs(p)<small) step=0;
            else if(p<0)       step=min2(1,(tb[0]-t0)/p);
            else               step=min2(1,(tb[1]-t0)/p);

            if(icycle==0 && step!=1 && step!=0) step*=0.99; /* avoid border */
            for (i=0; step>small; i++,step/=4) {
               t=t0+step*p;
               if(!lklSiteClass) lfunt(t, a,b,x, &l, space);
               else              lfunt_SiteClass(t, a,b,x, &l,space);
               if(l<l0) break;
            }
            if(step<=small) { t=t0; l=l0; break; }
            if(fabs(t-t0)<e*fabs(1+t) && fabs(l-l0)<e) break;
            t0=t; l0=l;
         }
         x[ib]=t; nodes[a].branch=t;

         if(fabs(t-t00)>diffb) diffb=fabs(t-t00);
      }   /* for (ib) */
      *lnL=l;
      if(noisy>2) printf("%14.6f\n",-l);
      if(fabs(*lnL-lnL0)<e && diffb<e*.1) break;
      lnL0= *lnL;
   }  /* for (icycle) */
   ReRootTree(oldroot);  /* did not update lkl */
   FOR(i,tree.nnode) _oldlkl[i]=0;

   return(0);
}





int lfunt(double t, int a,int b,double x[],double *l, double space[])
{
/* See notes for lfunt_dd and minbranches
*/
   int i,j,k, h,ig,n=com.ncode, nroot=n;
   int n1=(com.cleandata&&b<com.ns?1:n);
   double expt,uexpt=0,multiply;
   double *P=space, piqi,pqj, fh;

#if (CODEML)
   if (com.seqtype==CODONseq && com.model==2) {
      if(N_OMEGA<=5)
         { U=_UU[nodes[a].label]; V=_VV[nodes[a].label]; Root=_Root[nodes[a].label]; }
      else 
         EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, com.kappa,nodes[a].omega,PMat);
   }
#endif

#if (BASEML)
   if (com.nhomo==2)
      EigenTN93(com.model,nodes[a].kappa,1,com.pi,&nR,Root,Cijk);
   nroot=nR;
#endif

   *l=0;
   for (ig=0; ig<com.ngene; ig++) {
      if(com.Mgene>1) SetPGene(ig,1,1,0,x+tree.nbranch);
      FOR(i,n*n) P[i]=0;

      for(k=0,expt=1; k<nroot; k++) {
         multiply=com.rgene[ig]*Root[k];
         if(k) expt=exp(t*multiply);

#if (CODEML)  /* uses U & V */
         FOR(i,n) for(j=0,uexpt=U[i*n+k]*expt; j<n; j++)
            P[i*n+j]+=uexpt*V[k*n+j];
#elif (BASEML) /* uses Cijk */
         FOR(i,n) FOR(j,n)
            P[i*n+j]+=Cijk[i*n*nroot+j*nroot+k]*expt;
#endif
      }

      for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
         for(i=0,fh=0; i<n1; i++) {
            if(n1==1) piqi=com.pi[i=com.z[b][h]];
            else      piqi=com.pi[i]*nodes[b].lkl[h*n+i];

            for(j=0,pqj=0; j<n; j++)
               pqj+=P[i*n+j]*nodes[a].lkl[h*n+j];
            fh+=piqi*pqj;
         }
         if(fh<1e-250) 
            printf("a bit too small: fh = %10.6e\n",fh);
         *l-=log(fh)*com.fpatt[h];
         FOR(i,_nnodeScale) *l-=_nodeScaleF[i*com.npatt+h]*com.fpatt[h];
      }
   }
   return(0);
}


int lfuntdd(double t, int a,int b,double x[],double *l,double*dl,double*ddl,
    double space[])
{
/* Calculates lnL for branch length t for branch b->a.
   See notes in minbranches().
   partial likelihood updated correctly already.

   i for b, j for a?
*/
   int i,j,k, h,ig,n=com.ncode, nroot=n;
   int n1=(com.cleandata&&b<com.ns?1:n);
   double expt,uexpt=0,multiply;
   double *P=space, *dP=P+n*n,*ddP=dP+n*n, piqi,pqj,dpqj,ddpqj, fh,dfh,ddfh;

#if(CODEML)
   if (com.seqtype==CODONseq && com.model) {
      if(com.model==2 && N_OMEGA<=5)
         { U=_UU[nodes[a].label]; V=_VV[nodes[a].label]; Root=_Root[nodes[a].label]; }
      else 
         EigenQc(0,-1,NULL,NULL,NULL,Root,U,V, com.kappa,nodes[a].omega,PMat);
   }
#endif

#if (BASEML)
   if (com.nhomo==2)
      EigenTN93(com.model,nodes[a].kappa,1,com.pi,&nR,Root,Cijk);
   nroot=nR;
#endif
   *l=*dl=*ddl=0;
   for (ig=0; ig<com.ngene; ig++) {
      if(com.Mgene>1) SetPGene(ig,1,1,0,x+tree.nbranch);
      FOR(i,n*n) P[i]=dP[i]=ddP[i]=0;

      for(k=0,expt=1; k<nroot; k++) {
         multiply=com.rgene[ig]*Root[k];
         if(k) expt=exp(t*multiply);

#if (CODEML)  /* uses U & V */
         FOR(i,n) for(j=0,uexpt=U[i*n+k]*expt; j<n; j++) {
            P[i*n+j]+=uexpt*V[k*n+j];
            if(k) {
               dP[i*n+j]+=uexpt*V[k*n+j]*multiply;
               ddP[i*n+j]+=uexpt*V[k*n+j]*multiply*multiply;
            }
         }
#elif (BASEML) /* uses Cijk */
         FOR(i,n) FOR(j,n) {
            P[i*n+j]+=Cijk[i*n*nroot+j*nroot+k]*expt;
            if(k) {
               dP[i*n+j]+=Cijk[i*n*nroot+j*nroot+k]*expt*multiply;
               ddP[i*n+j]+=Cijk[i*n*nroot+j*nroot+k]*expt*multiply*multiply;
            }
         }
#endif
      }

      for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
         for(i=0,fh=dfh=ddfh=0; i<n1; i++) {
            if(n1==1) piqi=com.pi[i=com.z[b][h]];
            else      piqi=com.pi[i]*nodes[b].lkl[h*n+i];

            for(j=0,pqj=dpqj=ddpqj=0; j<n; j++) {
               pqj+=P[i*n+j]*nodes[a].lkl[h*n+j];
               dpqj+=dP[i*n+j]*nodes[a].lkl[h*n+j];
               ddpqj+=ddP[i*n+j]*nodes[a].lkl[h*n+j];
            }
            fh+=piqi*pqj;
            dfh+=piqi*dpqj;
            ddfh+=piqi*ddpqj;
         }
         if(fh<1e-250) {
            printf("too small: fh[%d] = %10.6e\n",h,fh);
            OutaTreeN(F0,0,1);
            exit(0);
         }
         *l-=log(fh)*com.fpatt[h];
         FOR(i,_nnodeScale) *l-=_nodeScaleF[i*com.npatt+h]*com.fpatt[h];
         *dl-=dfh/fh*com.fpatt[h];
         *ddl-=(-dfh*dfh+fh*ddfh)/(fh*fh) * com.fpatt[h];
      }
   }  /* for(ig) */
   return(0);
}


int lfunt_SiteClass(double t, int a,int b,double x[],double *l,double space[])
{
/* see notes in lfuntdd_SiteClass
   For branch&site models, look at the notes in GetPMatBranch()
*/
   int i,j,k, h,ig,ir,it, n=com.ncode, nroot=n;
   int n1=(com.cleandata&&b<com.ns?1:n);  /* one state for a tip? */
   double y,expt,uexpt=0,multiply, piqi,pqj;
   double *P=space, *fh=P+n*n;
   double *Sh=fh+com.npatt;  /* scale factor for each site pattern*/
   double *pK=com.fhK;  /* proportion for each site class after scaling */

   int iUVR=0,iUVR0[4][2]={{0,1}, {2,3}, {0,4}, {2,4}};

#if (BASEML)
   if (com.nhomo==2)
      EigenTN93(com.model,nodes[a].kappa,1,com.pi,&nR,Root,Cijk);
   nroot=nR;
#endif

   if(_nnodeScale==0) 
      FOR(ir,com.ncatG) FOR(h,com.npatt)  pK[ir*com.npatt+h]=com.freqK[ir];
   else {
      FOR(h,com.npatt) {
         for(ir=0,it=0; ir<com.ncatG; ir++) {
            for(k=0,y=0; k<_nnodeScale; k++)
               y+=_nodeScaleF[ir*_nnodeScale*com.npatt + k*com.npatt+h];
            if((pK[ir*com.npatt+h]=y)>pK[it*com.npatt+h]) it=ir;
         }
         Sh[h]=pK[it*com.npatt+h];
         FOR(ir,com.ncatG)
            pK[ir*com.npatt+h]=com.freqK[ir]*exp(pK[ir*com.npatt+h]-Sh[h]);
      }
   }

   FOR(h,com.npatt) fh[h]=0;
   FOR(ir,com.ncatG) {
      SetPSiteClass(ir, x+tree.nbranch);  /* com.ntime=0, x is original */

#if CODEML  /* branch b->a */
      if(com.seqtype==CODONseq && com.NSsites && com.model) { /* branch&site models */
         iUVR=iUVR0[ir][nodes[a].label];
         U=_UU[iUVR]; V=_VV[iUVR]; Root=_Root[iUVR]; 
      }
#endif

      if(ir) {
         for(i=com.ns;i<tree.nnode;i++)
            nodes[i].lkl+=(tree.nnode-com.ns)*n*com.npatt;
      }
      for (ig=0; ig<com.ngene; ig++) {
         if(com.Mgene>1)
            SetPGene(ig,com.Mgene>1,com.Mgene>1,com.nalpha>1,x+tree.nbranch);
         if(com.nalpha>1) SetPSiteClass(ir, x+tree.nbranch);

         FOR(i,n*n) P[i]=0;
         for(k=0,expt=1; k<nroot; k++) {
            multiply=com.rgene[ig]*Root[k]*_rateSite;
            if(k) expt=exp(t*multiply);

#if (CODEML)  /* uses U & V */
            FOR(i,n) for(j=0,uexpt=U[i*n+k]*expt; j<n; j++)
               P[i*n+j]+=uexpt*V[k*n+j];
#elif (BASEML) /* uses Cijk */
            FOR(i,n) FOR(j,n)
               P[i*n+j]+=Cijk[i*n*nroot+j*nroot+k]*expt;
#endif
         }  /* for (k), look through eigenroots */
         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            for(i=0; i<n1; i++) {
               if(n1==1) piqi=pK[ir*com.npatt+h]*com.pi[i=com.z[b][h]];
               else      piqi=pK[ir*com.npatt+h]*com.pi[i]*nodes[b].lkl[h*n+i];

               for(j=0,pqj=0; j<n; j++)
                  pqj+=P[i*n+j]*nodes[a].lkl[h*n+j];
               fh[h]+=piqi*pqj;
            }
         }  /* for (h) */
      }     /* for (ig) */
   }        /* for(ir) */

   for(i=com.ns;i<tree.nnode;i++)
      nodes[i].lkl-=(com.ncatG-1)*(tree.nnode-com.ns)*n*com.npatt;
   for(h=0,*l=0; h<com.npatt; h++) {
      
      if(fh[h]<1e-250) 
         printf("a bit too small: fh = %10.6e\n",fh[h]);

      *l-=log(fh[h])*com.fpatt[h];
      if(_nnodeScale) *l-=Sh[h]*com.fpatt[h];
   }
   return(0);
}


int lfuntdd_SiteClass(double t, int a,int b,double x[],
    double *l,double*dl,double*ddl,double space[])
{
/* dt and ddt for site-class models, modified from lfuntdd()
   nodes[].lkl (and _nodeScaleF if scaling is used) is shifted for ir, 
   and moved back to the rootal place at the end of the routine.

   At the start of this routine, nodes[].lkl has the conditional probabilties 
   for each node, each site pattern, for each site class (ir).  
   Scaling: When scaling is used, scale factors 
   _nodeScaleF[ir*_nnodeScale*com.npatt + k*com.npatt+h] for all nodes 
   are collected into Sh[h], after adjusting for rate classes, since the 
   sum is taken over ir.  Sh[h] and pK[ir*com.npatt+h] together store the 
   scale factors and proportions for site classes.  com.freqK[ir] is not 
   used in this routine beyond this point.
   if(com.Malpha), com.freqK[]=1/com.ncatG and does not change with ig, 
   and so the collection of Sh for sites at the start of the routine is o.k.

   The space for com.fhK[] is used.
   space[2*ncode*ncode + 4*npatt]:
    dP[ncode*ncode],ddP[ncode*ncode],fh[npatt],dfh[npatt],ddfh[npatt],Sh[npatt]
    pK[ncatG*npatt]=com.fhK[]
*/
   int i,j,k, h,ig,ir,it, n=com.ncode, nroot=n;
   int n1=(com.cleandata&&b<com.ns?1:n);  /* one state for a tip? */
   double y,expt,uexpt=0,multiply, piqi,pqj,dpqj,ddpqj;
   double *P=PMat, *dP=space,*ddP=dP+n*n;
   double *fh=ddP+n*n, *dfh=fh+com.npatt, *ddfh=dfh+com.npatt;
   double *Sh=ddfh+com.npatt;  /* scale factor for each site pattern */
   double *pK=com.fhK;  /* proportion for each site class after scaling */
   int iUVR=0,iUVR0[4][2]={{0,1}, {2,3}, {0,4}, {2,4}};

#if (BASEML)
   if (com.nhomo==2)
      EigenTN93(com.model,nodes[a].kappa,1,com.pi,&nR,Root,Cijk);
   nroot=nR;
#endif
   if(com.sspace<(int)(2*n*n+4*com.npatt)*(int)sizeof(double))
      error2("lfuntdd_SiteClass: not enough momery is allocated.\nLet me know");
   if(_nnodeScale==0) 
      FOR(ir,com.ncatG) FOR(h,com.npatt)  pK[ir*com.npatt+h]=com.freqK[ir];
   else {
      FOR(h,com.npatt) {
         for(ir=0,it=0; ir<com.ncatG; ir++) {
            for(k=0,y=0; k<_nnodeScale; k++)
               y+=_nodeScaleF[ir*_nnodeScale*com.npatt + k*com.npatt+h];
            if((pK[ir*com.npatt+h]=y)>pK[it*com.npatt+h]) it=ir;
         }
         Sh[h]=pK[it*com.npatt+h];
         FOR(ir,com.ncatG)
            pK[ir*com.npatt+h]=com.freqK[ir]*exp(pK[ir*com.npatt+h]-Sh[h]);
      }
   }

   FOR(h,com.npatt) fh[h]=dfh[h]=ddfh[h]=0;
   FOR(ir,com.ncatG) {
      SetPSiteClass(ir, x+tree.nbranch);

#if CODEML  /* branch b->a */
      if(com.seqtype==CODONseq && com.NSsites && com.model) { /* branch&site models */
         iUVR=iUVR0[ir][nodes[a].label];
         U=_UU[iUVR]; V=_VV[iUVR]; Root=_Root[iUVR]; 
      }
#endif

      if(ir) {
         for(i=com.ns;i<tree.nnode;i++)
            nodes[i].lkl+=(tree.nnode-com.ns)*n*com.npatt;
      }
      for (ig=0; ig<com.ngene; ig++) {
         if(com.Mgene>1)
            SetPGene(ig,com.Mgene>1,com.Mgene>1,com.nalpha>1,x+tree.nbranch);
         if(com.nalpha>1) SetPSiteClass(ir, x+tree.nbranch);

         FOR(i,n*n) P[i]=dP[i]=ddP[i]=0;
         for(k=0,expt=1; k<nroot; k++) {
            multiply=com.rgene[ig]*Root[k]*_rateSite;
            if(k) expt=exp(t*multiply);

#if (CODEML)  /* uses U & V */
            FOR(i,n) for(j=0,uexpt=U[i*n+k]*expt; j<n; j++) {
               P[i*n+j]+=uexpt*V[k*n+j];
               if(k) {
                  dP[i*n+j]+=uexpt*V[k*n+j]*multiply;
                  ddP[i*n+j]+=uexpt*V[k*n+j]*multiply*multiply;
               }
            }
#elif (BASEML) /* uses Cijk */
            FOR(i,n) FOR(j,n) {
               P[i*n+j]+=Cijk[i*n*nroot+j*nroot+k]*expt;
               if(k) {
                  dP[i*n+j]+=Cijk[i*n*nroot+j*nroot+k]*expt*multiply;
                  ddP[i*n+j]+=Cijk[i*n*nroot+j*nroot+k]*expt*multiply*multiply;
               }
            }
#endif
         }  /* for (k), look through eigenroots */
         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            for(i=0; i<n1; i++) {
               if(n1==1) piqi=pK[ir*com.npatt+h]*com.pi[i=com.z[b][h]];
               else      piqi=pK[ir*com.npatt+h]*com.pi[i]*nodes[b].lkl[h*n+i];

               for(j=0,pqj=dpqj=ddpqj=0; j<n; j++) {
                  pqj+=P[i*n+j]*nodes[a].lkl[h*n+j];
                  dpqj+=dP[i*n+j]*nodes[a].lkl[h*n+j];
                  ddpqj+=ddP[i*n+j]*nodes[a].lkl[h*n+j];
               }
               fh[h]+=piqi*pqj;
               dfh[h]+=piqi*dpqj;
               ddfh[h]+=piqi*ddpqj;
            }
         }  /* for (h) */
      }     /* for (ig) */
   }        /* for(ir) */

   for(i=com.ns;i<tree.nnode;i++)
      nodes[i].lkl-=(com.ncatG-1)*(tree.nnode-com.ns)*n*com.npatt;
   for(h=0,*l=*dl=*ddl=0; h<com.npatt; h++) {
      if(fh[h]<1e-250) 
         printf("a bit small: fh = %10.6e\n",fh[h]);

      *l-=log(fh[h])*com.fpatt[h];
      if(_nnodeScale) *l-=Sh[h]*com.fpatt[h];
      *dl-=dfh[h]/fh[h] * com.fpatt[h];
      *ddl-=(-dfh[h]*dfh[h]+fh[h]*ddfh[h])/(fh[h]*fh[h]) * com.fpatt[h];
   }

   return(0);
}




#endif         /* #ifdef LFUNCTIONS */

#ifdef BIRTHDEATH

void BranchLengthBD(int rooted, double birth, double death, double sample, 
     double mut)
{
/* Generate random branch lengths (nodes[].branch) using the birth and
   death process with species sampling, or the Yule (coalescent?) process
   if sample=0, when only parameter mut is used.
   Note: older interior nodes have larger node numbers, so root is at
   node com.ns*2-2 with time t[ns-2], while the youngest node is at 
   node com.ns with time t[0].  When unrooted=0, the root is removed with
   branch lengths adjusted.
   This works with the tree generated from RandomLHistory().
*/
   int i,j, it, imin,fixt0=1;
   double la=birth, mu=death, rho=sample, tmin, r, t[NS-1];
   double phi, eml, y;

   if (sample==0)  /* coalescent model */
      for (i=com.ns,y=0; i>1; i--) 
          nodes[com.ns*2-i].divtime=y += -log(rndu())/(i*(i-1.)/2.)*mut/2;
   else  {         /* BD with sampling */
      if (fixt0) t[com.ns-2]=1;
      if (fabs(la-mu)>1e-6) {
         eml=exp(mu-la);  phi=(rho*la*(eml-1)+(mu-la)*eml)/(eml-1);
         for (i=0; i<com.ns-1-(fixt0); i++) {
           r=rndu(); t[i]=log((phi-r*rho*la)/(phi-r*rho*la+r*(la-mu)))/(mu-la);
       }
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
      for (i=com.ns; i>1; i--) nodes[com.ns*2-i].divtime=t[com.ns-i]*mut;
   }
   FOR (i,com.ns) nodes[i].divtime=0;
   for (i=0; i<tree.nnode; i++) 
      if (i!=tree.root) 
         nodes[i].branch=nodes[nodes[i].father].divtime-nodes[i].divtime;
   if (!rooted) {
      it=nodes[tree.root].sons[2];
      nodes[it].branch =
      2*nodes[2*com.ns-2].divtime-nodes[tree.root].divtime-nodes[it].divtime;
   }
}

#endif


#ifdef NODESTRUCTURE

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
      j=(int)(i*rndu()); 
      nodes[nodea[j]].father=it; nodes[it].sons[0]=nodea[j];
      nodea[j]=nodea[i-1];
      j=(int)((i-1)*rndu()); 
      nodes[nodea[j]].father=it; nodes[it].sons[1]=nodea[j];
      nodea[j]=it;
      if (!rooted && i==3) {
         nodes[it].nson++; 
         nodes[nodea[1-j]].father=it; nodes[it].sons[2]=nodea[1-j];
      }
   }
   tree.root=it;  tree.nnode=ns*2-1-!rooted;
   NodeToBranch();
   return (0);
}


#endif  /* NODESTRUCTURE */


#ifdef EVOLVE

static double PMat[16], Qfactor, kappa1, kappa2;

int GenerateSeq (void)
{
/* makes sequence com.z[tree.root] with com.pi[] and evolves it
   along the tree, using nodes[].branch, com.ns, com.ls, com.ncode, com.model,
   com.kappa, com.alpha, com.rates (if com.alpha>0)
   z[ns*2-1] has 0,1,2,3 first and decoded before printing out.
*/
   int i,h;
   double T=com.pi[0],C=com.pi[1],A=com.pi[2],G=com.pi[3],Y=T+C,R=A+G;

   if (com.ncode!=4) error2 ("ncode");
   if (com.ns<2 || com.ls<1) error2 ("err GenerateSeq 1");

   if (com.model==F84) { kappa1=1+com.kappa/Y; kappa2=1+com.kappa/R; }
   else                  kappa1=kappa2=com.kappa;
   Qfactor=1/(2*T*C*kappa1+2*A*G*kappa2 + 2*Y*R);
   if (com.model>K80) 
      dnamaker (com.z[tree.root], com.ls, com.pi);
   else
      FOR(i,com.ls) com.z[tree.root][i]=(char)(rndu()*4.);

   if (com.model==JC69) EvolveJC (tree.root);
   else                 Evolve   (tree.root);
   com.npatt=0;
   FOR(i,com.ns) FOR(h,com.ls) {
      if(com.z[i][h]<0 || com.z[i][h]>3) error2("err GenerateSeq 2.");
#ifndef SIMULATE
      com.z[i][h]=BASEs[com.z[i][h]];
#endif
   }
   /* Ancestral sequences are not encoded here. */
   return (0);
}

void Evolve (int inode)
{
/* evolve sequence com.z[tree.root] along the tree to generate com.z[], 
   using nodes[].branch and  com.model
   Needs com.z[nnode], while com.z[0] -- com.z[ns-1] constitute the data.
   Nucleotides only.
*/
   int is, h,i,j, ison, ib, n=com.ncode;
   double t, r;
   
   for (is=0; is<nodes[inode].nson; is++) {
      ison=nodes[inode].sons[is];
      FOR(h,com.ls) com.z[ison][h]=com.z[inode][h];
      if (com.alpha==0) {
         if (com.model<=K80) PMatK80(PMat,nodes[ison].branch,com.kappa);
         else {
            t=nodes[ison].branch*Qfactor;
            PMatTN93 (PMat, t*kappa1, t*kappa2, t, com.pi);
         }
         FOR (i, n) for(j=1; j<n; j++)  PMat[i*n+j]+=PMat[i*n+j-1];
         FOR (h, com.ls) {
            for (j=0,ib=com.z[ison][h],r=rndu(); j<n; j++)
               if (r<PMat[ib*n+j]) break;
            if (j!=ib) com.z[ison][h]=j;
         }
      }
      else {
         FOR (h, com.ls) {
            if (h==0 || com.rates[h]!=com.rates[h-1]) {
               if (com.model<=K80) 
                  PMatK80 (PMat, nodes[ison].branch*com.rates[h], com.kappa);
               else {
                  t=nodes[ison].branch*Qfactor*com.rates[h];
                  PMatTN93 (PMat, t*kappa1, t*kappa2, t, com.pi);
               }
               FOR (i, n) for(j=1; j<n; j++)  PMat[i*n+j]+=PMat[i*n+j-1];
            }
            for (j=0,ib=com.z[ison][h],r=rndu(); j<n; j++)
               if (r<PMat[ib*n+j]) break;
            if (j-ib) com.z[ison][h]=j;
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
   When alpha>0, com.rates[] are the accumulative probabilities.
*/
   int is,j, h, nmut, imut, ison;
   double r;

   if (com.alpha && fabs(com.rates[com.ls-1]-1)>1e-4) 
      printf ("rates c.d.f.: 1 = %.6f?\n", com.rates[com.ls-1]);
   for (is=0; is<nodes[inode].nson; is++) {
      ison=nodes[inode].sons[is];
      FOR(h,com.ls) com.z[ison][h]=com.z[inode][h];
      nmut = rndpoisson (nodes[ison].branch*com.ls);
      for (imut=0; imut<nmut; imut++) {
         if (com.alpha==0) h=(int)(rndu()*com.ls);
         else 
            for (h=0,r=rndu(); h<com.ls; h++) if (r<com.rates[h]) break;
         j=(int)(rndu()*(com.ncode-1));
         if (j>=com.z[ison][h]) j++;
         com.z[ison][h]=j;
      }
      if (nodes[ison].nson) EvolveJC(ison);
   }
}

int PatternWeightSimple (int CollapsJC, double space[])
{
/* Counts site patterns, changes .z[], .fpatt[], .npatt, etc.
   This is modified from PatternWeight(), and does not deal with
   multiple genes.
*/   
   int  h, ht, b,j,k, same=0;
   char *zt[NS], zh[NS];

if(!com.cleandata) 
   puts("\nPatternWeightSimple: missing data not dealt with yet");

   FOR (j,com.ns) zt[j]=(char*)space+j*com.ls;
   FOR (h,com.ls) com.fpatt[h]=0; 

   for (h=0,com.npatt=0; h<com.ls; h++) {
      if (CollapsJC) {
         zh[0]=b=0;
         for (j=1; j<com.ns; j++) {
            for (k=0; k<j; k++) if (com.z[j][h]==com.z[k][h]) break;
            zh[j]=(k<j?zh[k]:++b);
         }
      }
      else 
         for (j=0; j<com.ns; j++) zh[j]=com.z[j][h];
      for (ht=0,same=0; ht<com.npatt; ht++) {
         for (j=0,same=1; j<com.ns; j++)
          if (zt[j][ht]!=zh[j]) { same=0; break; }
         if (same)  break; 
      }
      if (same)   com.fpatt[ht]++;
      else {
         FOR (j, com.ns) zt[j][com.npatt]=zh[j];
         com.fpatt[com.npatt++]=1;
      }
   }
   FOR (j,com.ns) {
      memcpy (com.z[j], zt[j], com.npatt);
      com.z[j][com.npatt]=0;
   }
   return (0);
}

#endif
