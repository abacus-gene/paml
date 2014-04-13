/* TREESUB.c
   subroutines that operates on trees, inserted into other programs 
   such as baseml, basemlg, codeml, and pamp etc.
*/

#include "tools.h"

extern char BASEs[], nBASEs[], *EquateNUC[], AAs[], BINs[];
extern int noisy;

#ifdef  BASEML
#define EIGEN
#define REALSEQUENCE
#define NODESTRUCTURE
#define TREESEARCH
#define PURTURBATION
#define STARDECOMPOSITION
#define STEPWISEADDITION
#define STEPWISEADDITION_MP
#define LSDISTANCE
#define LFUNCTIONS
#define RECONSTRUCTION
#endif 

#ifdef  CODEML
#define EIGEN
#define REALSEQUENCE
#define NODESTRUCTURE
#define TREESEARCH
#define PURTURBATION
#define STARDECOMPOSITION
#define STEPWISEADDITION
#define STEPWISEADDITION_MP
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
   puts("\nseq file looks strange.  Trying to read it as a paup data file.");

   for ( ; ; ) {
      if(fgets(line,lline,fseq)==NULL) error("seq err1: EOF");
      strcase(line,0);
      if(strstr(line,paupstart)) break;
   }
   for ( ; ; ) {
      if(fgets(line,lline,fseq)==NULL) error("seq err2: EOF");
      strcase(line,0);
      if((p=strstr(line,ntax))!=NULL) {
         while (*p != '=') { if(*p==0) error("seq err"); p++; }
         sscanf(p+1,"%d", &com.ns);
         if((p=strstr(line,nchar))==NULL) error("err: expect nchar");
         while (*p != '=') { if(*p==0) error("err"); p++; }
         sscanf(p+1,"%d", &com.ls);
         break;
      } 
   }
   printf("\nns: %d\tls: %d\n", com.ns, com.ls); 
   for ( ; ; ) {
      if(fgets(line,lline,fseq)==NULL) error("seq err1: EOF");
      strcase(line,0);
      if (strstr(line,paupend)) break;
   }
   return(0);
}

int ReadSeq (FILE *fout, FILE *fseq, int seqtype)
{
/* read in sequence, translate into protein (CODON2AAseq), and 
   count genes (ngene, lgene[]).
   seqtype: 0=nucleotides; 1=codons; 2:AAs; 3:CODON2AAs; 4:BINs
   com.pose[] used to store gene marks.  ls/3 gene marks for codon sequences
   char opt_c[]="AKGI";
      A:alpha given. K:kappa given
      G:many genes,  I:interlaved format
   lcode: # of characters in the final sequence 
*/
   char *line, *p,*p1, eq='.';
   int i,j,k, ch, noptline=0, lspname=LSPNAME, miss=0, simple=1, paupseq;
   int lline=255,lt[NS], igroup, lcode, Sequential=1,basecoding=0;
   int n31=(seqtype==CODONseq||seqtype==CODON2AAseq?3:1), gap=(n31==3?3:10);
   int nchar=(com.seqtype==AAseq?20:4);
   char *pch=((seqtype<=1||seqtype==CODON2AAseq)?BASEs:(seqtype==2?AAs:BINs));

   if (com.seqtype==4) error("seqtype==BINs, check with author");
   if (noisy && (com.seqtype<=CODONseq||com.seqtype<=CODON2AAseq)) {
      puts("\n\nAmbiguity character definition table:\n");
      FOR (i,strlen(BASEs)) {
         printf("%c (%d): ", BASEs[i],nBASEs[i]);
         FOR(j,nBASEs[i])  printf("%c ", EquateNUC[i][j]);
         FPN(F0);
      }
   }
   GetSeqFileType(fseq, &paupseq);

   if (com.ns>NS) error ("too many sequences.. raise NS?");
   if (com.ls%n31!=0) {
      printf ("\n%d nucleotides", com.ls); error ("seq. len. err.");
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
   if(com.z[com.ns-1]==NULL || com.pose==NULL||com.fpatt==NULL) error("oom");
   com.rgene[0]=1;   com.ngene=1;  

   lline=max(lline, com.ls+lspname+10);
   if((line=(char*)malloc(lline*sizeof(char)))==NULL) error("oom");
   FOR (j,com.ls) com.pose[j]=0;      /* gene #1, default */
   if(paupseq) goto readseq;

   /* first line */
   if (!fgets (line, lline, fseq)) error ("err ReadSeq: first line");
   for (j=0; j<lline && line[j] && line[j]!='\n'; j++) {
      if (!isalnum(line[j])) continue;
      line[j]=(char)toupper(line[j]);
      switch (line[j]) {
         case 'G': noptline++;   break;
         case 'C': basecoding=1; break;
         case 'S': Sequential=1; break;
         case 'I': Sequential=0; break;
         default : printf (".Bad option %c\n", line[j]);  exit (-1);
      }
   }
   if (strchr(line,'C')) {   /* protein-coding DNA sequences */
      if (seqtype!=0 || com.ls%3!=0 || noptline<1)
         error("Format for protein-coding sequences?");
      com.ngene=3; FOR(i,3) com.lgene[i]=com.ls/3;
      for (i=0;i<com.ls;i++) com.pose[i]=i%3;
      noptline--;
   }

   /* option lines */
   FOR (j, noptline) {
      line[0]=0;
      while (!isalnum(line[0]))  line[0]=(char)fgetc(fseq); 
      line[0]=(char)toupper(line[0]);
      switch (line[0]) {
      case ('G') :
         if(basecoding) error("incorrect option format, use GC?\n");
         if (fscanf(fseq,"%d",&com.ngene)!=1) error ("expecting #gene here..");
         if (com.ngene>NGENE) error ("raise NGENE?");
         if (com.fix_rgene) {    /* specified in the main program? */
            puts ("reading rates for genes from the seq. file, correct?");
            FOR (k, com.ngene) if (fscanf (fseq, "%lf", &com.rgene[k]) != 1)
               error ("err: fix_rgene?");
         }
         fgets (line, lline, fseq);
         if (!blankline(line)) {      /* #sites in each gene on the 2nd line */
            for (i=0,p=line; i<com.ngene; i++) {
               while (*p && !isalnum(*p)) p++;
               if (sscanf(p,"%d",&com.lgene[i])!=1) break;
               while (*p && isalnum(*p)) p++;
            }
            for (; i<com.ngene; i++)
               if (fscanf(fseq,"%d", &com.lgene[i])!=1) error("EOF at lgene");

            for (i=0,k=0; i<com.ngene; k+=com.lgene[i],i++)
               FOR (j, com.lgene[i]) com.pose[k+j]=i;
            for (i=0,k=0; i<com.ngene; i++) k+=com.lgene[i];
            if (k!=lcode) {
               matIout(F0, com.lgene, 1, com.ngene);
               printf("\n%6d != %d", lcode, k);  error("ls for genes");
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
            if(!fgets(line,lline,fseq)) error("err: option lines");
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
         if (!fgets (line, lline, fseq)) error ("err: EOF?");
         if (blankline(line)) {
            if (!fgets (line, lline, fseq)) error ("err: EOF");
            if (!hasbase(line))
               { printf("err: empty line (seq %d)\n",j+1); exit(-1); }
         }
         p=line+(line[0]=='=' || line[0]=='>') ;
         while(isspace(*p)) p++;
         if ((ch=strstr(p,"  ")-p)<lspname && ch>0) lspname=ch;
         strncpy (com.spname[j], p, lspname);
         p+=lspname;

         for (k=lspname; k>0; k--)
            if (!isgraph(com.spname[j][k]))   com.spname[j][k]=0;
            else    break;
         if (noisy) printf ("Reading seq #%2d: %s\n", j+1, com.spname[j]);
         for (k=0; k<com.ls; p++) {
            while (*p=='\n' || *p=='\0') {
               p=fgets(line, lline, fseq);
               if(p==NULL)
                  { printf("\nEOF at site %d, seq %d\n", k+1,j+1); exit(-1); }
            }
            *p=(char)toupper(*p);  
            p1=strchr(pch,*p);
            if (p1 && p1-pch>=nchar)  miss=1;
            if (*p==eq) {
               if (j==0) error ("err: . in 1st seq.?");
               com.z[j][k]=com.z[0][k];  k++;
            }
            else if (p1) 
               com.z[j][k++]=*p;
            else if (isalpha(*p)) 
               { printf("\nerr: %c at %d seq %d.",*p,k+1,j+1); exit(0); }
            else if (*p==EOF) error ("err: EOF?");
         }           /* for(k) */
      }                 /* for (j,com.ns) */
   }
   else   {              /* interlaved */
      if (noisy) printf ("\nInterlaved format..\n");
      FOR (j, com.ns) lt[j]=0;  /* temporary seq length */
      for (igroup=0; ; igroup++) {
         /*
         printf ("\nreading block %d ", igroup+1);  matIout(F0,lt,1,com.ns);*/

         FOR (j, com.ns) if (lt[j]<com.ls) break;
         if (j==com.ns) break;
         FOR (j,com.ns) {
            if (!fgets (line, lline, fseq)) {
               printf("\nerr reading site %d, seq %d group %d",
                  lt[j]+1,j+1,igroup+1);
               error ("err: EOF?");
            }
            if (!hasbase(line)) {
               if (j)  {
                  printf ("\n%d, seq %d group %d", lt[j]+1, j+1, igroup+1);
                  error ("err: empty line.");
               }
               else 
                  if (PopEmptyLines (fseq, lline, line)==-1) {
                     printf ("\n%d, seq %d group %d", lt[j]+1, j+1, igroup+1);
                     error ("err: EOF?");
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

               for (k=lspname; k>0; k--)   /* trim species names */
                  if (!isgraph(com.spname[j][k]))  com.spname[j][k]=0;
                  else   break;
               if(noisy) printf("Reading seq #%2d: %s\n",j+1,com.spname[j]);
            }
            for (; *p && *p!='\n'; p++) {
               if (lt[j]==com.ls) break;
               *p=(char)toupper(*p);
               p1=strchr(pch,*p);
               if (p1 && p1-pch>=nchar) miss=1;
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
               else if (*p==EOF) error ("err: EOF");
            }         /* for (*p) */
         }            /* for (j,com.ns) */
      }               /* for (igroup) */
   }
   free (line);
   if(!miss) com.cleandata=1;
   else if (com.cleandata) {  /* remove ambiguity characters */
      if(noisy)  puts("\nSites with gaps or missing data are removed.");
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
   IdenticalSeqs();

   if(n31==3) com.ls=lcode;

#ifdef CODEML
   if(com.seqtype==CODON2AAseq) {
      if (noisy>2) puts("\nTranslating into AA Sequences\n");
      FOR(i,com.ns) DNA2protein(com.z[i], com.z[i], com.ls,com.icode);
      com.seqtype=AAseq;
      if(fout) printsma(fout,com.spname,com.z,com.ns,com.ls,com.ls,10,0,NULL);
   }
#endif

   return (0);
}

int IdenticalSeqs(void)
{
/* This checks for identical sequences and create a data set of unique 
   sequences.  The file name is <SeqDataFile>.tmp.  This is casually 
   written and need more test.
   called right after the sequence data are read.
*/
   char tmpf[96], keep[NS];
   FILE *ftmp;
   int is,js,h, same, nkept=com.ns;

   FOR(is,com.ns) keep[is]=1;
   FOR(is,com.ns)  { 
      if(!keep[is]) continue;
      FOR(js,is)  {
         for(h=0,same=1; h<com.ls; h++)
            if(com.z[is][h]!=com.z[js][h]) break;
         if(h==com.ls) {
            printf("Seqs. %3d & %3d (%s & %s) seem identical!\n",
               js+1,is+1,com.spname[js],com.spname[is]);
            keep[is]=0;
         }
      }
   }
   FOR(is,com.ns) if(!keep[is]) nkept--;
   if(nkept<com.ns) {
      strcpy(tmpf,com.seqf);  strcat(tmpf,".tmp");
      if((ftmp=fopen(tmpf,"w"))==NULL)
         error("err IdenticalSeqs: file creation error");
      printSeqs(ftmp,NULL,keep,1);
      fclose(ftmp);
      printf("\nUniqe sequences collected in %s (Enter to continue).\n",tmpf);
      getchar();
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
*/
   int ncode=com.ncode, h,ht,j,k,newpatt,ig, *poset=(int*)space;
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
            if(j==com.ns&&k==nb) { puts("Do empty data now!"); getchar();}

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
   Ambiguity characters in nucleotide sequences are resolved by iteration. 
   For base frequencies in species, they are resolved within that sequence.
   For average base frequencies among species, they are resolved over all 
   species.

   This routine is called by baseml and aaml.  codonml uses another
   routine InitializeCode()
*/
   char *pch=(seqtype==0?BASEs:(seqtype==2?AAs:BINs)), indel[]="-?";
   int wname=30, h,j,k, ig, nconstp, nc=com.ncode;
   int irf,nrf=(com.seqtype==AAseq?1:4), miss=0, b,nb,ib[4];
   double pi0[20],t,lmax=0;

   PatternWeight(fout, space);

   if (fout)  fprintf (fout,"\nFrequencies..");
   for(h=0,nconstp=0; h<com.npatt; h++) {
      for (j=1;j<com.ns;j++)  if(com.z[j][h]!=com.z[0][h]) break;
      if (j==com.ns && com.z[0][h]!=indel[0] && com.z[0][h]!=indel[1])  
         nconstp+=(int)com.fpatt[h];
   }
   for (ig=0,zero(com.pi,nc); ig<com.ngene; ig++) {
      if (com.ngene>1)  fprintf (fout,"\n\nGene # %2d (len %4d)",
                             ig+1, com.lgene[ig]-(ig==0?0:com.lgene[ig-1]));
      fprintf(fout,"\n%*s",wname, ""); FOR(k,nc) fprintf(fout,"%7c",pch[k]);
      /* freqs in species first, com.piG[] are recalculated later. */
      for(j=0,zero(com.piG[ig],nc); j<com.ns; j++) {
         FOR(k,nc) { pi0[k]=1./nc; com.pi[k]=0; }
         for(irf=0; irf<nrf; irf++) {
            for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
               b=strchr(pch,com.z[j][h])-pch;
               if(b<0) error("Initialize: strange char");
               else if (b<nc)
                  com.pi[b]+=com.fpatt[h];
               else if(com.seqtype==BASEseq) {
                  NucListall(com.z[j][h], &nb, ib);
                  for(k=0,t=0; k<nb; k++) t+=pi0[ib[k]];
                  FOR(k,nb) com.pi[ib[k]]+=pi0[ib[k]]/t*com.fpatt[h];
               }
            }
            abyx(1/sum(com.pi,nc), com.pi, nc);
            if(distance(com.pi,pi0,nc)<1e-6) break;
            xtoy(com.pi,pi0,nc);
         }   /* for(irf) */
         fprintf(fout,"\n%-*s",wname,com.spname[j]);
         FOR(k,nc) fprintf(fout,"%7.4f",com.pi[k]);
         FOR(k,nc) com.piG[ig][k]+=com.pi[k]/com.ns;
      }    /* for(j) */
/*
      fprintf(fout,"\n\n%-*s",wname,"Mean");
      FOR(k,nc) fprintf(fout,"%7.4f",com.piG[ig][k]);
*/

      /* average freqs over species */
      for(irf=0,zero(com.piG[ig],nc); irf<nrf; irf++) {
         FOR (j,com.ns) {
            for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
               b=strchr(pch,com.z[j][h])-pch;
               if (b<nc)
                  com.piG[ig][b]+=com.fpatt[h];
               else if(com.seqtype==BASEseq) {
                  NucListall(com.z[j][h], &nb, ib);
                  for(k=0,t=0,miss=1; k<nb; k++) t+=pi0[ib[k]];
                  FOR(k,nb) com.piG[ig][ib[k]]+=pi0[ib[k]]/t*com.fpatt[h];
               }
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
   FOR (j,nc) fprintf (fout,"%7.4f", com.pi[j]);
   if(miss) ;

   fprintf (fout,"\n\n# constant sites: %6d (%6.2f%%)",
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
     FOR(j,com.ns) if(transform(com.z[j],com.npatt,1,com.seqtype)) error("uh?");

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
   int  h,k, is,lnew,nindel, n31=(com.seqtype==1||com.seqtype==3?3:1);
   int nchar=(com.seqtype==AAseq?20:4);
   char b, *p1,*pch=((com.seqtype<=1)?BASEs:(com.seqtype==2?AAs:BINs));
   char *miss;  /* miss[h]=1 if site (codon) h is missing, 0 otherwise */

   if (com.ls%n31) error ("err in RemoveIndel.");
   if((miss=(char*)malloc(com.ls/n31 *sizeof(char)))==NULL) 
      error("RemoveIndel: oom");
   FOR (h,com.ls/n31) miss[h]=0;
   for (is=0; is<com.ns; is++) {
      for (h=0,nindel=0; h<com.ls/n31; h++) {
         for (k=0; k<n31; k++) {
            b=(char)toupper(com.z[is][h*n31+k]);
            p1=strchr(pch,b);
            if (p1 && p1-pch>=nchar) { miss[h]=1; nindel++; }
         }
      }
      if (nindel) printf ("\n%6d characters removed in seq. %d", nindel, is+1);
   }
   for (h=0,lnew=0; h<com.ls/n31; h++)  {
      if(miss[h]) continue;
      for (is=0; is<com.ns; is++) {
         for (k=0; k<n31; k++) 
            com.z[is][lnew*n31+k]=com.z[is][h*n31+k];
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
   if (finf==NULL || fninf==NULL) error("MPInformSites: file creation error");

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

puts("PatternJC69like: does not work with missing data\n");
return(0);
   for (h=0,com.npatt=0,ig=-1; h<npatt0; h++) {
      if (ig<com.ngene-1 && h==com.posG[ig+1]) { com.posG[++ig]=com.npatt; }
      zh[0]=b=1;
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
   FOR (j, com.ns) com.z[j][com.npatt]=0;

   if (noisy>=2) printf ("\nnew no. site patterns:%7d\n", com.npatt);
   if (fout) fprintf (fout, "\nno. site patterns:%7d\n", com.npatt);

   if (fout) {
      FOR (h,com.npatt) {
         fprintf (fout," %3.0f", com.fpatt[h]);
         if ((h+1)%15==0) FPN (fout);
      }
      printsma (fout,com.spname,com.z,com.ns,com.npatt,60,10,1,NULL);
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
   int n=com.ncode, nb=(com.seqtype==CODONseq?3:1);

   FOR(h,ls) {
      hp=(pose?pose[h]:h);
      if(encoded) {
         if(com.seqtype!=CODONseq) fputc(pch[z[hp]],fout);
         else {
#ifdef CODEML
            fprintf(fout,"%s",getcodon(str,FROM61[z[hp]])); /* 0,1,...,60 */
#else
            fprintf(fout,"%s",getcodon(str,z[hp]));         /* 0,1,...,63 */
#endif
         }
      }
      else  /* raw data, not coded */
         FOR(i,nb) fputc(z[hp*nb+i],fout);
      if(gap && (com.seqtype==CODONseq || (h+1)%gap==0)) fputc(' ',fout);
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
   This uses com.seqtype.
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
      fprintf(fout,"%8d%8d\n",nskept,ls1);

   for(j=0; j<com.ns; j++,FPN(fout)) {
      if(keep && !keep[j]) continue;
      fprintf(fout,"%s%-*s  ", (format?"      ":""),wname,com.spname[j]);
      print1seq(fout, com.z[j],com.ls, com.cleandata, com.pose);
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

   if (testXMat(x)) error ("X err..");
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
/*
      if (1-fd/b<=0) return (-1);
*/
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
      /* *kappa = min(*kappa, largek);   */    *kappa = max(*kappa, -.5);
      return  4*bt*(tc*(1+ *kappa/Y)+ag*(1+ *kappa/R)+Y*R);
   case (HKY85):         /* TN93  */
      a1=1-Y*P1/(2*tc)-Q/(2*Y);  a2=1-R*P2/(2*ag)-Q/(2*R);   b=1-Q/(2*Y*R);
/*      if (a1<=0 || a2<=0 || b<=0) return (-1); */
      if (a1<=0 || a2<=0 || b<=0) return (fd);
      if (alpha<=0) { a1=-log(a1); a2=-log(a2); b=-log(b); }
      else   { a1=-gammap(a1,alpha); a2=-gammap(a2,alpha); b=-gammap(b,alpha);}
      a1t=.5/Y*(a1-R*b);  a2t=.5/R*(a2-Y*b);  bt=.5*b;
      *kappa = largek;
      if (bt>0) *kappa = min((a1t+a2t)/2/bt, largek);
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
   char *models[]={"JC69", "K80", "F81", "F84", "TN93"}, distf[32]="Distance.nuc";
   int is,js,k1,k2, h, miss=0, status=0, nc=4;
   double x[16], kappat=0, t,bigD=5;
   FILE *fdist=(FILE*)fopen(distf,"w");
   
   if(fdist==NULL) puts("DistnaceMatNuc: file creation error");
   else            fprintf(fdist,"%6d\n", com.ns);
   if (model>=HKY85) model=4; /* TN93 here */
   if (fout) {
      fprintf(fout,"\nDistances:%5s", models[model]);
      if (model!=JC69 && model!=F81) fprintf (fout, "(kappa) ");
      fprintf(fout," (alpha set at %.2f)\n", alpha);
      fprintf(fout,"\nThis matrix is not used in later m.l. analysis.\n");
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
            { puts("distance formula inapplicable..");   t=-1;  status=-1; }
         SeqDistance[is*(is-1)/2+js] = (t>=0?t:bigD);
         if(fdist) fprintf(fdist," %7.4f", t);
         if (fout) fprintf(fout,"%8.4f", t);
         if (fout && (model==K80 || model==F84 || model==HKY85))
            fprintf(fout,"(%7.4f)", kappat);
/*       fprintf(frst,"%5d%5d %8.4f %8.4f\n", is+1,js+1, t, kappat);
*/
       }
       if(fdist) FPN(fdist);
   }
   if(miss) fputs("\n(Ambiguity characters are not used in the above.)",fout);
   FPN(fout);
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
         error ("err:model");
      }
   }
#ifdef BASEMLG
   FOR (i,64) CijkIs0[i] = (Cijk[i]==0);
#endif
   return(0);
}

#ifdef EIGEN

int EigenREV (FILE* fout, double kappa[], double pi[], 
              int *nR, double Root[], double Cijk[])
{
/* pi[] is constant
*/
   int i,j,k;
   double Q[16], U[16], V[16], T1[16], T2[16], mr;

   *nR=4;
   zero (Q, 16);
   for (i=0,k=0; i<3; i++) for (j=i+1; j<4; j++)
      if (i*4+j!=11) Q[i*4+j]=Q[j*4+i]=kappa[k++];
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
      error ("\ncomplex roots in EigenREV?");
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

double *SeqDistance=NULL; 
int *ancestor=NULL;

int fun_LS (double x[], double diff[], int np, int npair);

int fun_LS (double x[], double diff[], int np, int npair)
{
   int i,j, aa, it;
   double dexp;

   if (SetBranch(x) && noisy>2) puts ("branch len.");
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
           testx, com.ns*(com.ns-1)/2, 1e-3, 1e-3);
   return (i);
}

#endif 

#ifdef NODESTRUCTURE

void BranchToNode (void)
{
/* tree.origin need to be specified before calling this
*/
   int i,j, from, to;
   
   tree.nnode=tree.nbranch+1;
/*
   if (tree.origin<0 || tree.origin>com.ns*2-2) 
      { printf ("root at %d", tree.origin+1); error ("tree origin"); }
*/
   FOR (j,tree.nnode) ClearNode (j);   /* Is this necessary? */
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
   NodeToBranchSub (tree.origin);
   if (tree.nnode!=tree.nbranch+1)  puts ("nd!=nb+1");
}

void ClearNode (int inode)
{
   nodes[inode].father=nodes[inode].ibranch=-1;
   nodes[inode].nson=0;
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
            error("err ReadaTreeB.");
      }
      if (tree.branches[j][0]<com.ns) YoungAncestor=1;
/*
      printf ("\nBranch #%3d: %3d -> %3d",
         j+1, tree.branches[j][0]+1,tree.branches[j][1]+1);
*/
   }
   if(popline) fgets(line, 254, ftree);
   tree.origin=tree.branches[0][0];
   if (YoungAncestor) error("Ancestors in the data?  To be fixed later.");

   com.ntime = com.clock ? tree.nnode-com.ns+(tree.origin<com.ns)
                         : tree.nbranch;

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
   int i, lline=10000;
   char line[10000];
   char *paupstart="begin tree", *paupend="translate", paupch=';',ch;

   fscanf (ftree, "%d%d", &i, ntree);
   if (i==com.ns) { *pauptree=0; return(0); }
   puts("\n# seqs in tree file does not match.  Read as a paup tree file.");
   for ( ; ; ) {
      if(fgets(line,lline,ftree)==NULL) error("tree err1: EOF");
      strcase(line,0);
      if (strstr(line,paupstart)) { *pauptree=1; *ntree=-1; break; }
   }
   if(shortform) return(0);
   for ( ; ; ) {
      if(fgets(line,lline,ftree)==NULL) error("tree err2: EOF");
      strcase(line,0);
      if (strstr(line,paupend)) break;
   }
   for ( ; ; ) {
      if((ch=fgetc(ftree))==EOF) error("tree err3: EOF");
      if (ch==paupch) break;
   }
   if(fgets(line,lline,ftree)==NULL) error("tree err4: EOF");

   return(0);
}

int PaupTreeRubbish(FILE *ftree);
int PaupTreeRubbish(FILE *ftree)
{
/* This pops out the rubbish, typically "tree PAUP_1 = [&U]" with
   "[&U]" optional, before each tree representation in a paup tree file
*/
   int lline=10000;
   char line[10000], *paup1="tree", treename[64],ch;

   fscanf (ftree, "%s", line);
   strcase(line,0);
   if (strstr(line,paup1)==NULL) error("err or end of treefile: expect tree");
   fscanf(ftree, "%s%s%c", treename, line, &ch);

   for (; ;) {
      if((ch=fgetc(ftree))==EOF) { 
         puts("err or end of treefile: expect ]");
         return(-1);
      }
      else if (ch==']') break;
      else if (ch=='(') error("err treefile: strange");
   }
   return(0);
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
   char check[NS], line[255], delimiters[]="(),:#";

   *hasbranch=0;
   while(isspace(ch)) ch=fgetc(ftree);
   ungetc(ch, ftree);
   if (isdigit(ch)) { ReadaTreeB (ftree, popline); return (0); }

   tree.nnode=com.ns;  tree.nbranch=0;
   FOR (i, 2*com.ns-1) ClearNode (i);
   FOR (i, com.ns) check[i]=0;
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
            tree.origin=cnode;
         cfather=cnode;
      }
      else if (ch==')')
         { level--;  inodeb=cfather; cfather=nodes[cfather].father; }
      else if (ch==':') {
         *hasbranch=1; fscanf(ftree,"%lf",&nodes[inodeb].branch);
      }
#ifdef BASEML
      else if (ch=='#') {
         fscanf(ftree,"%lf",&nodes[inodeb].parmark);
      }
#endif
      else if (ch==',') ;
      else { /* read species name or number */
         line[0]=(char)ch;  line[1]=(char)fgetc(ftree);
         if (com.ns<10 && isdigit(line[0]) && isdigit(line[1])) 
           { ungetc(line[1], ftree); line[1]=0; }
         else 
            for (i=1; ; )  { 
               if (strchr(delimiters,line[i]))
                  { ungetc(line[i], ftree); line[i]=0; break; }
               line[++i]=(char)fgetc(ftree);
         }
         for(j=i-1; j>0; j--)    if (!isgraph(line[j])) line[j]=0;
         for(i=0,hasname=0; line[i]; i++)  if (!isdigit(line[i])) hasname=1;

         if (hasname) {   /* name */
            for (i=0; i<com.ns; i++) if (!strcmp(line,com.spname[i])) break;
            if ((cnode=i)==com.ns) printf("\nSpecies %s?\n", line);
         }
         else {           /* number */
            sscanf(line, "%d", &cnode);   cnode--;
            if (cnode<0 || cnode>=com.ns)
               { printf("\nspecies # %d in ReadaTreeN\n", cnode+1); exit(-1); }
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

   FOR(i,com.ns)  if(check[i]!=1) return(-1);
   if(tree.nbranch>2*com.ns-2) { 
      printf("nbranch %d", tree.nbranch); error("too many branches in tree?");
   }

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
      if (branchlen)fprintf(fout,":%.6f", nodes[nodes[inode].sons[i]].branch);

      if (com.ns>9 || spnames || branchlen)
         if (i<nodes[inode].nson-1) fprintf(fout, ", ");
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
/* This points the nodes[com.ns+inode].lkl to the right space in com.chunk.
   The space is different depending on com.cleandata (0 or 1)
*/
   int inode, start=(com.cleandata?0:com.ns);

   for (inode=0; inode<tree.nnode-com.ns; inode++)
      nodes[com.ns+inode].lkl=com.chunk+com.ncode*com.npatt*(inode+start);
}

void SetDivTime(int inode, double x[]);
void SetDivTime(int inode, double x[])
{
   int i,ison;

   FOR (i,nodes[inode].nson) {
      nodes[ ison=nodes[inode].sons[i] ].divtime=0;
      if(nodes[ison].nson) {
         nodes[ison].divtime=nodes[inode].divtime*x[ison-com.ns];
         SetDivTime(ison, x);
      }
   }
}

int SetBranch (double x[])
{
   int i, status=0;
   double small=-1e-6;

   if (com.clock) {
      nodes[tree.origin].divtime=x[0];
      SetDivTime(tree.origin, x);
      FOR (i,tree.nnode)
         if (i!=tree.origin) {
            nodes[i].branch=nodes[nodes[i].father].divtime-nodes[i].divtime;
            if(nodes[i].branch<small)  status=-1;
         }

#if (BASEML || CODEML)
  /*  for clock=2 for baseml, disabled right now for debugging listtree.c */

      if(N_rateBranch>1) FOR(i,tree.nbranch) if(rateBranch[i]) 
       nodes[tree.branches[i][1]].branch*=x[tree.nnode-com.ns+rateBranch[i]-1];
/*
OutaTreeN(F0,0,1); FPN(F0);
matIout(F0,rateBranch,1,tree.nbranch);
matout(F0,x,1,com.np);
*/
#endif


   }
   else
      FOR(i,tree.nnode) {
         if(i!=tree.origin) 
            if ((nodes[i].branch=x[nodes[i].ibranch])<small)  status=-1;
      }

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
   if (nib!=tree.nbranch-com.ns) error("err BranchPartition"); 
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
      it=tree.origin;  if(tree.origin>=is) it+=2;
      FOR(i,2) tree.branches[tree.nbranch+i][0]=tree.origin=is+1;
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
      if (tree.origin>=is) tree.origin+=2;
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
#endif

#ifdef PURTURBATION

int Perturbation(FILE* fout, int initialMP, double space[]);

int Perturbation(FILE* fout, int initialMP, double space[])
{
/* heuristic tree search by the NNI tree perturbation algorithm.  
   Some trees are evaluated multiple times as no trees are kept.
   This needs more work.
*/
   int step=0, ntree=0, nmove=0, improve=0, ineighb, i;
   double *x=treestar.x;
   FILE *ftree;

if(com.clock) error("\n\aerr: pertubation does not work with a clock yet.\n");
if(initialMP) error("\ncannot get initial parsimony tree yet.");

   fprintf(fout, "\n\nHeuristic tree search by NNI perturbation\n");
   if (initialMP) {
      if (noisy) printf("\nInitial tree from stepwise addition with MP:\n");
      fprintf(fout, "\nInitial tree from stepwise addition with MP:\n");
      StepwiseAdditionMP (space);
   }
   else {
      if (noisy) printf ("\nInitial tree read from file %s:\n", com.treef);
      fprintf(fout, "\nInitial tree read from file.\n");
      if ((ftree=fopen (com.treef,"r"))==NULL) error ("treefile not exist?");
      fscanf (ftree, "%d%d", &i, &ntree);
      if (i!=com.ns) error ("ns in the tree file");
      if(ReadaTreeN (ftree, &i, 1)) error ("err tree..");
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
   treebest.tree=tree;  memcpy(treebest.nodes, nodes, sizeof(nodes));

   for (step=0; ; step++) {
      for (ineighb=0,improve=0; ineighb<(tree.nbranch-com.ns)*2; ineighb++) {
         tree=treebest.tree; memcpy (nodes, treebest.nodes, sizeof(nodes));
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
            treebest.tree=tree;  memcpy (treebest.nodes, nodes, sizeof(nodes));
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
   tree=treebest.tree;  memcpy (nodes, treebest.nodes, sizeof(nodes));
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

#endif /* PURTURBATION */

#ifdef STEPWISEADDITION_MP

static int *_U0, *_step0, _mnnode;
/* up pass characters and changes for the star tree: each of size npatt*nnode*/

int StepwiseAdditionMP (double space[])
{
/* tree search by species addition.
*/
   char *z0[NS];
   int  ns0=com.ns, is, i,j,h, tiestep=0,tie,bestbranch=0;
   double bestscore=0,score;

   _mnnode=com.ns*2-1;
   _U0=(int*)malloc(com.npatt*_mnnode*sizeof(int));
   _step0=(int*)malloc(com.npatt*_mnnode*sizeof(int));
   if (noisy>2) 
     printf("\n%ld bytes for MP (U0 & N0)\n", 2*com.npatt*_mnnode*sizeof(int));
   if (_U0==NULL || _step0==NULL) error ("oom in StepwiseAdditionMP()");

   FOR (i,ns0)  z0[i]=com.z[i];
   tree.nbranch=tree.origin=com.ns=3;
   FOR (i, tree.nbranch) { tree.branches[i][0]=com.ns; tree.branches[i][1]=i; }
   BranchToNode ();
   FOR (h, com.npatt)
      FOR (i,com.ns)
        { _U0[h*_mnnode+i]=1<<(com.z[i][h]-1); _step0[h*_mnnode+i]=0; }
   for (is=com.ns,tie=0; is<ns0; is++) {
      treestar.tree=tree;  memcpy (treestar.nodes, nodes, sizeof(nodes));

      for (j=0; j<treestar.tree.nbranch; j++,com.ns--) {
         tree=treestar.tree;  memcpy (nodes, treestar.nodes, sizeof(nodes));
         com.ns++;
         AddSpecies (is, j);
         score=MPScoreStepwiseAddition (is, space, 0);
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
      tree=treestar.tree;  memcpy (nodes, treestar.nodes, sizeof(nodes));
      com.ns=is+1;
      AddSpecies (is, bestbranch);
      score=MPScoreStepwiseAddition (is, space, 1);

      if (noisy)
       { printf("\r  Added %d [%5.0f steps]",is+1,bestscore); fflush(F0);}
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
      for (ist=is; (father=nodes[ist].father)!=tree.origin; ist=father) {
         if ((son2=nodes[father].sons[0])==ist)  son2=nodes[father].sons[1];
         N[father]=N[ist]+N[son2];
         if ((U[father]=U[ist]&U[son2])==0)
            { U[father]=U[ist]|U[son2];  N[father]++; }
      }
      FOR (i,3) U3[i]=U[nodes[tree.origin].sons[i]];
      N[tree.origin]=2;
      if (U3[0]&U3[1]&U3[2]) N[tree.origin]=0;
      else if (U3[0]&U3[1] || U3[1]&U3[2] || U3[0]&U3[2]) N[tree.origin]=1;
      FOR(i,3) N[tree.origin]+=N[nodes[tree.origin].sons[i]];

      if (save) {
         memcpy (pU0, U, tree.nnode*sizeof(int));
         memcpy (pN0, N, tree.nnode*sizeof(int));
      }
      score+=N[tree.origin]*com.fpatt[h];
   }
   return (score);
}

#endif

#ifdef STEPWISEADDITION

int StepwiseAddition (FILE* fout, double space[])
{
/* heuristic tree search by species addition.  Species are added in the order 
   of occurrence in the data.
   No attempts were made to get good initial values, although this should 
   be possible.
*/
   char *z0[NS];
   int  ns0=com.ns, is, i,j, bestbranch=0;
   double bestscore=0,score, *x=treestar.x;

   if (noisy) printf("\n\nHeuristic tree search by stepwise addition\n");
   if (fout) fprintf(fout, "\n\nHeuristic tree search by stepwise addition\n");
   FOR (i,ns0)  z0[i]=com.z[i];
   tree.nbranch=tree.origin=com.ns=(com.clock?2:3);  

   FOR (i, tree.nbranch) { tree.branches[i][0]=com.ns; tree.branches[i][1]=i; }
   BranchToNode ();
   for (is=com.ns; is<ns0; is++) {                  /* add the is_th species */
      treestar.tree=tree;  memcpy (treestar.nodes, nodes, sizeof(nodes));

      for (j=0; j<treestar.tree.nbranch+(com.clock>0); j++,com.ns--) { 
         tree=treestar.tree;  memcpy (nodes, treestar.nodes, sizeof(nodes));
         com.ns++;
         AddSpecies (is, j);
         score=TreeScore(x, space);
         if (noisy>1)
            { printf("\n "); OutaTreeN(F0, 0, 0); printf("%12.2f",-score); }

         if (j==0 || score<bestscore || (score==bestscore&&rndu()<.2)) {
            treebest.tree=tree;  memcpy(treebest.nodes, nodes, sizeof(nodes));
            xtoy (x, treebest.x, com.np);
            bestscore=score; bestbranch=j;
         }
      }
      tree=treebest.tree;  memcpy (nodes, treebest.nodes, sizeof(nodes));
      xtoy (treebest.x, x, com.np);
      com.ns=is+1;

      if (noisy) {
         printf("\n\nAdded sp. %d, %s [%.2f]\n",is+1,com.spname[is],bestscore);
         OutaTreeN(F0,1,1);  puts(";");
         if (com.np>com.ntime) {
            printf("\tparameters:"); 
            for(i=com.ntime; i<com.np; i++) printf("%9.5f", x[i]);
            FPN(F0);
         }
      }
      if (fout) {
         fprintf(fout,"\n\nAdded sp. %d, %s [%.2f]\n",
                 is+1, com.spname[is], bestscore);
         OutaTreeN(fout,0,0); fputs(";\n",fout);
         OutaTreeN(fout,1,1); fputs(";\n",fout);
         if (com.np>com.ntime) {
            fprintf(fout, "\tparameters:"); 
            for(i=com.ntime; i<com.np; i++) fprintf(fout, "%9.5f", x[i]);
            FPN(fout);
         }
         fflush(fout);
      }
   }  
   tree.lnL=bestscore;
   return (0);
}

#endif

#ifdef STARDECOMPOSITION

int DecompTree (int inode, int ison1, int ison2);
#define hdID(i,j) (max(i,j)*(max(i,j)-1)/2+min(i,j))

int StarDecomposition (FILE *fout, double space[])
{
/* automatic tree search by star decomposition, nhomo<=1
   returns (0,1,2,3) for the 4s problem.
*/
   int status=0,stage=0, i,j, itree,ntree=0,ntreet,best=0,improve=1,collaps=0;
   int inode, nson=0, ison1,ison2, son1, son2;
   double x[NP],xb[NP][2], e=1e-7;
   FILE *ftree, *fsum=frst;
/*
fsum=NULL;
*/
   if (com.runmode==1) {   /* read the star-like tree from trees.*s */
      if ((ftree=fopen (com.treef,"r"))==NULL) error ("no treefile");
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
   GetInitials (x);   NFunCall=0;
   SetxBound (com.np, xb);
   if(!com.cleandata) InitPartialLikelihood ();
   status=ming2(frub,&tree.lnL,com.plfun,NULL,x,xb, space,e,com.np);

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
      treestar=treebest;  memcpy(nodes,treestar.nodes,sizeof(nodes));

      if (collaps && stage) { 
         printf ("\ncollapsing nodes\n");
         OutaTreeN (F0, 1, 1);  FPN(F0);

         tree=treestar.tree;  memcpy(nodes, treestar.nodes, sizeof(nodes));
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
            for (i=1,x[0]=max(x[0],.01); i<com.ntime; i++)  x[i]=.5;
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
         SetxBound (com.np, xb);
         if(!com.cleandata) InitPartialLikelihood ();

         status=ming2(frub,&tree.lnL,com.plfun,NULL,x,xb, space,e,com.np);

         if (tree.lnL<treebest.tree.lnL) {
            treebest.tree=tree;  memcpy (treebest.nodes, nodes, sizeof(nodes));
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
   int ig, j, ngene0, npatt0, lgene0[NGENE], posG0[NGENE+1];

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
      FOR (j,com.ns) com.z[j]+=posG0[ig];   com.fpatt+=posG0[ig];
      xtoy (com.piG[ig], com.pi, com.ncode);

      printf ("\n\nGene #%4d  ls:%4d  npatt:%4d\n",ig+1,com.ls,com.npatt);
      fprintf(fout,"\nGene #%4d  ls:%4d  npatt:%4d\n",ig+1,com.ls,com.npatt);

      fprintf(frst,"\n%2d %4d %4d ",ig+1,com.ls,com.npatt);

      if (com.runmode==0)  Forestry(fout, space);
#ifdef CODEML
      if (com.runmode==-2) {
         if(com.seqtype==CODONseq) PairwiseCodon(fout, space);
         else                      PairwiseAA(fout);
      }
#endif
      else                 StepwiseAddition(fout, space);

      FOR (j,com.ns) com.z[j]-=posG0[ig];
      com.fpatt-=posG0[ig];
   }
   com.ngene=ngene0;  com.npatt=npatt0;  com.ls=lgene0[ngene0-1];
   FOR (ig, ngene0)   com.lgene[ig]=lgene0[ig];
   FOR (ig, ngene0+1) com.posG[ig]=posG0[ig];

   return (0);
}

int printSeqsMgenes (void)
{
/* separate sites from different partitions (genes) into different files.
   called before sequences are recoded
*/
   FILE *fseq;
   char seqf[20];
   int ig, lg, i,j,h, n13=1;

   printf ("Separating sites in genes into different files.\n");
   if (com.seqtype==CODONseq||com.seqtype==CODON2AAseq) {
      n13=3;
      puts ("Check the resulting files carefully before use (Press Enter).");
      getchar ();
   }
   for (ig=0, FPN(F0); ig<com.ngene; ig++) {
      for (h=0,lg=0; h<com.ls/n13; h++)  if (com.pose[h]==ig) lg++;
      sprintf (seqf, "Gene%d.seq\0", ig+1);
      if ((fseq=fopen(seqf,"w"))==NULL) error("file creation err.");
      printf ("%d sites in gene %d go to file %s\n", lg, ig+1,seqf);

      fprintf (fseq, "%8d%8d\n", com.ns, lg*n13);
      for (j=0; j<com.ns; FPN(fseq),j++) {
         fprintf(fseq,"%-20s  ", com.spname[j]);
         if (n13==1)  {       /* nucleotide or aa sequences */
            FOR (h,com.ls)   
               if (com.pose[h]==ig) fprintf (fseq, "%c", com.z[j][h]);
         }
         else {               /* codon sequences */
            FOR (h,com.ls/n13)
               if (com.pose[h]==ig) {
                  FOR (i,3) fprintf(fseq,"%c", com.z[j][h*3+i]);
                  fputc(' ', fseq);
               }
         }
      }
      fclose (fseq);
   }
   return (0);
}

#endif   /* ifdef REALSEQUENCE */
#endif   /* ifdef STARDECOMPOSITION */
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
*/
   int h, i, bintree,U[3],change;
   double score;

   Nsteps=(int*)space;
   bintree=(tree.nnode==2*com.ns-1 - (nodes[tree.origin].nson==3));
   if (bintree)  chUB=Nsteps+tree.nnode;
   else {
      chU=(char*)(Nsteps+tree.nnode);
      NchU=chU+tree.nnode*com.ncode;  Kspace=NchU+tree.nnode;
   }
   for (h=0,score=0; h<com.npatt; h++) {
      FOR (i,tree.nnode) Nsteps[i]=0;
      if (bintree) {           /* binary trees, use bit operation */
         FOR (i,com.ns)  chUB[i]=1<<(com.z[i][h]);
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
         FOR(i,com.ns)
            {chU[i*com.ncode]=(char)(com.z[i][h]); NchU[i]=(char)1; }
         for (i=com.ns; i<tree.nnode; i++)  NchU[i]=0;
         UpPassScoreOnly (tree.origin);
      }
      score+=Nsteps[tree.origin]*com.fpatt[h];
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
   UpPass (tree.origin);
   *nchange=Nsteps[tree.origin];
   if (job==0) return (0);
   FOR (i,n) chMark[tree.origin*n+i]=chMarkU[tree.origin*n+i];
   DownPass (tree.origin);
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

   for (j=0,visit[i=0]=(char)(tree.origin-com.ns); j<tree.nbranch; j++) 
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
      ICharaCur[j=tree.origin-com.ns]=0;  PATHWay[j]=CharaCur[j*n+0];
      FOR (j,nid) Equivoc[j]=(char)(NCharaCur[j]>1);
      DownStates (tree.origin);

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
                     if (it==tree.origin || nodes[it].father==i) break;
                  if (it==tree.origin)
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



#ifdef LFUNCTIONS
#ifdef RECONSTRUCTION


int ProbSitePattern (double x[], double *lnL, double *fhs);
int AncestralMarginal (FILE *fout, double x[], double fhs[]);
int AncestralJoint (FILE *fout, double x[], double fhs[], double space[]);

int ProbSitePattern (double x[], double *lnL, double *fhs)
{
/* allocate memory for fhs[] and calculate probabilities for 
   observed site patterns.  This routine works with constant or gamma rates.
*/
   int ig, i,h, ir;
   double fh;

   if (SetParameters(x)) puts ("par err.");
   zero (fhs, com.npatt);
   if (com.alpha==0) {
      for (ig=0,*lnL=0; ig<com.ngene; ig++) {
         if (com.Mgene>1) SetPGene(ig, 1, 1, 0, x);
         PartialLikelihood (tree.origin, ig);
         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            for (i=0; i<com.ncode; i++) 
               fhs[h] += com.pi[i]*nodes[tree.origin].lkl[h*com.ncode+i];
            *lnL -= log(fhs[h])*com.fpatt[h];
         }
      }
   }
   else {
      for (ig=0; ig<com.ngene; ig++) {
         SetPGene(ig, com.Mgene>1, com.Mgene>1, com.nalpha>1, x);
         for (ir=0; ir<com.ncatG; ir++) {
            FOR (i, tree.nnode)
               nodes[i].branch *= (ir==0?com.rK[ir]:com.rK[ir]/com.rK[ir-1]);
            PartialLikelihood (tree.origin, ig);
            for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
               for (i=0,fh=0; i<com.ncode; i++)
                  fh += com.pi[i]*nodes[tree.origin].lkl[h*com.ncode+i];
               fhs[h] += com.freqK[ir]*fh;
            }
         }
         FOR(i,tree.nnode)  nodes[i].branch/=com.rK[com.ncatG-1];
      }
      for (h=0,*lnL=0; h<com.npatt; h++)
         *lnL-=log(fhs[h])*com.fpatt[h];
   }
   return (0);
}


int AncestralMarginal (FILE *fout, double x[], double fhs[])
{
/* Ancestral reconstruction for each interior node.  This works under both 
   the one rate and gamma rates models.
   pnode[npatt*nid] stores the prob for the best chara at a node and site.
   The best character is kept in za[0,nnode-ns][], coded as 0,...,n-1.
   The data may be coded (com.cleandata==1) or not (com.cleandata=0).
   Call ProbSitePatt() before running this routine.
   pMAPnode[NS-1], pMAPnodeA[] stores the MAP probabilities for a site and 
      for the entire sequence.
   pChar1node[npatt*ncode] stores prob for each char at each pattern at 
   the current node (in).
*/
   char *pch=(com.seqtype==0?BASEs:(com.seqtype==2?AAs:BINs)), *za[NNODE];
   char str[4]="",aa3[4]="";
   int n=com.ncode, inode,ir, ic,i,j,h,hp,ig, best, oldroot=tree.origin;
   int nid=tree.nnode-com.ns, nb=(com.seqtype==CODONseq?3:1), wname=30;
   double lnL=0, fh, y, pbest, *pChar1node, *pnode;
   double pMAPnode[NS-1], pMAPnodeA[NS-1], smallp=0.0001;

   if(noisy) puts("\nMarginal reconstruction of ancestral sequences\n");
   if(com.ns>100||(com.seqtype==CODONseq&&com.ns>60)) 
      puts("\n(need test for large data sets)\n");

   fprintf (fout,"\n(1) Marginal reconstruction of ancestral sequences ");
   fprintf (fout,"(eq. 4 in Yang et al. 1995).\n    This algorithm is OK.\n");
   pChar1node=(double*)malloc(com.npatt*n*sizeof(double));
   pnode=(double*)malloc((nid*com.npatt+1)*(sizeof(double)+sizeof(char)));
   if (pnode==NULL||pChar1node==NULL) error ("oom");
   FOR (j,nid) za[j]=(char*)(pnode+nid*com.npatt)+j*com.npatt;

   if (SetParameters(x)) puts ("par err.");
   if (com.verbose>1) 
      fprintf(fout,"\nProb distribs at nodes, those with p<%.4f not listed\n",
         smallp);
   for (inode=com.ns; inode<tree.nnode; inode++) {
      /* This loop reroots the tree at inode, to reconstruct seq at inode */
      zero (pChar1node,com.npatt*n);
      if(noisy) printf ("Node %2d: ", inode+1);
      ReRootTree (inode);
      if (com.alpha==0) {   /* constant rate for sites (ngene>1 OK) */
         for (ig=0,lnL=0; ig<com.ngene; ig++) {
            if (com.Mgene>1) SetPGene(ig, 1, 1, 0, x);
            PartialLikelihood (tree.origin, ig);
            for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
               for (i=0,fh=0,pbest=0,best=-1; i<n; i++) {
                  y=com.pi[i]*nodes[tree.origin].lkl[h*n+i];
                  fh +=  y;
                  if ((y=y/fhs[h])>pbest) { pbest=y; best=i; }
                  pChar1node[h*n+i]=y;
               }
               if (fabs(fh-fhs[h])>min(1e-6,fh*.001)) error("err: fh 1?");
               lnL -= log(fh)*(double)com.fpatt[h];
               za[inode-com.ns][h]=(char)best;
               pnode[(inode-com.ns)*com.npatt+h]=pbest;
            }
         }
      }
      else {                /* gamma rates among sites */
         /* sum over the rate distribution (ir) to calculate pChar1node[] */
         zero(com.fhK,com.npatt); /* calculates lnL for error checking */
         for (ig=0; ig<com.ngene; ig++) {
            SetPGene(ig, com.Mgene>1, com.Mgene>1, com.nalpha>1, x);
            for (ir=0; ir<com.ncatG; ir++) {
               FOR (i, tree.nnode)
                 nodes[i].branch *= (ir==0?com.rK[ir]:com.rK[ir]/com.rK[ir-1]);
               PartialLikelihood (tree.origin, ig);
               for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
                  for (i=0,fh=0; i<n; i++) {
                     fh += (y=com.pi[i]*nodes[tree.origin].lkl[h*n+i]);
                     pChar1node[h*n+i] += com.freqK[ir]*y;
                  }
                  com.fhK[h] += com.freqK[ir]*fh;
               }
            }
            FOR(i,tree.nnode)  nodes[i].branch/=com.rK[com.ncatG-1];
         }
         for (h=0,lnL=0; h<com.npatt; h++) {
            FOR (i,n) pChar1node[h*n+i]/=fhs[h];
            if (fabs(1-sum(pChar1node+h*n,n))>1e-5) error("Err: sum!=1");

            for (i=0,best=-1,pbest=-1; i<n; i++)
               if (pChar1node[h*n+i]>pbest) {best=i; pbest=pChar1node[h*n+i];}
            za[inode-com.ns][h]=(char)best;
            pnode[(inode-com.ns)*com.npatt+h]=pbest;
            if (fabs(com.fhK[h]-fhs[h])>1e-6) error("err: fh 2?");
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
            FOR(j,n)  if(pChar1node[hp*n+j]>smallp) {
               if (com.seqtype!=CODONseq) { str[0]=pch[j]; str[1]=0; }
               else {
#ifdef CODEML
                  getcodon(str,FROM61[j]);
#endif
               }
               fprintf(fout,"%s(%5.3f) ", str,pChar1node[hp*n+j]);
            }
         }
      }     /* if (verbose) */
   }        /* for (inode) */

   zero(pMAPnode,nid);  fillxc(pMAPnodeA, 1., nid);
   for (inode=0; inode<tree.nnode-com.ns; inode++) {
      FOR (h,com.npatt) {
         pMAPnode[inode]+=com.fpatt[h]*pnode[inode*com.npatt+h]/com.ls;
         pMAPnodeA[inode]*=pow(pnode[inode*com.npatt+h], com.fpatt[h]);
      }
   }

   fputs("\n\nProb of best character at each node, listed by site ",fout);
   if (com.verbose==0) fputs("(constant sites not listed)",fout);
   if (com.ngene>1) fputs("\n\nSite (g) Freq  Data: \n",fout);
   else             fputs("\n\nSite   Freq   Data: \n",fout);

   FOR (h, com.ls) {
      hp=com.pose[h];
      fprintf (fout,"\n%4d ", h+1);
      if (com.ngene>1) {  /* which gene the site is from */
         for(ig=1;ig<com.ngene;ig++) if (hp<com.posG[ig]) break;
         fprintf(fout,"(%d)", ig);
      }
      fprintf(fout," %5.0f   ", com.fpatt[hp]);
      print1site(fout,hp);   fputs(": ",fout);
      FOR(j,nid) {
         if (com.seqtype!=CODONseq)
            fprintf(fout,"%c(%5.3f) ",pch[za[j][hp]],pnode[j*com.npatt+hp]);
         else {
#ifdef CODEML
            i = GenetCode[com.icode][ic=FROM61[za[j][hp]]] - 1;
            fprintf(fout,"%s(%1c) (%5.3f) ",
               getcodon(str,ic),AAs[i],pnode[j*com.npatt+hp]);
#endif
         }
      }
   }

   fputs("\n\n\nList of extant and reconstructed sequences\n\n",fout);
   /* extant seqs may and may not be coded. */
   for(j=0;j<tree.nnode;j++,FPN(fout)) {
      if(j<com.ns) fprintf(fout,"%-*s ", wname,com.spname[j]);
      else         fprintf(fout,"node #%-*d", wname-5,j+1);
      i=(com.cleandata||j>=com.ns);  /* coded or not */
      print1seq(fout, (j<com.ns?com.z[j]:za[j-com.ns]),com.ls,i, com.pose);
   }

   fprintf(fout,"\n\nOverall accuracy of the %d ancestral sequences:", nid);
   matout2 (fout,pMAPnode, 1, nid, 9,5);  fputs("for a site.\n",fout);
   matout2 (fout,pMAPnodeA, 1, nid, 9,5); fputs("for the sequence.\n", fout);
   ReRootTree (oldroot);
   free (pnode);  free(pChar1node);
   return (0);
}


int AncestralJoint (FILE *fout, double x[], double fhs[], double space[])
{
/* This does joint ancestral reconstruction, taking character states for 
   all nodes at a site as one entity.  The algorithm selects some candadate
   reconstructions to be evaluated, and may miss important reconstructions.
   The algorithm needs improvement.
   29 Sept 1998: Made changes so that outputs are by site rather than
   by pattern.  Calculations are duplicated if the two sites share the
   same data.

   fhs[]: fh[] for each site pattern
   pnode[nid*npatt]: prob at each node, for all patterns
   pChar1site[nid][ncode]: prob for chara at nodes, for one site pattern
   nid: Number of internal nodes
   evaluates user-specified reconstructions if com.print>=2.  
*/
   char *pch=(com.seqtype==0?BASEs:(com.seqtype==2?AAs:BINs));
   char NChara[NS-1], Chara[(NS-1)*NCODE], line[250];
   int n=com.ncode,j,h,hp,ig,ig0=-1, npathMP,nid=tree.nbranch-com.ns+1, it;
   int ipath,npath=0;
   int nodeb[NNODE], markb[NCODE], nchange,nchange0, user=(com.print>1);
   double lnL=0,y, *Ptb, ppath=0, ppathMP,pmin=.01, *pChar1site;

if(!com.cleandata) 
{ puts("\nJoint reconstruction not implemented for unclean data."); return(0);}

   if(noisy) printf("\nJoint reconstruction of ancestral sequences\n");
   fputs("\n(2) Joint reconstruction of ancestral sequences ", fout);
   fputs("(eqn. 2 in Yang et al. 1995).\nThe algorithm may miss ",fout);
   fputs("reconstructions but the probabilities are correct (see DOC).\n",fout);

   Ptb=(double*)malloc(tree.nbranch*n*n*sizeof(double));
   pChar1site=(double*)malloc(nid*com.ncode*sizeof(double));
   if (Ptb==NULL||pChar1site==NULL) error("oom");
   fputs("\nListed by site, reconstruction (prob.)\n",fout);
   fputs("parsimony reconstruction marked with *, |# paths(# changes)\n",fout);
   if (com.ngene>1) fputs("\n\nSite (g) Freq  Data: \n",fout);
   else             fputs("\n\nSite   Freq   Data: \n",fout);

   /* reconstruct sequences, joint.  One rate for all sites */
   for (h=0; h<com.ls; h++) {
      if (user) {
         printf("\nsite (-1 to exit) & reconstruction? ");
         scanf("%d",&h);  if(--h<0||h>com.ls-1) return(0);
         scanf ("%s", line);
         for(j=0;j<nid;j++) { 
            line[j]=(char)toupper(line[j]); 
            Chara[j*n]=CodeChara(line[j],com.seqtype);
         }
         npath=1; FOR(j,nid) NChara[j]=1;
      }
      hp=com.pose[h];
      for(ig=1;ig<com.ngene;ig++) if (hp<com.posG[ig]) break;
      ig--;

      if (user) {
         printf(" Site %d ", h+1);
         if (com.ngene>1) printf("is in gene (partition) %d, ", ig+1);
         printf("has freq %.0f\n ", com.fpatt[hp]);
         print1site(F0,hp);  printf(": ");
      }
      if (com.Mgene>1) SetPGene(ig, 1, 1, 0, x);
      if(ig!=ig0) {
         for (j=0,y=com.rgene[ig]; j<tree.nbranch; j++) 
#ifdef CODEML
          PMatUVRoot(Ptb+j*n*n,nodes[tree.branches[j][1]].branch*y,n,U,V,Root);
#elif BASEML
          PMatCijk(Ptb+j*n*n, nodes[tree.branches[j][1]].branch*y);
#endif
         ig0=ig;
      }
      InteriorStatesMP (0, hp, &nchange, NULL, NULL, space); 
      FOR(j,nid*n) pChar1site[j]=0;
      FOR(j,n) markb[j]=0;
      FOR(j,com.ns) { nodeb[j]=com.z[j][hp]; markb[nodeb[j]]++; }
      if (!user) {
         /* Decide on the number of pathways. */
         if (((com.seqtype==BASEseq && nid<15) || nid<7)) { /*small data set */
            if ((com.seqtype==BASEseq && nid<10) || nid<5)     /* all paths */
               FOR (j,nid) { NChara[j]=(char)n;  FOR(it,n) Chara[j*n+it]=(char)it; }
            else {
               FOR (j,nid) NChara[j]=0;
               FOR (j,n) 
                  if(markb[j])  FOR(it,nid) Chara[it*n+NChara[it]++]=(char)j;
            }
         }
         else 
            InteriorStatesMP (1, hp, &nchange, NChara, Chara, space);
         for (j=0,npath=1; j<nid; j++)  npath*=NChara[j]; 
      }
      if (noisy && !user && ((h+1)%10==0 || h+1==com.ls)) 
         printf("\nSite %3d/%3d: %3d paths evaluated", h+1,com.ls,npath);

      fprintf (fout,"\n%4d ", h+1);
      if (com.ngene>1) fprintf(fout,"(%d)", ig+1);
      fprintf(fout," %5.0f   ", com.fpatt[hp]);

      print1site(fout,hp);
      fprintf(fout, ":  ");
      for (ipath=0,ppath=ppathMP=0,npathMP=0; ipath<npath; ipath++) {
         for (j=nid-1,it=ipath; j>=0; it/=NChara[j--])
            nodeb[com.ns+j]=Chara[j*n+it%NChara[j]];
         for (j=0,nchange0=0; j<tree.nbranch; j++) 
            nchange0+=(nodeb[tree.branches[j][0]]!=nodeb[tree.branches[j][1]]);
         /* if (ipath && nchange0>nchange+4) continue;  */

         for (j=0,y=com.pi[nodeb[tree.origin]]; j<tree.nbranch; j++) 
         y*=Ptb[j*n*n+nodeb[tree.branches[j][0]]*n+nodeb[tree.branches[j][1]]];

         FOR (j,nid) pChar1site[j*n+nodeb[j+com.ns]]+=y;
         if (nchange0==nchange)  { npathMP++; ppathMP+=y; }
         if (y>ppath) ppath=y;
         if (!user && (y/fhs[hp]>pmin || nchange0==nchange)) {
            FOR (j, nid)  fprintf (fout, "%c", pch[nodeb[com.ns+j]]);
            fprintf (fout, nchange0==nchange?"* ":"  ");
            fprintf (fout, "(%.3f) ", y/fhs[hp]);
         }
         else if (user) {
            FOR(j,nid)  printf("%c",pch[nodeb[com.ns+j]]);
            printf (" prob = %.4f\n", y/fhs[hp]);
         }
      }  /* for (ipath) */
      if (!user) fprintf (fout, "|%3d (%d)", npathMP, nchange);
   }  /* for (h) */
   if (noisy)  printf("\n");
   free (Ptb); free(pChar1site);
   return (0);
}


int AncestralSeqs (FILE *fout, double x[], double space[])
{
/* Ancestral sequence reconstruction using likelihood (Yang et al. 1995).
   Marginal works with constant rate and variable rates among sites.
   Joint works only with constant rate among sites.
*/
   double lnL=0, *fhs;

   if(_nScaleF) { puts("AncestralSeqs: scaling not fixed yet."); return(-1);} 
   if (tree.nnode==com.ns) 
      { puts("\nNo ancestral nodes to reconstruct..\n");  return(0); }
   if (noisy) printf ("\nReconstructed ancestral states go into file rst.\n");
   fprintf(fout, "\nAncestral reconstruction by %sML.\n",
          (com.seqtype==0?"BASE":(com.seqtype==1?"CODON":"AA")));
   FPN(fout);  OutaTreeN(fout,1,1);  FPN(fout);  FPN(fout);
   OutaTreeN(fout,0,0);  FPN(fout);  FPN(fout);
   OutaTreeB (fout);     FPN(fout);
   fprintf (fout, "\nNodes %d to %d are ancestral\n", com.ns+1,tree.nnode);
   if ((fhs=(double*)malloc(com.npatt*sizeof(double)))==NULL) error("oom");
   
   ProbSitePattern (x, &lnL, fhs);
   if (com.alpha) {
      puts("\nRates are variable among sites.");
      puts("Only marginal reconstructions are implemented.");
   }
   if(com.seqtype==CODONseq)
      puts("\nOnly marginal reconstructions are available for codon models.");
   if(com.verbose==0) fputs("Constant sites not listed for verbose=0\n",fout);

   AncestralMarginal (fout, x, fhs);  fflush(fout);
   if (com.alpha==0 && com.seqtype!=CODONseq && tree.nnode>com.ns+1)
      AncestralJoint (fout, x, fhs, space);
   FPN(fout);
   free (fhs);
   return (0);
}

#endif

int SetnScale(int inode);
extern char *_nScale;  /* _nScale[inode]=1 if the node is used for scaling */

int SetnScale(int inode)
{
/* This marks nodes for applying scaling factors when calculating f[h].
*/
   int i,d=0;
   FOR(i,nodes[inode].nson)
      d+=(nodes[nodes[inode].sons[i]].nson?SetnScale(nodes[inode].sons[i]):1);
   if(d>30 && inode!=tree.origin)
      { _nScale[inode]=1; d=1; }

   return (d);
}

int logfx_r(double x[], int np);
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

   if (noisy) printf("\nEstimated rates for sites go into file %s\n",ratef);
   if (SetParameters(x)) puts ("par err. lfunAdG_rate");

   fprintf(fout, "\nEstimated rates for sites from %sML.\n",
          (com.seqtype==0?"BASE":(com.seqtype==1?"CODON":"AA")));
   OutaTreeN(fout,1,1); FPN(fout);
   fprintf (fout,"\nFrequencies and rates for categories (K=%d)", com.ncatG);
   matout (fout, com.freqK, 1, com.ncatG);
   matout (fout, com.rK, 1, com.ncatG);
   if (com.rho) {
      fprintf(fout, "\nTransition prob matrix over sites");
      matout2(fout, com.MK, com.ncatG, com.ncatG, 8, 4);
   }

   if((fhs=(double*)malloc(com.npatt*sizeof(double)))==NULL) error("rate oom");
   if(_nScaleF) 
      for(h=0,logfx_r(x,np); h<com.npatt; h++) {
         t=com.fhK[0*com.npatt+h];  lnL-=t*com.fpatt[h];
         for(ir=1,fh=com.freqK[ir]*(com.fhK[h]=1); ir<com.ncatG; ir++)
            fh+=com.freqK[ir]*
               (com.fhK[ir*com.npatt+h]=exp(com.fhK[ir*com.npatt+h]-t));
         fhs[h]=fh*exp(t);
      }
   else
      fx_r(x,np);
   if (com.rho==0) {     /* dG model */
      fputs("\nSite Freq  Data    ln(f)    Nexp.   Rates\n\n",fout);
      for (h=0,mre=vre=0; h<com.npatt; h++) {
         for (ir=0,fh=0,it=0,t=re=0; ir<com.ncatG; ir++) {
            if((fh1=com.freqK[ir]*com.fhK[ir*com.npatt+h])>t) 
               { t=fh1; it=ir; }
            fh+=fh1;  re+=fh1*com.rK[ir];
         }
         if(!_nScaleF) fhs[h]=fh;
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
*/
   int  nscale=1, h,il, ir, j, FPE=0;
   int  direction=-1;  /* 1: n->1;  -1: 1->n */
   double lnL=0, b1[NCATG], b2[NCATG], fh;

   NFunCall++;
   if(!_nScaleF)  
      fx_r(x,np);
   else   {  /* scaling */
      logfx_r(x,np);
      FOR(h, com.npatt) {
         lnL-=(fh=com.fhK[0*com.npatt+h])*com.fpatt[h];
         for(ir=1,com.fhK[h]=1; ir<com.ncatG; ir++) 
            com.fhK[ir*com.npatt+h]=exp(com.fhK[ir*com.npatt+h]-fh);
      }
   }

   h = (direction==1?com.ls-1:0);
   for (il=0; il<com.ls; h-=direction,il++) {
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
      if((il+1)%nscale==0) {
         fh=sum(b1,com.ncatG);
         if (fh<1e-40) {
            printf ("h,fh%6d %12.4e\n", h+1,fh);
            print1site(F0,h);  FPN(F0);  getchar();
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
/* discrete gamma rates for sites
   This calculates log{f(x_h|r_ir)} for the purpose of scaling.
*/
   int  h,ir, FPE=0;
   double lnL=0, fh=0;

   NFunCall++;
   if(_nScaleF)
      for(h=0,logfx_r(x,np); h<com.npatt; h++) {
         lnL-=(fh=com.fhK[0*com.npatt+h])*com.fpatt[h];
         for(ir=1,com.fhK[h]=1; ir<com.ncatG; ir++) 
            com.fhK[ir*com.npatt+h]=exp(com.fhK[ir*com.npatt+h]-fh);
      }
   else 
      fx_r(x,np);
   FOR(h,com.npatt) {
      if (com.fpatt[h]==0) continue;
      for(ir=0,fh=0; ir<com.ncatG;ir++) 
         fh+=com.freqK[ir]*com.fhK[ir*com.npatt+h];

      if (fh<=0 && FPE==0) {
         FPE=1;  matout(F0,x,1,np);
         printf("\nlfundG: h=%4d  fhK=%9.5f", h, fh);
         print1site(F0,h);  FPN(F0);  getchar();
      }
      lnL-=log(fh)*com.fpatt[h];

      if (com.print<0){
         fprintf(flfh,"\n%6d%6.0f%16.10f%9.2f ", 
            h+1,com.fpatt[h],log(fh),com.ls*fh);
         print1site(flfh,h);
      }
   }
   return (lnL);
}



int logfx_r(double x[], int np)
{
/* This calculates log{f(x|r)} where f is the probability of observing 
   data x at a site, given the rate r for the site.
   The results are stored in com.fhK[com.ncatG*com.npatt].
   This deals with underflows with large trees and uses global variables 
   _nScale and _nScaleF[com.npatt].
*/
   int  h, ir, i, ig, FPE=0;
   double fh;

   if (SetParameters(x)) puts ("\npar err..");
   for (ig=0; ig<com.ngene; ig++) { /* alpha may differ over ig */
      SetPGene(ig, com.Mgene>1, com.Mgene>1, com.nalpha>1, x);
      for (ir=0; ir<com.ncatG; ir++) {
         FOR(i,tree.nnode)  nodes[i].branch*=com.rK[ir];
         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) _nScaleF[h]=0;
         PartialLikelihood (tree.origin, ig);
         FOR(i,tree.nnode)  nodes[i].branch/=com.rK[ir];
         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            if (com.fpatt[h]==0) continue;
            for (i=0,fh=0; i<com.ncode; i++)
               fh += com.pi[i]*nodes[tree.origin].lkl[h*com.ncode+i];

            if (FPE==0 && fh<0) {
               FPE=1;  matout(F0,x,1,np);
               printf ("\nlogfx_r: h=%3d  fhK=%9.5e\nData: ", h+1,fh);
               print1site(F0,h);  FPN(F0);  getchar();
            }

            com.fhK[ir*com.npatt+h] = _nScaleF[h]+log(fh);
/*
            printf("\nir & h:%4d%4d %9.5f ",ir+1,h+1,com.fhK[ir*com.npatt+h]);
*/
         }
      }
   }  /* for(ig) */
   return(0);
}


int fx_r(double x[], int np)
{
/* This calculates f(x|r), the probability of observing data x at a site, 
   given the rate r for the site.
   The results are stored in com.fhK[com.ncatG*com.npatt].
   This does not deal with underflows with large trees
*/
   int  h, ir, i, ig, FPE=0;
   double fh;

   if (SetParameters (x)) puts ("\npar err..");
   zero (com.fhK, com.npatt);
   for (ig=0; ig<com.ngene; ig++) {
      SetPGene(ig, com.Mgene>1, com.Mgene>1, com.nalpha>1, x);
      for (ir=0; ir<com.ncatG; ir++) {
         FOR(i,tree.nnode) nodes[i].branch*=com.rK[ir];
         PartialLikelihood (tree.origin, ig);
         FOR(i,tree.nnode) nodes[i].branch /= com.rK[ir];
         for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
            if (com.fpatt[h]==0) continue;
            for (i=0,fh=0; i<com.ncode; i++)
               fh += com.pi[i]*nodes[tree.origin].lkl[h*com.ncode+i];
            com.fhK[ir*com.npatt+h] = max(fh,0);

            if (FPE==0 && fh<-1e-10) {
               FPE=1; matout(F0,x,1,np);
               printf ("\nlogfx_r: h=%3d ir=%2d fhK=%9.5e\nData: ", h+1,ir,fh);
               print1site(F0,h);  FPN(F0);  getchar();
            }
         }
      }
   }
   return(0);
}


double lfun (double x[], int np)
{
   int  h, i, ig,FPE=0;
   double lnL=0, fh;

   NFunCall++;
   if(_nScaleF) FOR(h,com.npatt) _nScaleF[h]=0;
   if (SetParameters(x)) puts ("\npar err..");
   for (ig=0; ig<com.ngene; ig++) {
      if (com.Mgene>1) SetPGene(ig, 1, 1, 0, x);
      PartialLikelihood (tree.origin, ig);
      for (h=com.posG[ig]; h<com.posG[ig+1]; h++) {
         if (com.fpatt[h]==0) continue;
         for (i=0,fh=0; i<com.ncode; i++) 
            fh += com.pi[i]*nodes[tree.origin].lkl[h*com.ncode+i];
         if (fh<=0 && FPE==0) {
            FPE=1; matout(F0,x,1,np);
            printf("\alfun: h=%4d  fh=%9.5e (hit Enter) \nData: ", h+1,fh);
            print1site(F0,h);  FPN(F0);  getchar();
         }
         fh = log(fh) + (_nScaleF?_nScaleF[h]:0);
         lnL-=fh*com.fpatt[h];
         if (com.print<0) {
            fprintf(flfh,"\n%6d%6.0f%16.10f%9.2f ",
               h+1,com.fpatt[h],fh,com.ls*exp(fh));
            print1site(flfh,h);
         }
      }  /* for(h) */
   }     /* for(ig) */
   return (lnL);
}



int print1site (FILE*fout, int h)
{
/* This print out one site in the data com.z[].  It may be the h-th 
   site in the original data file or the h-th pattern.
   The data are coded if com.cleandata==1 or not if otherwise.
*/
   char str[4]="";
   char *pch=(com.seqtype==0?BASEs:(com.seqtype==2?AAs:BINs));
   int i, ib,nb=(com.seqtype==CODONseq?3:1), iaa,ic;

   if(!com.cleandata) 
      FOR(i,com.ns) {
         FOR(ib,nb) fprintf(fout,"%c", com.z[i][h*nb+ib]);
#ifdef CODEML
         if (nb==3)  {
            Codon2AA(com.z[i]+h*nb, str, com.icode, &iaa);
            fprintf(fout," (%c) ",AAs[iaa]);
		    }
#endif
      }
   else {
      FOR(i,com.ns) {
         if(nb==1) fprintf(fout,"%c", pch[com.z[i][h]]);
#ifdef CODEML
         else {
           ic=FROM61[com.z[i][h]];  iaa=GenetCode[com.icode][ic]-1;
           fprintf(fout,"%s (%c) ", getcodon(str,ic), AAs[iaa]);
	 }
#endif
      }
   }
   return(0);
}
   

int InitPartialLikelihood (void)
{
/* set partial likelihood at tips of the tree, considering missing data.
   This need to be modified when more proper algorithm for dealing with
   missing data is worked out.
   Need testing if sequences in the data are ancestors.
*/
   int n=com.ncode, is,j,k,h;
   int i, nb[3],ib[3][4], i0,i1,i2,ic,nsense;
   char *pch=(com.seqtype==0?BASEs:(com.seqtype==2?AAs:BINs));

   FOR(is,com.ns) nodes[is].lkl=com.chunk+n*com.npatt*is;
   for (is=0;is<com.ns;is++) {
      zero(nodes[is].lkl, com.npatt*n);
      for (h=0;h<com.npatt;h++) {
         if(com.seqtype==CODONseq) {
#ifdef CODEML
            for (h=0;h<com.npatt;h++) {
               FOR(i,3)  NucListall(com.z[is][h*3+i],&nb[i],ib[i]);
               for(i0=0,nsense=0; i0<nb[0]; i0++) FOR(i1,nb[1]) FOR(i2,nb[2]) {
                  ic=ib[0][i0]*16+ib[1][i1]*4+ib[2][i2];         
                  ic=FROM64[ic];
                  if(ic>-1) { nsense++;  nodes[is].lkl[h*n+ic]=1; }
               }
            }
            if(nsense==0) error("Initlkl: stop codon");
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
      if (i!=tree.origin) 
         nodes[i].branch=nodes[nodes[i].father].divtime-nodes[i].divtime;
   if (!rooted) {
      it=nodes[tree.origin].sons[2];
      nodes[it].branch =
      2*nodes[2*com.ns-2].divtime-nodes[tree.origin].divtime-nodes[it].divtime;
   }
}

#endif




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
   tree.origin=it;  tree.nnode=ns*2-1-!rooted;
   NodeToBranch();
   return (0);
}


#ifdef EVOLVE

static double PMat[16], Qfactor, kappa1, kappa2;

int GenerateSeq (void)
{
/* makes sequence com.z[tree.origin] with com.pi[] and evolves it
   along the tree, using nodes[].branch, com.ns, com.ls, com.ncode, com.model,
   com.kappa, com.alpha, com.rates (if com.alpha>0)
   z[ns*2-1] has 0,1,2,3 first and decoded before printing out.
*/
   int i,h;
   double T=com.pi[0],C=com.pi[1],A=com.pi[2],G=com.pi[3],Y=T+C,R=A+G;

   if (com.ncode!=4) error ("ncode");
   if (com.ns<2 || com.ls<1) error ("err GenerateSeq()");

   if (com.model==F84) { kappa1=1+com.kappa/Y; kappa2=1+com.kappa/R; }
   else                  kappa1=kappa2=com.kappa;
   Qfactor=1/(2*T*C*kappa1+2*A*G*kappa2 + 2*Y*R);
   if (com.model>K80) 
      dnamaker (com.z[tree.origin], com.ls, com.pi);
   else
      FOR(i,com.ls) com.z[tree.origin][i]=(char)(rndu()*4.);

   if (com.model==JC69) EvolveJC (tree.origin);
   else                 Evolve   (tree.origin);
   com.npatt=0;
   FOR(i,com.ns) FOR(h,com.ls) {
      if(com.z[i][h]<0 || com.z[i][h]>3) error("err GenerateSeq.");
      com.z[i][h]=BASEs[com.z[i][h]];
   }
   /* Ancestral sequences are not encoded here. */
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
            if (j-ib) com.z[ison][h]=j;
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
      if (nodes[ison].nson) EvolveJC (ison); 
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

puts("\nPatternWeightSimple: missing data not dealt with yet");

   FOR (j,com.ns) zt[j]=(char*)space+j*com.ls;
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
