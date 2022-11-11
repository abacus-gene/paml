/* treesub.c
 subroutines that operates on trees, inserted into other programs
 such as baseml, basemlg, codeml, and pamp.
 */

extern char BASEs[], * EquateBASE[], AAs[], BINs[], CODONs[][4], nChara[], CharaMap[][64];

extern int noisy;

#ifdef  BASEML
#define REALSEQUENCE
#define NODESTRUCTURE
#define TREESEARCH
#define LSDISTANCE
#define LFUNCTIONS
#define RECONSTRUCTION
#define MINIMIZATION
#endif

#ifdef  CODEML
#define REALSEQUENCE
#define NODESTRUCTURE
#define TREESEARCH
#define LSDISTANCE
#define LFUNCTIONS
#define RECONSTRUCTION
#define MINIMIZATION
#endif

#ifdef  BASEMLG
#define REALSEQUENCE
#define NODESTRUCTURE
#define LSDISTANCE
#endif

#ifdef  RECONSTRUCTION
#define PARSIMONY
#endif

#ifdef  MCMCTREE
#define REALSEQUENCE
#define NODESTRUCTURE
#define LFUNCTIONS
#endif

#if(defined CODEML || defined YN00)
double SS, NN, Sd, Nd;  /* kostas, # of syn. sites,# of non syn. sites,# of syn. subst.,# of non syn. subst. */
#endif



#if(defined BASEML || defined BASEMLG || defined CODEML || defined BPP || defined EVOLVER || defined MCMCTREE || defined THREES)


int PatternWeightSimple(void)
{
   /* This is modified from PatternWeight(), and does not deal with multiple genes in
      the same alignment (com.ngene, com.lgene[], com.posG[], etc.)
      Binary search is used to sort site patterns, with patterns represented using 0-ended strings.
      The routine works with nucleotide, amino acid, or codon sequences.
      This should work with both encoded and un-encoded sequences.
      If com.pose is not NULL, this generates the site-to-pattern map in com.pose[].
      com.fpatt holds site-pattern counts.
      Sequences z[ns][ls] are copied into patterns zt[ls*lpatt], and bsearch is used
      twice to avoid excessive copying, first to count npatt and identify the site patterns and
      second to generate fpatt[], pose[] etc.
   */
   int maxnpatt = com.ls, h, l, u, ip, j, k, same;
   int n31 = (com.seqtype == CODONseq ? 3 : 1);
   int lpatt = com.ns * n31 + 1;
   int* p2s;  /* point patterns to sites in zt */
   char timestr[36];
   char* p, * zt;
   double nc = (com.seqtype == 1 ? 64 : com.ncode) + !com.cleandata + 1;
   int debug = 0;

   /* (A) Collect and sort patterns.  Get com.npatt.
          Move sequences com.z[ns][ls] into sites zt[ls*lpatt].
          Use p2s to map patterns to sites in zt to avoid copying.
    */

   if ((com.seqtype == 1 && com.ns < 5) || (com.seqtype != 1 && com.ns < 7))
      maxnpatt = (int)(pow(nc, (double)com.ns) + 0.5);
   if (maxnpatt > com.ls)   maxnpatt = com.ls;
   p2s = (int*)malloc(maxnpatt * sizeof(int));
   zt = (char*)malloc(com.ls * lpatt * sizeof(char));
   if (p2s == NULL || zt == NULL)  error2("oom p2s or zt");
   memset(zt, 0, com.ls * lpatt * sizeof(char));
   for (j = 0; j < com.ns; j++)
      for (h = 0; h < com.ls; h++)
         for (k = 0; k < n31; k++)
            zt[h * lpatt + j * n31 + k] = (char)(com.z[j][h * n31 + k] + 1);

   com.npatt = l = u = ip = 0;
   for (h = 0; h < com.ls; h++) {
      if (debug) printf("\nh %3d %s", h, zt + h * lpatt);
      /* bsearch in existing patterns.  Knuth 1998 Vol3 Ed2 p.410
         ip is the loc for match or insertion.  [l,u] is the search interval.
       */
      same = 0;
      if (h != 0) {  /* not 1st pattern? */
         for (l = 0, u = com.npatt - 1;;) {
            if (u < l) break;
            ip = (l + u) / 2;
            k = strcmp(zt + h * lpatt, zt + p2s[ip] * lpatt);
            if (k < 0)        u = ip - 1;
            else if (k > 0)   l = ip + 1;
            else { same = 1;  break; }
         }
      }
      if (!same) {
         if (com.npatt > maxnpatt)
            error2("npatt > maxnpatt");
         if (l > ip) ip++;        /* last comparison in bsearch had k > 0. */
         /* Insert new pattern at ip.  This is the expensive step. */

         if (ip < com.npatt)
            memmove(p2s + ip + 1, p2s + ip, (com.npatt - ip) * sizeof(int));
         p2s[ip] = h;
         com.npatt++;
      }

      if (debug) {
         printf(": %3d (%c ilu %3d%3d%3d) ", com.npatt, (same ? '.' : '+'), ip, l, u);
         for (j = 0; j < com.npatt; j++)
            printf(" %s", zt + p2s[j] * lpatt);
      }
      if (noisy > 2 && ((h + 1) % 10000 == 0 || h + 1 == com.ls))
         printf("\r%12d patterns at %8d / %8d sites (%.1f%%), %s",
            com.npatt, h + 1, com.ls, (h + 1.) * 100 / com.ls, printtime(timestr));
   }     /* for (h)  */
   if (noisy > 2) printf("\n%d patterns\n", com.npatt);

   /* (B) count pattern frequencies */
   com.fpatt = (double*)realloc(com.fpatt, com.npatt * sizeof(double));
   if (com.fpatt == NULL) error2("oom fpatt");
   memset(com.fpatt, 0, com.npatt * sizeof(double));

   for (h = 0; h < com.ls; h++) {
      for (same = 0, l = 0, u = com.npatt - 1;;) {
         if (u < l) break;
         ip = (l + u) / 2;
         k = strcmp(zt + h * lpatt, zt + p2s[ip] * lpatt);
         if (k < 0)        u = ip - 1;
         else if (k > 0)   l = ip + 1;
         else { same = 1;  break; }
      }
      com.fpatt[ip]++;
      if (com.pose) com.pose[h] = ip;
   }     /* for (h)  */
   for (j = 0; j < com.ns; j++) {
      com.z[j] = p = (char*)realloc(com.z[j], com.npatt * sizeof(char));
      for (ip = 0; ip < com.npatt; ip++)
         for (k = 0; k < n31; k++)
            *p++ = (char)(zt[p2s[ip] * lpatt + j * n31 + k] - 1);
   }
   free(p2s);  free(zt);

   return (0);
}

int ConvertSiteJC69like(char* z[], int ns, int h, char zh[])
{
   /* This converts a site (or pattern) of nucleotides or amino acids for JC69-like models.
      Sequence alignments in com.z[] are already encoded into 0, 1, ...
      If the data have no ambiguities (com.cleandata=1), the routine converts,
      for example, a site 1120 (CCAT) into 0012 (TTCA) before checking against old
      patterns already found.  If a site contain non-gap ambiguities, it is not
      converted.  For every site, the routine changes ? or N into - first, and then
      convert the site iff there are no non-gap ambiguities.  Thus CC?T will be
      converted into CC-T and then into TT-C.  A site with CCRT will not be
      convertd.  In theory such sites may be compressed as well, but the effort is
      perhaps not worthwhile.
   */
   char b, gap;
   char* pch = (com.seqtype == 0 ? BASEs : (com.seqtype == 2 ? AAs : BINs));
   int j, k, convert;

   gap = (char)(strchr(pch, (int)'-') - pch);

   if (com.cleandata) { /* clean data, always convert */
      zh[0] = b = 0;
      b++;
      for (j = 1; j < com.ns; j++) {
         for (k = 0; k < j; k++)
            if (z[j][h] == z[k][h]) break;
         zh[j] = (k < j ? zh[k] : b++);
      }
   }
   else { /* convert only if there are no non-gap ambiguity characters */
      for (j = 0; j < ns; j++)
         zh[j] = z[j][h];

      /* After this loop, convert = 0 or 1 decides whether to convert. */
      for (j = 0, convert = 1; j < ns; j++) {
         if (zh[j] < com.ncode)
            continue;
         if (nChara[(int)zh[j]] == com.ncode) {
            zh[j] = gap;
            continue;
         }
         convert = 0;
         break;
      }
      if (convert) {
         b = 0;
         if (zh[0] != gap)
            zh[0] = b++;
         for (j = 1; j < ns; j++) {
            if (zh[j] != gap) {
               for (k = 0; k < j; k++)
                  if (zh[j] == z[k][h]) break;
               if (k < j) zh[j] = zh[k];
               else    zh[j] = b++;
            }
         }
      }
   }
   for (j = 0; j < ns; j++)  zh[j] ++;  /* change 0 to 1. */
   return(0);
}


int PatternWeightJC69like(void)
{
   /*  This collaps site patterns further for JC69-like models, called after PatternWeight().
       This is used for JC and poisson amino acid models.  The routine could be merged into
       PatternWeight(), which should lead to faster computation, but this is not done because
       right now InitializeBaseAA() prints out base or amino acid frequencies after PatternWeight()
       and before this routine.

       If com.pose is not NULL, the routine also updates com.pose.  This allows the program
       to work if com.readpattern==1.
       This works for nucleotide and amino acid models, but not codon models.
       This routine is nearly identical to PatternWeight.
       fpatt0 stores the old com.fpatt info, while com.fpatt is shrunk.
       Think about merging them (encode sequences first and compress sites).
   */
   int npatt0 = com.npatt, lpatt = com.ns + 1, h, l, u, ip, j, k, same;
   int* p2s;  /* point patterns to sites in zt */
   int* pose0 = NULL;
   char timestr[36];
   char* p, * zt;
   double* fpatt0;
   int debug = 0;

   /* (A) Collect and sort patterns.  Get com.npatt.
          Move sequences com.z[ns][ls] into sites zt[ls*lpatt].
          Use p2s to map patterns to sites in zt to avoid copying.
   */
   if (noisy > 2) printf("Counting site patterns again, for JC69.. %s\n", printtime(timestr));
   if (com.seqtype == 1) error2("PatternWeightJC69like does not work for codon seqs");
   if (com.ngene > 1) error2("PatternWeightJC69like does not work when ngene > 1");

   p2s = (int*)malloc(npatt0 * sizeof(int));
   zt = (char*)malloc(npatt0 * lpatt * sizeof(char));
   fpatt0 = (double*)malloc(npatt0 * sizeof(double));
   if (p2s == NULL || zt == NULL || fpatt0 == NULL)  error2("oom p2s or zt or fpatt0");
   memset(zt, 0, npatt0 * lpatt * sizeof(char));
   memmove(fpatt0, com.fpatt, npatt0 * sizeof(double));
   if (com.pose && com.ls > npatt0) {
      if ((pose0 = (int*)malloc(com.ls * sizeof(int))) == NULL)  error2("oom pose0");
      memmove(pose0, com.pose, com.ls * sizeof(int));
   }
   for (h = 0; h < npatt0; h++)
      ConvertSiteJC69like(com.z, com.ns, h, zt + h * lpatt);

   l = u = ip = com.npatt = 0;
   for (h = 0; h < npatt0; h++) {
      if (debug) printf("\nh %3d %s", h, zt + h * lpatt);

      /* bsearch in existing patterns.  Knuth 1998 Vol3 Ed2 p.410
      ip is the loc for match or insertion.  [l,u] is the search interval.
      */
      same = 0;
      if (h != 0) {  /* not 1st pattern? */
         for (l = 0, u = com.npatt - 1; ; ) {
            if (u < l) break;
            ip = (l + u) / 2;
            k = strcmp(zt + h * lpatt, zt + p2s[ip] * lpatt);
            if (k < 0)        u = ip - 1;
            else if (k > 0)   l = ip + 1;
            else { same = 1;  break; }
         }
      }
      if (!same) {
         if (l > ip) ip++;        /* last comparison in bsearch had k > 0. */
         /* Insert new pattern at ip.  This is the expensive step. */
         if (ip < com.npatt)
            memmove(p2s + ip + 1, p2s + ip, (com.npatt - ip) * sizeof(int));
         p2s[ip] = h;
         com.npatt++;
      }
      if (debug) {
         printf(": %3d (%c ilu %3d%3d%3d) ", com.npatt, (same ? '.' : '+'), ip, l, u);
         for (j = 0; j < com.npatt; j++)
            printf(" %s", zt + p2s[j] * lpatt);
      }
      if (noisy > 2 && ((h + 1) % 10000 == 0 || h + 1 == npatt0))
         printf("\rCompressing, %6d patterns at %6d / %6d sites (%.1f%%), %s",
            com.npatt, h + 1, npatt0, (h + 1.) * 100 / npatt0, printtime(timestr));
   }     /* for (h)  */
   if (noisy > 2) printf("\n");

   /* (B) count pattern frequencies and collect pose[] */
   com.fpatt = (double*)realloc(com.fpatt, com.npatt * sizeof(double));
   memset(com.fpatt, 0, com.npatt * sizeof(double));

   for (h = 0; h < npatt0; h++) {
      for (same = 0, l = 0, u = com.npatt - 1; ; ) {
         if (u < l) break;
         ip = (l + u) / 2;
         k = strcmp(zt + h * lpatt, zt + p2s[ip] * lpatt);
         if (k < 0)        u = ip - 1;
         else if (k > 0)   l = ip + 1;
         else { same = 1;  break; }
      }
      com.fpatt[ip] += fpatt0[h];
      if (com.pose) {
         if (com.ls < npatt0)
            com.pose[h] = ip;
         else
            for (j = 0; j < com.ls; j++)
               if (pose0[j] == h)  /* In original alignment, site j has old pattern h. */
                  com.pose[j] = ip;
      }
      if (noisy > 2 && ((h + 1) % 10000 == 0 || h + 1 == npatt0))
         printf("\rCollecting fpatt[] & pose[], %6d patterns at %6d / %6d sites (%.1f%%), %s",
            com.npatt, h + 1, npatt0, (h + 1.) * 100 / npatt0, printtime(timestr));
   }     /* for (h)  */
   if (noisy > 2) printf("\n");

   for (j = 0; j < com.ns; j++) {
      com.z[j] = (char*)realloc(com.z[j], com.npatt * sizeof(char));
      for (ip = 0, p = com.z[j]; ip < com.npatt; ip++)
         *p++ = (char)(zt[p2s[ip] * lpatt + j] - 1);
   }
   free(p2s);  free(zt);  free(fpatt0);
   if (pose0) free(pose0);
   if (com.ngene > 1) error2("rewrite PatternWeightJC69like to deal with partitioned data");
   com.posG[1] = com.npatt;
   return (0);
}


#endif




#ifdef REALSEQUENCE

int hasbase(char* str)
{
   char* p = str, * eqdel = ".-?";
   while (*p)
      if (*p == eqdel[0] || *p == eqdel[1] || *p == eqdel[2] || isalpha(*p++))
         return(1);
   return(0);
}


int GetSeqFileType(FILE* fseq, int* NEXUSseq);
int IdenticalSeqs(void);
void RemoveEmptySequences(void);

int GetSeqFileType(FILE* fseq, int* format)
{
   /* NEXUSstart="begin data" and NEXUSdata="matrix" identify nexus file format.
    Modify if necessary.
    format: 0: alignment; 1: fasta; 2: nexus.
    */
   int  lline = 1000, ch, aligned;
   char fastastarter = '>';
   char line[1000], * NEXUSstart = "begin data", * NEXUSdata = "matrix", * p;
   char* ntax = "ntax", * nchar = "nchar";

   while (isspace(ch = fgetc(fseq)))
      ;
   ungetc(ch, fseq);
   if (ch == fastastarter) {
      *format = 1;
      ScanFastaFile(fseq, &com.ns, &com.ls, &aligned);
      if (aligned)
         return(0);
      else
         error2("The seq file appears to be in fasta format, but not aligned?");
   }
   if (fscanf(fseq, "%d%d", &com.ns, &com.ls) == 2) {
      *format = 0; return(0);
   }
   *format = 2;
   printf("\nseq file is not paml/phylip format.  Trying nexus format.");

   for (; ; ) {
      if (fgets(line, lline, fseq) == NULL) error2("seq err1: EOF");
      strcase(line, 0);
      if (strstr(line, NEXUSstart)) break;
   }
   for (; ; ) {
      if (fgets(line, lline, fseq) == NULL) error2("seq err2: EOF");
      strcase(line, 0);
      if ((p = strstr(line, ntax)) != NULL) {
         while (*p != '=') { if (*p == 0) error2("seq err"); p++; }
         sscanf(p + 1, "%d", &com.ns);
         if ((p = strstr(line, nchar)) == NULL) error2("expect nchar");
         while (*p != '=') { if (*p == 0) error2("expect ="); p++; }
         sscanf(p + 1, "%d", &com.ls);
         break;
      }
   }
   /* printf("\nns: %d\tls: %d\n", com.ns, com.ls);  */
   for (; ; ) {
      if (fgets(line, lline, fseq) == NULL) error2("seq err1: EOF");
      strcase(line, 0);
      if (strstr(line, NEXUSdata)) break;
   }
   return(0);
}

int PopupNEXUSComment(FILE* fseq)
{
   int ch, comment1 = ']';
   for (; ; ) {
      ch = fgetc(fseq);
      if (ch == EOF) error2("expecting ]");
      if (ch == comment1) break;
      if (noisy) putchar(ch);
   }
   return(0);
}


#if(MCMCTREE)

int ReadMorphology(FILE* fout, FILE* fin, int locus)
{
   int i, j;

   if ((data.zmorph[locus][0] = (double*)malloc((com.ns * 2 - 1) * com.ls * sizeof(double))) == NULL)
      error2("oom zmorph");
   if ((data.Rmorph[locus] = (double*)malloc(com.ls * com.ls * sizeof(double))) == NULL)
      error2("oom Rmorph");

   printf("Locus %d has morphological alignment\n", locus + 1);
   for (i = 1; i < com.ns * 2 - 1; i++) {
      data.zmorph[locus][i] = data.zmorph[locus][0] + i * com.ls;
   }
   for (i = 0; i < com.ns; i++) {
      fscanf(fin, "%s", com.spname[i]);
      printf("Reading data for species #%2d: %s     \r", i + 1, com.spname[i]);
      for (j = 0; j < com.ls; j++)
         fscanf(fin, "%lf", &data.zmorph[locus][i][j]);
   }

   for (i = 0; i < com.ns; i++) {
      fprintf(fout, "%-10s ", com.spname[i]);
      for (j = 0; j < com.ls; j++)
         fprintf(fout, " %8.5f", data.zmorph[locus][i][j]);
      fprintf(fout, "\n");
   }

#if(0)
   char line[1024], str[64];
   fscanf(fin, "%s", str);
   fgets(line, 1024, fin);
   i = j = -1;
   if (strstr("Correlation", str)) {
      for (i = 0; i < com.ls; i++) {
         for (j = 0; j < com.ls; j++)
            if (fscanf(fin, "%lf", &data.Rmorph[locus][i * com.ls + j]) != 1) break;
         if (j < com.ls) break;
      }
   }
   if (i != com.ls || j != com.ls) {
      printf("\ndid not find a good R matrix.  Setting it to identity matrix I.\n");
      for (i = 0; i < com.ls; i++)
         for (j = 0; j < com.ls; j++)
            data.Rmorph[locus][i * com.ls + j] = (i == j);
   }
#endif
   return(0);
}

#endif

int ReadSeq(FILE* fout, FILE* fseq, int cleandata, int locus, int read_seq_only)
{
   /* read in sequence, translate into protein (CODON2AAseq), and
    This counts ngene but does not initialize lgene[].
    It also codes (transforms) the sequences.
    com.seqtype: 0=nucleotides; 1=codons; 2:AAs; 3:CODON2AAs; 4:BINs
    com.pose[] is used to store gene or site-partition labels.
    ls/3 gene marks for codon sequences.
    char opt_c[]="GIPM";
    G:many genes;  I:interlaved format;  P:patterns;  M:morphological characters

    Use cleandata=1 to clean up ambiguities.  In return, com.cleandata=1 if the
    data are clean or are cleaned, and com.cleandata=0 is the data are unclean.
    */
   char* p, * p1, eq = '.', comment0 = '[', * line;
   int format = 0;  /* 0: paml/phylip, 1: fasta; 2: NEXUS/nexus */
   int i, j, k, ch, noptline = 0, lspname = LSPNAME, miss = 0, nb;
   int lline = 10000, lt[NS], igroup, Sequential = 1, basecoding = 0;
   int n31 = (com.seqtype == CODONseq || com.seqtype == CODON2AAseq ? 3 : 1);
   int gap = (n31 == 3 ? 3 : 10), nchar = (com.seqtype == AAseq ? 20 : 4);
   int h;
   char* pch = ((com.seqtype <= 1 || com.seqtype == CODON2AAseq) ? BASEs : (com.seqtype == 2 ? AAs : BINs));
   char* NEXUSend = "end;";
   double lst;

#if(MCMCTREE)
   data.datatype[locus] = com.seqtype;
#endif
   h = -1; /* avoid warning */
   com.readpattern = 0;
   if (com.seqtype == 4) error2("seqtype==BINs, check with author");
   if (noisy >= 9 && (com.seqtype <= CODONseq || com.seqtype == CODON2AAseq)) {
      puts("\n\nAmbiguity character definition table:\n");
      for (i = 0; i < (int)strlen(BASEs); i++) {
         nb = (int)strlen(EquateBASE[i]);
         printf("%c (%d): ", BASEs[i], nb);
         for (j = 0; j < nb; j++)  printf("%c ", EquateBASE[i][j]);
         printf("\n");
      }
   }
   GetSeqFileType(fseq, &format);

   if (com.ns > NS) error2("too many sequences.. raise NS?");
   if (com.ls % n31 != 0) {
      printf("\n%d nucleotides, not a multiple of 3!", com.ls); exit(-1);
   }
   if (noisy) printf("ns = %d  \tls = %d\n", com.ns, com.ls);

   for (j = 0; j < com.ns; j++) {
      if (com.spname[j]) free(com.spname[j]);
      com.spname[j] = (char*)malloc((lspname + 1) * sizeof(char));
      for (i = 0; i < lspname + 1; i++) com.spname[j][i] = 0;
      if ((com.z[j] = (char*)realloc(com.z[j], com.ls * sizeof(char))) == NULL)
         error2("oom z");
   }
   com.rgene[0] = 1;   com.ngene = 1;
   lline = max2(lline, com.ls / n31 * (n31 + 1) + lspname + 50);
   if ((line = (char*)malloc(lline * sizeof(char))) == NULL) error2("oom line");

   /* first line */
   if (format == 0) {
      if (!fgets(line, lline, fseq)) error2("ReadSeq: first line");
      com.readpattern = (strchr(line, 'P') || strchr(line, 'p'));
#if(MCMCTREE)
      char tmp[32];
      if (strchr(line, 'M') || strchr(line, 'm')) {
         data.datatype[locus] = MORPHC;
         data.zpopvar[locus] = 0; data.ldetRm[locus] = 0;                /* MdR */
         sscanf(line, "%s %lf %lf", tmp, &data.zpopvar[locus], &data.ldetRm[locus]);  /* MdR */
         printf("MdR: got %f %f, locus: %d\n", data.zpopvar[locus], data.ldetRm[locus], locus);
         /* TODO: Do some error checking. MdR. */
      }
#endif
   }
#if(MCMCTREE)
   if (data.datatype[locus] == MORPHC) { /* morhpological data */
      ReadMorphology(fout, fseq, locus);
      return(0);
   }
   else
#endif
      if (!com.readpattern) {
         if (com.pose) free(com.pose);
         if ((com.pose = (int*)malloc(com.ls / n31 * sizeof(int))) == NULL)
            error2("oom pose");
         for (j = 0; j < com.ls / n31; j++) com.pose[j] = 0;      /* gene #1, default */
      }
      else {
         if (com.pose) free(com.pose);
         com.pose = NULL;
      }
   if (format)
      goto readseq;

   for (j = 0; j < lline && line[j] && line[j] != '\n'; j++) {
      if (!isalnum(line[j])) continue;
      line[j] = (char)toupper(line[j]);
      switch (line[j]) {
      case 'G': noptline++;   break;
      case 'C': basecoding = 1; break;
      case 'S': Sequential = 1; break;
      case 'I': Sequential = 0; break;
      case 'P':               break;  /* already dealt with. */
      default:
         printf("\nBad option '%c' in first line of seqfile\n", line[j]);
         exit(-1);
      }
   }
   if (strchr(line, 'C')) {   /* protein-coding DNA sequences */
      if (com.seqtype == 2) error2("option C?");
      if (com.seqtype == 0) {
         if (com.ls % 3 != 0 || noptline < 1)  error2("option C?");
         com.ngene = 3;
         for (i = 0; i < 3; i++) com.lgene[i] = com.ls / 3;
#if(defined(BASEML) || defined(BASEMLG))
         com.coding = 1;
         if (com.readpattern)
            error2("partterns for coding sequences (G C P) not implemented.");
         else
            for (i = 0; i < com.ls; i++) com.pose[i] = (char)(i % 3);

#endif
      }
      noptline--;
   }

   /* option lines */
   for (j = 0; j < noptline; j++) {
      for (ch = 0; ; ) {
         ch = (char)fgetc(fseq);
         if (ch == comment0)
            PopupNEXUSComment(fseq);
         if (isalnum(ch)) break;
      }

      ch = (char)toupper(ch);
      switch (ch) {
      case ('G'):
         if (basecoding) error2("sequence data file: incorrect option format, use GC?\n");
#if(defined MCMCTREE || defined BPP)
         error2("sequence data file: option G should not be used.");
#endif
         if (fscanf(fseq, "%d", &com.ngene) != 1) error2("expecting #gene here..");
         if (com.ngene > NGENE) error2("raise NGENE?");

         fgets(line, lline, fseq);
         if (!blankline(line)) {    /* #sites in genes on the 2nd line */
            for (i = 0, p = line; i < com.ngene; i++) {
               while (*p && !isalnum(*p)) p++;
               if (sscanf(p, "%d", &com.lgene[i]) != 1) break;
               while (*p && isalnum(*p)) p++;
            }
            /* if ngene is large and some lgene is on the next line */
            for (; i < com.ngene; i++)
               if (fscanf(fseq, "%d", &com.lgene[i]) != 1) error2("EOF at lgene");

            for (i = 0, k = 0; i < com.ngene; i++)
               k += com.lgene[i];
            if (k != com.ls / n31) {
               matIout(F0, com.lgene, 1, com.ngene);
               printf("\n%6d != %d", com.ls / n31, k);
               puts("\nOption G: total length over genes is not correct");
               if (com.seqtype == 1) {
                  puts("Note: gene length is in number of codons.");
               }
               puts("Sequence length in number of nucleotides.");
               exit(-1);
            }
            if (!com.readpattern)
               for (i = 0, k = 0; i < com.ngene; k += com.lgene[i], i++)
                  for (j = 0; j < com.lgene[i]; j++)
                     com.pose[k + j] = i;

         }
         else {                   /* site marks on later line(s)  */
            if (com.readpattern)
               error2("option PG: use number of patterns in each gene and not site marks");
            for (k = 0; k < com.ls / n31; ) {
               if (com.ngene > 9)  fscanf(fseq, "%d", &ch);
               else {
                  do ch = fgetc(fseq); while (!isdigit(ch));
                  ch = ch - (int)'1' + 1;  /* assumes 1,2,...,9 are consecutive */
               }
               if (ch<1 || ch>com.ngene)
               {
                  printf("\ngene mark %d at %d?\n", ch, k + 1);  exit(-1);
               }
               com.pose[k++] = ch - 1;
            }
            if (!fgets(line, lline, fseq)) error2("sequence file, gene marks");
         }
         break;
      default:
         printf("Bad option '%c' in option lines in seqfile\n", line[0]);
         exit(-1);
      }
   }

readseq:
   /* read sequence */
   if (Sequential) {    /* sequential */
      if (noisy) printf("Reading sequences, sequential format..\n");
      for (j = 0; j < com.ns; j++) {
         lspname = LSPNAME;
         for (i = 0; i < 2 * lspname; i++) line[i] = '\0';
         if (!fgets(line, lline, fseq)) error2("EOF?");
         if (blankline(line)) {
            if (PopEmptyLines(fseq, lline, line))
            {
               printf("error in sequence data file: empty line (seq %d)\n", j + 1); exit(-1);
            }
         }
         p = line + (line[0] == '=' || line[0] == '>');
         while (isspace(*p)) p++;
         if ((ch = (int)(strstr(p, "  ") - p)) < lspname && ch > 0) lspname = ch;
         strncpy(com.spname[j], p, lspname);
         k = (int)strlen(com.spname[j]);
         p += (k < lspname ? k : lspname);

         for (; k > 0; k--) {/* trim spaces */
            if (!isgraph(com.spname[j][k]))   com.spname[j][k] = 0;
            else    break;
         }

         if (noisy >= 2) printf("Reading seq #%2d: %s       %s", j + 1, com.spname[j], (noisy > 3 ? "\n" : "\r"));
         for (k = 0; k < com.ls; p++) {
            while (*p == '\n' || *p == '\0') {
               p = fgets(line, lline, fseq);
               if (p == NULL)
               {
                  printf("\nEOF at site %d, seq %d\n", k + 1, j + 1); exit(-1);
               }
            }
            *p = (char)toupper(*p);
            if ((com.seqtype == BASEseq || com.seqtype == CODONseq) && *p == 'U')
               *p = 'T';
            p1 = strchr(pch, *p);
            if (p1 && p1 - pch >= nchar)
               miss = 1;
            if (*p == eq) {
               if (j == 0) error2("Error in sequence data file: . in 1st seq.?");
               com.z[j][k] = com.z[0][k];  k++;
            }
            else if (p1)
               com.z[j][k++] = *p;
            else if (isalpha(*p)) {
               printf("\nError in sequence data file: %c at %d seq %d.\n", *p, k + 1, j + 1);
               puts("Make sure to separate the sequence from its name by 2 or more spaces.");
               exit(0);
            }
            else if (*p == (char)EOF) error2("EOF?");
         }           /* for(k) */
         if (strchr(p, '\n') == NULL) /* pop up line return */
            while ((ch = fgetc(fseq)) != '\n' && ch != EOF);
      }   /* for (j,com.ns) */
   }
   else { /* interlaved */
      if (noisy) printf("Reading sequences, interlaved format..\n");
      for (j = 0; j < com.ns; j++) lt[j] = 0;  /* temporary seq length */
      for (igroup = 0; ; igroup++) {
         /*
          printf ("\nreading block %d ", igroup+1);  matIout(F0,lt,1,com.ns);*/

         for (j = 0; j < com.ns; j++)
            if (lt[j] < com.ls) break;
         if (j == com.ns) break;
         for (j = 0; j < com.ns; j++) {
            if (!fgets(line, lline, fseq)) {
               printf("\nerr reading site %d, seq %d group %d\nsites read in each seq:",
                  lt[j] + 1, j + 1, igroup + 1);
               error2("EOF?");
            }
            if (!hasbase(line)) {
               if (j) {
                  printf("\n%d, seq %d group %d", lt[j] + 1, j + 1, igroup + 1);
                  error2("empty line.");
               }
               else
                  if (PopEmptyLines(fseq, lline, line) == -1) {
                     printf("\n%d, seq %d group %d", lt[j] + 1, j + 1, igroup + 1);
                     error2("EOF?");
                  }
            }
            p = line;
            if (igroup == 0) {
               lspname = LSPNAME;
               while (isspace(*p)) p++;
               if ((ch = (int)(strstr(p, "  ") - p)) < lspname && ch > 0)
                  lspname = ch;
               strncpy(com.spname[j], p, lspname);
               k = (int)strlen(com.spname[j]);
               p += (k < lspname ? k : lspname);

               for (; k > 0; k--)   /* trim spaces */
                  if (!isgraph(com.spname[j][k]))
                     com.spname[j][k] = 0;
                  else
                     break;
               if (noisy >= 2) printf("Reading seq #%2d: %s     \r", j + 1, com.spname[j]);
            }
            for (; *p && *p != '\n'; p++) {
               if (lt[j] == com.ls) break;
               *p = (char)toupper(*p);
               if ((com.seqtype == BASEseq || com.seqtype == CODONseq) && *p == 'U')
                  *p = 'T';
               p1 = strchr(pch, *p);
               if (p1 && p1 - pch >= nchar)
                  miss = 1;
               if (*p == eq) {
                  if (j == 0) {
                     printf("err: . in 1st seq, group %d.\n", igroup);
                     exit(-1);
                  }
                  com.z[j][lt[j]] = com.z[0][lt[j]];
                  lt[j]++;
               }
               else if (p1)
                  com.z[j][lt[j]++] = *p;
               else if (isalpha(*p)) {
                  printf("\nerr: unrecognised character %c at %d seq %d block %d.",
                     *p, lt[j] + 1, j + 1, igroup + 1);
                  exit(-1);
               }
               else if (*p == (char)EOF) error2("EOF");
            }         /* for (*p) */
         }            /* for (j,com.ns) */

         if (noisy > 2) {
            printf("\nblock %3d:", igroup + 1);
            for (j = 0; j < com.ns; j++) printf(" %6d", lt[j]);
         }

      }               /* for (igroup) */
   }
   if (format == 2) {  /* NEXUS format: pop up extra lines until "end;"  */
      for (; ; ) {
         if (fgets(line, lline, fseq) == NULL) error2("seq err1: EOF");
         strcase(line, 0);
         if (strstr(line, NEXUSend)) break;
      }
   }
   free(line);

   if (read_seq_only) return (0);

   /*** delete empty sequences ******************/
#if(0)
   { int ns1 = com.ns, del[100] = { 0 };
   for (i = 0; i < com.ns; i++) {
      for (h = 0; h < com.ls; h++)   if (com.z[i][h] != '?') break;
      if (h == com.ls) {
         del[i] = 1;
         ns1--;
      }
   }
   fprintf(fnew, "\n\n%4d %6d\n", ns1, com.ls);
   for (i = 0; i < com.ns; i++) {
      if (!del[i]) {
         fprintf(fnew, "\n%-40s  ", com.spname[i]);
         for (h = 0; h < com.ls; h++) fprintf(fnew, "%c", com.z[i][h]);
      }
   }
   fflush(fnew);
   }
#endif


#ifdef CODEML
   /* mask stop codons as ???.  */
   if (com.seqtype == 1 && MarkStopCodons())
      miss = 1;
#endif

   if (!miss)
      com.cleandata = 1;
   else if (cleandata) {  /* forced removal of ambiguity characters */
      if (noisy > 2)  puts("\nSites with gaps or missing data are removed.");
      if (fout) {
         fprintf(fout, "\nBefore deleting alignment gaps\n");
         fprintf(fout, " %6d %6d\n", com.ns, com.ls);
         printsma(fout, com.spname, com.z, com.ns, com.ls, com.ls, gap, com.seqtype, 0, 0, NULL);
      }
      RemoveIndel();
      if (fout) fprintf(fout, "\nAfter deleting gaps. %d sites\n", com.ls);
   }

   if (fout && !com.readpattern) {/* verbose=1, listing sequences again */
      fprintf(fout, " %6d %6d\n", com.ns, com.ls);
      printsma(fout, com.spname, com.z, com.ns, com.ls, com.ls, gap, com.seqtype, 0, 0, NULL);
   }

   if (n31 == 3) com.ls /= n31;

   /* IdenticalSeqs(); */

#ifdef CODEML
   if (com.seqtype == 1 && com.verbose) Get4foldSites();

   if (com.seqtype == CODON2AAseq) {
      if (noisy > 2) puts("\nTranslating into AA sequences\n");
      for (j = 0; j < com.ns; j++) {
         if (noisy > 2) printf("Translating sequence %d\n", j + 1);
         DNA2protein(com.z[j], com.z[j], com.ls, com.icode);
      }
      com.seqtype = AAseq;

      if (fout) {
         fputs("\nTranslated AA Sequences\n", fout);
         fprintf(fout, "%4d  %6d", com.ns, com.ls);
         printsma(fout, com.spname, com.z, com.ns, com.ls, com.ls, 10, com.seqtype, 0, 0, NULL);
      }
   }
#endif

#if (defined CODEML || defined BASEML)
   if (com.ngene == 1 && com.Mgene == 1) com.Mgene = 0;
   if (com.ngene > 1 && com.Mgene == 1 && com.verbose)  printSeqsMgenes();

   if (com.bootstrap) { BootstrapSeq("boot.txt");  exit(0); }
#endif


#if (0 && defined CODEML)
   /* list sites with 2 types of serine codons: TC? and TCY.  19 March 2014, Ziheng. */
   if (com.seqtype == 1) {
      char codon[4] = "";
      int nbox0, nbox1, present = 0;
      for (h = 0; h < com.ls; h++) {
         for (i = 0, nbox0 = nbox1 = 0; i < com.ns; i++) {
            codon[0] = com.z[i][h * 3 + 0];
            codon[1] = com.z[i][h * 3 + 1];
            codon[2] = com.z[i][h * 3 + 2];
            if (codon[0] == 'T' && codon[1] == 'C') nbox0++;
            else if (codon[0] == 'A' && codon[1] == 'G' && (codon[2] == 'T' || codon[2] == 'C')) nbox1++;
         }
         if (nbox0 && nbox1 && nbox0 + nbox1 == com.ns) {
            present = 1;
            printf("\ncodon site %6d: ", h + 1);
            for (i = 0; i < com.ns; i++)
               printf("%c%c%c ", com.z[i][h * 3 + 0], com.z[i][h * 3 + 1], com.z[i][h * 3 + 2]);
         }
      }
      if (present) printf("\nAbove are 'synonymous' sites with 2 types of serine codons: TC? and TCY.\n");
   }
#endif

   if (noisy >= 2) printf("\nSequences read..\n");
   if (com.ls == 0) {
      puts("no sites. Got nothing to do");
      return(1);
   }

#if (defined MCMCTREE)
   /* Check and remove empty sequences.  */
   if (com.cleandata == 1)
      RemoveEmptySequences();

#endif


   /***** GIBBON, 2017.10.8 **/
   /* return(0); */

   if (!com.readpattern) {
      PatternWeight();
   }
   else {  /*  read pattern counts */
      com.npatt = com.ls;
      if ((com.fpatt = (double*)realloc(com.fpatt, com.npatt * sizeof(double))) == NULL)
         error2("oom fpatt");
      for (h = 0, lst = 0; h < com.npatt; h++) {
         fscanf(fseq, "%lf", &com.fpatt[h]);
         lst += com.fpatt[h];
         if (com.fpatt[h] < 0 || com.fpatt[h]>1e66)
            printf("fpatth[%d] = %.6g\n", h + 1, com.fpatt[h]);
      }
      if (lst > 1.00001) {
         com.ls = (int)lst;
         if (noisy) printf("\n%d site patterns read, %d sites\n", com.npatt, com.ls);
      }
      if (com.ngene == 1) {
         com.lgene[0] = com.ls;
         com.posG[0] = 0;
         com.posG[1] = com.npatt;
      }
      else {
         for (j = 0, com.posG[0] = 0; j < com.ngene; j++)
            com.posG[j + 1] = com.posG[j] + com.lgene[j];

         for (j = 0; j < com.ngene; j++) {
            com.lgene[j] = (j == 0 ? 0 : com.lgene[j - 1]);
            for (h = com.posG[j]; h < com.posG[j + 1]; h++)
               com.lgene[j] += (int)com.fpatt[h];
         }
      }
   }

   EncodeSeqs();
   if (fout && !com.readpattern) {
      fprintf(fout, "\nPrinting out site pattern counts\n\n");
      printPatterns(fout);

      if (com.verbose >= 2) {
         fprintf(fout, "\nSite-to-pattern map: ");
         for (i = 0; i < com.ls; i++)
            fprintf(fout, " %2d", com.pose[i] + 1);
         fprintf(fout, "\n");
      }
   }

   return (0);
}


#if (defined CODEML)

int MarkStopCodons(void)
{
   /* this converts the whole column into ??? if there is a stop codon in one sequence.
      Data in com.z[] are just read in and not encoded yet.
   */
   int i, j, h, NColumnEdited = 0, nstops = 0;
   char codon[4] = "", stops[6][4] = { "", "", "" };

   if (com.seqtype != 1) error2("should not be here");

   for (i = 0; i < 64; i++)
      if (GeneticCode[com.icode][i] == -1)
         getcodon(stops[nstops++], i);

   for (h = 0; h < com.ls / 3; h++) {
      for (i = 0; i < com.ns; i++) {
         codon[0] = com.z[i][h * 3 + 0];
         codon[1] = com.z[i][h * 3 + 1];
         codon[2] = com.z[i][h * 3 + 2];
         for (j = 0; j < nstops; j++)
            if (strcmp(codon, stops[j]) == 0) {
               printf("\nstop codon %s in seq. # %3d (%s), nucleotide site %d", codon, i + 1, com.spname[i], h * 3 + 1);
               break;
            }
         if (j < nstops) break;
      }
      if (i < com.ns) {
         for (i = 0; i < com.ns; i++)
            com.z[i][h * 3 + 0] = com.z[i][h * 3 + 1] = com.z[i][h * 3 + 2] = '?';
         NColumnEdited++;
      }
   }
   if (NColumnEdited) {
      printf("\n%2d columns are converted into ??? because of stop codons\nPress Enter to continue", NColumnEdited);
      getchar();
   }
   return(NColumnEdited);
}

#endif


void RemoveEmptySequences(void)
{
   /* this removes empty sequences (? or - only) and adjust com.ns
    */
   int j, h, nsnew;
   char emptyseq[NS];

   for (j = 0; j < com.ns; j++) {
      emptyseq[j] = 1;
      for (h = 0; h < com.ls * (com.seqtype == 1 ? 3 : 1); h++)
         if (com.z[j][h] != '?' && com.z[j][h] != '-') {
            emptyseq[j] = 0;
            break;
         }
   }
   for (j = 0, nsnew = 0; j < com.ns; j++) {
      if (emptyseq[j]) {
         printf("seq #%3d: %-30s is removed\n", j + 1, com.spname[j]);
         free(com.z[j]);
         free(com.spname[j]);
         continue;
      }
      com.z[nsnew] = com.z[j];
      com.spname[nsnew] = com.spname[j];
      nsnew++;
   }
   for (j = nsnew; j < com.ns; j++) {
      com.z[j] = NULL;
      com.spname[j] = NULL;
   }
   com.ns = nsnew;
   if (com.ns <= 0)
      error2("\nThis locus has nothing left after empty sequences are removed.  Use cleandata = 0.\n");
}


int printPatterns(FILE* fout)
{
   int j, h, n31 = (com.seqtype == CODONseq || com.seqtype == CODON2AAseq ? 3 : 1);
   int gap = (n31 == 3 ? 3 : 10); //, n=(com.seqtype==AAseq?20:4);

   fprintf(fout, "\n%10d %10d  P", com.ns, com.npatt * n31);
   if (com.ngene > 1) {
      fprintf(fout, " G\nG %d  ", com.ngene);
      for (j = 0; j < com.ngene; j++)
         fprintf(fout, "%7d", com.posG[j + 1] - com.posG[j]);
   }
   fprintf(fout, "\n");

   if (com.seqtype == 1 && com.cleandata) {
      ; /* nothing is printed out for yn00, as the coding is different. */
#if(defined CODEML || defined YN00)
      printsmaCodon(fout, com.z, com.ns, com.npatt, com.npatt, 1);
#endif
   }
   else
      printsma(fout, com.spname, com.z, com.ns, com.npatt, com.npatt, gap, com.seqtype, 1, 0, NULL);
   if (com.ls > 1.0001) {
      fprintf(fout, "\n");
      for (h = 0; h < com.npatt; h++) {
         fprintf(fout, " %4.0f", com.fpatt[h]);
         if ((h + 1) % 15 == 0) fprintf(fout, "\n");
      }
      fprintf(fout, "\n\n");
   }
   return(0);
}



void EncodeSeqs(void)
{
   /* This encodes sequences and set up com.TipMap[][], called after sites are collapsed
      into patterns.
   */
   int is, h, k, ch;
   char* pch = ((com.seqtype == 0 || com.seqtype == 1) ? BASEs : (com.seqtype == 2 ? AAs : BINs));

   if (com.seqtype != 1) {
      for (is = 0; is < com.ns; is++) {
         for (h = 0; h < com.npatt; h++) {
            ch = com.z[is][h];
            com.z[is][h] = (char)(k = (int)(strchr(pch, ch) - pch));
            if (k < 0) {
               printf("strange character %c in seq %d site %d\n", ch, is + 1, h + 1);
               exit(-1);
            }
         }
      }
   }
#if (defined CODEML || defined YN00)
   else if (com.seqtype == 1) {
      int n = com.ncode, i, j, nA;
      char b[3], ic, c[4] = "";

      /* collect all observed codons into CODONs, with a maximum of 256 distinct codons. */
      memset(&CODONs[0][0], 0, 256 * 4 * sizeof(char));
      for (nA = 0; nA < n; nA++) {
         ic = FROM61[nA]; b[0] = ic / 16; b[1] = (ic / 4) % 4; b[2] = ic % 4;
         for (i = 0; i < 3; i++) CODONs[nA][i] = BASEs[(int)b[i]];
      }
      for (j = 0, nA = n; j < com.ns; j++) {
         for (h = 0; h < com.npatt; h++) {
            for (k = 0; k < 3; k++) {
               c[k] = com.z[j][h * 3 + k];
               b[k] = (int)(strchr(BASEs, c[k]) - BASEs);
               if (b[k] < 0) printf("strange nucleotide %c in seq %d\n", c[k], j + 1);
            }
            if (b[0] < 4 && b[1] < 4 && b[2] < 4) {
               k = FROM64[b[0] * 16 + b[1] * 4 + b[2]];
               if (k < 0) {
                  printf("\nstop codon %s in seq #%2d: %s\n", c, j + 1, com.spname[j]);
                  printf("\ncodons in other sequences are\n");
                  for (i = 0; i < com.ns; i++) {
                     for (k = 0; k < 3; k++) c[k] = com.z[i][h * 3 + k];
                     printf("seq #%2d %-30s %s\n", i + 1, com.spname[i], c);
                  }
                  exit(-1);
               }
            }
            else {  /* an ambiguous codon */
               for (k = n; k < nA; k++)
                  if (strcmp(CODONs[k], c) == 0) break;
            }
            if (k == nA) {
               if (++nA > 256)
                  error2("too many ambiguity codons in the data.  Contact author");
               strcpy(CODONs[nA - 1], c);
            }
            com.z[j][h] = (char)k;
         }
         com.z[j] = (char*)realloc(com.z[j], com.npatt);
      }
      if (nA > n) {
         printf("%d ambiguous codons are seen in the data:\n", nA - n);
         for (k = n; k < nA; k++)  printf("%4s", CODONs[k]);
         printf("\n");
      }
   }
#endif
}


int print1site(FILE* fout, int h)
{
   /* This print out one site in the sequence data, com.z[].  It may be the h-th
   site in the original data file or the h-th pattern.  The data are coded.
   naa > 1 if the codon codes for more than one amino acid.
   */
   char* pch = (com.seqtype == 0 ? BASEs : (com.seqtype == 2 ? AAs : BINs));
   int i, b;

   for (i = 0; i < com.ns; i++) {
      b = com.z[i][h];
      if (com.seqtype == 0 || com.seqtype == 2)
         fprintf(fout, "%c", pch[b]);
#if defined(CODEML)
      else if (com.seqtype == 1) {
         if (com.model == FromCodon0)
            fprintf(fout, "%s ", CODONs[b]);
         else {
            int aa = GetAASiteSpecies(i, h);
            fprintf(fout, "%s (%c) ", CODONs[b], aa);
         }
      }
#endif
   }
   return(0);
}


void SetMapAmbiguity(int seqtype, int ModelAA2Codon)
{
   /* This sets up CharaMap, the map from the ambiguity characters to resolved characters.
   */
   int n = com.ncode, i, j;
   char* pch = (seqtype == 0 ? BASEs : (seqtype == 2 ? AAs : BINs));
   char* pbases = (seqtype == 0 ? BASEs : NULL);
   char** pEquateBASE = (seqtype == 0 ? EquateBASE : NULL);
   char debug = 0;

   for (j = 0; j < n; j++) {  /* basic characters, coded according to the definition in pch. */
      nChara[j] = (char)1;
      CharaMap[j][0] = (char)j;
   }
   if (seqtype != 1 && ModelAA2Codon == 0) {  /* nucleotides or amino acids */
      for (j = n, pch += n; *pch; j++, pch++) {
         if (com.seqtype == 0 || com.seqtype == 5) {  /* ambiguities are allowed for those 2 types */
            nChara[j] = (char)strlen(pEquateBASE[j]);
            for (i = 0; i < nChara[j]; i++)
               CharaMap[j][i] = (char)(strchr(pbases, pEquateBASE[j][i]) - pbases);
         }
         else {  /* for non-nucleotide characters, ambiguity characters must be ? or -. */
            nChara[j] = (char)n;
            for (i = 0; i < n; i++)
               CharaMap[j][i] = (char)i;
         }
         if (debug) {
            printf("character %c (%d): ", pbases[j], nChara[j]);
            for (i = 0; i < nChara[j]; i++)
               printf("%c", pbases[(int)CharaMap[j][i]]);
            printf("\n");
         }
      }
   }
#ifdef CODEML
   else if (seqtype == 1 && ModelAA2Codon == 0) {    /* codons */
      int i0, i1, i2, nb[3], ib[3][4], ic;

      for (j = n; j < 256 && CODONs[j][0]; j++) {
         nChara[j] = (char)0;
         for (i = 0; i < 3; i++)
            NucListall(CODONs[j][i], &nb[i], ib[i]);
         for (i0 = 0; i0 < nb[0]; i0++) {
            for (i1 = 0; i1 < nb[1]; i1++)
               for (i2 = 0; i2 < nb[2]; i2++) {
                  ic = ib[0][i0] * 16 + ib[1][i1] * 4 + ib[2][i2];
                  if (GeneticCode[com.icode][ic] != -1)
                     CharaMap[j][(int)nChara[j]++] = FROM64[ic];
               }
         }
         if (nChara[j] == 0) {
            printf("\ncodon %s is stop codon", CODONs[j]);
            exit(-1);
         }
      }
   }
   else if (seqtype == 2 && ModelAA2Codon) {  /* code the 20 amino acids as ambiguous codons.  X represents all codons */
      int iaa, ic;
      for (j = 64; j < 64 + 20; j++) {
         iaa = j - 64;
         for (ic = 0; ic < Nsensecodon; ic++) {
            if (iaa == GeneticCode[com.icode][FROM61[ic]])
               CharaMap[j][(int)nChara[j]++] = ic;
         }
      }
   }

#endif
}


int IdenticalSeqs(void)
{
   /* This checks for identical sequences and create a data set of unique
    sequences.  The file name is <SeqDataFile.unique.  This is casually
    written and need more testing.
    The routine is called right after the sequence data are read.
    For codon sequences, com.ls has the number of codons, which are NOT
    coded.
    */
   char tmpf[96], keep[NS];
   FILE* ftmp;
   int is, js, h, nkept = com.ns;
   int ls1 = com.ls * (com.seqtype == CODONseq || com.seqtype == CODON2AAseq ? 3 : 1);

   puts("\nIdenticalSeqs: not tested\a");
   for (is = 0; is < com.ns; is++)
      keep[is] = 1;
   for (is = 0; is < com.ns; is++) {
      if (!keep[is]) continue;
      for (js = 0; js < is; js++) {
         for (h = 0; h < ls1; h++)
            if (com.z[is][h] != com.z[js][h]) break;
         if (h == ls1) {
            printf("Seqs. %3d & %3d (%s & %s) are identical!\n",
               js + 1, is + 1, com.spname[js], com.spname[is]);
            keep[is] = 0;
         }
      }
   }
   for (is = 0; is < com.ns; is++)
      if (!keep[is]) nkept--;
   if (nkept < com.ns) {
      strcpy(tmpf, com.seqf);
      strcat(tmpf, ".unique");
      if ((ftmp = fopen(tmpf, "w")) == NULL) error2("IdenticalSeqs: file error");
      printSeqs(ftmp, com.z, com.spname, com.ns, com.ls, com.npatt, com.fpatt, NULL, keep, 1);
      fclose(ftmp);
      printf("\nUnique sequences collected in %s.\n", tmpf);
   }
   return(0);
}


void AllPatterns(FILE* fout)
{
   /* This prints out an alignment containting all possible site patterns, and then exits.
    This alignment may be useful to generate a dataset of infinitely long sequences,
    summarized in the site pattern probabilities.
    Because the PatternWeight() function changes the order of patterns, this routine
    prints out the alignment as one of patterns, with lots of 1's below it, to avoid
    baseml or codeml calling that routine to collaps sites.
    You then replace those 1'with the calculated pattern probabilities for further
    analysis.
    */
   int j, h, it, ic;
   // char codon[4]="   ";
   int n31 = (com.seqtype == CODONseq || com.seqtype == CODON2AAseq ? 3 : 1);
   int gap = (n31 == 3 ? 3 : 10);

   com.ns = 3;
   for (j = 0, com.npatt = 1; j < com.ns; j++) com.npatt *= com.ncode;
   printf("%3d species, %d site patterns\n", com.ns, com.npatt);
   com.cleandata = 1;
   for (j = 0; j < com.ns; j++) {
      com.spname[j] = (char*)realloc(com.spname[j], 11 * sizeof(char));
      sprintf(com.spname[j], "%c ", 'a' + j);
   }
   for (j = 0; j < com.ns; j++)
      if ((com.z[j] = (char*)malloc(com.npatt * sizeof(char))) == NULL)
         error2("oom in AllPatterns");
   for (h = 0; h < com.npatt; h++) {
      for (j = 0, it = h; j < com.ns; j++) {
         ic = it % com.ncode;
         it /= com.ncode;
         com.z[com.ns - 1 - j][h] = (char)ic;
      }
   }
   com.ls = com.npatt;

   fprintf(fout, " %6d %6d  P\n", com.ns, com.ls * n31);
   if (com.seqtype == 1) {
#if(defined CODEML || defined YN00)
      printsmaCodon(fout, com.z, com.ns, com.ls, com.ls, 0);
#endif
   }
   else
      printsma(fout, com.spname, com.z, com.ns, com.ls, com.ls, gap, com.seqtype, 1, 0, NULL);

   for (h = 0; h < com.npatt; h++) {
      fprintf(fout, " 1");
      if ((h + 1) % 40 == 0) fprintf(fout, "\n");
   }
   fprintf(fout, "\n");
   exit(0);
}


int PatternWeight(void)
{
   /* This collaps sites into patterns, for nucleotide, amino acid, or codon sequences.
      Site patterns are represented by 0-ended character strings.  The code adds 1 to the character
      during transposing so that it works for both unoded and coded data.
      com.pose[i] has labels for genes as input and maps sites to patterns in return.
      com.fpatt has site-pattern counts.
      This deals with multiple genes/partitions, and uses com.ngene, com.lgene[], com.posG[] etc.
      Sequences z[ns][ls] are copied into patterns zt[ls*lpatt], and bsearch is used
      twice to avoid excessive copying of com.fpatt[] and com.pose[], the first round to count
      npatt and generate the site patterns and the second round to generate fpatt[] & com.pose[].
   */
   int maxnpatt = com.ls, h, ip = -1, l, u, j, k, same, ig, * poset;
   // int gap = (com.seqtype==CODONseq ? 3 : 10);
   int n31 = (com.seqtype == CODONseq ? 3 : 1);
   int lpatt = com.ns * n31 + 1;   /* extra 0 used for easy debugging, can be voided */
   int* p2s;  /* point patterns to sites in zt */
   char timestr[36];
   char* p, * zt;
   double nc = (com.seqtype == 1 ? 64 : com.ncode) + !com.cleandata + 1;
   int debug = 0;

   /* (A) Collect and sort patterns.  Get com.npatt, com.lgene, com.posG.
   Move sequences com.z[ns][ls] into sites zt[ls*lpatt].
   Use p2s to map patterns to sites in zt to avoid copying.
   */
   if (noisy) printf("Counting site patterns.. %s\n", printtime(timestr));

   if ((com.seqtype == 1 && com.ns < 5) || (com.seqtype != 1 && com.ns < 7))
      maxnpatt = (int)(pow(nc, (double)com.ns) + 0.5) * com.ngene;
   if (maxnpatt > com.ls) maxnpatt = com.ls;
   p2s = (int*)malloc(maxnpatt * sizeof(int));
   zt = (char*)malloc(com.ls * lpatt * sizeof(char));
   if (p2s == NULL || zt == NULL)  error2("oom p2s or zt");
   memset(zt, 0, com.ls * lpatt * sizeof(char));
   for (j = 0; j < com.ns; j++)
      for (h = 0; h < com.ls; h++)
         for (k = 0; k < n31; k++)
            zt[h * lpatt + j * n31 + k] = (char)(com.z[j][h * n31 + k] + 1);

   for (ig = 0; ig < com.ngene; ig++) com.lgene[ig] = 0;
   for (ig = 0, com.npatt = 0; ig < com.ngene; ig++) {
      com.posG[ig] = l = u = ip = com.npatt;
      for (h = 0; h < com.ls; h++) {
         if (com.pose[h] != ig) continue;
         if (debug) printf("\nh %3d %s", h, zt + h * lpatt);

         /* bsearch in existing patterns.  Knuth 1998 Vol3 Ed2 p.410
         ip is the loc for match or insertion.  [l,u] is the search interval.
         */
         same = 0;
         if (com.lgene[ig]++ != 0) {  /* not 1st pattern? */
            for (l = com.posG[ig], u = com.npatt - 1; ; ) {
               if (u < l) break;
               ip = (l + u) / 2;
               k = strcmp(zt + h * lpatt, zt + p2s[ip] * lpatt);
               if (k < 0)        u = ip - 1;
               else if (k > 0)   l = ip + 1;
               else { same = 1;  break; }
            }
         }
         if (!same) {
            if (com.npatt > maxnpatt)
               error2("npatt > maxnpatt");
            if (l > ip) ip++;        /* last comparison in bsearch had k > 0. */
                                     /* Insert new pattern at ip.  This is the expensive step. */
            if (ip < com.npatt)
               memmove(p2s + ip + 1, p2s + ip, (com.npatt - ip) * sizeof(int));
            p2s[ip] = h;
            com.npatt++;
         }

         if (debug) {
            printf(": %3d (%c ilu %3d%3d%3d) ", com.npatt, (same ? '.' : '+'), ip, l, u);
            for (j = 0; j < com.npatt; j++)
               printf(" %s", zt + p2s[j] * lpatt);
         }
         if (noisy && ((h + 1) % 10000 == 0 || h + 1 == com.ls))
            printf("\rCompressing, %6d patterns at %6d / %6d sites (%.1f%%), %s",
               com.npatt, h + 1, com.ls, (h + 1.) * 100 / com.ls, printtime(timestr));

      }     /* for (h)  */
      if (noisy) printf("\n");
   }        /* for (ig) */

            /* (B) count pattern frequencies and collect pose[] */
   com.posG[com.ngene] = com.npatt;
   for (j = 0; j < com.ngene; j++)
      if (com.lgene[j] == 0)
         error2("some genes do not have any sites?");
   for (j = 1; j < com.ngene; j++)
      com.lgene[j] += com.lgene[j - 1];

   com.fpatt = (double*)realloc(com.fpatt, com.npatt * sizeof(double));
   poset = (int*)malloc(com.ls * sizeof(int));
   if (com.fpatt == NULL || poset == NULL) error2("oom poset");
   memset(com.fpatt, 0, com.npatt * sizeof(double));

   for (ig = 0; ig < com.ngene; ig++) {
      for (h = 0; h < com.ls; h++) {
         if (com.pose[h] != ig) continue;
         for (same = 0, l = com.posG[ig], u = com.posG[ig + 1] - 1; ; ) {
            if (u < l) break;
            ip = (l + u) / 2;
            k = strcmp(zt + h * lpatt, zt + p2s[ip] * lpatt);
            if (k < 0)        u = ip - 1;
            else if (k > 0)   l = ip + 1;
            else { same = 1;  break; }
         }
         com.fpatt[ip]++;
         poset[h] = ip;
         if (noisy && ((h + 1) % 10000 == 0 || h + 1 == com.ls))
            printf("\rCollecting fpatt[] & pose[], %6d patterns at %6d / %6d sites (%.1f%%), %s",
               com.npatt, h + 1, com.ls, (h + 1.) * 100 / com.ls, printtime(timestr));
      }     /* for (h)  */
      if (noisy) printf("\n");
   }        /* for (ig) */

   if (com.seqtype == CODONseq && com.ngene == 3 && com.lgene[0] == com.ls / 3)
      puts("\nCheck option G in data file?\n");

   for (j = 0; j < com.ns; j++) {
      com.z[j] = (char*)realloc(com.z[j], com.npatt * n31 * sizeof(char));
      for (ip = 0, p = com.z[j]; ip < com.npatt; ip++)
         for (k = 0; k < n31; k++)
            *p++ = (char)(zt[p2s[ip] * lpatt + j * n31 + k] - 1);
   }
   memcpy(com.pose, poset, com.ls * sizeof(int));
   free(poset);  free(p2s);  free(zt);
   return (0);
}



void AddFreqSeqGene(int js, int ig, double pi0[], double pi[]);


void Chi2FreqHomo(double f[], int ns, int nc, double X2G[2])
{
   /* This calculates a chi-square like statistic for testing that the base
    or amino acid frequencies are identical among sequences.
    f[ns*nc] where ns is #sequences (rows) and nc is #states (columns).
    */
   int i, j;
   double mf[64] = { 0 }, smallv = 1e-50;

   X2G[0] = X2G[1] = 0;
   for (i = 0; i < ns; i++)
      for (j = 0; j < nc; j++)
         mf[j] += f[i * nc + j] / ns;

   for (i = 0; i < ns; i++) {
      for (j = 0; j < nc; j++) {
         if (mf[j] > smallv) {
            X2G[0] += square(f[i * nc + j] - mf[j]) / mf[j];
            if (f[i * nc + j])
               X2G[1] += 2 * f[i * nc + j] * log(f[i * nc + j] / mf[j]);
         }
      }
   }
}

int InitializeBaseAA(FILE* fout)
{
   /* Count site patterns (com.fpatt) and calculate base or amino acid frequencies
    in genes and species.  This works on raw (uncoded) data.
    Ambiguity characters in sequences are resolved by iteration.
    For frequencies in each species, they are resolved within that sequence.
    For average base frequencies among species, they are resolved over all
    species.

    This routine is called by baseml and aaml.  codonml uses another
    routine InitializeCodon()
    */
   char* pch = (com.seqtype == 0 ? BASEs : (com.seqtype == 2 ? AAs : BINs));
   char indel[] = "-?";
   int wname = 30, h, js, k, ig, nconstp, n = com.ncode;
   int irf, nrf = 20;
   double pi0[20], t, lmax = 0, X2G[2], * pisg;  /* freq for species & gene, for X2 & G */

   if (noisy) printf("Counting frequencies..");
   if (fout)  fprintf(fout, "\nFrequencies..");
   if ((pisg = (double*)malloc(com.ns * n * sizeof(double))) == NULL)
      error2("oom pisg");
   for (h = 0, nconstp = 0; h < com.npatt; h++) {
      for (js = 1; js < com.ns; js++)
         if (com.z[js][h] != com.z[0][h])  break;
      if (js == com.ns && com.z[0][h] != indel[0] && com.z[0][h] != indel[1])
         nconstp += (int)com.fpatt[h];
   }
   for (ig = 0, zero(com.pi, n); ig < com.ngene; ig++) {
      if (com.ngene > 1)
         fprintf(fout, "\n\nGene %2d (len %4d)", ig + 1, com.lgene[ig] - (ig == 0 ? 0 : com.lgene[ig - 1]));
      fprintf(fout, "\n%*s", wname, "");
      for (k = 0; k < n; k++) fprintf(fout, "%7c", pch[k]);

      /* The following block calculates freqs in each species for each gene.
       Ambiguities are resolved in each species.  com.pi and com.piG are
       used for output only, and are not be used later with missing data.
       */
      zero(com.piG[ig], n);
      zero(pisg, com.ns * n);
      for (js = 0; js < com.ns; js++) {
         fillxc(pi0, 1.0 / n, n);
         for (irf = 0; irf < nrf; irf++) {
            zero(com.pi, n);
            AddFreqSeqGene(js, ig, pi0, com.pi);
            t = sum(com.pi, n);
            if (t < 1e-10) {
               printf("Some sequences are empty.\n");
               fillxc(com.pi, 1.0 / n, n);
            }
            else
               abyx(1 / t, com.pi, n);
            if (com.cleandata || com.cleandata || (t = distance(com.pi, pi0, n)) < 1e-8)
               break;
            xtoy(com.pi, pi0, n);
         }   /* for(irf) */
         fprintf(fout, "\n%-*s", wname, com.spname[js]);
         for (k = 0; k < n; k++) fprintf(fout, "%8.5f", com.pi[k]);
         if (com.ncode == 4 && com.ngene == 1) fprintf(fout, " GC = %5.3f", com.pi[1] + com.pi[3]);
         for (k = 0; k < n; k++) com.piG[ig][k] += com.pi[k] / com.ns;
         xtoy(com.pi, pisg + js * n, n);
      }    /* for(js,ns) */
      if (com.ngene > 1) {
         fprintf(fout, "\n\n%-*s", wname, "Mean");
         for (k = 0; k < n; k++) fprintf(fout, "%7.4f", com.piG[ig][k]);
      }

      Chi2FreqHomo(pisg, com.ns, n, X2G);

      fprintf(fout, "\n\nHomogeneity statistic: X2 = %.5f G = %.5f ", X2G[0], X2G[1]);

      /* fprintf(frst1,"\t%.5f", X2G[1]); */

   }  /* for(ig) */
   if (noisy) printf("\n");

   /* If there are missing data, the following block calculates freqs
    in each gene (com.piG[]), as well as com.pi[] for the entire sequence.
    Ambiguities are resolved over entire data sets across species (within
    each gene for com.piG[]).  These are used in ML calculation later.
    */
   if (com.cleandata) {
      for (ig = 0, zero(com.pi, n); ig < com.ngene; ig++) {
         t = (ig == 0 ? com.lgene[0] : com.lgene[ig] - com.lgene[ig - 1]) / (double)com.ls;
         for (k = 0; k < n; k++)  com.pi[k] += com.piG[ig][k] * t;
      }
   }
   else {
      for (ig = 0; ig < com.ngene; ig++) {
         xtoy(com.piG[ig], pi0, n);
         for (irf = 0; irf < nrf; irf++) {  /* com.piG[] */
            zero(com.piG[ig], n);
            for (js = 0; js < com.ns; js++)
               AddFreqSeqGene(js, ig, pi0, com.piG[ig]);
            t = sum(com.piG[ig], n);
            if (t < 1e-10) {
               puts("empty sequences?");
               for (k = 0; k < n; k++)  com.piG[ig][k] = 1.0 / n;
               t = 1;
            }
            abyx(1 / t, com.piG[ig], n);
            if (distance(com.piG[ig], pi0, n) < 1e-8) break;
            xtoy(com.piG[ig], pi0, n);
         }         /* for(irf) */
      }            /* for(ig) */
      zero(pi0, n);
      for (k = 0; k < n; k++) for (ig = 0; ig < com.ngene; ig++)
         pi0[k] += com.piG[ig][k] / com.ngene;
      for (irf = 0; irf < nrf; irf++) {  /* com.pi[] */
         zero(com.pi, n);
         for (ig = 0; ig < com.ngene; ig++)  for (js = 0; js < com.ns; js++)
            AddFreqSeqGene(js, ig, pi0, com.pi);

         t = sum(com.pi, n);
         if (t < 1e-10) {
            for (k = 0; k < n; k++)  com.pi[k] = 1.0 / n;
            t = 1;
         }
         abyx(1 / sum(com.pi, n), com.pi, n);
         if (distance(com.pi, pi0, n) < 1e-8) break;
         xtoy(com.pi, pi0, n);
      }            /* for(ig) */
   }
   fprintf(fout, "\n\n%-*s", wname, "Average");
   for (k = 0; k < n; k++) fprintf(fout, "%8.5f", com.pi[k]);
   if (!com.cleandata) fputs("\n(Ambiguity characters are used to calculate freqs.)\n", fout);

   fprintf(fout, "\n\n# constant sites: %6d (%.2f%%)",
      nconstp, (double)nconstp * 100. / com.ls);

   if (com.model == 0 || (com.seqtype == BASEseq && com.model == 1)) {
      fillxc(com.pi, 1. / n, n);
      for (ig = 0; ig < com.ngene; ig++)
         xtoy(com.pi, com.piG[ig], n);
   }
   if (com.seqtype == BASEseq && com.model == 5) { /* T92 model */
      com.pi[0] = com.pi[2] = (com.pi[0] + com.pi[2]) / 2;
      com.pi[1] = com.pi[3] = (com.pi[1] + com.pi[3]) / 2;
      for (ig = 0; ig < com.ngene; ig++) {
         com.piG[ig][0] = com.piG[ig][2] = (com.piG[ig][0] + com.piG[ig][2]) / 2;
         com.piG[ig][1] = com.piG[ig][3] = (com.piG[ig][1] + com.piG[ig][3]) / 2;
      }
   }

   /* this is used only for REV & REVu in baseml and model==3 in aaml */
   if (com.seqtype == AAseq) {
      for (k = 0, t = 0; k < n; k++)  t += (com.pi[k] > 0);
      if (t <= 4)
         puts("\n\a\t\tAre these a.a. sequences?");
   }
   if (com.cleandata && com.ngene == 1) {
      for (h = 0, lmax = -(double)com.ls * log((double)com.ls); h < com.npatt; h++)
         if (com.fpatt[h] > 1) lmax += com.fpatt[h] * log((double)com.fpatt[h]);
   }
   if (fout) {
      if (lmax) fprintf(fout, "\nln Lmax (unconstrained) = %.6f\n", lmax);
      fflush(fout);
   }

   free(pisg);
   return(0);
}


void AddFreqSeqGene(int js, int ig, double pi0[], double pi[])
{
   /* This adds the character counts in sequence js in gene ig to pi,
    using pi0, by resolving ambiguities.  The data are coded.  com.cleandata==1 or 0.
    This is for nucleotide and amino acid sequences only.
    */
    //char *pch = (com.seqtype==0 ? BASEs : (com.seqtype==2 ? AAs : BINs));
   int k, h, b, n = com.ncode;
   double t;

   if (com.cleandata) {
      for (h = com.posG[ig]; h < com.posG[ig + 1]; h++)
         pi[(int)com.z[js][h]] += com.fpatt[h];
   }
   else {
      for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
         b = com.z[js][h];
         if (b < n)
            pi[b] += com.fpatt[h];
         else {
            /*
             if(com.seqtype==BASEseq) {
             NucListall(BASEs[b], &nb, ib);
             for(k=0,t=0; k<nb; k++) t += pi0[ib[k]];
             for(k=0; k<nb; k++)
             pi[ib[k]] += pi0[ib[k]]/t * com.fpatt[h];
             }
             */
            if (com.seqtype == BASEseq) {
               for (k = 0, t = 0; k < nChara[b]; k++)
                  t += pi0[(int)CharaMap[b][k]];
               for (k = 0; k < nChara[b]; k++)
                  pi[(int)CharaMap[b][k]] += pi0[(int)CharaMap[b][k]] / t * com.fpatt[h];
            }
            else if (com.seqtype == AAseq)  /* unrecognized AAs are treated as "?". */
               for (k = 0; k < n; k++) pi[k] += pi0[k] * com.fpatt[h];
         }
      }
   }
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
   int  n = com.ncode, h, k, j, js, lnew, nindel, n31 = 1;
   char b, * miss;  /* miss[h]=1 if site (codon) h is missing, 0 otherwise */
   char* pch = ((com.seqtype <= 1 || com.seqtype == CODON2AAseq) ? BASEs : (com.seqtype == 2 ? AAs : BINs));

   if (com.seqtype == CODONseq || com.seqtype == CODON2AAseq) {
      n31 = 3; n = 4;
   }

   if (com.ls % n31) error2("ls in RemoveIndel.");
   if ((miss = (char*)malloc(com.ls / n31 * sizeof(char))) == NULL)
      error2("oom miss");
   for (h = 0; h < com.ls / n31; h++)
      miss[h] = 0;
   for (js = 0; js < com.ns; js++) {
      for (h = 0, nindel = 0; h < com.ls / n31; h++) {
         for (k = 0; k < n31; k++) {
            b = (char)toupper(com.z[js][h * n31 + k]);
            for (j = 0; j < n; j++)
               if (b == pch[j]) break;
            if (j == n) {
               miss[h] = 1; nindel++;
            }
         }
      }
      if (noisy > 2 && nindel)
         printf("\n%6d ambiguity characters in seq. %d", nindel, js + 1);
   }
   if (noisy > 2) {
      for (h = 0, k = 0; h < com.ls / n31; h++)  if (miss[h]) k++;
      printf("\n%d sites are removed. ", k);
      if (k < 1000)
         for (h = 0; h < com.ls / n31; h++)  if (miss[h]) printf(" %2d", h + 1);
   }

   for (h = 0, lnew = 0; h < com.ls / n31; h++) {
      if (miss[h]) continue;
      for (js = 0; js < com.ns; js++) {
         for (k = 0; k < n31; k++)
            com.z[js][lnew * n31 + k] = com.z[js][h * n31 + k];
      }
      com.pose[lnew] = com.pose[h];
      lnew++;
   }
   com.ls = lnew * n31;
   free(miss);
   return (0);
}



int MPInformSites(void)
{
   /* Outputs parsimony informative and noninformative sites into
    two files named MPinf.seq and MPninf.seq
    Uses transformed sequences.
    Not used for a long time.  Does not work if com.pose is NULL.
    */
   char* imark;
   char* pch = (com.seqtype == 0 ? BASEs : (com.seqtype == 2 ? AAs : BINs));
   int h, i, markb[NS], inf, lsinf;
   FILE* finf, * fninf;

   puts("\nMPInformSites: missing data not dealt with yet?\n");

   finf = fopen("MPinf.seq", "w");
   fninf = fopen("MPninf.seq", "w");
   if (finf == NULL || fninf == NULL) error2("MPInformSites: file creation error");

   puts("\nSorting parsimony-informative sites: MPinf.seq & MPninf.seq");
   if ((imark = (char*)malloc(com.ls * sizeof(char))) == NULL) error2("oom imark");
   for (h = 0, lsinf = 0; h < com.ls; h++) {
      for (i = 0; i < com.ns; i++) markb[i] = 0;
      for (i = 0; i < com.ns; i++) markb[(int)com.z[i][com.pose[h]]]++;

      for (i = 0, inf = 0; i < com.ncode; i++)  if (markb[i] >= 2)  inf++;
      if (inf >= 2) { imark[h] = 1; lsinf++; }
      else imark[h] = 0;
   }
   fprintf(finf, "%6d%6d\n", com.ns, lsinf);
   fprintf(fninf, "%6d%6d\n", com.ns, com.ls - lsinf);
   for (i = 0; i < com.ns; i++) {
      fprintf(finf, "\n%s\n", com.spname[i]);
      fprintf(fninf, "\n%s\n", com.spname[i]);
      for (h = 0; h < com.ls; h++)
         fprintf((imark[h] ? finf : fninf), "%c", pch[(int)com.z[i][com.pose[h]]]);
      fprintf(finf, "\n\n");
   }
   free(imark);
   fclose(finf);  fclose(fninf);
   return (0);
}



int Site2Pattern(FILE* fout)
{
   int h;
   fprintf(fout, "\n\nMapping site to pattern (i.e. site %d has pattern %d):\n",
      com.ls - 1, com.pose[com.ls - 2] + 1);
   for (h = 0; h < com.ls; h++) {
      fprintf(fout, "%6d", com.pose[h] + 1);
      if ((h + 1) % 10 == 0)  fprintf(fout, "\n");
   }
   fprintf(fout, "\n");
   return (0);
}


#endif



int print1seq(FILE* fout, char* z, int ls, int pose[])
{
   /* This prints out one sequence, and the sequences are encoded.
    z[] contains patterns if (pose!=NULL)
    This uses com.seqtype.
    */
   int h, hp, gap = 10;
   char* pch = (com.seqtype == 0 ? BASEs : (com.seqtype == 2 ? AAs : BINs));
   // char str[4]="";
   // int nb = (com.seqtype==CODONseq?3:1);

   for (h = 0; h < ls; h++) {
      hp = (pose ? pose[h] : h);
      if (com.seqtype != CODONseq) {
         fprintf(fout, "%c", pch[(int)z[hp]]);
         if ((h + 1) % gap == 0) fputc(' ', fout);
      }
      else
         fprintf(fout, "%s ", CODONs[(int)z[hp]]);
   }
   return(0);
}

void printSeqs(FILE* fout, char* z[], char* spnames[], int ns, int ls, int npatt, double fpatt[], int* pose, char keep[], int format)
{
   /* Print sequences into fout, using paml (format=0 or 1) or NEXUS (format=2)
      formats.
      Use pose=NULL if called before site patterns are collapsed.
      keep[] marks the sequences to be printed.  Use NULL for keep if all sequences
      are to be printed.
      Sequences may (com.cleandata==1) and may not (com.cleandata=0) be coded.
      z[] has site patterns if pose!=NULL.
      This uses com.seqtype, and ls is the number of codons for codon seqs.
      See notes in print1seq()

      format = 0,1: PAML sites or patterns
      2: NEXUS Nexus format.
      3: NEXUS Nexus JC69 format

      This is used by evolver.  Check and merge with printsma().
   */
   int h, i, j, ls1, n31 = (com.seqtype == 1 ? 3 : 1), nskept = ns, wname = 10;
   char* dt = (com.seqtype == AAseq ? "protein" : "dna");
   char* pch = (com.seqtype == 0 ? BASEs : AAs);

   ls1 = (format == 1 ? npatt : ls);
   if (keep)
      for (j = 0; j < ns; j++) nskept -= !keep[j];
   if (format == 0 || format == 1)
      fprintf(fout, "\n\n%d %d %s\n\n", nskept, ls1 * n31, (format == 1 ? " P" : ""));
   else if (format == 2 || format == 3) {  /* NEXUS format */
      fprintf(fout, "\nbegin data;\n");
      fprintf(fout, "   dimensions ntax=%d nchar=%d;\n", nskept, ls1 * n31);
      fprintf(fout, "   format datatype=%s missing=? gap=-;\n   matrix\n", dt);
   }

   for (j = 0; j < ns; j++) {
      if (keep && !keep[j]) continue;
      fprintf(fout, "%s%-*s  ", (format == 2 || format == 3 ? "      " : ""), wname, spnames[j]);
      if (format == 3) {     /* NEXUS Nexus JC69 format */
         for (h = 0, ls1 = 0; h < com.npatt; h++)
            for (i = 0; i < (int)com.fpatt[h]; i++) {
               fprintf(fout, "%c", pch[(int)z[j][h]]);
               if (++ls1 % 10 == 0) fprintf(fout, " ");
            }
      }
      else
         print1seq(fout, z[j], (format == 1 ? npatt : ls), pose);
      fprintf(fout, "\n");
   }
   if (format == 2 || format == 3) fprintf(fout, "   ;\nend;");
   else if (format == 1) {
      fprintf(fout, "\n");
      for (h = 0; h < com.npatt; h++) {
         /* fprintf(fout," %12.8f", com.fpatt[h]/(double)ls); */
         fprintf(fout, " %4.0f", com.fpatt[h]);
         if ((h + 1) % 15 == 0) fprintf(fout, "\n");
      }
   }
   fprintf(fout, "\n\n");
   fflush(fout);
}


#define gammap(x,alpha) (alpha*(1-pow(x,-1.0/alpha)))
/* DistanceREV () used to be here, moved to pamp.
 */

#if (defined BASEML || defined BASEMLG || defined MCMCTREE || defined PROBTREE || defined YULETREE)

double SeqDivergence(double x[], int model, double alpha, double* kappa)
{
   /* alpha=0 if no gamma
    return -1 if in error.
    Check DistanceF84() if variances are wanted.
    */
   int i, j;
   double p[4], Y, R, a1, a2, b, P1, P2, Q, fd, tc, ag, GC;
   double d, smallv = 1e-10 / com.ls, largek = 999, larged = 9;

   if (testXMat(x)) {
      matout(F0, x, 4, 4);
      printf("\nfrequency matrix error, setting distance to large d");
      return(larged);
   }
   for (i = 0, fd = 1, zero(p, 4); i < 4; i++) {
      fd -= x[i * 4 + i];
      for (j = 0; j < 4; j++)
      {
         p[i] += x[i * 4 + j] / 2;  p[j] += x[i * 4 + j] / 2;
      }
   }
   P1 = x[0 * 4 + 1] + x[1 * 4 + 0];
   P2 = x[2 * 4 + 3] + x[3 * 4 + 2];
   Q = x[0 * 4 + 2] + x[0 * 4 + 3] + x[1 * 4 + 2] + x[1 * 4 + 3] + x[2 * 4 + 0] + x[2 * 4 + 1] + x[3 * 4 + 0] + x[3 * 4 + 1];
   if (fd < smallv)
      return(0);
   if (P1 < smallv) P1 = 0;
   if (P2 < smallv) P2 = 0;
   if (Q < smallv) Q = 0;
   Y = p[0] + p[1];    R = p[2] + p[3];  tc = p[0] * p[1]; ag = p[2] * p[3];

   switch (model) {
   case (JC69):
      for (i = 0; i < 4; i++) p[i] = .25;
   case (F81):
      for (i = 0, b = 0; i < 4; i++)  b += p[i] * (1 - p[i]);
      if (1 - fd / b <= 0) return (larged);

      if (alpha <= 0) return (-b * log(1 - fd / b));
      else return  (-b * gammap(1 - fd / b, alpha));
   case (K80):
      a1 = 1 - 2 * (P1 + P2) - Q;   b = 1 - 2 * Q;
      /*      if (a1<=0 || b<=0) return (-1); */
      if (a1 <= 0 || b <= 0) return (larged);
      if (alpha <= 0) { a1 = -log(a1);  b = -log(b); }
      else { a1 = -gammap(a1, alpha); b = -gammap(b, alpha); }
      a1 = .5 * a1 - .25 * b;  b = .25 * b;
      if (b > smallv) *kappa = a1 / b; else *kappa = largek;
      return (a1 + 2 * b);
   case (F84):
      if (Y < smallv || R < smallv) {
         *kappa = -1;  d = larged;
      }
      else {
         a1 = (2 * (tc + ag) + 2 * (tc * R / Y + ag * Y / R) * (1 - Q / (2 * Y * R)) - P1 - P2) / (2 * tc / Y + 2 * ag / R);
         b = 1 - Q / (2 * Y * R);
         /*      if (a1<=0 || b<=0) return (-1); */
         if (a1 <= 0 || b <= 0) return (larged);
         if (alpha <= 0) { a1 = -log(a1); b = -log(b); }
         else { a1 = -gammap(a1, alpha); b = -gammap(b, alpha); }
         a1 = .5 * a1;  b = .5 * b;
         *kappa = a1 / b - 1;
         *kappa = max2(*kappa, -.5);
         d = 4 * b * (tc * (1 + *kappa / Y) + ag * (1 + *kappa / R) + Y * R);
      }
      return d;
   case (HKY85):         /* HKY85, from Rzhetsky & Nei (1995 MBE 12, 131-51) */
      if ((0 && Y < smallv) || R < smallv) {
         *kappa = -1;  d = larged;
      }
      else {
         *kappa = largek;
         a1 = 1 - Y * P1 / (2 * tc) - Q / (2 * Y);
         a2 = 1 - R * P2 / (2 * ag) - Q / (2 * R);
         b = 1 - Q / (2 * Y * R);
         if (a1 <= 0 || a2 <= 0 || b <= 0) return (larged);
         if (alpha <= 0) { a1 = -log(a1); a2 = -log(a2); b = -log(b); }
         else { a1 = -gammap(a1, alpha); a2 = -gammap(a2, alpha); b = -gammap(b, alpha); }
         a1 = -R / Y * b + a1 / Y;
         a2 = -Y / R * b + a2 / R;
         if (b > 0) *kappa = min2((a1 + a2) / (2 * b), largek);
         d = 2 * (p[0] * p[1] + p[2] * p[3]) * (a1 + a2) / 2 + 2 * Y * R * b;
      }
      return d;
   case (T92):
      *kappa = largek;
      GC = p[1] + p[3];
      a1 = 1 - Q - (P1 + P2) / (2 * GC * (1 - GC));   b = 1 - 2 * Q;
      if (a1 <= 0 || b <= 0) return (larged);
      if (alpha <= 0) { a1 = -log(a1); b = -log(b); }
      else { a1 = -gammap(a1, alpha); b = -gammap(b, alpha); }
      if (Q > 0) *kappa = 2 * a1 / b - 1;
      return 2 * GC * (1 - GC) * a1 + (1 - 2 * GC * (1 - GC)) / 2 * b;
   case (TN93):         /* TN93  */
      if (Y < smallv || R < smallv || tc < smallv || ag < smallv) {
         *kappa = -1;  d = larged;
      }
      else {
         a1 = 1 - Y * P1 / (2 * tc) - Q / (2 * Y);
         a2 = 1 - R * P2 / (2 * ag) - Q / (2 * R);
         b = 1 - Q / (2 * Y * R);
         /*      if (a1<=0 || a2<=0 || b<=0) return (-1); */
         if (a1 <= 0 || a2 <= 0 || b <= 0) return (larged);
         if (alpha <= 0) { a1 = -log(a1); a2 = -log(a2); b = -log(b); }
         else { a1 = -gammap(a1, alpha); a2 = -gammap(a2, alpha); b = -gammap(b, alpha); }
         a1 = .5 / Y * (a1 - R * b);  a2 = .5 / R * (a2 - Y * b);  b = .5 * b;
         *kappa = largek;
         /*
         printf("\nk1&k2 = %.6f %.6f\n", a1/b,a2/b);
         */
         if (b > 0) *kappa = min2((a1 + a2) / (2 * b), largek);
         d = 4 * p[0] * p[1] * a1 + 4 * p[2] * p[3] * a2 + 4 * Y * R * b;
      }
      return d;
   }
   return (-1);
}


double DistanceIJ(int is, int js, int model, double alpha, double* kappa)
{
   /* Distance between sequences is and js.
    See DistanceMatNuc() for more details.
    */
   char b0, b1;
   int h, n = 4;
   double x[16], sumx, larged = 9;

   zero(x, 16);
   if (com.cleandata && com.seqtype == 0) {
      for (h = 0; h < com.npatt; h++)
         x[com.z[is][h] * n + com.z[js][h]] += com.fpatt[h];
   }
   else {
      for (h = 0; h < com.npatt; h++) {
         b0 = com.z[is][h];
         b1 = com.z[js][h];
         if (b0 < n && b1 < n)
            x[b0 * n + b1] += com.fpatt[h];
         /*
         else
            missing = 1;
         */
      }
   }
   sumx = sum(x, 16);

   if (sumx <= 0) return(larged);    /* questionable??? */
   abyx(1. / sum(x, 16), x, 16);
   return SeqDivergence(x, model, alpha, kappa);
}


#if (defined LSDISTANCE && defined REALSEQUENCE)

extern double* SeqDistance;

int DistanceMatNuc(FILE* fout, FILE* f2base, int model, double alpha)
{
   /* This calculates pairwise distances.  The data may be clean and coded
    (com.cleandata==1) or not.  In the latter case, ambiguity sites are not
    used (pairwise deletion).  Site patterns are used.
    */
   int is, js, status = 0;
   double kappat = 0, t, bigD = 9;

   if (f2base) fprintf(f2base, "%6d\n", com.ns);
   if (model >= REV) model = TN93; /* TN93 here */
   if (fout) {
      fprintf(fout, "\nDistances:%5s", models[model]);
      if (model != JC69 && model != F81) fprintf(fout, " (kappa) ");
      fprintf(fout, " (alpha set at %.2f)\n", alpha);
      fprintf(fout, "This matrix is not used in later m.l. analysis.\n");
      if (!com.cleandata) fprintf(fout, "\n(Pairwise deletion.)");
   }
   for (is = 0; is < com.ns; is++) {
      if (fout) fprintf(fout, "\n%-15s  ", com.spname[is]);
      if (f2base) fprintf(f2base, "%-15s   ", com.spname[is]);
      for (js = 0; js < is; js++) {
         t = DistanceIJ(is, js, model, alpha, &kappat);
         if (t < 0) { t = bigD; status = -1; }
         SeqDistance[is * (is - 1) / 2 + js] = t;
         if (f2base) fprintf(f2base, " %8.6f", t);
         if (fout) fprintf(fout, " %9.6f", t);
         if (fout && (model == K80 || model >= F84))
            fprintf(fout, " (%7.4f)", kappat);
      }
      if (f2base) fprintf(f2base, "\n");
   }
   if (fout) fprintf(fout, "\n");
   if (status) puts("\ndistance formula sometimes inapplicable..");
   return(status);
}



#endif


#ifdef BASEMLG
extern int CijkIs0[];
#endif

extern int nR;
extern double Cijk[], Root[];

int QTN93(int model, double Q[], double kappa1, double kappa2, double pi[])
{
   int i, j;
   double T = pi[0], C = pi[1], A = pi[2], G = pi[3], Y = T + C, R = A + G, scalefactor;

   if (model == JC69 || model == F81) kappa1 = kappa2 = com.kappa = 1;
   else if (com.model < TN93)       kappa2 = kappa1;
   if (model == F84) { kappa2 = 1 + kappa1 / R; kappa1 = 1 + kappa1 / Y; }
   scalefactor = 1 / (2 * T * C * kappa1 + 2 * A * G * kappa2 + 2 * Y * R);

   for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) Q[i * 4 + j] = (i == j ? 0 : 1);
   Q[0 * 4 + 1] = Q[1 * 4 + 0] = kappa1;
   Q[2 * 4 + 3] = Q[3 * 4 + 2] = kappa2;
   for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) Q[i * 4 + j] *= pi[j] * scalefactor;
   for (i = 0; i < 4; i++) { Q[i * 4 + i] = 0;  Q[i * 4 + i] = -sum(Q + i * 4, 4); }

   return (0);
}

int RootTN93(int model, double kappa1, double kappa2, double pi[],
   double* scalefactor, double Root[])
{
   double T = pi[0], C = pi[1], A = pi[2], G = pi[3], Y = T + C, R = A + G;

   if (model == JC69 || model == F81) kappa1 = kappa2 = com.kappa = 1;
   else if (com.model < TN93)       kappa2 = kappa1;
   if (model == F84) { kappa2 = 1 + kappa1 / R; kappa1 = 1 + kappa1 / Y; }

   *scalefactor = 1 / (2 * T * C * kappa1 + 2 * A * G * kappa2 + 2 * Y * R);

   Root[0] = 0;
   Root[1] = -(*scalefactor);
   Root[2] = -(Y + R * kappa2) * (*scalefactor);
   Root[3] = -(Y * kappa1 + R) * (*scalefactor);
   return (0);
}


int eigenTN93(int model, double kappa1, double kappa2, double pi[],
   int* nR, double Root[], double Cijk[])
{
   /* initialize Cijk[] & Root[], which are the only part to be changed
    for a new substitution model
    for JC69, K80, F81, F84, HKY85, TN93
    Root: real Root divided by v, the number of nucleotide substitutions.
    */
   int i, j, k, nr;
   double scalefactor, U[16], V[16], t;
   double T = pi[0], C = pi[1], A = pi[2], G = pi[3], Y = T + C, R = A + G;

   if (model == JC69 || model == F81) kappa1 = kappa2 = com.kappa = 1;
   else if (com.model < TN93)       kappa2 = kappa1;
   RootTN93(model, kappa1, kappa2, pi, &scalefactor, Root);

   *nR = nr = 2 + (model == K80 || model >= F84) + (model >= HKY85);
   U[0 * 4 + 0] = U[1 * 4 + 0] = U[2 * 4 + 0] = U[3 * 4 + 0] = 1;
   U[0 * 4 + 1] = U[1 * 4 + 1] = 1 / Y;   U[2 * 4 + 1] = U[3 * 4 + 1] = -1 / R;
   U[0 * 4 + 2] = U[1 * 4 + 2] = 0;  U[2 * 4 + 2] = G / R;  U[3 * 4 + 2] = -A / R;
   U[2 * 4 + 3] = U[3 * 4 + 3] = 0;  U[0 * 4 + 3] = C / Y;  U[1 * 4 + 3] = -T / Y;

   xtoy(pi, V, 4);
   V[1 * 4 + 0] = R * T;   V[1 * 4 + 1] = R * C;
   V[1 * 4 + 2] = -Y * A;  V[1 * 4 + 3] = -Y * G;
   V[2 * 4 + 0] = V[2 * 4 + 1] = 0;  V[2 * 4 + 2] = 1;   V[2 * 4 + 3] = -1;
   V[3 * 4 + 0] = 1;  V[3 * 4 + 1] = -1;   V[3 * 4 + 2] = V[3 * 4 + 3] = 0;

   for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
      Cijk[i * 4 * nr + j * nr + 0] = U[i * 4 + 0] * V[0 * 4 + j];
      switch (model) {
      case JC69:
      case F81:
         for (k = 1, t = 0; k < 4; k++) t += U[i * 4 + k] * V[k * 4 + j];
         Cijk[i * 4 * nr + j * nr + 1] = t;
         break;
      case K80:
      case F84:
         Cijk[i * 4 * nr + j * nr + 1] = U[i * 4 + 1] * V[1 * 4 + j];
         for (k = 2, t = 0; k < 4; k++) t += U[i * 4 + k] * V[k * 4 + j];
         Cijk[i * 4 * nr + j * nr + 2] = t;
         break;
      case HKY85:   case T92:   case TN93:
         for (k = 1; k < 4; k++)  Cijk[i * 4 * nr + j * nr + k] = U[i * 4 + k] * V[k * 4 + j];
         break;
      default:
         error2("model in eigenTN93");
      }
   }
#ifdef BASEMLG
   for (i = 0; i < 64; i++)  CijkIs0[i] = (Cijk[i] == 0);
#endif
   return(0);
}


#endif



#if (defined(CODEML) || defined(YN00))

int printfcode(FILE* fout, double fb61[], double space[])
{
   /* space[64*2]
    */
   int i, n = Nsensecodon;

   fprintf(fout, "\nCodon freq.,  x 10000\n");
   zero(space, 64);
   for (i = 0; i < n; i++) space[FROM61[i]] = fb61[i] * 10000;
   printcu(fout, space, com.icode);
   return(0);
}


int printsmaCodon(FILE* fout, char* z[], int ns, int ls, int lline, int simple)
{
   /* print, in blocks, multiple aligned and transformed codon sequences.
    indels removed.
    This is needed as codons are coded 0,1, 2, ..., 60, and
    printsma won't work.
    */
   int ig, ngroup, lt, il, is, i, b, lspname = 30;
   char equal = '.', * pz, c0[4], c[4];

   if (ls == 0) return(1);
   ngroup = (ls - 1) / lline + 1;
   fprintf(fout, "\n");
   for (ig = 0; ig < ngroup; ig++) {
      /* fprintf (fout,"%-8d\n", ig*lline+1); */
      for (is = 0; is < ns; is++) {
         fprintf(fout, "%-*s  ", lspname, com.spname[is]);
         lt = 0;
         for (il = ig * lline, pz = z[is] + il; lt < lline && il < ls; il++, lt++, pz++) {
            b = *pz;
            b = FROM61[b];
            c[0] = (char)(b / 16);
            c[1] = (char)((b % 16) / 4);
            c[2] = (char)(b % 4);
            c[3] = 0;
            for (i = 0; i < 3; i++)
               c[i] = BASEs[(int)c[i]];
            if (is && simple) {
               b = z[0][il];
               b = FROM61[b];
               c0[0] = (char)(b / 16); c0[1] = (char)((b % 16) / 4); c0[2] = (char)(b % 4);
               for (i = 0; i < 3; i++)
                  if (c[i] == BASEs[(int)c0[i]]) c[i] = equal;
            }
            fprintf(fout, "%3s ", c);
         }
         fprintf(fout, "\n");
      }
   }
   return (0);
}


int setmark_61_64(void)
{
   /* This sets two matrices FROM61[], and FROM64[], which translate between two
    codings of codons.  In one coding, codons go from 0, 1, ..., 63 while in
    the other codons range from 0, 1, ..., 61 with the three stop codons removed.
    FROM61[] translates from the 61-state coding to the 64-state coding, while
    FROM64[] translates from the 64-state coding to the 61-state coding.

    This routine also sets up FourFold[4][4], which defines the 4-fold codon
    boxes.
    */
   int i, j, k, * code = GeneticCode[com.icode];
   int c[3], aa0, aa1, by[3] = { 16,4,1 };
   double nSilent, nStop, nRepl;

   Nsensecodon = 0;
   for (i = 0; i < 64; i++) {
      if (code[i] == -1)  FROM64[i] = -1;
      else { FROM61[Nsensecodon] = i; FROM64[i] = Nsensecodon++; }
   }
   com.ncode = Nsensecodon;

   for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) {
      k = i * 16 + j * 4;
      FourFold[i][j] = (code[k] == code[k + 1] && code[k] == code[k + 2] && code[k] == code[k + 3]);
   }

   for (i = 0, nSilent = nStop = nRepl = 0; i < 64; i++) {
      c[0] = i / 16; c[1] = (i / 4) % 4; c[2] = i % 4;
      if ((aa0 = code[i]) == -1) continue;
      for (j = 0; j < 3; j++) for (k = 0; k < 3; k++) {
         aa1 = code[i + ((c[j] + k + 1) % 4 - c[j]) * by[j]];
         if (aa1 == -1)        nStop++;
         else if (aa0 == aa1)  nSilent++;
         else                  nRepl++;
      }
   }
   /*
    printf("\ncode Stop Silent Replace\n");
    printf("%3d (%d)  %6.0f%6.0f%6.0f  %12.6f%12.6f\n",
    com.icode, 64-com.ncode, nStop,nSilent,nRepl,nStop*3/(com.ncode*9),nSilent*3/(com.ncode*9));
    */
   return (0);
}

int DistanceMatNG86(FILE* fout, FILE* fds, FILE* fdn, FILE* ft, double alpha)
{
/* Estimation of dS and dN by the method of Nei & Gojobori (1986)
   This works with both coded (com.cleandata==1) and uncoded data.
   In the latter case (com.cleandata==0), the method does pairwise delection.

   alpha for gamma rates is used for dN only.
   */
   char* codon[2];
   int is, js, k, h, wname = 20, status = 0, ndiff, nsd[4];
   double ns, na, nst = -1, nat = -1, S = -1, N = -1, St, Nt, dS, dN, dN_dS, y, bigD = 3, lst;
   double SEds=-1, SEdn=-1, p;

   if (fout) {
      fputs("\n\n\nNei & Gojobori 1986. dN/dS (dN, dS)", fout);
      if (com.cleandata == 0) fputs("\n(Pairwise deletion)", fout);
      fputs("\n(Note: This matrix is not used in later ML. analysis.\n", fout);
      fputs("Use runmode = -2 for ML pairwise comparison.)\n", fout);
   }

   if (fds) {
      fprintf(fds, "%6d\n", com.ns);
      fprintf(fdn, "%6d\n", com.ns);
      fprintf(ft, "%6d\n", com.ns);
   }
   if (noisy > 1 && com.ns > 10)  puts("NG distances for seqs.:");
   for (is = 0; is < com.ns; is++) {
      if (fout)
         fprintf(fout, "\n%-*s", wname, com.spname[is]);
      if (fds) {
         fprintf(fds, "%-*s ", wname, com.spname[is]);
         fprintf(fdn, "%-*s ", wname, com.spname[is]);
         fprintf(ft, "%-*s ", wname, com.spname[is]);
      }
      for (js = 0; js < is; js++) {
         for (k = 0; k < 4; k++) nsd[k] = 0;
         for (h = 0, lst = 0, nst = nat = S = N = 0; h < com.npatt; h++) {
            if (com.z[is][h] >= com.ncode || com.z[js][h] >= com.ncode)
               continue;
            codon[0] = CODONs[(int)com.z[is][h]];
            codon[1] = CODONs[(int)com.z[js][h]];
            lst += com.fpatt[h];
            ndiff = difcodonNG(codon[0], codon[1], &St, &Nt, &ns, &na, 0, com.icode);
            nsd[ndiff] += (int)com.fpatt[h];
            S += St * com.fpatt[h];
            N += Nt * com.fpatt[h];
            nst += ns * com.fpatt[h];
            nat += na * com.fpatt[h];
         }  /* for(h) */
         if (S <= 0 || N <= 0)
            y = 0;
         else {       /* rescale for stop codons */
            y = lst * 3. / (S + N);
            S *= y;
            N *= y;
         }
         if (noisy >= 9)
            printf("\n%3d %3d:Sites %7.1f +%7.1f =%7.1f\tDiffs %7.1f +%7.1f =%7.1f",
               is + 1, js + 1, S, N, S + N, nst, nat, nst + nat);

         dS = (S <= 0 ? 0 : 1 - 4. / 3 * nst / S);
         dN = (N <= 0 ? 0 : 1 - 4. / 3 * nat / N);
         if (noisy >= 9 && (dS <= 0 || dN <= 0))
         {
            puts("\nNG86 unusable."); status = -1;
         }
         if (dS == 1) dS = 0;
         else      dS = (dS <= 0 ? -1 : 3. / 4 * (-log(dS)));
         if (dN == 1) dN = 0;
         else      dN = (dN <= 0 ? -1 : 3. / 4 * (alpha == 0 ? -log(dN) : alpha * (pow(dN, -1 / alpha) - 1)));

         dN_dS = (dS > 0 && dN > 0 ? dN / dS : -1);
         if (fout) fprintf(fout, "%7.4f (%5.4f %5.4f)", dN_dS, dN, dS);

         if (N > 0 && dN < 0)  dN = bigD;
         if (S > 0 && dS < 0)  dS = bigD;

#ifdef CODEML
         SeqDistance[is * (is - 1) / 2 + js] = (S <= 0 || N <= 0 ? 0 : (S * dS + N * dN) * 3 / (S + N));
#endif

         if (fds) {
            fprintf(fds, " %7.4f", dS);
            fprintf(fdn, " %7.4f", dN);
            fprintf(ft, " %7.4f", (S * dS + N * dN) * 3 / (S + N));
         }
         if (alpha == 0 && dS < bigD) {
            p = nst / S;  SEds = sqrt(9 * p * (1 - p) / (square(3 - 4 * p) * S));
         }
         if (alpha == 0 && dN < bigD) {
            p = nat / N;  SEdn = sqrt(9 * p * (1 - p) / (square(3 - 4 * p) * N));
         }
      }    /* for(js) */
      if (fds) {
         fprintf(fds, "\n"); fprintf(fdn, "\n"); fprintf(ft, "\n");
      }
      if (noisy > 1 && com.ns > 10)  printf(" %3d", is + 1);
   }    /* for(is) */
   printf("\n");
   if (status) fprintf(fout, "NOTE: -1 means that NG86 is inapplicable.\n");

   SS = S, NN = N, Sd = nst, Nd = nat;  /* kostas */

   return (0);
}


#endif



#ifdef BASEML

int eigenQREVbase(FILE* fout, double Q[NCODE * NCODE], double kappa[], double pi[], int* nR, double Root[], double Cijk[])
{
   /* This returns the Q matrix in Q.  pi[] is constant.
   */
   int n = com.ncode, i, j, k;
   int nr = (com.ngene > 1 && com.Mgene >= 3 ? com.nrate / com.ngene : com.nrate);
   double Q0[NCODE * NCODE], U[NCODE * NCODE], V[NCODE * NCODE], mr, space_pisqrt[NCODE * NCODE];

   NPMatUVRoot = 0;
   *nR = n;
   zero(Q, n * n);
   if (com.model == REV) {
      if (n != 4) error2("ncode != 4 for REV");
      Q[3 * n + 2] = Q[2 * n + 3] = 1;  /* r_AG = r_GA = 1. */
      for (i = 0, k = 0; i < n - 1; i++) for (j = i + 1; j < n; j++)
         if (i * n + j != 2 * n + 3)
            Q[i * n + j] = Q[j * n + i] = kappa[k++];
   }
   else       /* (model==REVu) */
      for (i = 0; i < n - 1; i++) for (j = i + 1; j < n; j++)
         Q[i * n + j] = Q[j * n + i] = (StepMatrix[i * n + j] ? kappa[StepMatrix[i * n + j] - 1] : 1);

   for (i = 0; i < n; i++) for (j = 0; j < n; j++)
      Q[i * n + j] *= pi[j];

   for (i = 0, mr = 0; i < n; i++) {
      Q[i * n + i] = 0;
      Q[i * n + i] = -sum(Q + i * n, n);
      mr -= pi[i] * Q[i * n + i];
   }
   abyx(1 / mr, Q, n * n);

   if (fout) {
      mr = 2 * pi[0] * Q[0 * n + 1] + 2 * pi[2] * Q[2 * n + 3];
      if (com.nhomo == 0) {
         fprintf(fout, "\nRate parameters:  ");
         for (j = 0; j < nr; j++)
            fprintf(fout, " %8.5f", kappa[j]);
         fprintf(fout, "\nBase frequencies: ");
         for (j = 0; j < n; j++)
            fprintf(fout, " %8.5f", pi[j]);
      }
      fprintf(fout, "\nRate matrix Q, Average Ts/Tv =%9.4f", mr / (1 - mr));
      matout(fout, Q, n, n);
   }
   else {
      xtoy(Q, Q0, n * n);
      eigenQREV(Q0, pi, n, Root, U, V, space_pisqrt);
      for (i = 0; i < n; i++) for (j = 0; j < n; j++) for (k = 0; k < n; k++)
         Cijk[i * n * n + j * n + k] = U[i * n + k] * V[k * n + j];
   }
   return (0);
}


int QUNREST(FILE* fout, double Q[], double rate[], double pi[])
{
   /* This constructs the rate matrix Q for the unrestricted model.
    pi[] is changed in the routine.
    */
   int n = com.ncode, i, j, k;
   double mr, ts, space[20];

   if (com.model == UNREST) {
      if (n != 4) error2("ncode != 4 for UNREST");
      for (i = 0, k = 0, Q[14] = 1; i < n; i++) for (j = 0; j < n; j++)
         if (i != j && i * n + j != 14)  Q[i * n + j] = rate[k++];
   }
   else  /* (model==UNRESTu) */
      for (i = 0; i < n; i++) for (j = 0; j < n; j++)
         if (i != j)
            Q[i * n + j] = (StepMatrix[i * n + j] ? rate[StepMatrix[i * n + j] - 1] : 1);

   for (i = 0; i < n; i++) {
      Q[i * n + i] = 0;
      Q[i * n + i] = -sum(Q + i * n, n);
   }

   /* get pi */
   QtoPi(Q, com.pi, n, space);

   for (i = 0, mr = 0; i < n; i++)  mr -= pi[i] * Q[i * n + i];
   for (i = 0; i < n * n; i++)  Q[i] /= mr;

   if (fout) {
      ts = pi[0] * Q[0 * n + 1] + pi[1] * Q[1 * n + 0] + pi[2] * Q[2 * n + 3] + pi[3] * Q[3 * n + 2];

      fprintf(fout, "Rate parameters:  ");
      for (j = 0; j < com.nrate; j++)  fprintf(fout, " %8.5f", rate[j]);
      fprintf(fout, "\nBase frequencies: ");
      for (j = 0; j < n; j++) fprintf(fout, " %8.5f", pi[j]);
      if (n == 4)
         fprintf(fout, "\nrate matrix Q, Average Ts/Tv (similar to kappa/2) =%9.4f", ts / (1 - ts));
      else
         fprintf(fout, "\nrate matrix Q");
      matout(fout, Q, n, n);
   }
   return (0);
}

#endif


#ifdef LSDISTANCE

double* SeqDistance = NULL;
int* ancestor = NULL;

int SetAncestor()
{
   /* This finds the most recent common ancestor of species is and js.
    */
   int is, js, it, a1, a2;

   PrintTree(0);

   for (is = 0; is < com.ns; is++) for (js = 0; js < is; js++) {
      it = is * (is - 1) / 2 + js;
      ancestor[it] = -1;
      for (a1 = is; a1 != -1; a1 = nodes[a1].father) {
         for (a2 = js; a2 != -1; a2 = nodes[a2].father)
            if (a1 == a2) { ancestor[it] = a1; break; }
         if (ancestor[it] != -1) break;
      }
      if (ancestor[it] == -1) error2("no ancestor");
   }
   return(0);
}

int fun_LS(double x[], double diff[], int np, int npair);

int fun_LS(double x[], double diff[], int np, int npair)
{
   int i, j, aa, it = -np;
   double dexp;

   if (SetBranch(x) && noisy > 2) puts("branch len.");
   if (npair != com.ns * (com.ns - 1) / 2) error2("# seq pairs err.");
   for (i = 0; i < com.ns; i++) for (j = 0; j < i; j++) {
      it = i * (i - 1) / 2 + j;
      for (aa = i, dexp = 0; aa != ancestor[it]; aa = nodes[aa].father)
         dexp += nodes[aa].branch;
      for (aa = j; aa != ancestor[it]; aa = nodes[aa].father)
         dexp += nodes[aa].branch;
      diff[it] = SeqDistance[it] - dexp;

      if (fabs(diff[it]) > 1000) {
         printf("\ndistances very different: diff = %12.6f ", diff[it]);
      }

   }
   return(0);
}

int LSDistance(double* ss, double x[], int(*testx)(double x[], int np))
{
   /* get Least Squares estimates of branch lengths for a given tree topology
    This uses nls2, a general least squares algorithm for nonlinear programming
    to estimate branch lengths, and it thus inefficient.
    */
   int i;

   if ((*testx)(x, com.ntime)) {
      matout(F0, x, 1, com.ntime);
      puts("initial err in LSDistance()");
   }
   SetAncestor();
   i = nls2((com.ntime > 20 && noisy >= 3 ? F0 : NULL),
      ss, x, com.ntime, fun_LS, NULL, testx, com.ns * (com.ns - 1) / 2, 1e-6);

   return (i);
}

double PairDistanceML(int is, int js)
{
   /* This calculates the ML distance between is and js, the sum of ML branch
    lengths along the path between is and js.
    LSdistance() has to be called once to set ancestor before calling this
    routine.
    */
   int it, a;
   double dij = 0;

   if (is == js) return(0);
   if (is < js) { it = is; is = js; js = it; }

   it = is * (is - 1) / 2 + js;
   for (a = is; a != ancestor[it]; a = nodes[a].father)
      dij += nodes[a].branch;
   for (a = js; a != ancestor[it]; a = nodes[a].father)
      dij += nodes[a].branch;
   return(dij);
}


int GroupDistances()
{
   /* This calculates average group distances by summing over the ML
    branch lengths */
   int newancestor = 0, i, j, ig, jg;
   /*   int ngroup=2, Ningroup[10], group[200]={1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}; */ /* dloop for HC200.paup */
   int ngroup = 10, Ningroup[10], group[115] = {
       10, 9, 9, 9, 9, 9, 9, 9, 9, 10,
       9, 9, 3, 3, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 2, 2, 2, 2, 2, 2, 4, 4, 4,
       4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
       4, 4, 4, 4, 4, 4, 5, 5, 5, 5,
       5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
       5, 5, 5, 5, 6, 6, 6, 6, 6, 6,
       6, 7, 7, 7, 7, 7, 7, 7, 7, 7,
       8, 8, 8, 8, 8 };  /* dloop data for Anne Yoder, ns=115 */
   double dgroup, npairused;

   /* ngroup=2; for(j=0;j<com.ns; j++) group[j]=1+(group[j]>2); */

   for (j = 0; j < ngroup; j++) Ningroup[j] = 0;
   for (j = 0; j < com.ns; j++) Ningroup[group[j] - 1]++;
   printf("\n# sequences in group:");
   matIout(F0, Ningroup, 1, ngroup);
   if (ancestor == NULL) {
      newancestor = 1;
      ancestor = (int*)realloc(ancestor, com.ns * (com.ns - 1) / 2 * sizeof(int));
      if (ancestor == NULL) error2("oom ancestor");
   }
   SetAncestor();

   for (ig = 0; ig < ngroup; ig++) {
      printf("\ngroup %2d", ig + 1);
      for (jg = 0; jg < ig + 1; jg++) {
         dgroup = 0;  npairused = 0;
         for (i = 0; i < com.ns; i++) for (j = 0; j < com.ns; j++) {
            if (i != j && group[i] == ig + 1 && group[j] == jg + 1) {
               dgroup += PairDistanceML(i, j);
               npairused++;
            }
         }
         dgroup /= npairused;
         printf("%9.4f", dgroup);

         /* printf("%6.1f", dgroup/0.2604*5); */ /* 0.2604, 0.5611 */
      }
   }
   if (newancestor == 1)  free(ancestor);
   return(0);
}

#endif

#ifdef NODESTRUCTURE


double CountTrees(int ns, int rooted)
{
   double i, ntree = 1;

   if (ns > 70) error2("ns too large in NumberofTrees().");
   for (i = 4; i <= ns + rooted; i++)
      ntree *= 2 * i - 5;
   return (ntree);
}

double CountLHs(int ns)
{
   double i, nLH = 1;

   if (ns > 70) error2("ns too large in NumberofLHs().");
   for (i = 3; i <= ns; i++)
      nLH *= i * (i - 1) / 2;
   return (nLH);
}



void BranchToNode(void)
{
   /* tree.root need to be specified before calling this
    */
   int i, from, to;

   tree.nnode = tree.nbranch + 1;
   for (i = 0; i < tree.nnode; i++) {
      nodes[i].father = nodes[i].ibranch = -1;  
      nodes[i].nson = 0;
      /* nodes[i].label = -1; */
   }
   for (i = 0; i < tree.nbranch; i++) {
      from = tree.branches[i][0];
      to = tree.branches[i][1];
      nodes[from].sons[nodes[from].nson++] = to;
      nodes[to].father = from;
      nodes[to].ibranch = i;
   }
   /*  nodes[tree.root].branch=0;  this breaks method=1 */
}

void NodeToBranchSub(int inode);

void NodeToBranchSub(int inode)
{
   int i, ison;

   for (i = 0; i < nodes[inode].nson; i++) {
      tree.branches[tree.nbranch][0] = inode;
      tree.branches[tree.nbranch][1] = ison = nodes[inode].sons[i];
      nodes[ison].ibranch = tree.nbranch++;
      if (nodes[ison].nson > 0)  NodeToBranchSub(ison);
   }
}

void NodeToBranch(void)
{
   tree.nbranch = 0;
   NodeToBranchSub(tree.root);
   if (tree.nnode != tree.nbranch + 1)
      error2("nnode != nbranch + 1?");
}


void ClearNode(int inode)
{
   /* a source of confusion. Try not to use this routine.
    */
   nodes[inode].father = nodes[inode].ibranch = -1;
   nodes[inode].nson = 0;
   nodes[inode].branch = nodes[inode].age = 0;
   /* nodes[inode].label = -1; */
   /* nodes[inode].branch = 0; clear node structure only, not branch lengths */
   /* for(i=0; i<com.ns; i++) nodes[inode].sons[i]=-1; */
}

int ReadTreeB(FILE* ftree, int popline)
{
   char line[255];
   int nodemark[NS * 2 - 1] = { 0 }; /* 0: absent; 1: father only (root); 2: son */
   int i, j, state = 0, YoungAncestor = 0;

   if (com.clock) {
      puts("\nbranch representation of tree might not work with clock model.");
      getchar();
   }

   fscanf(ftree, "%d", &tree.nbranch);
   for (j = 0; j < tree.nbranch; j++) {
      for (i = 0; i < 2; i++) {
         if (fscanf(ftree, "%d", &tree.branches[j][i]) != 1) state = -1;
         tree.branches[j][i]--;
         if (tree.branches[j][i]<0 || tree.branches[j][i]>com.ns * 2 - 1)
            error2("ReadTreeB: node numbers out of range");
      }
      nodemark[tree.branches[j][1]] = 2;
      if (nodemark[tree.branches[j][0]] != 2) nodemark[tree.branches[j][0]] = 1;
      if (tree.branches[j][0] < com.ns)  YoungAncestor = 1;

      printf("\nBranch #%3d: %3d -> %3d", j + 1, tree.branches[j][0] + 1, tree.branches[j][1] + 1);

   }
   if (popline) fgets(line, 254, ftree);
   for (i = 0, tree.root = -1; i < tree.nbranch; i++)
      if (nodemark[tree.branches[i][0]] != 2) tree.root = tree.branches[i][0];
   if (tree.root == -1) error2("root err");
   for (i = 0; i < com.ns; i++)
      if (nodemark[i] == 0) {
         matIout(F0, nodemark, 1, com.ns);
         error2("branch specification of tree");
      }

   if (YoungAncestor) {
      puts("\nAncestors in the data?  Take care.");
      if (!com.cleandata) {
         puts("This kind of tree does not work with unclean data.");
         getchar();
      }
   }

   /*
    com.ntime = com.clock ? (tree.nbranch+1)-com.ns+(tree.root<com.ns)
    : tree.nbranch;
    */

   BranchToNode();
   return (state);
}


int OutTreeB(FILE* fout)
{
   int j;
   char* fmt[] = { " %3d..%-3d", " %2d..%-2d" };
   for (j = 0; j < tree.nbranch; j++)
      fprintf(fout, fmt[0], tree.branches[j][0] + 1, tree.branches[j][1] + 1);
   return (0);
}

int GetTreeFileType(FILE* ftree, int* ntree, int* pauptree, int shortform);

int GetTreeFileType(FILE* ftree, int* ntree, int* pauptree, int shortform)
{
   /* paupstart="begin trees" and paupend="translate" identify paup tree files.
    paupch=";" will be the last character before the list of trees.
    Modify if necessary.
    */
   int i, k, lline = 32000, ch = 0, paupch = ';';
   char line[32000];
   char* paupstart = "begin tree", * paupend = "translate";

   *pauptree = 0;
   k = fscanf(ftree, "%d%d", &i, ntree);
   if (k == 2) {
      if (i == com.ns)  return(0);                 /* old paml style */
      else              error2("Number of sequences different in tree and seq files.");
   }
   else if (k == 1) { *ntree = i; return(0); }           /* phylip & molphy style */
   while (ch != '(' && !isalnum(ch) && ch != EOF)
      ch = fgetc(ftree);  /* treeview style */
   if (ch == '(') {
      *ntree = -1;  ungetc(ch, ftree);
      return(0);
   }

   puts("\n# seqs in tree file does not match.  Read as the nexus format.");
   for ( ; ; ) {
      if (fgets(line, lline, ftree) == NULL) error2("tree err1: EOF");
      strcase(line, 0);
      if (strstr(line, paupstart)) { *pauptree = 1; *ntree = -1; break; }
   }
   if (shortform) return(0);
   for ( ; ; ) {
      if (fgets(line, lline, ftree) == NULL) error2("tree err2: EOF");
      strcase(line, 0);
      if (strstr(line, paupend)) break;
   }
   for ( ; ; ) {
      if ((ch = fgetc(ftree)) == EOF) error2("tree err3: EOF");
      if (ch == paupch) break;
   }
   if (fgets(line, lline, ftree) == NULL) error2("tree err4: EOF");

   return(0);
}

int PopPaupTree(FILE* ftree);
int PopPaupTree(FILE* ftree)
{
   /* This reads out the string in front of the tree in the nexus format,
      typically "tree PAUP_1 = [&U]" with "[&U]" optional
   */
   int ch;

   for (; ;) {
      ch = fgetc(ftree);
      if (ch == '(')
      {
         ungetc(ch, ftree); return(0);
      }
      else if (ch == EOF || ch == '/')
         return(-1);
   }
   return(0);
}


void DownTreeCladeLabel(int inode, int label)
{
   /* This goes down the tree to change $ labels (stored in nodes[].label2) into
      # labels (stored in nodes[].label).  To deal with nested clade labels,
      branches within a clade are labeled by negative numbers initially, and
      converted to positive labels at the end of the algorithm.
      nodes[].label and nodes[].label2 are initialized to -1 before this routine is called.
   */
   int i;

   if (nodes[inode].label2 != -1)
      label = (int)nodes[inode].label2;
   if (inode != tree.root && nodes[inode].label == -1)
      nodes[inode].label = label;
   for (i = 0; i < nodes[inode].nson; i++)
      DownTreeCladeLabel(nodes[inode].sons[i], label);
}

int IsNameNumber(char line[])
{
   /* returns 0 if line has species number; 1 if it has name.
    It uses com.ns.
    */
   int isname = 1, alldigits = 1, n;
   char* p = line;

   while (*p)
      if (!isdigit(*p++)) { alldigits = 0; break; }
   if (alldigits) {
      n = atoi(line);
      if (n >= 1 && n <= com.ns) isname = 0;
   }
   return(isname);
}


int ReadUntil(FILE* fin, char line[], char delimiters[], int maxline)
{
   /* read into line[] from fin until any character in delimiters[] or maxline,
      & trim spaces at the end.
   */
   int i;

   for (i = 0; i < maxline; i++) {  /* read notes into line[] */
      line[i] = (char)fgetc(fin);
      if (line[i] == EOF)
         error2("EOF when reading node label");
      if (strchr(delimiters, (int)line[i])) {
         ungetc(line[i], fin);
         line[i] = '\0';
         break;
      }
   }
   for (; i > 0; i--) { /* trim spaces*/
      if (isgraph(line[i]))
         break;
      else
         line[i] = 0;
   }
   return(i + 1);
}



/* Read a newick tree from ftree.
Branch lengths are read in nodes[].branch, and branch (node) labels
(integers) are preceeded by # and read in nodes[].label.  If the clade label
$ is used, the label is read nodes[].label2 first and then moved into
nodes[].label in the routine DownTreeCladeLabel().

Calibration information for mcmctree may be read into nodes[].branch and nodes[].label,
as well as nodes[].NodeStr, and is processed inside mcmctree.
*haslength is set to 1 (branch lengths), 2 (calibration info) or 3 (both).
However, the bit for calibrations is set only if the symbols > < exist and not for
calibrations specified using L, U, G, etc, which will be stored in nodes[].annotation
and processed using ProcessizenodeAnnotation() in mcmctree.

This assumes that com.ns is known.  Names are considered case-sensitive, with trailing spaces ignored.

copyname = 0: species numbers and names are both accepted, but names have to match the names
in com.spname[], which are from the sequence data file.
Used by baseml and codeml, for example.
1: species names are copied into com.spname[], but species numbers are accepted.
Used by evolver for simulation, in which case no species names were read before.
2: the tree must have species names, which are copied into com.spname[].

isname = 0:   species number; 1: species name;
*/
int ReadTreeN(FILE* ftree, int* haslength, int copyname, int popline)
{
   int cnode, cfather = -1;  /* current node and father */
   int inode = 0;  /* node number that will have the next branch length */
   int i, k, level = 0, isname, ch = ' ', tip = 0;
   /* delimiters for names & labels (labels if not in quotes) */
   char check[NS], delimitersN[] = "(),:;#[", delimitersL[] = "(),:;", quote[2][4] = { "\"\'[", "\"\']" };
   int lline = 32000, namelabel = -1;  /* namelabel = 0 for name 1 for label */
   char line[32000];
   char debug = 0;

   /*
   (((A1, B1) ' #1 B(0.1, 0.3)', (A2, B2) #1), C) 'B(0.9, 1.1)';
   ((A)H#H[&gamma = 0:3];H#H)R;
   ((A : 0:02; (B : 0:01)H#H1[&gamma=0:3] : 0:01)S : 0:03; (H#H1 : 0:02;C : 0:03)T : 0:02)R :0:03;
   */
   if (com.ns <= 0)  error2("you should specify # seqs in the tree file.");

   for (i = 0; i < com.ns; i++) check[i] = 0;
   *haslength = 0;
   tree.nnode = com.ns;  tree.nbranch = 0;
   for (i = 0; i < 2 * com.ns - 1; i++) {
      nodes[i].father = nodes[i].ibranch = -1;
      nodes[i].nson = 0;  nodes[i].label = nodes[i].label2 = -1;  nodes[i].branch = 0;
      nodes[i].age = 0;     /* For TipDate models this is set later for each tree. */
      nodes[i].annotation = NULL;
#if (defined(BASEML) || defined(CODEML))
      nodes[i].fossil = 0;
#endif
   }
   while (isspace(ch))
      ch = fgetc(ftree);  /* skip spaces */
   ungetc(ch, ftree);
   if (isdigit(ch)) {
      ReadTreeB(ftree, popline); return(0);
   }
   if (PopPaupTree(ftree) == -1) return(-1);

   for (; ; ) {
      ch = fgetc(ftree);
      if (ch == EOF) return(-1);
      else if (ch == ';') {
         if (level != 0)
            error2("; in treefile");
         else
            break;
      }
      else if (ch == ',') namelabel = 0;
      else if (!isgraph(ch))
         continue;
      else if (ch == '(') {       /* left (  */
         level++;
         namelabel = 0;
         cnode = tree.nnode++;
         if (tree.nnode > 2 * com.ns - 1)
            error2("check #seqs and tree: perhaps too many '('?");
         if (cfather >= 0) {
            if (nodes[cfather].nson >= MAXNSONS) {
               printf("there are at least %d daughter nodes, raise MAXNSONS?", nodes[cfather].nson);
               exit(-1);
            }
            nodes[cfather].sons[nodes[cfather].nson++] = cnode;
            nodes[cnode].father = cfather;
            tree.branches[tree.nbranch][0] = cfather;
            tree.branches[tree.nbranch][1] = cnode;
            nodes[cnode].ibranch = tree.nbranch++;
         }
         else
            tree.root = cnode;
         cfather = cnode;
      }
      /* treating : and > in the same way is risky. */
      else if (ch == ')') {
         level--;  namelabel = 1;  inode = cfather;  cfather = nodes[cfather].father;
      }
      else if (ch == ':') {
         *haslength = 1;
         fscanf(ftree, "%lf", &nodes[inode].branch);
      }
      else if (namelabel == 0) { /* read species name or number for tip */
         if (level <= 0)
            error2("expecting ; in the tree file");
         ungetc(ch, ftree);
         ReadUntil(ftree, line, delimitersN, lline);
         isname = IsNameNumber(line);
         if (isname == 0) {     /* number */
            if (copyname == 2) error2("Use names in tree.");
            sscanf(line, "%d", &cnode);
            cnode--;
         }
         else {                 /* name */
            if (copyname == 0) {
               for (i = 0; i < com.ns; i++)
                  if (!strcmp(line, com.spname[i])) break;
               if ((cnode = i) == com.ns) {
                  printf("\nError: species %s in tree file not found in the sequence alignment file.\n", line);
                  exit(-1);
               }
            }
            else {
               if (tip > com.ns - 1)
                  error2("error in tree: too many species in tree");
               strcpy(com.spname[cnode = tip++], line);
            }
         }
         nodes[cnode].father = cfather;
         if (nodes[cfather].nson >= MAXNSONS)
            error2("too many daughter nodes, raise MAXNSONS");
         nodes[cfather].sons[nodes[cfather].nson++] = cnode;
         tree.branches[tree.nbranch][0] = cfather;
         tree.branches[tree.nbranch][1] = cnode;
         nodes[cnode].ibranch = tree.nbranch++;
         inode = cnode;
         check[cnode]++;
         namelabel = 1;  /* it is possile that a name is followed by a label */
      }
      else {     /* label or annotation (node can be tip or internal node) */
         if (strchr(quote[0], ch)) {
            k = ReadUntil(ftree, line, quote[1], lline);
            ch = fgetc(ftree);   /* pop off ' "  */
         }
         else {
            ungetc(ch, ftree);
            k = ReadUntil(ftree, line, delimitersL, lline);
         }
         nodes[inode].annotation = (char*)malloc((k + 1) * sizeof(char));
         if (nodes[inode].annotation == NULL) error2("oom annotation");
         strcpy(nodes[inode].annotation, line);
         /* This line is used in MCcoal.  ProcessNodeAnnotation is used to process node annotation elsewhere. */
         if (line[0] == '#')        /* branch labels used in baseml and codeml */
            sscanf(line + 1, "%lf", &nodes[inode].label);
         else if (line[0] == '$')   /* clade labels used to label branches in baseml and codeml */
            sscanf(line + 1, "%lf", &nodes[inode].label2);

         else if (line[0] == '@' || line[0] == '=') {
            sscanf(line + 1, "%lf", &nodes[inode].age);
#if (defined(BASEML) || defined(CODEML))
            if (com.clock) nodes[inode].fossil = 1;
#if (defined(CODEML))
            nodes[inode].omega = 0;
#endif
#endif
         }
      }
   }   /* for( ; ; ) */

   if (popline)
      fgets(line, lline, ftree);
   for (i = 0; i < com.ns; i++) {
      if (check[i] > 1) {
         printf("\nSeq/species #%d (%s) occurs more than once in the tree\n", i + 1, com.spname[i]); exit(-1);
      }
      else if (check[i] < 1) {
         printf("\nSeq #%d (%s) is missing in the tree\n", i + 1, com.spname[i]);
         exit(-1);
      }
   }
   if (tree.nbranch > 2 * com.ns - 2) {
      printf("nbranch %d", tree.nbranch); 
      error2("too many branches in tree?");
   }
   if (tree.nnode != tree.nbranch + 1) {
      printf("\nnnode%6d != nbranch%6d + 1\n", tree.nnode, tree.nbranch);
      exit(-1);
   }

   if (debug) PrintTree(1);
   return (0);
}

int OutSubTreeN(FILE* fout, int inode, int spnames, int printopt);

int OutSubTreeN(FILE* fout, int inode, int spnames, int printopt)
{
   int i, dad, nsib;

   if (inode > com.ns * 2 - 1)
      error2("inode large?");
   dad = nodes[inode].father;
   nsib = (inode == tree.root ? 0 : nodes[dad].nson);

   if (inode != tree.root && inode == nodes[dad].sons[0])
      fputc('(', fout);

   for (i = 0; i < nodes[inode].nson; i++)
      OutSubTreeN(fout, nodes[inode].sons[i], spnames, printopt);

   if (nodes[inode].nson == 0) { /* inode is tip */
      if (spnames) {
         if (printopt & PrNodeNum) fprintf(fout, "%d_", inode + 1);
         fprintf(fout, "%s", com.spname[inode]);
      }
      else
         fprintf(fout, "%d", inode + 1);
   }
   if ((printopt & PrNodeNum) && nodes[inode].nson)
      fprintf(fout, " %d ", inode + 1);

#if (defined BPP)
   if ((printopt & PrNodeStr) && inode < com.ns && nodes[inode].label > 0)
      fprintf(fout, " [&theta=%.6f]", nodes[inode].label);
#elif (defined CODEML || defined BASEML)
   if ((printopt & PrLabel) && nodes[inode].label > 0)
      fprintf(fout, " #%.0f", nodes[inode].label);
#else
   if ((printopt & PrLabel) && nodes[inode].label > 0)
      fprintf(fout, " #%.6f", nodes[inode].label);
   /* fprintf(fout, " [&label=%.6f]", nodes[inode].label); */
#endif

   if ((printopt & PrAge) && nodes[inode].age)
      fprintf(fout, " @%.6f", nodes[inode].age);

   /*  Add branch labels to be read by Rod Page's TreeView.  */
#if (defined CODEML)
   if ((printopt & PrOmega) && inode != tree.root)
      fprintf(fout, " #%.6g ", nodes[inode].omega);
#elif (defined (EVOLVER) || defined (MCMCTREE) || defined (BPP))
   if ((printopt & PrNodeStr) && inode >= com.ns && nodes[inode].annotation)
      fprintf(fout, " %s", nodes[inode].annotation);
#endif

   if ((printopt & PrBranch) && (inode != tree.root || nodes[inode].branch > 0))
      fprintf(fout, ": %.6f", nodes[inode].branch);
   if (nsib == 0)            /* root */
      fputc(';', fout);
   else if (inode == nodes[dad].sons[nsib - 1])  /* last sib */
      fputc(')', fout);
   else                     /* not last sib */
      fprintf(fout, ", ");

   return (0);
}


int OutTreeN(FILE* fout, int spnames, int printopt)
{
   OutSubTreeN(fout, tree.root, spnames, printopt);
   return(0);
}


int DeRoot(void)
{
   /* This cnages the bifurcation at the root into a trifurcation, but setting one of
    the sons to be the new root.  The new root is the first son that is not a tip.
    tree.nnode is updated, but the routine does not re-number the nodes, so the new
    node labels do not go from ns, ns + 1, ..., as they normally should.
    */
   int i, ison, sib, root = tree.root;

   if (nodes[root].nson != 2) error2("in DeRoot?");

   ison = nodes[root].sons[i = 0];
   if (nodes[ison].nson == 0)
      ison = nodes[root].sons[i = 1];
   sib = nodes[root].sons[1 - i];
   nodes[sib].branch += nodes[ison].branch;
   nodes[sib].father = tree.root = ison;
   nodes[tree.root].father = -1;
   nodes[tree.root].sons[nodes[tree.root].nson++] = sib;  /* sib added as the last child of the new root */
   nodes[tree.root].branch = 0;
   tree.nnode--;  /* added 2007/4/9 */
   return(0);
}

int PruneSubTreeN(int inode, int keep[])
{
/* This prunes tips from the tree, using keep[com.ns].  
   Removed nodes in the big tree has nodes[].father=-1 and nodes[].nson=0.
   If merged branches have the same label, the label is kept.  Otherwise 
   the node label is set to -1.
   The code changes nodes[inode].nson and nodes[inode].sons[] only after 
   all inode's descendent nodes are processed. 
   */
   int i, j, ison, father = nodes[inode].father, nson0 = nodes[inode].nson;

   for (i = 0; i < nson0; i++) {
      PruneSubTreeN(nodes[inode].sons[i], keep);
   }
   /* remove inode because of no descendents.
    Note this does not touch the father node */
   if (inode < com.ns && keep[inode] == 0)  /* inode is tip */
      nodes[inode].father = -1;
   else if (inode >= com.ns) {              /* inode is inner node */
      for (i = 0, nodes[inode].nson = 0; i < nson0; i++) {
         ison = nodes[inode].sons[i];
         if (nodes[ison].father != -1)
            nodes[inode].sons[nodes[inode].nson++] = ison;
      }
      if (nodes[inode].nson == 0)
         nodes[inode].father = -1;
   }

   /* remove inode if it has a single descendent ison */
   if (inode >= com.ns && nodes[inode].nson == 1 && inode != tree.root) {
      ison = nodes[inode].sons[0];
      nodes[ison].father = father;
      nodes[ison].branch += nodes[inode].branch;
      if (nodes[ison].label != nodes[inode].label)
         nodes[ison].label = -1;
      for (j = 0; j < nodes[father].nson; j++) {
         if (nodes[father].sons[j] == inode) {
            nodes[father].sons[j] = ison;  break;
         }
      }
      nodes[inode].nson = 0;
      nodes[inode].father = -1;
   }
   else if (nodes[inode].nson == 1 && inode == tree.root) { /* move down root if root has 1 descendent */
      nodes[inode].father = -1;
      nodes[inode].nson = 0;
      ison = nodes[tree.root].sons[0];
      tree.root = ison;
      nodes[tree.root].father = -1;
      nodes[tree.root].branch = 0;
   }

   /*
    printf("\nVisiting inode %d\n", inode);
    for(i=0; i<tree.nnode; i++) printf(" %2d", i);
    printf("\n);
    */
   return(0);
}


int GetSubTreeN(int keep[], int newnodeNO[])
{
   /* This removes some tips to generate the subtree.  Branch lengths are
    preserved by summing them up when some nodes are removed.
    The algorithm use post-order tree traversal to remove tips and nodes.  It
    then switches to the branch representation to renumber nodes.
    space[] can be NULL.  If not, it returns newnodeNO[], which holds the
    new node numbers; for exmaple, newnodeNO[12]=5 means that old node 12 now
    becomes node 5.

    The routine does not change com.ns or com.spname[], which have to be updated
    outside.

    CHANGE OF ROOT happens if the root in the old tree had >=3 sons, but has 2
    sons in the new tree and if (!com.clock).  In that case, the tree is derooted.

    This routine does not work if a current seq is ancestral to some others
    and if that sequence is removed. (***check this comment ***)

    Different formats for keep[] are used.  Suppose the current tree is for
    nine species: a b c d e f g h i.

    (A) keep[]={1,0,1,1,1,0,0,1,0} means that a c d e h are kept in the tree.
    The old tip numbers are not changed, so that OutTreeN(?,1,?) gives the
    correct species names or OutTreeN(?,0,?) gives the old species numbers.

    (B) keep[]={1,0,2,3,4,0,0,5,0} means that a c d e h are kept in the tree, and
    they are renumbered 0 1 2 3 4 and all the internal nodes are renumbered
    as well to be consecutive.  The positive numbers are consecutive natural numbers.

    keep[]={5,0,2,1,4,0,0,3,0} means that a c d e h are kept in the tree.
    However, the order of the sequences are changed to d c h e a, so that the
    numbers are now 0 1 2 3 4 for d c h e a.  This is useful when the subtree
    is extracted from a big tree for a subset of the sequence data, while the
    species are odered d c h e a in the sequence data file.
    This option can be used to renumber the tips in the complete tree.
    */
   int nsnew, i, j, k, nnode0 = tree.nnode, sumnumber = 0, ison, sib;
   int unrooted; /* tree is rooted is the root in main tree has 2 sons or if nsnew=2.*/
   struct TREEN* nodestmp;
   char* spnames0[NS];
   int debug = 0;

   if (debug) {
      for (i = 0; i < com.ns; i++)
         printf("%-30s %2d\n", com.spname[i], keep[i]);
   }
   for (i = 0; i < com.ns; i++)  spnames0[i] = com.spname[i];

   for (i = 0, nsnew = 0; i < com.ns; i++)
      if (keep[i]) { nsnew++; sumnumber += keep[i]; }
   if (nsnew < 2)  return(-1);

   unrooted = (nsnew > 2 && nodes[tree.root].nson >= 3);
   /* mark removed nodes in the big tree by father=-1 && nson=0.
    nodes[].label records the number of nodes collapsed. 
    */
   PruneSubTreeN(tree.root, keep);
   /* If unrooted tree has a bifurcation at the new root, collapse root.  */
   if (unrooted && nodes[tree.root].nson == 2) {
      /* find ison that is inner node, and merge branches ison and sib */
      ison = nodes[tree.root].sons[i = 0];
      if (nodes[ison].nson == 0)
         ison = nodes[tree.root].sons[i = 1];
      sib = nodes[tree.root].sons[1 - i];

      nodes[sib].branch += nodes[ison].branch;
      if(nodes[sib].label != nodes[ison].label) 
         nodes[sib].label = -1;
      nodes[sib].father = tree.root = ison;
      nodes[tree.root].father = -1;
      nodes[tree.root].sons[nodes[tree.root].nson++] = sib;  /* sib added as the last child of the new root */
      nodes[tree.root].branch = 0;
   }
   if (debug) PrintTree(1);

   for (i = 0, k = 1; i < tree.nnode; i++) if (nodes[i].father != -1) k++;
   tree.nnode = k;

   nodestmp = (struct TREEN*)malloc(nnode0 * sizeof(struct TREEN));
   if (nodestmp == NULL) error2("oom GetSubTreeN");

   if (sumnumber > nsnew) {
      /* to renumber the nodes.  k is current node number */
      if (sumnumber != nsnew * (nsnew + 1) / 2)
         error2("keep[] not right in GetSubTreeN");
      memmove(nodestmp, nodes, nnode0 * sizeof(struct TREEN));

      for (i = 0; i < nnode0; i++) newnodeNO[i] = -1;
      for (i = 0; i < com.ns; i++) if (keep[i]) newnodeNO[i] = keep[i] - 1;
      newnodeNO[tree.root] = k = nsnew;
      tree.root = k++;
      for (; i < nnode0; i++)
         if (nodes[i].father > -1) newnodeNO[i] = k++;
     
      /* construct nodes using newnodeNO: old node i becomes new node k */
      for (i = 0; i < nnode0; i++) {
         k = newnodeNO[i];  /* k < i? */
         if (k == -1) continue;
         nodes[k].father = (nodestmp[i].father > 0 ? newnodeNO[nodestmp[i].father] : -1);
         nodes[k].nson = nodestmp[i].nson;
         for(j=0; j<nodes[k].nson; j++)
            nodes[k].sons[j] = newnodeNO[ nodestmp[i].sons[j] ];
         nodes[k].branch = nodestmp[i].branch;
         nodes[k].label = nodestmp[i].label;
         nodes[k].annotation = nodestmp[i].annotation;
         nodes[k].age = nodestmp[i].age;
         nodes[k].name = nodestmp[i].name;
         com.spname[k] = spnames0[i];
         nodes[k].ibranch = -1;
      }
      if (debug) OutTreeN(stdout, 1, PrLabel);
   }

   free(nodestmp);
   return (0);
}


void PrintTree(int timebranches)
{
   int i, j;

   printf("\nns = %d  nnode = %d", com.ns, tree.nnode);
   printf("\n%6s%6s", "dad", "node");
   if (timebranches)  printf("%10s%10s%10s", "age", "branch", "label");
   printf(" %7s%7s", "nson:", "sons");
   for (i = 0; i < tree.nnode; i++) {
      printf("\n%6d%6d", nodes[i].father, i);
      if (timebranches)
         printf(" %9.6f %9.6f %9.0f", nodes[i].age, nodes[i].branch, nodes[i].label);

      printf("%7d: ", nodes[i].nson);
      for (j = 0; j < nodes[i].nson; j++)  printf(" %2d", nodes[i].sons[j]);
      if (nodes[i].annotation) printf(" annotation: |%s| ", nodes[i].annotation);
   }
   printf("\n");
   OutTreeN(F0, 0, 0); printf("\n");
   OutTreeN(F0, 1, 0); printf("\n");
   OutTreeN(F0, 1, 1); printf("\n");
}


void PointconPnodes(void)
{
   /* This points the nodes[com.ns+inode].conP to the right space in com.conP.
      The space is different depending on com.cleandata (0 or 1)
      This routine updates internal nodes com.conP only.
      End nodes (com.conP0) are updated in InitConditionalPNode().
   */
   int nintern = 0, i;

   for (i = 0; i < tree.nbranch + 1; i++)
      if (nodes[i].nson > 0)  /* more thinking */
         nodes[i].conP = com.conP + (size_t)com.ncode * com.npatt * nintern++;
}


int SetxInitials(int np, double x[], double xb[][2])
{
   /* This forces initial values into the boundary of the space
    */
   int i;

   for (i = com.ntime; i < np; i++) {
      if (x[i] < xb[i][0] * 1.005) x[i] = xb[i][0] * 1.05;
      if (x[i] > xb[i][1] / 1.005) x[i] = xb[i][1] / 1.05;
   }
   for (i = 0; i < com.np; i++) {
      if (x[i] < xb[i][0]) x[i] = xb[i][0] * 1.2;
      if (x[i] > xb[i][1]) x[i] = xb[i][1] * 0.8;
   }
   return(0);
}


#if(defined(BASEML) || defined(CODEML) || defined(MCMCTREE))
int GetTipDate(void)
{
   /* This scans species/sequence names to collect the sampling dates.  The last field of
      the name contains the date.
      Divergence times are rescaled by using TipDate_TimeUnit.
   */
   int i, j, indate, ndates = 0;
   double young = -1, old = -1;     /* sample dates are between young and old */
   char* p;
   struct tm sampletime = { 0 };    /* number of seconds from 1970-01-01 */
   int year, month, day;  /* yyyy-mm-dd, example: 2020-02-10 */

   tipdate.youngest = 0;
   for (i = 0, ndates = 0; i < stree.nspecies; i++) {
      stree.nodes[i].age = 0;
      j = (int)strlen(com.spname[i]);
      for (indate = 0, p = com.spname[i] + j - 1; j >= 0; j--, p--) {
         if (isdigit(*p) || *p == '.' || *p == '-') indate = 1;
         else if (indate)
            break;
      }
      if (strchr(p, '-')) {
         if (sscanf(p + 1, "%d-%d-%d", &year, &month, &day) == EOF) {
            printf("sequence name: %s\n", com.spname[i]);
            error2("date format wrong for tipdate analysis?");
         }
         tipdate.ymd = 1;
         sampletime.tm_year = year - 1900;
         sampletime.tm_mon = (month > 0 ? month - 1 : 0);
         sampletime.tm_mday = (day > 0 ? day : 15);           /* mid of month if day is missing */
         stree.nodes[i].age = mktime(&sampletime) / (24 * 3600.0);  /* days since 1970-01-01 */
      }
      else
         sscanf(p + 1, "%lf", &stree.nodes[i].age);
      if (stree.nodes[i].age <= 0)
         error2("Tip date <= 0: somehow i am using positive numbers only for tip date");
      else
         ndates++;

      if (i == 0)
         young = old = stree.nodes[i].age;
      else {
         old = min2(old, stree.nodes[i].age);
         young = max2(young, stree.nodes[i].age);
      }
   }
   if (ndates == 0) {
      if (tipdate.timeunit == -1) tipdate.timeunit = 1;
      return(0);
   }
   if (ndates != stree.nspecies)
      printf("TipDate model requires sampling date for each sequence.");

   // if (ymd) *TipDate = 2;  /* flag to change time back to yyyy-mm-dd */
   tipdate.youngest = young;
   if (tipdate.timeunit <= 0)
      tipdate.timeunit = (young - old) * 2.5;
   if (young - old < 1e-100)
      error2("TipDate: all sequences are of the same age?");
   for (i = 0; i < stree.nspecies * 2 - 1; i++) {
      if (i < stree.nspecies || nodes[i].fossil) {
         stree.nodes[i].age = (young - stree.nodes[i].age) / tipdate.timeunit;
         if (stree.nodes[i].age < 1e-100) stree.nodes[i].age = 0;
      }
   }

   if (noisy)
      printf("\nTipDate model\nDate range: (%.6f, %.6f) => (0, %.6f). TimeUnit = %.6f.\n",
         young, old, (young - old) / tipdate.timeunit, tipdate.timeunit);

   return(0);
}

#endif


#if(defined(BASEML) || defined(CODEML))

double* AgeLow = NULL;
int NFossils = 0, AbsoluteRate = 0;
/* TipDate models:
 MutationRate = mut/TipDate_TimeUnit;
 age = age*TipDate_TimeUnit
 */

void SetAge(int inode, double x[]);
void GetAgeLow(int inode);
/* number of internal node times, usd to deal with known ancestors.  Broken? */
static int innode_time = 0;

/* Ziheng Yang, 25 January 2003
 The following routines deal with clock and local clock models, including
 the TipDate model (Rambaut 2000 Bioinformatics 16:395-399; Yoder & Yang 2000
 Mol Biol Evol 17:1081-1090; Yang & Yoder 2003 Syst Biol).
 The tree is rooted.  The routine SetAge assumes that ancestral nodes are
 arranged in the increasing order and so works only if the input tree uses
 the parenthesis notation and not the branch notation.  The option of known
 ancestors is probably broken.

 The flag AbsoluteRate=1 if(TipDate || NFossils).  This could be removed
 as the flags TipDate and NFossils are sufficient.
 clock = 1: global clock, deals with TipDate with no or many fossils,
            ignores branch rates (#) in tree if any.
       = 2: local clock models, as above, but requires branch rates # in tree.
       = 3: as 2, but requires Mgene and option G in sequence file.

 Order of variables in x[]: divergence times, rates for branches, rgene, ...
 In the following ngene=4, com.nbtype=3, with r_ij to be the rate of gene i and branch class j.

 clock=1 or 2:
       [times, r00(if absolute) r01 r02  rgene1 rgene2 rgene3]
       NOTE: rgene[] has relative rates
 clock=3:
       [times, r00(if absolute) r01 r02  r11 r12  r21 r22 r31 r32 rgene1 rgene2 rgene3]
       NOTE: rgene1=r10, rgene2=r20, rgene3=r30

 If(nodes[tree.root].fossil==0) x[0] has absolute time for the root.
 Otherwise x[0] has proportional ages.
 */


double GetBranchRate(int igene, int ibrate, double x[], int* ix)
{
   /* This finds the right branch rate in x[].  The rate is absolute if AbsoluteRate.
    ibrate=0,1,..., indicates the branch rate class.
    This routine is used in the likeihood calculation and in formatting output.
    ix (k) has the position in x[] for the branch rate if the rate is a parameter.
    and is -1 if the rate is not a parameter in the ML iteration.  This is
    for printing SEs.
    */
   int nage = tree.nnode - com.ns - NFossils, k = nage + AbsoluteRate;
   double rate00 = (AbsoluteRate ? x[nage] : 1), brate = rate00;

   if (igene == 0 && ibrate == 0)
      k = (AbsoluteRate ? nage : -1);
   else if (com.clock == GlobalClock) {
      brate = x[k = com.ntime + igene - 1];  /* igene>0, rgene[] has absolute rates */
   }
   else if (com.clock == LocalClock) {  /* rgene[] has relative rates */
      if (igene == 0 && ibrate) { brate = x[k += ibrate - 1]; }
      else if (igene && ibrate == 0) { brate = rate00 * x[com.ntime + igene - 1]; k = -1; }
      else if (igene && ibrate) { brate = x[k + ibrate - 1] * x[com.ntime + igene - 1]; k = -1; }
   }
   else if (com.clock == ClockCombined) {
      if (ibrate == 0 && igene)  brate = x[k = com.ntime + igene - 1];
      else                    brate = x[k += ibrate - 1 + igene * (com.nbtype - 1)]; /* ibrate>0 */
   }

   if (ix) *ix = k;
   return(brate);
}


void SetAge(int inode, double x[])
{
   /* This is called from SetBranch(), to set up age for nodes under clock
    models (clock=1,2,3).
    if(TipDate||NFossil), that is, if(AbsoluteRate), this routine sets up
    times (nodes[].age) and then SetBranch() sets up branch lengths by
    multiplying times with rate:
    [].age[i] = AgeLow[i] + ([father].age - AgeLow[i])*x[i]

    The routine assumes that times are arranged in the order of node numbers,
    and should work if parenthesis notation of tree is used in the tree file,
    but not if the branch notation is used.
    */
   int i, ison;

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      if (nodes[ison].nson) {
         if (AbsoluteRate) {
            if (!nodes[ison].fossil)
               nodes[ison].age = AgeLow[ison]
               + (nodes[inode].age - AgeLow[ison]) * x[innode_time++];
         }
         else
            nodes[ison].age = nodes[inode].age * x[innode_time++];
         SetAge(ison, x);
      }
   }
}

void GetAgeLow(int inode)
{
   /* This sets AgeLow[], the minimum age of each node.  It moves down the tree to
    scan [].age, which has tip dates and fossil dates.  It is needed if(AbsoluteRate)
    and is called by GetInitialsTimes().
    */
   int i, ison;
   double tlow = 0;

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      if (nodes[ison].nson)
         GetAgeLow(ison);
      tlow = max2(tlow, nodes[ison].age);
   }
   if (nodes[inode].fossil) {
      if (nodes[inode].age < tlow)
         error2("age in tree is in conflict.");
      AgeLow[inode] = nodes[inode].age;
   }
   else
      AgeLow[inode] = nodes[inode].age = tlow;
}



int SetBranch(double x[])
{
   /* if(AbsoluteRate), mutation rate is not multiplied here, but during the
    likelihood calculation.  It is copied into com.rgene[0].
    */
   int i, status = 0;
   double smallv = -1e-5;

   if (com.fix_blength == 3) {  /* proportional branch lengths */
      for (i = 0; i < tree.nnode; i++) {
         if (i != tree.root)
            nodes[i].branch = com.blengths0[nodes[i].ibranch] * x[0];
      }
      return(0);
   }

   if (com.clock == 0) {
      for (i = 0; i < tree.nnode; i++) {
         if (i != tree.root)
            if ((nodes[i].branch = x[nodes[i].ibranch]) < smallv)  status = -1;
      }
      return(status);
   }
   innode_time = 0;
   if (!LASTROUND) { /* transformed variables (proportions) are used */
      if (!nodes[tree.root].fossil) /* note order of times in x[] */
         nodes[tree.root].age = x[innode_time++];
      SetAge(tree.root, x);
   }
   else {           /* times are used */
      for (i = com.ns; i < tree.nnode; i++)
         if (!nodes[i].fossil) nodes[i].age = x[innode_time++];
   }

   for (i = 0; i < tree.nnode; i++) {  /* [].age to [].branch */
      if (i == tree.root) continue;
      nodes[i].branch = nodes[nodes[i].father].age - nodes[i].age;
      if (nodes[i].branch < smallv)
         status = -1;
   }
   return(status);
}


int GetInitialsTimes(double x[])
{
   /* this counts com.ntime and initializes x[] under clock and local clock models,
    including TipDate and ClockCombined models.  See above for notes.
    Under local clock models, com.ntime includes both times and rates for
    lineages.
    A recursive algorithm is used to specify initials if(TipDate||NFossil).
    */
   int i, j, k;
   double maxage, t;

   if (com.fix_blength == 2) {
      com.ntime = 0; com.method = 0; return(0);
   }
   if (com.fix_blength == 3) {
      com.ntime = 1; com.method = 0;
      for (i = 0; i < tree.nnode; i++) {
         if (i != tree.root)
            com.blengths0[nodes[i].ibranch] = nodes[i].branch;
      }
      x[0] = 0.5 + rndu();
      return(0);
   }

   /* no clock */
   if (com.clock == 0) {
      com.ntime = tree.nbranch;
      if (com.fix_blength == 1)  return(0);
      for (i = 0; i < com.ntime; i++)
         x[i] = rndu() * 0.1 + 0.01;

      if (com.fix_blength == 0 && com.clock < 5 && ancestor && com.ntime < 100)
         LSDistance(&t, x, testx);

      return(0);
   }

   /* clock models: check branch rate labels and fossil dates first */
   if (com.clock < 5) {
      com.nbtype = 1;
      if (com.clock == 1)
         for (i = 0; i < tree.nnode; i++) nodes[i].label = 0;
      else {
         for (i = 0; i < tree.nnode; i++) {
            if (i != tree.root && (j = (int)nodes[i].label + 1) > com.nbtype) {
               com.nbtype = j;
               if (j<0 || j>tree.nbranch - 1) error2("branch label in the tree.");
            }
         }
         for (j = 0; j < com.nbtype; j++) {
            for (i = 0; i < tree.nnode; i++)
               if (i != tree.root && j == (int)nodes[i].label) break;
            if (i == tree.nnode)
               printf("\nNot all branch labels (0, ..., %d) are found on tree?", com.nbtype - 1);
         }
         if (noisy) printf("\nfound %d branch rates in tree.\n", com.nbtype);
         if (com.nbtype <= 1) error2("use clock = 1 or add branch rate labels in tree");

         for (i = 0; i < tree.nbranch; i++)
            printf("%3.0f", nodes[tree.branches[i][1]].label);
         printf("\n");
      }
   }
   for (i = 0, NFossils = 0, maxage = 0; i < tree.nnode; i++) {
      if (nodes[i].nson && nodes[i].fossil) {
         NFossils++;
         maxage = max2(maxage, nodes[i].age);
      }
   }
   if (NFossils && maxage > 10)
      error2("Change time unit so that fossil dates fall in (0.00001, 10).");

   if (tipdate.flag)
      GetTipDate();

   AbsoluteRate = (tipdate.flag || NFossils);
   if (com.clock >= 5 && AbsoluteRate == 0)
      error2("needs fossil calibrations");

   com.ntime = AbsoluteRate + (tree.nnode - com.ns - NFossils) + (com.nbtype - 1);
   if (com.clock == ClockCombined)
      com.ntime += (com.ngene - 1) * (com.nbtype - 1);
   com.ntime += (tree.root < com.ns); /* root is a known sequence. Broken? */

   /* DANGER! AgeLow is not freed in the program. Fix this? */
   k = 0;
   if (AbsoluteRate) {
      AgeLow = (double*)realloc(AgeLow, tree.nnode * sizeof(double));
      GetAgeLow(tree.root);
   }
   if (!nodes[tree.root].fossil)
      x[k++] = (AbsoluteRate ? nodes[tree.root].age * (1.2 + rndu()) : rndu() * .5 + .1);  /* root age */
   for (; k < tree.nnode - com.ns - NFossils; k++)   /* relative times */
      x[k] = 0.4 + .5 * rndu();
   if (com.clock != 6)                           /* branch rates */
      for (; k < com.ntime; k++)
         x[k] = 0.1 * (.5 + rndu());
   else
      for (j = 0, k = com.ntime - 1; j < data.ngene; j++, k++)
         x[k] = 0.1 * (.5 + rndu());
   return(0);
}

int OutputTimesRates(FILE* fout, double x[], double var[])
{
   /* SetBranch() has been called before calling this, so that [].age is up
    to date.
    */
   int i, j, k = AbsoluteRate + tree.nnode - com.ns - NFossils;
   double scale = (tipdate.flag ? tipdate.timeunit : 1);

   /* rates */
   if (AbsoluteRate && com.clock < 5) {
      fputs("\nSubstitution rate is per time unit\n", fout);
      if (com.nbtype > 1) fprintf(fout, "Rates for branch groups\n");
      fprintf(fout, "\n");
      for (i = 0; i < com.ngene; i++) {
         if (com.ngene > 1) fprintf(fout, "Gene %2d: ", i + 1);
         for (j = 0; j < com.nbtype; j++) {
            fprintf(fout, "%12.6f", GetBranchRate(i, j, x, &k));
            if (i == 0 && j == 0 && !AbsoluteRate) continue;
            if ((com.clock != LocalClock || com.ngene == 1) && com.getSE) {
               if (k == -1) error2("we are in trouble. k should not be -1 here.");
               fprintf(fout, " +- %8.6f", sqrt(var[k * com.np + k]));
            }
         }
      }
   }
   else
      if (com.clock == 2) {
         fprintf(fout, "rates for branches:    1");
         for (k = tree.nnode - com.ns; k < com.ntime; k++) fprintf(fout, " %8.5f", x[k]);
      }

   /* times */
   if (AbsoluteRate) {
      fputs("\nNodes and Times\n\n", fout);
   }
   if (tipdate.flag) { /* DANGER! SE not printed if(TipDate && NFossil). */
      struct tm* sampletime; /* yyyy-mm-dd, example: 2020-02-10 */
      time_t seconds = 0;    /* seconds since 1970 */
      fprintf(fout, "\n");
      for (i = 0, k = 0; i < tree.nnode; i++) {
         fprintf(fout, "Node %3d Time %7.2f ", i + 1, tipdate.youngest - nodes[i].age * scale);
         if (tipdate.ymd) {  /* original sequence names have yyyy-mm-dd format */
            nodes[i].age = (tipdate.youngest - nodes[i].age * tipdate.timeunit);
            seconds = (int)(nodes[i].age * 24 * 3600);
            sampletime = gmtime(&seconds);
            fprintf(fout, " (%04d-%02d-%02d)", sampletime->tm_year + 1900, sampletime->tm_mon + 1, sampletime->tm_mday);
            if (i < com.ns) fprintf(fout, "  %-20s", com.spname[i]);
         }
         if (com.getSE && i >= com.ns && !nodes[i].fossil) {
            fprintf(fout, " +- %6.2f", sqrt(var[k * com.np + k]) * scale);
            k++;
         }
      }
   }
   else if (AbsoluteRate) {
      for (i = com.ns, k = 0; i < tree.nnode; i++) {
         fprintf(fout, "Node %3d Time %9.5f", i + 1, nodes[i].age);
         if (com.getSE && i >= com.ns && !nodes[i].fossil) {
            fprintf(fout, " +- %7.5f", sqrt(var[k * com.np + k]));
            if (fabs(nodes[i].age - x[k]) > 1e-5) error2("node order wrong.");
            k++;
            fprintf(fout, "\n");
         }
      }
   }
   return(0);
}


int SetxBoundTimes(double xb[][2])
{
   /* This sets bounds for times (or branch lengths) and branch rates
    */
   int i = -1, j, k;
   double tb[] = { 4e-6,50 }, rateb[] = { 1e-5, 999 }, pb[] = { .000001, .999999 };

   if (com.fix_blength == 3) {
      xb[0][0] = 1e-4;
      xb[0][1] = 9999;
   }
   if (com.clock == 0) {
      for (i = 0; i < com.ntime; i++) {
         xb[i][0] = tb[0];
         xb[i][1] = tb[1];
      }
   }
   else {
      k = 0;  xb[0][0] = tb[0];  xb[0][1] = tb[1];
      if (!nodes[tree.root].fossil) {
         if (AbsoluteRate)  xb[0][0] = AgeLow[tree.root];
         k = 1;
      }
      for (; k < tree.nnode - com.ns - NFossils; k++)  /* proportional ages */
      {
         xb[k][0] = pb[0]; xb[k][1] = pb[1];
      }
      for (; k < com.ntime; k++)                    /* rate and branch rates */
         for (j = 0; j < 2; j++) xb[k][j] = rateb[j];
   }
   return(0);
}

#endif


#if(defined(BASEML) || defined(BASEMLG) || defined(CODEML))


int readx(double x[], int* fromfile)
{
   /* this reads parameters from file, used as initial values
    if(runmode>0), this reads common substitution parameters only into x[], which
    should be copied into another place before heuristic tree search.  This is broken
    right now.  Ziheng, 9 July 2003.
    fromfile=0: if nothing read from file, 1: read from file, -1:fix parameters
    */
   static int times = 0;
   int i, npin;
   double* xin;

   times++;  *fromfile = 0;
   if (finitials == NULL || (com.runmode > 0 && times > 1)) return(0);
   if (com.runmode <= 0) { npin = com.np; xin = x; }
   else { npin = com.np - com.ntime; xin = x + com.ntime; }

   if (npin <= 0) return(0);
   if (com.runmode > 0 && com.seqtype == 1 && com.model) error2("option or in.codeml");
   printf("\nReading initials/paras from file (np=%d). Stop if wrong.\n", npin);
   fscanf(finitials, "%lf", &xin[i = 0]);
   *fromfile = 1;
   if (xin[0] == -1) { *fromfile = -1; LASTROUND = 1; }
   else           i++;
   for (; i < npin; i++)
      if (fscanf(finitials, "%lf", &xin[i]) != 1) break;
   if (i < npin)
   {
      printf("err at #%d. Edit or remove it.\n", i + 1); exit(-1);
   }
   if (com.runmode > 0) {
      matout(F0, xin, 1, npin);
      puts("Those are fixed for tree search.  Stop if wrong.");
   }
   return(0);
}

#endif

#if(defined(BASEML) || defined(CODEML))

int CollapsNode(int inode, double x[])
{
   /* Merge inode to its father. Update the first com.ntime elments of
    x[] only if (x!=NULL), by using either x[] if clock=1 or
    nodes[].branch if clock=0.  So when clock=0, the routine works
    properly only if SetBranch() is called before this routine, which
    is true if m.l. or l.s. has been used to estimate branch lengths.
    */
   int i, j, ifather, ibranch, ison;

   if (inode == tree.root || inode < com.ns) error2("err CollapsNode");
   ibranch = nodes[inode].ibranch;   ifather = nodes[inode].father;
   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      tree.branches[nodes[ison].ibranch][0] = ifather;
   }
   for (i = ibranch + 1; i < tree.nbranch; i++)
      for (j = 0; j < 2; j++) tree.branches[i - 1][j] = tree.branches[i][j];
   tree.nbranch--; com.ntime--;
   for (i = 0; i < tree.nbranch; i++)  for (j = 0; j < 2; j++)
      if (tree.branches[i][j] > inode)  tree.branches[i][j]--;
   BranchToNode();

   if (x) {
      if (com.clock)
         for (i = inode + 1; i < tree.nnode + 1; i++) x[i - 1 - com.ns] = x[i - com.ns];
      else {
         for (i = ibranch + 1; i < tree.nbranch + 1; i++)  x[i - 1] = x[i];
         SetBranch(x);
      }
   }
   return (0);
}

#endif


#if(defined(BPP) || defined(EVOLVER))

void Tree2PartitionDescentTree(int inode, char split[])
{
   int i, ison;

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      if (ison < com.ns)
         split[ison] = '1';
      else
         Tree2PartitionDescentTree(ison, split);
   }
}

void Tree2Partition(char splits[])
{
   /* This generates branch partitions in splits[nib*(com.ns+1)].
    splits[0,...,ns-1] is the first split, splits[ns,...,2*ns-1] the second, and so on.
    For large trees, the algorithm is inefficient.
    The root node has 2 sons if the tree is rooted and >=3 sons if the tree
    is unrooted.  For unrooted tree, the mark for the first species is set to 0.
    For rooted trees, the mark for the first species can be either 0 or 1.
    */
   int unrooted = (nodes[tree.root].nson >= 3);
   int s = com.ns, lsplit = s + 1, nsplit = tree.nnode - s - 1, i, j, k;
   char* split;

   if (s <= 2) return;
   memset(splits, 0, nsplit * lsplit * sizeof(char));
   for (i = com.ns, k = 0; i < tree.nnode; i++) {
      if (i == tree.root) continue;
      split = splits + k * lsplit;
      for (j = 0; j < s; j++) split[j] = '0';
      Tree2PartitionDescentTree(i, split);
      /* If unrooted tree, set first species to 0 if tree is unrooted */
      if (unrooted && split[0] == '1')
         for (j = 0; j < s; j++) split[j] = '0' + '1' - split[j];
      k++;
   }
}

int Partition2Tree(char splits[], int lsplit, int ns, int nsplit, double label[])
{
   /* This generates the tree in nodes[] using splits or branch partitions.
      It constructs pptable[(2*s-1)*(2*s-1)] to generate the tree.
      Split i corresponds to node ns+i, while the root is ns + nsplit.
      label[] has labels for splits and become labels for nodes on the tree.
      This works even if ns=1 and ns=2.
   */
   int i, j, k, s21 = ns * 2 - 1, a, minndesc, ndesc[NS] = { 0 };  /* clade size */
   char debug = 0, * p, * pptable;

   if (nsplit > ns - 2) error2("too many splits for ns");
   if (nsplit < 0) nsplit = 0;
   if ((pptable = (char*)malloc((s21 * s21 + 1) * sizeof(char))) == NULL)
      error2("oom in Partition2Tree");
   memset(pptable, 0, s21 * s21 * sizeof(char));

   /* initialize tree */
   tree.nnode = ns + nsplit + 1;
   tree.root = ns + nsplit;
   tree.nbranch = 0;
   for (i = 0; i < tree.nnode; i++) {
      nodes[i].father = nodes[i].ibranch = -1;
      nodes[i].nson = 0;  nodes[i].label = -1;  nodes[i].branch = nodes[i].age = 0;
   }

   /* set up pptable */
   for (i = 0, p = splits, ndesc[tree.root - ns] = ns; i < nsplit; i++, p += lsplit) {
      for (j = 0; j < ns; j++) {
         if (p[j] == '1') {  /* clade (split) i includes tip j */
            pptable[j * s21 + ns + i] = 1;
            ndesc[i]++;
         }
      }
   }
   for (i = 0; i < tree.nnode - 1; i++) pptable[i * s21 + tree.root] = 1;
   for (i = 0; i < nsplit; i++) {
      for (j = 0; j < i; j++) {
         if (pptable[(ns + i) * s21 + ns + j] || pptable[(ns + j) * s21 + ns + i] || ndesc[i] == ndesc[j])
            continue;
         for (k = 0; k < ns; k++)
            if (pptable[k * s21 + ns + i] == 1 && pptable[k * s21 + ns + j] == 1) break;
         if (k < ns) {  /* i and j are ancestral to k, and are ancestral to each other. */
            if (ndesc[i] < ndesc[j])   pptable[(ns + i) * s21 + ns + j] = 1;
            else                       pptable[(ns + j) * s21 + ns + i] = 1;
         }
      }
   }
   if (debug) {
      printf("\npptable\n");
      for (i = 0; i < s21; i++) {
         for (j = 0; j < s21; j++)
            printf(" %1d", (int)pptable[i * s21 + j]);
         printf("\n");
      }
      printf("ndesc: ");
      for (i = 0; i < nsplit; i++) printf(" %2d", ndesc[i]);
      printf("\n");
   }

   /* generate tree nodes and labels.  For each nonroot node, youngest ancestor is dad. */
   for (i = 0; i < tree.nnode - 1; i++) {
      minndesc = ns + 1;  a = -1;
      for (j = ns; j < tree.nnode; j++) {
         if (pptable[i * s21 + j] == 1 && minndesc > ndesc[j - ns]) {
            minndesc = ndesc[j - ns];
            a = j;
         }
      }
      if (a < 0)
         error2("jgl");
      nodes[i].father = a;
      if (nodes[a].nson > MAXNSONS)
         error2("Partition2Tree: raise MAXNSONS?");
      nodes[a].sons[nodes[a].nson++] = i;
      if (a != tree.root && label) nodes[a].label = label[a - ns];
   }
   if (debug) {
      PrintTree(1);
      OutTreeN(F0, 1, 0);   printf("\n");
   }
   free(pptable);
   return(0);
}


int GetNSfromTreeFile(FILE* ftree, int* ns, int* ntree)
{
   /* This gets the number of species/sequence from the tree file.
   */
   char separators[] = "(,):#";
   int inname = 0, k, c;
   double y;

   *ns = *ntree = -1;
   k = fscanf(ftree, "%d%d", ns, ntree);
   if (k == 1)
   {
      *ntree = *ns; *ns = -1;
   }
   else if (k == 0) {
      *ns = 0;
      while ((c = fgetc(ftree)) != ';') {
         if (c == EOF) return(-1);
         if (strchr(separators, c)) {
            if (c == ':' || c == '#')
               fscanf(ftree, "%lf", &y);
            if (inname)
            {
               inname = 0; (*ns)++;
            }
         }
         else if (isgraph(c))
            inname = 1;
      }
      rewind(ftree);
   }
   return(0);
}

void CladeSupport(FILE* fout, char treef[], int getSnames, char maintreef[], int pick1tree)
{
   /* This reads all bootstrap or Bayesian trees from treef to identify the best trees
    (trees with the highest posterior probabilities), and to construct the majority-rule
    consensus tree.  The set of the best trees constitute the 95% or 99% credibility set
    of trees.
    A tree (ptree) is represented by its splits, ordered lexicographically, and concatenated.
    It can also read a main tree file and goes through main trees and attach support
    values on splits on each main tree.
    split1 if for one tree, and split50 is for the majority-rule consensus tree.
    */
   int i, j, k, i1, ntreeM, ntree, itreeM, sizetree, found, lline = 1024;
   int* index, * indexspace, maxnsplits, nsplits = 0, s, nsplit1=-1, nsplit50, lsplit, same;
   int maxnptree, nptree = 0, sizeptree;
   char* split1, * splits = NULL, * splitM = NULL, * split50 = NULL, * pM, * ptree, * ptrees = NULL, * ptree0, * pc;
   double* countsplits = NULL, * countptree = NULL, * Psplit50, * Psame, y, cdf;
   double Eff_MAPtree, Pjump = 0, * MAPtreeID, mMAPtree, vMAPtree, rho1MAPtree;
   char pick1treef[32] = "pick1tree.tre", line[1024];
   struct TREEN* besttree;
   int besttreeroot = -1;
   FILE* ft, * fM = NULL, * f1tree = NULL;
   int debug = 0;

   /* Count trees and splits */
   printf("\nRead tree sample, count trees & splits \n");
   ft = gfopen(treef, "r");

   if (getSnames)      /* species ordered as in first tree in file */
      GetNSfromTreeFile(ft, &com.ns, &k);
   if (com.ns < 3) error2("need >=3 species to justify this much effort.");
   s = com.ns;  lsplit = s + 1;  maxnsplits = maxnptree = s; sizeptree = (s - 2) * lsplit;
   if ((split1 = (char*)malloc(4 * (s - 2) * lsplit * sizeof(char))) == NULL)
      error2("oom splits");
   ptree = split1 + (s - 2) * lsplit;
   ptree0 = ptree + (s - 2) * lsplit;
   split50 = ptree0 + (s - 2) * lsplit;
   memset(split1, 0, 3 * (s - 2) * lsplit * sizeof(char));
   if ((Psplit50 = (double*)malloc(s * sizeof(double))) == NULL)
      error2("oom Psplit50");

   sizetree = (s * 2 - 1) * sizeof(struct TREEN);
   if ((nodes = (struct TREEN*)malloc(sizetree * 2)) == NULL) error2("oom");
   for (i = 0; i < s * 2 - 1; i++) nodes[i].annotation = NULL;
   besttree = nodes + s * 2 - 1;

   /* read trees in treef and collect splits in a sorted list */
   for (ntree = 0; ; ntree++) {
      if (nptree + s >= maxnptree) {
         maxnptree = (int)(maxnptree * (ntree < 1000 ? 2 : 1.2));
         ptrees = (char*)realloc(ptrees, maxnptree * sizeptree);
         countptree = (double*)realloc(countptree, maxnptree * sizeof(double));
         if (ptrees == NULL || countptree == NULL) error2("oom ptrees || countptree");
         memset(ptrees + nptree * sizeptree, 0, (maxnptree - nptree) * sizeptree);
         memset(countptree + nptree, 0, (maxnptree - nptree) * sizeof(double));
      }
      if (nsplits + s >= maxnsplits) {
         maxnsplits *= 2;
         splits = (char*)realloc(splits, maxnsplits * lsplit * sizeof(char));
         countsplits = (double*)realloc(countsplits, maxnsplits * sizeof(double));
         if (splits == NULL || countsplits == NULL) error2("oom splits realloc");
      }

      /* if (getSnames), copy species/sequence names from first tree in file.  */
      if (ReadTreeN(ft, &i, (getSnames && ntree == 0), 1)) break;
      if (debug || (ntree + 1) % 5000 == 0) {
         printf("\rtree %5d  ", ntree + 1);
         if (s < 15) OutTreeN(F0, 1, 0);
      }
      Tree2Partition(split1);
      nsplit1 = tree.nnode - s - 1;

      /* Process the tree */
      qsort(split1, nsplit1, lsplit, (int(*)(const void*, const void*))strcmp);
      if (debug) {
         for (i = 0; i < nsplit1; i++) printf(" %s", split1 + i * lsplit);
         printf("\n");
      }

      /* ptree is the same as split1, but with trailing 0s removed. */
      for (i = 0, pc = ptree; i < nsplit1; i++)
         for (j = 0; j < s; j++) *pc++ = split1[i * lsplit + j];
      j = binarysearch(ptree, ptrees, nptree, sizeptree, (int(*)(const void*, const void*))strcmp, &found);
      if (found)
         countptree[j]++;
      else {
         if (j < nptree) {
            memmove(ptrees + (j + 1) * sizeptree, ptrees + j * sizeptree, (nptree - j) * sizeptree);
            memmove(countptree + j + 1, countptree + j, (nptree - j) * sizeof(double));
         }
         memmove(ptrees + j * sizeptree, ptree, sizeptree);
         nptree++;
         countptree[j] = 1;
      }

      if (ntree > 0 && strcmp(ptree0, ptree))
         Pjump++;
      memmove(ptree0, ptree, (s - 2) * lsplit * sizeof(char));

      /* Process the splits in the tree */
      for (i = 0; i < nsplit1; i++) {  /* going through splits in current tree */
         j = binarysearch(split1 + i * lsplit, splits, nsplits, lsplit, (int(*)(const void*, const void*))strcmp, &found);
         if (found)  /* found */
            countsplits[j]++;
         else {
            if (j < nsplits) {  /* check the size of the moved block here */
               memmove(splits + (j + 1) * lsplit, splits + j * lsplit, (nsplits - j) * lsplit * sizeof(char));
               memmove(countsplits + (j + 1), countsplits + j, (nsplits - j) * sizeof(double));
            }
            memcpy(splits + j * lsplit, split1 + i * lsplit, lsplit * sizeof(char));
            nsplits++;
            countsplits[j] = 1;
         }
      }
      if (debug) {
         printf("%4d splits: ", nsplits);
         for (k = 0; k < nsplits; k++) printf(" %s (%.0f)", splits + k * lsplit, countsplits[k]);
         printf("\n");
      }
   }
   printf("\n%6d trees read, %d distinct trees.\n", ntree, nptree);

   k = max2(nsplits, nptree);
   if ((index = (int*)malloc(k * 2 * sizeof(int))) == NULL) error2("oom index");
   indexspace = index + k;

   printf("\nSpecies in order:\n");
   for (i = 0; i < s; i++) printf("%2d. %s\n", i + 1, com.spname[i]);
   printf("\n(A) Best trees in the sample (%d distinct trees in all)\n", nptree);
   fprintf(fout, "\n\nSpecies in order:\n");
   for (i = 0; i < s; i++) fprintf(fout, "%2d. %s\n", i + 1, com.spname[i]);
   fprintf(fout, "\n(A) Best trees in the sample (%d distinct trees in all)\n", nptree);

   indexing(countptree, nptree, index, 1, indexspace);

   for (k = 0, cdf = 0; k < nptree; k++) {
      j = index[k];  y = countptree[j];
      for (i = 0, pc = split1; i < nsplit1; i++, *pc++ = '\0')
         for (i1 = 0; i1 < s; i1++)
            *pc++ = ptrees[j * sizeptree + i * s + i1];
      Partition2Tree(split1, lsplit, s, nsplit1, NULL);
      printf(" %6.0f %8.6f %8.6f  ", y, y / ntree, (cdf += y / ntree));
      OutTreeN(F0, 1, 0);
      /* for(i=0; i<nsplit1; i++) printf("%s", split1+i*lsplit);  */
      printf("\n");

      fprintf(fout, " %6.0f %8.6f %8.6f  ", y, y / ntree, cdf);
      OutTreeN(fout, 1, 0);
      /* for(i=0; i<nsplit1; i++) fprintf(fout, "%s", split1+i*lsplit); */
      fprintf(fout, "\n");
      if (k == 0) {
         memmove(besttree, nodes, sizetree);
         besttreeroot = tree.root;
      }
      if (k > 9 && cdf > 0.99 && y / ntree < 0.001) break;
   }

   printf("\n(B) Best splits in the sample of trees (%d splits in all)\n", nsplits);
   indexing(countsplits, nsplits, index, 1, indexspace);
   for (k = 0; k < nsplits; k++) {
      j = index[k];  y = countsplits[j];
      printf(" %6.0f %9.6f  %s\n", y, y / ntree, splits + j * lsplit);
      if (k > 9 && y / ntree < 0.001) break;
   }
   fprintf(fout, "\n(B) Best splits in the sample of trees (%d splits in all)\n", nsplits);
   for (k = 0; k < nsplits; k++) {
      j = index[k];  y = countsplits[j];
      fprintf(fout, " %6.0f %9.6f  %s\n", y, y / ntree, splits + j * lsplit);
      if (k > 9 && y / ntree < 0.001) break;
   }

   /* Majority-rule consensus tree */
   for (k = 0, nsplit50 = 0; k < nsplits; k++)
      if (countsplits[k] / ntree >= 0.5) nsplit50++;
   for (k = 0, nsplit50 = 0; k < nsplits; k++) {
      if (countsplits[k] / ntree > 0.5) {
         memmove(split50 + nsplit50 * lsplit, splits + k * lsplit, lsplit);
         Psplit50[nsplit50++] = countsplits[k] / ntree;
      }
   }
   Partition2Tree(split50, lsplit, s, nsplit50, Psplit50);
   printf("\n(C) Majority-rule consensus tree\n");
   OutTreeN(F0, 1, PrLabel);    printf("\n");
   fprintf(fout, "\n(C) Majority-rule consensus tree\n");
   OutTreeN(fout, 1, PrLabel);        fprintf(fout, "\n");

   /* Probabilities of trees in the main tree file */
   if (maintreef) fM = fopen(maintreef, "r");
   if (fM == NULL) {
      fM = fopen("besttree", "w+");
      fprintf(fM, "%6d  %6d\n\n", com.ns, 1);
      memmove(nodes, besttree, sizetree);  tree.root = besttreeroot;
      OutTreeN(fM, 1, 0);  fprintf(fM, "\n");
      rewind(fM);
   }

   fscanf(fM, "%d%d", &i, &ntreeM);
   if (i != s || ntreeM < 1) error2("<ns> <ntree> on the first line in main tree.");

   splitM = (char*)malloc(ntreeM * (s - 2) * lsplit * sizeof(char));
   Psame = (double*)malloc(ntreeM * sizeof(double));
   if (splitM == NULL || Psame == NULL) error2("oom splitM");
   zero(Psame, ntreeM);
   if (pick1tree >= 1 && pick1tree <= ntreeM && (f1tree = (FILE*)fopen(pick1treef, "w")) == NULL)
      error2("oom");
   for (itreeM = 0, pM = splitM; itreeM < ntreeM; itreeM++, pM += (s - 2) * lsplit) {
      if (ReadTreeN(fM, &i, 0, 1)) break;
      if (tree.nnode<s * 2 - 2 || tree.nnode>s * 2 - 1) error2("main trees have to be binary");
      Tree2Partition(pM);
      qsort(pM, tree.nnode - s - 1, lsplit, (int(*)(const void*, const void*))strcmp);
      if (debug) {
         printf("\nmain tree %2d: ", itreeM + 1);
         OutTreeN(F0, 1, 0);
         for (i = 0; i < tree.nnode - s - 1; i++) printf(" %s", pM + i * lsplit);
      }
   }
   if ((MAPtreeID = (double*)malloc(ntree * sizeof(double))) == NULL)
      error2("oom MAPtreeID");
   memset(MAPtreeID, 0, ntree * sizeof(double));

   /* read the tree sample again */
   rewind(ft);
   for (ntree = 0; ; ntree++) {
      if (ReadTreeN(ft, &i, 0, 0)) break;
      fgets(line, lline, ft);
      Tree2Partition(split1);
      for (itreeM = 0, pM = splitM; itreeM < ntreeM; itreeM++, pM += (s - 2) * lsplit) {
         for (i = 0, same = 1; i < tree.nnode - s - 1; i++) {
            if (bsearch(split1 + i * lsplit, pM, tree.nnode - s - 1, lsplit, (int(*)(const void*, const void*))strcmp) == NULL) {
               same = 0; break;
            }
         }
         if (same) {
            Psame[itreeM] ++;
            if (itreeM == 0) MAPtreeID[ntree] = 1;
            if (pick1tree - 1 == itreeM) {
               OutTreeN(f1tree, 1, 1);   fprintf(f1tree, "%s", line);
            }
         }
      }
   }

   printf("\n(D) Best tree (or trees from the maintree file) with support values\n");
   fprintf(fout, "\n(D) Best tree (or trees from the maintree file) with support values\n");
   /* read the main trees another round just for printing. */
   rewind(fM);
   fscanf(fM, "%d%d", &s, &ntreeM);
   for (itreeM = 0; itreeM < ntreeM; itreeM++) {
      if (ReadTreeN(fM, &i, 0, 1)) break;
      Tree2Partition(split1);
      for (i = s, k = 0; i < tree.nnode; i++) {
         if (i == tree.root) continue;
         j = binarysearch(split1 + k * lsplit, splits, nsplits, lsplit, (int(*)(const void*, const void*))strcmp, &found);
         if (found) nodes[i].label = countsplits[j] / ntree;
         k++;
      }
      OutTreeN(F0, 1, PrLabel);
      printf("  [P = %8.6f]\n", Psame[itreeM] / ntree);
      OutTreeN(fout, 1, PrLabel);
      fprintf(fout, "  [P = %8.6f]\n", Psame[itreeM] / ntree);
   }
   if (pick1tree) printf("\ntree #%d collected into %s\n", pick1tree, pick1treef);

#if(1)
   Eff_MAPtree = Eff_IntegratedCorrelationTime(MAPtreeID, ntree, &mMAPtree, &vMAPtree, &rho1MAPtree);
#else
   Eff_MAPtree = Eff_IntegratedCorrelationTime2(MAPtreeID, ntree, (int)sqrt(ntree) / 2, &mMAPtree, &vMAPtree);
#endif

   if (fM) {
      free(splitM);  free(Psame);
      fclose(fM);    if (f1tree) fclose(f1tree);
   }

   printf("\nPjump =       %12.9f\nEff_MAPtree = %12.9f\n", Pjump / (ntree - 1), Eff_MAPtree);

   free(split1);  free(splits);  free(countsplits);  free(Psplit50);
   free(ptrees);  free(countptree);  free(index);
   free(nodes);   free(MAPtreeID);
   fclose(ft);
   exit(0);
}

#endif


int NSameBranch(char partition1[], char partition2[], int nib1, int nib2, int IBsame[])
{
   /* counts the number of correct (identical) bipartitions.
    nib1 and nib2 are the numbers of interior branches in the two trees
    correctIB[0,...,(correctbranch-1)] lists the correct interior branches,
    that is, interior branches in tree 1 that is also in tree 2.
    IBsame[i]=1 if interior branch i is correct.
    */
   int i, j, nsamebranch = 0;

#if(1)
   for (i = 0, nsamebranch = 0; i < nib1; i++)
      for (j = 0, IBsame[i] = 0; j < nib2; j++) {
         if (strcmp(partition1 + i * (com.ns + 1), partition2 + j * (com.ns + 1)) == 0) {
            nsamebranch++;  IBsame[i] = 1;  break;
         }
      }
#else
   int k = 0;
   for (i = 0; i < nib1; i++)
      for (j = 0, IBsame[i] = 0; j < nib2; j++) {
         for (k = 0; k < com.ns; k++)
            if (partition1[i * (com.ns + 1) + k] != partition2[j * (com.ns + 1) + k]) break;
         if (k == com.ns) {
            nsamebranch++;  IBsame[i] = 1;  break;
         }
      }
#endif
   return (nsamebranch);
}


int AddSpecies(int is, int ib)
{
   /* Add species (is) to tree at branch ib.  The tree currently has
    is+1-1 species.  Interior node numbers are increased by 2 to make
    room for the new nodes.
    if(com.clock && ib==tree.nbranch), the new species is added as an
    outgroup to the rooted tree.
    */
   int i, j, it;

   if (ib > tree.nbranch + 1 || (ib == tree.nbranch && !com.clock)) return(-1);

   if (ib == tree.nbranch && com.clock) {
      for (i = 0; i < tree.nbranch; i++) for (j = 0; j < 2; j++)
         if (tree.branches[i][j] >= is) tree.branches[i][j] += 2;
      it = tree.root;  if (tree.root >= is) it += 2;
      for (i = 0; i < 2; i++) tree.branches[tree.nbranch + i][0] = tree.root = is + 1;
      tree.branches[tree.nbranch++][1] = it;
      tree.branches[tree.nbranch++][1] = is;
   }
   else {
      for (i = 0; i < tree.nbranch; i++) for (j = 0; j < 2; j++)
         if (tree.branches[i][j] >= is) tree.branches[i][j] += 2;
      it = tree.branches[ib][1];
      tree.branches[ib][1] = is + 1;
      tree.branches[tree.nbranch][0] = is + 1;
      tree.branches[tree.nbranch++][1] = it;
      tree.branches[tree.nbranch][0] = is + 1;
      tree.branches[tree.nbranch++][1] = is;
      if (tree.root >= is) tree.root += 2;
   }
   BranchToNode();
   return (0);
}


#ifdef TREESEARCH

static struct TREE
{
   struct TREEB tree; struct TREEN nodes[2 * NS - 1]; double x[NP];
}
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
   int step = 0, ntree = 0, nmove = 0, improve = 0, ineighb, i;
   int sizetree = (2 * com.ns - 1) * sizeof(struct TREEN);
   double* x = treestar.x;
   FILE* ftree;

   if (com.clock) error2("\n\aerr: pertubation does not work with a clock yet.\n");
   if (initialMP && !com.cleandata)
      error2("\ncannot get initial parsimony tree for gapped data yet.");

   fprintf(fout, "\n\nHeuristic tree search by NNI perturbation\n");
   if (initialMP) {
      if (noisy) printf("\nInitial tree from stepwise addition with MP:\n");
      fprintf(fout, "\nInitial tree from stepwise addition with MP:\n");
      StepwiseAdditionMP(space);
   }
   else {
      if (noisy) printf("\nInitial tree read from file %s:\n", com.treef);
      fprintf(fout, "\nInitial tree read from file.\n");
      if ((ftree = fopen(com.treef, "r")) == NULL) error2("treefile not exist?");
      fscanf(ftree, "%d%d", &i, &ntree);
      if (i != com.ns) error2("ns in the tree file");
      if (ReadTreeN(ftree, &i, 0, 1)) error2("err tree..");
      fclose(ftree);
   }
   if (noisy) {
      printf("\n");   OutTreeN(F0, 0, 0);  printf("\n");
   }
   tree.lnL = TreeScore(x, space);
   if (noisy) { OutTreeN(F0, 0, 1);  printf("\n lnL = %.4f\n", -tree.lnL); }
   OutTreeN(fout, 1, 1);  fprintf(fout, "\n lnL = %.4f\n", -tree.lnL);
   if (com.np > com.ntime) {
      fprintf(fout, "\tparameters:");
      for (i = com.ntime; i < com.np; i++) fprintf(fout, "%9.5f", x[i]);
      fprintf(fout, "\n");
   }
   fflush(fout);
   treebest.tree = tree;  memcpy(treebest.nodes, nodes, sizetree);

   for (step = 0; ; step++) {
      for (ineighb = 0, improve = 0; ineighb < (tree.nbranch - com.ns) * 2; ineighb++) {
         tree = treebest.tree; memcpy(nodes, treebest.nodes, sizetree);
         NeighborNNI(ineighb);
         if (noisy) {
            printf("\nTrying tree # %d (%d move[s]) \n", ++ntree, nmove);
            OutTreeN(F0, 0, 0);  printf("\n");
         }
         tree.lnL = TreeScore(x, space);
         if (noisy) { OutTreeN(F0, 1, 1); printf("\n lnL = %.4f\n", -tree.lnL); }
         if (noisy && com.np > com.ntime) {
            printf("\tparameters:");
            for (i = com.ntime; i < com.np; i++) printf("%9.5f", x[i]);
            printf("\n");
         }
         if (tree.lnL <= treebest.tree.lnL) {
            treebest.tree = tree;  memcpy(treebest.nodes, nodes, sizetree);
            improve = 1; nmove++;
            if (noisy) printf(" moving to this tree\n");
            if (fout) {
               fprintf(fout, "\nA better tree:\n");
               OutTreeN(fout, 0, 0); fprintf(fout, "\n");
               OutTreeN(fout, 1, 1); fprintf(fout, "\n");
               fprintf(fout, "\nlnL = %.4f\n", tree.lnL);
               if (com.np > com.ntime) {
                  fprintf(fout, "\tparameters:");
                  for (i = com.ntime; i < com.np; i++) fprintf(fout, "%9.5f", x[i]);
                  fprintf(fout, "\n");
               }
               fflush(fout);
            }
         }
      }
      if (!improve) break;
   }
   tree = treebest.tree;  memcpy(nodes, treebest.nodes, sizetree);
   if (noisy) {
      printf("\n\nBest tree found:\n");
      OutTreeN(F0, 0, 0);  printf("\n");
      OutTreeN(F0, 1, 1);  printf("\n");
      printf("\nlnL = %.4f\n", tree.lnL);
   }
   if (fout) {
      fprintf(fout, "\n\nBest tree found:\n");
      OutTreeN(fout, 0, 0);  fprintf(fout, "\n");
      OutTreeN(fout, 1, 1);  fprintf(fout, "\n");
      fprintf(fout, "\nlnL = %.4f\n", tree.lnL);
   }
   return (0);
}


static int* _U0, * _step0, _mnnode;
/* up pass characters and changes for the star tree: each of size npatt*nnode*/

int StepwiseAdditionMP(double space[])
{
   /* tree search by species addition. */
   int  ns0 = com.ns, is, i, j, h, tiestep = 0, tie, bestbranch = 0;
   int sizetree = (2 * com.ns - 1) * sizeof(struct TREEN);
   double bestscore = 0, score;

   _mnnode = com.ns * 2 - 1;
   _U0 = (int*)malloc(com.npatt * _mnnode * sizeof(int));
   _step0 = (int*)malloc(com.npatt * _mnnode * sizeof(int));
   if (noisy > 2)
      printf("\n%9zd bytes for MP (U0 & N0)\n", 2 * com.npatt * _mnnode * sizeof(int));
   if (_U0 == NULL || _step0 == NULL) error2("oom U0&step0");

   tree.nbranch = tree.root = com.ns = 3;
   for (i = 0; i < tree.nbranch; i++) { tree.branches[i][0] = com.ns; tree.branches[i][1] = i; }
   BranchToNode();
   for (h = 0; h < com.npatt; h++)
      for (i = 0; i < com.ns; i++) {
         _U0[h * _mnnode + i] = 1 << (com.z[i][h] - 1); _step0[h * _mnnode + i] = 0;
      }
   for (is = com.ns, tie = 0; is < ns0; is++) {
      treestar.tree = tree;  memcpy(treestar.nodes, nodes, sizetree);

      for (j = 0; j < treestar.tree.nbranch; j++, com.ns--) {
         tree = treestar.tree;  memcpy(nodes, treestar.nodes, sizetree);
         com.ns++;
         AddSpecies(is, j);
         score = MPScoreStepwiseAddition(is, space, 0);
         /*
          OutTreeN(F0, 0, 0);
          printf(" Add sp %d (ns=%d) at branch %d, score %.0f\n", is+1,com.ns,j+1,score);
          */
         if (j && score == bestscore) tiestep = 1;
         if (j == 0 || score < bestscore || (score == bestscore && rndu() < .1)) {
            tiestep = 0;
            bestscore = score; bestbranch = j;
         }
      }
      tie += tiestep;
      tree = treestar.tree;  memcpy(nodes, treestar.nodes, sizetree);
      com.ns = is + 1;
      AddSpecies(is, bestbranch);
      score = MPScoreStepwiseAddition(is, space, 1);

      if (noisy)
      {
         printf("\r  Added %d [%5.0f steps]", is + 1, -bestscore); fflush(F0);
      }
   }
   if (noisy > 2) printf("  %d stages with ties, ", tie);
   tree.lnL = bestscore;
   free(_U0); free(_step0);
   return (0);
}

double MPScoreStepwiseAddition(int is, double space[], int save)
{
   /* this changes only the part of the tree affected by the newly added
    species is.
    save=1 for the best tree, so that _U0 & _step0 are updated
    */
   int* U, * N, U3[3], h, ist, i, father, son2, * pU0 = _U0, * pN0 = _step0;
   double score;

   U = (int*)space;  N = U + _mnnode;
   for (h = 0, score = 0; h < com.npatt; h++, pU0 += _mnnode, pN0 += _mnnode) {
      for (i = 0; i < tree.nnode; i++)
      {
         U[i] = pU0[i - 2 * (i >= is)]; N[i] = pN0[i - 2 * (i >= is)];
      }
      U[is] = 1 << (com.z[is][h] - 1);  N[is] = 0;
      for (ist = is; (father = nodes[ist].father) != tree.root; ist = father) {
         if ((son2 = nodes[father].sons[0]) == ist)  son2 = nodes[father].sons[1];
         N[father] = N[ist] + N[son2];
         if ((U[father] = U[ist] & U[son2]) == 0)
         {
            U[father] = U[ist] | U[son2];  N[father]++;
         }
      }
      for (i = 0; i < 3; i++) U3[i] = U[nodes[tree.root].sons[i]];
      N[tree.root] = 2;
      if (U3[0] & U3[1] & U3[2]) N[tree.root] = 0;
      else if (U3[0] & U3[1] || U3[1] & U3[2] || U3[0] & U3[2]) N[tree.root] = 1;
      for (i = 0; i < 3; i++) N[tree.root] += N[nodes[tree.root].sons[i]];

      if (save) {
         memcpy(pU0, U, tree.nnode * sizeof(int));
         memcpy(pN0, N, tree.nnode * sizeof(int));
      }
      score += N[tree.root] * com.fpatt[h];
   }
   return (score);
}


double TreeScore(double x[], double space[])
{
   static int fromfile = 0;
   int i;
   double xb[NP][2], e = 1e-9, lnL = 0;

   if (com.clock == 2) error2("local clock in TreeScore");
   com.ntime = com.clock ? tree.nnode - com.ns : tree.nbranch;

   GetInitials(x, &i);  /* this shoulbe be improved??? */
   if (i) fromfile = 1;
   PointconPnodes();

   if (com.method == 0 || !fromfile) SetxBound(com.np, xb);

   if (fromfile) {
      lnL = com.plfun(x, com.np);
      com.np = com.ntime;
   }
   NFunCall = 0;
   if (com.method == 0 || com.ntime == 0)
      ming2(NULL, &lnL, com.plfun, NULL, x, xb, space, e, com.np);
   else
      minB(NULL, &lnL, x, xb, e, space);

   return(lnL);
}


int StepwiseAddition(FILE* fout, double space[])
{
   /* heuristic tree search by species addition.  Species are added in the order
    of occurrence in the data.
    Try to get good initial values.
    */
   char* z0[NS], * spname0[NS];
   int ns0 = com.ns, is, i, j, randadd = 0, order[NS];
   int sizetree = (2 * com.ns - 1) * sizeof(struct TREEN);
   double bestscore = 0, score, * x = treestar.x;

   if (com.ns > 50) printf("if this crashes, increase com.sspace?");

   if (com.ns < 3) error2("2 sequences, no need for tree search");
   if (noisy) printf("\n\nHeuristic tree search by stepwise addition\n");
   if (fout) fprintf(fout, "\n\nHeuristic tree search by stepwise addition\n");
   for (i = 0; i < ns0; i++) { z0[i] = com.z[i]; spname0[i] = com.spname[i]; }
   tree.nbranch = tree.root = com.ns = (com.clock ? 2 : 3);

   for (i = 0; i < ns0; i++) order[i] = i;
   if (randadd) {
      for (i = 0; i < ns0; i++) {
         j = (int)(ns0 * rndu()); is = order[i]; order[i] = order[j]; order[j] = is;
      }
      if (noisy) for (i = 0; i < ns0; i++) printf(" %d", order[i] + 1);
      if (fout) {
         fputs("\nOrder of species addition:\n", fout);
         for (i = 0; i < ns0; i++) fprintf(fout, "%3d  %-s\n", order[i] + 1, com.spname[order[i]]);
      }
      for (i = 0; i < ns0; i++) {
         com.z[i] = z0[order[i]];
         com.spname[i] = spname0[order[i]];
      }
   }

   for (i = 0; i < tree.nbranch; i++) {
      tree.branches[i][0] = com.ns; tree.branches[i][1] = i;
   }
   BranchToNode();
   for (is = com.ns; is < ns0; is++) {                  /* add the is_th species */
      treestar.tree = tree;  memcpy(treestar.nodes, nodes, sizetree);

      for (j = 0; j < treestar.tree.nbranch + (com.clock > 0); j++, com.ns--) {
         tree = treestar.tree;  memcpy(nodes, treestar.nodes, sizetree);
         com.ns++;
         AddSpecies(is, j);
         score = TreeScore(x, space);
         if (noisy > 1)
         {
            printf("\n "); OutTreeN(F0, 0, 0); printf("%12.3f", -score);
         }

         if (j == 0 || score < bestscore || (score == bestscore && rndu() < .2)) {
            treebest.tree = tree;  memcpy(treebest.nodes, nodes, sizetree);
            xtoy(x, treebest.x, com.np);
            bestscore = score;
         }
      }
      tree = treebest.tree;  memcpy(nodes, treebest.nodes, sizetree);
      xtoy(treebest.x, x, com.np);
      com.ns = is + 1;

      if (noisy) {
         printf("\n\nAdded sp. %d, %s [%.3f]\n", is + 1, com.spname[is], -bestscore);
         OutTreeN(F0, 0, 0);  printf("\n");
         OutTreeN(F0, 1, 0);  printf("\n");
         if (com.np > com.ntime) {
            printf("\tparameters:");
            for (i = com.ntime; i < com.np; i++) printf("%9.5f", x[i]);
            printf("\n");
         }
      }
      if (fout) {
         fprintf(fout, "\n\nAdded sp. %d, %s [%.3f]\n",
            is + 1, com.spname[is], -bestscore);
         OutTreeN(fout, 0, 0); fprintf(fout, "\n");
         OutTreeN(fout, 1, 1); fprintf(fout, "\n");
         if (com.np > com.ntime) {
            fprintf(fout, "\tparameters:");
            for (i = com.ntime; i < com.np; i++) fprintf(fout, "%9.5f", x[i]);
            fprintf(fout, "\n");
         }
         fflush(fout);
      }
   }
   tree.lnL = bestscore;

   return (0);
}


int DecompTree(int inode, int ison1, int ison2);
#define hdID(i,j) (max2(i,j)*(max2(i,j)-1)/2+min2(i,j))

int StarDecomposition(FILE* fout, double space[])
{
   /* automatic tree search by star decomposition, nhomo<=1
    returns (0,1,2,3) for the 4s problem.
    */
   int status = 0, stage = 0, i, j, itree, ntree = 0, ntreet, best = 0, improve = 1, collaps = 0;
   int inode, nson = 0, ison1, ison2;
   int sizetree = (2 * com.ns - 1) * sizeof(struct TREEN);
   double x[NP];
   FILE* ftree, * fsum = frst;

   if (com.runmode == 1) {   /* read the star-like tree from tree file */
      if ((ftree = fopen(com.treef, "r")) == NULL)
         error2("no treefile");
      fscanf(ftree, "%d%d", &i, &ntree);
      if (ReadTreeN(ftree, &i, 0, 1)) error2("err tree file");
      fclose(ftree);
   }
   else {                  /* construct the star tree of ns species */
      tree.nnode = (tree.nbranch = tree.root = com.ns) + 1;
      for (i = 0; i < tree.nbranch; i++)
      {
         tree.branches[i][0] = com.ns; tree.branches[i][1] = i;
      }
      com.ntime = com.clock ? 1 : tree.nbranch;
      BranchToNode();
   }
   if (noisy) { printf("\n\nstage 0: ");       OutTreeN(F0, 0, 0); }
   if (fsum) { fprintf(fsum, "\n\nstage 0: ");  OutTreeN(fsum, 0, 0); }
   if (fout) { fprintf(fout, "\n\nstage 0: ");  OutTreeN(fout, 0, 0); }

   tree.lnL = TreeScore(x, space);

   if (noisy)  printf("\nlnL:%14.6f%6d", -tree.lnL, NFunCall);
   if (fsum) fprintf(fsum, "\nlnL:%14.6f%6d", -tree.lnL, NFunCall);
   if (fout) {
      fprintf(fout, "\nlnL(ntime:%3d  np:%3d):%14.6f\n",
         com.ntime, com.np, -tree.lnL);
      OutTreeB(fout);  fprintf(fout, "\n");
      for (i = 0; i < com.np; i++) fprintf(fout, "%9.5f", x[i]);
      fprintf(fout, "\n");
   }
   treebest.tree = tree;  memcpy(treebest.nodes, nodes, sizetree);
   for (i = 0; i < com.np; i++) treebest.x[i] = x[i];
   for (ntree = 0, stage = 1; ; stage++) {
      for (inode = treebest.tree.nnode - 1; inode >= 0; inode--) {
         nson = treebest.nodes[inode].nson;
         if (nson > 3) break;
         if (com.clock) { if (nson > 2) break; }
         else if (nson > 2 + (inode == treebest.tree.root)) break;
      }
      if (inode == -1 || /*stage>com.ns-3+com.clock ||*/ !improve) { /* end */
         tree = treebest.tree;  memcpy(nodes, treebest.nodes, sizetree);

         if (noisy) {
            printf("\n\nbest tree: ");  OutTreeN(F0, 0, 0);
            printf("   lnL:%14.6f\n", -tree.lnL);
         }
         if (fsum) {
            fprintf(fsum, "\n\nbest tree: ");  OutTreeN(fsum, 0, 0);
            fprintf(fsum, "   lnL:%14.6f\n", -tree.lnL);
         }
         if (fout) {
            fprintf(fout, "\n\nbest tree: ");  OutTreeN(fout, 0, 0);
            fprintf(fout, "   lnL:%14.6f\n", -tree.lnL);
            OutTreeN(fout, 1, 1);  fprintf(fout, "\n");
         }
         break;
      }
      treestar = treebest;  memcpy(nodes, treestar.nodes, sizetree);

      if (collaps && stage) {
         printf("\ncollapsing nodes\n");
         OutTreeN(F0, 1, 1); printf("\n");

         tree = treestar.tree;  memcpy(nodes, treestar.nodes, sizetree);
         for (i = com.ns, j = 0; i < tree.nnode; i++)
            if (i != tree.root && nodes[i].branch < 1e-7)
            {
               CollapsNode(i, treestar.x);  j++;
            }
         treestar.tree = tree;  memcpy(treestar.nodes, nodes, sizetree);

         if (j) {
            fprintf(fout, "\n%d node(s) collapsed\n", j);
            OutTreeN(fout, 1, 1);  fprintf(fout, "\n");
         }
         if (noisy) {
            printf("\n%d node(s) collapsed\n", j);
            OutTreeN(F0, 1, 1);  printf("\n");
            /*            if (j) getchar (); */
         }
      }

      ntreet = nson * (nson - 1) / 2;
      if (!com.clock && inode == treestar.tree.root && nson == 4)  ntreet = 3;
      com.ntime++;  com.np++;

      if (noisy) {
         printf("\n\nstage %d:%6d trees, ntime:%3d  np:%3d\nstar tree: ",
            stage, ntreet, com.ntime, com.np);
         OutTreeN(F0, 0, 0);
         printf("  lnL:%10.3f\n", -treestar.tree.lnL);
      }
      if (fsum) {
         fprintf(fsum, "\n\nstage %d:%6d trees, ntime:%3d  np:%3d\nstar tree: ",
            stage, ntreet, com.ntime, com.np);
         OutTreeN(fsum, 0, 0);
         fprintf(fsum, "  lnL:%10.6f\n", -treestar.tree.lnL);
      }
      if (fout) {
         fprintf(fout, "\n\nstage %d:%6d trees\nstar tree: ", stage, ntreet);
         OutTreeN(fout, 0, 0);
         fprintf(fout, " lnL:%14.6f\n", -treestar.tree.lnL);
         OutTreeN(fout, 1, 1);  fprintf(fout, "\n");
      }

      for (ison1 = 0, itree = improve = 0; ison1 < nson; ison1++) {
         for (ison2 = ison1 + 1; ison2 < nson && itree < ntreet; ison2++, itree++, ntree++) {
            DecompTree(inode, ison1, ison2);

            for (i = com.np - 1; i > 0; i--)  x[i] = treestar.x[i - 1];
            if (!com.clock)
               for (i = 0; i < tree.nbranch; i++)
                  x[i] = max2(nodes[tree.branches[i][1]].branch * 0.99, 0.0001);
            else
               for (i = 1, x[0] = max2(x[0], .01); i < com.ntime; i++)  x[i] = .5;

            if (noisy) {
               printf("\nS=%d:%3d/%d  T=%4d  ", stage, itree + 1, ntreet, ntree + 1);
               OutTreeN(F0, 0, 0);
            }
            if (fsum) {
               fprintf(fsum, "\nS=%d:%3d/%d  T=%4d  ", stage, itree + 1, ntreet, ntree + 1);
               OutTreeN(fsum, 0, 0);
            }
            if (fout) {
               fprintf(fout, "\nS=%d:%4d/%4d  T=%4d ", stage, itree + 1, ntreet, ntree + 1);
               OutTreeN(fout, 0, 0);
            }
            tree.lnL = TreeScore(x, space);

            if (tree.lnL < treebest.tree.lnL) {
               treebest.tree = tree;  memcpy(treebest.nodes, nodes, sizetree);
               for (i = 0; i < com.np; i++) treebest.x[i] = x[i];
               best = itree + 1;   improve = 1;
            }
            if (noisy)
               printf("%6d%2c %+8.6f", NFunCall, (status ? '?' : 'X'), treestar.tree.lnL - tree.lnL);
            if (fsum) {
               fprintf(fsum, "%6d%2c", NFunCall, (status ? '?' : 'X'));
               for (i = com.ntime; i < com.np; i++)  fprintf(fsum, "%7.3f", x[i]);
               fprintf(fsum, " %+8.6f", treestar.tree.lnL - tree.lnL);
               fflush(fsum);
            }
            if (fout) {
               fprintf(fout, "\nlnL(ntime:%3d  np:%3d):%14.6f\n",
                  com.ntime, com.np, -tree.lnL);
               OutTreeB(fout);   fprintf(fout, "\n");
               for (i = 0; i < com.np; i++) fprintf(fout, "%9.5f", x[i]);
               fprintf(fout, "\n");
               fflush(fout);
            }
         }  /* for (itree) */
      }
   }    /* for (stage) */

   if (com.ns <= 4 && !improve && best) error2("strange");

   if (com.ns <= 4) return (best);
   else return (0);
}

int DecompTree(int inode, int ison1, int ison2)
{
   /* decompose treestar at NODE inode into tree and nodes[]
    */
   int i, son1, son2;
   int sizetree = (2 * com.ns - 1) * sizeof(struct TREEN);
   double bt, fmid = 0.001, fclock = 0.0001;

   tree = treestar.tree;  memcpy(nodes, treestar.nodes, sizetree);
   for (i = 0, bt = 0; i < tree.nnode; i++)
      if (i != tree.root) bt += nodes[i].branch / tree.nbranch;

   nodes[tree.nnode].nson = 2;
   nodes[tree.nnode].sons[0] = son1 = nodes[inode].sons[ison1];
   nodes[tree.nnode].sons[1] = son2 = nodes[inode].sons[ison2];
   nodes[tree.nnode].father = inode;
   nodes[son1].father = nodes[son2].father = tree.nnode;

   nodes[inode].sons[ison1] = tree.nnode;
   for (i = ison2; i < nodes[inode].nson; i++)
      nodes[inode].sons[i] = nodes[inode].sons[i + 1];
   nodes[inode].nson--;

   tree.nnode++;
   NodeToBranch();
   if (!com.clock)
      nodes[tree.nnode - 1].branch = bt * fmid;
   else
      nodes[tree.nnode - 1].age = nodes[inode].age * (1 - fclock);

   return(0);
}


#ifdef REALSEQUENCE


int MultipleGenes(FILE* fout, FILE* ftree, FILE* fpair[], double space[])
{
   /* This does the separate analysis of multiple-gene data.
    Note that com.pose[] is not correct and so RateAncestor = 0 should be set
    in baseml and codeml.
    */
   int ig = 0, j, ngene0, npatt0, lgene0[NGENE] = { 0 }, posG0[NGENE + 1] = {0};
   int nb = ((com.seqtype == 1 && !com.cleandata) ? 3 : 1);

   if (com.ndata > 1) error2("multiple data sets & multiple genes?");

   ngene0 = com.ngene;  npatt0 = com.npatt;
   for (ig = 0; ig < ngene0; ig++)   lgene0[ig] = com.lgene[ig];
   for (ig = 0; ig < ngene0 + 1; ig++) posG0[ig] = com.posG[ig];

   ig = 0;
   /*
    printf("\nStart from gene (1-%d)? ", com.ngene);
    scanf("%d", &ig);
    ig--;
    */

   for (; ig < ngene0; ig++) {

      com.ngene = 1;
      com.ls = com.lgene[0] = ig == 0 ? lgene0[0] : lgene0[ig] - lgene0[ig - 1];
      com.npatt = ig == ngene0 - 1 ? npatt0 - posG0[ig] : posG0[ig + 1] - posG0[ig];
      com.posG[0] = 0;  com.posG[1] = com.npatt;
      for (j = 0; j < com.ns; j++)  com.z[j] += posG0[ig] * nb;
      com.fpatt += posG0[ig];
      xtoy(com.piG[ig], com.pi, com.ncode);

      printf("\n\nGene %2d  ls:%4d  npatt:%4d\n", ig + 1, com.ls, com.npatt);
      fprintf(fout, "\nGene %2d  ls:%4d  npatt:%4d\n", ig + 1, com.ls, com.npatt);
      fprintf(frst, "\nGene %2d  ls:%4d  npatt:%4d\n", ig + 1, com.ls, com.npatt);
      fprintf(frst1, "%d\t%d\t%d", ig + 1, com.ls, com.npatt);

#ifdef CODEML
      if (com.seqtype == CODONseq) {
         DistanceMatNG86(fout, fpair[0], fpair[1], fpair[2], 0);
         if (com.codonf >= F1x4MG) com.pf3x4 = com.f3x4[ig];
      }
#else
      if (com.fix_alpha)
         DistanceMatNuc(fout, fpair[0], com.model, com.alpha);
#endif

      if (com.runmode == 0)  Forestry(fout, ftree);
#ifdef CODEML
      else if (com.runmode == -2) {
         if (com.seqtype == CODONseq)
            PairwiseCodon(fout, fpair[3], fpair[4], fpair[5], space);
         else
            PairwiseAA(fout, fpair[0]);
      }
#endif
      else
         StepwiseAddition(fout, space);

      for (j = 0; j < com.ns; j++) com.z[j] -= posG0[ig] * nb;
      com.fpatt -= posG0[ig];
      fprintf(frst1, "\n");
   }
   com.ngene = ngene0;
   com.npatt = npatt0;
   com.ls = lgene0[ngene0 - 1];
   for (ig = 0; ig < ngene0; ig++)
      com.lgene[ig] = lgene0[ig];
   for (ig = 0; ig < ngene0 + 1; ig++)
      com.posG[ig] = posG0[ig];
   return (0);
}

void printSeqsMgenes(void)
{
   /* separate sites from different partitions (genes) into different files.
    called before sequences are coded.
    Note that this is called before PatternWeight and so posec or posei is used
    and com.pose is not yet allocated.
    In case of codons, com.ls is the number of codons.
    */
   FILE* fseq;
   char seqf[20];
   int ig, lg, i, j, h;
   int n31 = (com.seqtype == CODONseq || com.seqtype == CODON2AAseq ? 3 : 1);

   puts("Separating sites in genes into different files.\n");
   printf("\n");
   for (ig = 0; ig < com.ngene; ig++) {
      for (h = 0, lg = 0; h < com.ls; h++)
         if (com.pose[h] == ig)
            lg++;
      sprintf(seqf, "Gene%d.seq", ig + 1);
      if ((fseq = fopen(seqf, "w")) == NULL) error2("file creation err.");
      printf("%d sites in gene %d go to file %s\n", lg, ig + 1, seqf);

      fprintf(fseq, "%8d%8d\n", com.ns, lg * n31);
      for (j = 0; j < com.ns; j++) {

         /* fprintf(fseq,"*\n>\n%s\n", com.spname[j]); */
         fprintf(fseq, "%-20s  ", com.spname[j]);
         if (n31 == 1) {       /* nucleotide or aa sequences */
            for (h = 0; h < com.ls; h++) {
               if (com.pose[h] == ig)
                  fprintf(fseq, "%c", com.z[j][h]);
            }
         }
         else {               /* codon sequences */
            for (h = 0; h < com.ls; h++) {
               if (com.pose[h] == ig) {
                  for (i = 0; i < 3; i++) fprintf(fseq, "%c", com.z[j][h * 3 + i]);
                  fputc(' ', fseq);
               }
            }
         }
         fprintf(fseq, "\n");
      }
      fclose(fseq);
   }
   return;
}

void printSeqsMgenes2(void)
{
   /* This print sites from certain genes into one file.
    called before sequences are coded.
    In case of codons, com.ls is the number of codons.
    */
   FILE* fseq;
   char seqf[20] = "newseqs";
   int ig, lg, i, j, h;
   int n31 = (com.seqtype == CODONseq || com.seqtype == CODON2AAseq ? 3 : 1);

   int ngenekept = 0;
   char* genenames[44] = { "atpa", "atpb", "atpe", "atpf", "atph", "petb", "petg", "psaa",
       "psab", "psac", "psaj", "psba", "psbb", "psbc", "psbd", "psbe",
       "psbf", "psbh", "psbi", "psbj", "psbk", "psbl", "psbn", "psbt",
       "rl14", "rl16", "rl2", "rl20", "rl36", "rpob", "rpoc", "rpod", "rs11",
       "rs12", "rs14", "rs18", "rs19", "rs2", "rs3", "rs4", "rs7", "rs8",
       "ycf4", "ycf9" };
   int wantgene[44] = { 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0 };
   /*
    for(ig=0,lg=0; ig<com.ngene; ig++) wantgene[ig]=!wantgene[ig];
    */

   if (com.ngene != 44) error2("ngene!=44");
   for (h = 0; h < com.ls; h++) {
      printf("%3d", com.pose[h]);
      if ((h + 1) % 20 == 0) printf("\n");
      if ((h + 1) % 500 == 0) getchar();
   }
   matIout(F0, com.lgene, 1, com.ngene);
   matIout(F0, wantgene, 1, com.ngene);

   for (ig = 0, lg = 0; ig < com.ngene; ig++)
      if (wantgene[ig]) { ngenekept++; lg += com.lgene[ig]; }

   if ((fseq = fopen(seqf, "w")) == NULL) error2("file creation err.");
   fprintf(fseq, "%4d %4d  G\nG  %d  ", com.ns, lg * n31, ngenekept);
   for (ig = 0; ig < com.ngene; ig++)
      if (wantgene[ig]) fprintf(fseq, " %3d", com.lgene[ig]);
   fprintf(fseq, "\n");

   for (j = 0; j < com.ns; j++) {
      fprintf(fseq, "%-20s  ", com.spname[j]);
      if (n31 == 1) {       /* nucleotide or aa sequences */
         for (h = 0; h < com.ls; h++)
            if (wantgene[ig = com.pose[h]]) fprintf(fseq, "%c", com.z[j][h]);
      }
      else {               /* codon sequences */
         for (h = 0; h < com.ls; h++)
            if (wantgene[ig = com.pose[h]]) {
               for (i = 0; i < 3; i++) fprintf(fseq, "%c", com.z[j][h * 3 + i]);
               fputc(' ', fseq);
            }
      }
      fprintf(fseq, "\n");
   }
   fprintf(fseq, "\n");
   for (ig = 0; ig < com.ngene; ig++)
      if (wantgene[ig]) fprintf(fseq, " %s", genenames[ig]);
   fprintf(fseq, "\n");
   fclose(fseq);

   exit(0);
}

#endif   /* ifdef REALSEQUENCE */
#endif   /* ifdef TREESEARCH */
#endif   /* ifdef NODESTRUCTURE */



#ifdef PARSIMONY

void UpPassScoreOnly(int inode);
void UpPassScoreOnlyB(int inode);

static int* Nsteps, * chUB;   /* MM */
static char* Kspace, * chU, * NchU;
/* Elements of chU are character states (there are NchU of them).  This
 representation is used to speed up calculation for large trees.
 Bit operations on chUB are performed for binary trees
 */

void UpPassScoreOnly(int inode)
{
   /* => VU, VL, & MM, theorem 2 */
   int ison, i, j;
   char* K = Kspace, maxK;  /* chMark (VV) not used in up pass */

   for (i = 0; i < nodes[inode].nson; i++)
      if (nodes[nodes[inode].sons[i]].nson > 0)
         UpPassScoreOnly(nodes[inode].sons[i]);

   for (i = 0; i < com.ncode; i++) K[i] = 0;
   for (i = 0; i < nodes[inode].nson; i++)
      for (j = 0, ison = nodes[inode].sons[i]; j < NchU[ison]; j++)
         K[(int)chU[ison * com.ncode + j]]++;
   for (i = 0, maxK = 0; i < com.ncode; i++)  if (K[i] > maxK) maxK = K[i];
   for (i = 0, NchU[inode] = 0; i < com.ncode; i++)
      if (K[i] == maxK)  chU[inode * com.ncode + NchU[inode]++] = (char)i;
   Nsteps[inode] = nodes[inode].nson - maxK;
   for (i = 0; i < nodes[inode].nson; i++)
      Nsteps[inode] += Nsteps[nodes[inode].sons[i]];
}

void UpPassScoreOnlyB(int inode)
{
   /* uses bit operation, for binary trees only
    */
   int ison1, ison2, i, change = 0;

   for (i = 0; i < nodes[inode].nson; i++)
      if (nodes[nodes[inode].sons[i]].nson > 0)
         UpPassScoreOnlyB(nodes[inode].sons[i]);

   ison1 = nodes[inode].sons[0];  ison2 = nodes[inode].sons[1];
   if ((chUB[inode] = (chUB[ison1] & chUB[ison2])) == 0) {
      chUB[inode] = (chUB[ison1] | chUB[ison2]);  change = 1;
   }
   Nsteps[inode] = change + Nsteps[ison1] + Nsteps[ison2];
}


double MPScore(double space[])
{
   /* calculates MP score for a given tree using Hartigan's (1973) algorithm.
    sizeof(space) = nnode*sizeof(int)+(nnode+2)*ncode*sizeof(char).
    Uses Nsteps[nnode], chU[nnode*ncode], NchU[nnode].
    if(BitOperation), bit operations are used on binary trees.
    */
   int h, i, BitOperation, U[3], change;
   double score;

   Nsteps = (int*)space;
   BitOperation = (tree.nnode == 2 * com.ns - 1 - (nodes[tree.root].nson == 3));
   BitOperation = (BitOperation && com.ncode < 32);
   if (BitOperation)  chUB = Nsteps + tree.nnode;
   else {
      chU = (char*)(Nsteps + tree.nnode);
      NchU = chU + tree.nnode * com.ncode;  Kspace = NchU + tree.nnode;
   }
   for (h = 0, score = 0; h < com.npatt; h++) {
      for (i = 0; i < tree.nnode; i++) Nsteps[i] = 0;
      if (BitOperation) {
         for (i = 0; i < com.ns; i++)  chUB[i] = 1 << (com.z[i][h]);
         UpPassScoreOnlyB(tree.root);
         if (nodes[tree.root].nson > 2) {
            for (i = 0; i < 3; i++) U[i] = chUB[nodes[tree.root].sons[i]];
            change = 2;
            if (U[0] & U[1] & U[2]) change = 0;
            else if (U[0] & U[1] || U[1] & U[2] || U[0] & U[2]) change = 1;
            for (i = 0, Nsteps[tree.root] = change; i < 3; i++)
               Nsteps[tree.root] += Nsteps[nodes[tree.root].sons[i]];
         }
      }
      else {                   /* polytomies, use characters */
         for (i = 0; i < com.ns; i++) {
            chU[i * com.ncode] = (char)(com.z[i][h]); NchU[i] = (char)1;
         }
         for (i = com.ns; i < tree.nnode; i++)  NchU[i] = 0;
         UpPassScoreOnly(tree.root);
      }
      score += Nsteps[tree.root] * com.fpatt[h];
      /*
       printf("\nh %3d:    ", h+1);
       for (i = 0; i < com.ns; i++) printf("%2d  ", com.z[i][h]);
       printf(" %6d ", Nsteps[tree.root]);
       if((h+1)%10==0) exit(1);
       */
   }

   return (score);
}

double RemoveMPNinfSites(double* nsiteNinf)
{
   /* Removes parsimony-noninformative sites and return the number of changes
    at those sites.
    Changes .z[], .fpatt[], .npatt, etc.
    */
   int  h, j, it, npatt0 = com.npatt, markb[NCODE], gt2;
   double MPScoreNinf;

   for (h = 0, com.npatt = 0, MPScoreNinf = 0, *nsiteNinf = 0; h < npatt0; h++) {
      for (j = 0; j < com.ncode; j++) markb[j] = 0;
      for (j = 0; j < com.ns; j++)  markb[(int)com.z[j][h]]++;
      for (j = 0, it = gt2 = 0; j < com.ncode; j++)
         if (markb[j] >= 2) { it++; gt2 = 1; }
      if (it < 2) {                         /* non-informative */
         *nsiteNinf += com.fpatt[h];
         for (j = 0; j < com.ncode; j++)
            if (markb[j] == 1) MPScoreNinf += com.fpatt[h];
         if (!gt2) MPScoreNinf -= com.fpatt[h];
      }
      else {
         for (j = 0; j < com.ns; j++)
            com.z[j][com.npatt] = com.z[j][h];
         com.fpatt[com.npatt++] = com.fpatt[h];
      }
   }
   return (MPScoreNinf);
}

#endif


#ifdef RECONSTRUCTION

static char* chMark, * chMarkU, * chMarkL; /* VV, VU, VL */
/* chMark, chMarkU, chMarkL (VV, VU, VL) have elements 0 or 1, marking
 whether the character state is present in the set */
static char* PATHWay, * NCharaCur, * ICharaCur, * CharaCur;
/* PATHWay, NCharaCur, ICharaCur, CharaCur are for the current
 reconstruction.
 */

int UpPass(int inode);
int DownPass(int inode);

int UpPass(int inode)
{
   /* => VU, VL, & MM, theorem 2 */
   int n = com.ncode, i, j;
   char* K = chMark, maxK;   /* chMark (VV) not used in up pass */

   for (i = 0; i < nodes[inode].nson; i++)
      if (nodes[nodes[inode].sons[i]].nson > 0)
         UpPass(nodes[inode].sons[i]);

   for (i = 0; i < n; i++) K[i] = 0;
   for (i = 0; i < nodes[inode].nson; i++) {
      for (j = 0; j < n; j++) {
         if (chMarkU[nodes[inode].sons[i] * n + j])
            K[j]++;
      }
   }
   for (i = 0, maxK = 0; i < n; i++)
      if (K[i] > maxK) maxK = K[i];
   for (i = 0; i < n; i++) {
      if (K[i] == maxK) 
         chMarkU[inode * n + i] = 1;
      else if (K[i] == maxK - 1)
         chMarkL[inode * n + i] = 1;
   }
   Nsteps[inode] = nodes[inode].nson - maxK;
   for (i = 0; i < nodes[inode].nson; i++)
      Nsteps[inode] += Nsteps[nodes[inode].sons[i]];
   return (0);
}

int DownPass(int inode)
{
   /* VU, VL => VV, theorem 3 */
   int n = com.ncode, i, j, ison;

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      for (j = 0; j < n; j++)
         if (chMark[inode * n + j] > chMarkU[ison * n + j])
            break;
      if (j == n)
         for (j = 0; j < n; j++)
            chMark[ison * n + j] = chMark[inode * n + j];
      else
         for (j = 0; j < n; j++)
            chMark[ison * n + j] = (char)(chMarkU[ison * n + j] || (chMark[inode * n + j] && chMarkL[ison * n + j]));
   }
   for (i = 0; i < nodes[inode].nson; i++)
      if (nodes[nodes[inode].sons[i]].nson > 0)
         DownPass(nodes[inode].sons[i]);
   return (0);
}


int DownStates(int inode)
{
   /* VU, VL => NCharaCur, CharaCur, theorem 4 */
   int i;

   for (i = 0; i < nodes[inode].nson; i++)
      if (nodes[inode].sons[i] >= com.ns)
         DownStatesOneNode(nodes[inode].sons[i], inode);
   return (0);
}

int DownStatesOneNode(int ison, int father)
{
   /* States down inode, given father */
   char chi = PATHWay[father - com.ns];
   int n = com.ncode, j, in;

   if ((in = ison - com.ns) < 0) return (0);
   if (chMarkU[ison * n + chi]) {
      NCharaCur[in] = 1;
      CharaCur[in * n + 0] = chi;
   }
   else if (chMarkL[ison * n + chi]) {
      for (j = 0, NCharaCur[in] = 0; j < n; j++)
         if (chMarkU[ison * n + j] || j == chi)
            CharaCur[in * n + NCharaCur[in]++] = (char)j;
   }
   else {
      for (j = 0, NCharaCur[in] = 0; j < n; j++)
         if (chMarkU[ison * n + j])
            CharaCur[in * n + NCharaCur[in]++] = (char)j;
   }
   PATHWay[in] = CharaCur[in * n + (ICharaCur[in] = 0)];

   for (j = 0; j < nodes[ison].nson; j++)
      if (nodes[ison].sons[j] >= com.ns) break;
   if (j < nodes[ison].nson) DownStates(ison);

   return (0);
}

int InteriorStatesMP(int job, int h, int* nchange, char NChara[NS - 1],
   char Chara[(NS - 1) * NCODE], double space[]);

int InteriorStatesMP(int job, int h, int* nchange, char NChara[NS - 1],
   char Chara[(NS - 1) * NCODE], double space[])
{
   /* sizeof(space) = nnode*sizeof(int)+3*nnode*ncode*sizeof(char)
    job: 0=# of changes; 1:equivocal states
    */
   int n = com.ncode, i, j;

   Nsteps = (int*)space;            chMark = (char*)(Nsteps + tree.nnode);
   chMarkU = chMark + tree.nnode * n;   chMarkL = chMarkU + tree.nnode * n;
   for (i = 0; i < tree.nnode; i++) Nsteps[i] = 0;
   for (i = 0; i < 3 * n * tree.nnode; i++) chMark[i] = 0;
   for (i = 0; i < com.ns; i++)
      chMark[i * n + com.z[i][h]] = chMarkU[i * n + com.z[i][h]] = 1;
   UpPass(tree.root);
   *nchange = Nsteps[tree.root];
   if (job == 0) return (0);
   for (i = 0; i < n; i++) 
      chMark[tree.root * n + i] = chMarkU[tree.root * n + i];
   DownPass(tree.root);
   for (i = 0; i < tree.nnode - com.ns; i++) {
      for (j = 0, NChara[i] = 0; j < n; j++) {
         if (chMark[(i + com.ns) * n + j])
            Chara[i * n + NChara[i]++] = (char)j;
      }
   }
   return (0);
}


int PathwayMP(FILE* fout, double space[])
{
   /* Hartigan, JA.  1973.  Minimum mutation fits to a given tree.
    Biometrics, 29:53-65.
    */
   char* pch = (com.seqtype == 0 ? BASEs : AAs), visit[NS - 1];
   int n = com.ncode, nid = tree.nbranch - com.ns + 1, it, i, j, k, h, npath;
   int nchange, nchange0;
   char nodeb[NNODE], Equivoc[NS - 1];

   PATHWay = (char*)malloc(nid * (n + 3) * sizeof(char));
   NCharaCur = PATHWay + nid;  ICharaCur = NCharaCur + nid;  CharaCur = ICharaCur + nid;

   for (j = 0, visit[i = 0] = (char)(tree.root - com.ns); j < tree.nbranch; j++)
      if (tree.branches[j][1] >= com.ns)
         visit[++i] = (char)(tree.branches[j][1] - com.ns);
   /*
    printf ("\nOrder in nodes: ");
    for (j = 0; j < nid; j++) printf ("%4d", visit[j]+1+com.ns); printf("\n");
    */
   for (h = 0; h < com.npatt; h++) {
      fprintf(fout, "\n%4d%6.0f  ", h + 1, com.fpatt[h]);
      for (j = 0; j < com.ns; j++) fprintf(fout, "%c", pch[(int)com.z[j][h]]);
      fprintf(fout, ":  ");

      for (j = 0; j < com.ns; j++)
         nodeb[j] = (char)(com.z[j][h]);

      InteriorStatesMP(1, h, &nchange, NCharaCur, CharaCur, space);
      ICharaCur[j = tree.root - com.ns] = 0;  PATHWay[j] = CharaCur[j * n + 0];
      for (j = 0; j < nid; j++)
         Equivoc[j] = (char)(NCharaCur[j] > 1);
      DownStates(tree.root);

      for (npath = 0; ;) {
         for (j = 0, k = visit[nid - 1]; j < NCharaCur[k]; j++) {
            PATHWay[k] = CharaCur[k * n + j]; npath++;
            for (i = 0; i < nid; i++)
               fprintf(fout, "%c", pch[(int)PATHWay[i]]);
            fprintf(fout, "  ");

            for (i = 0; i < nid; i++)
               nodeb[i + com.ns] = PATHWay[i];
            for (i = 0, nchange0 = 0; i < tree.nbranch; i++)
               nchange0 += (nodeb[tree.branches[i][0]] != nodeb[tree.branches[i][1]]);
            if (nchange0 != nchange)
            {
               puts("\a\nerr:PathwayMP"); fprintf(fout, ".%d. ", nchange0);
            }

         }
         for (j = nid - 2; j >= 0; j--) {
            if (Equivoc[k = visit[j]] == 0) continue;
            if (ICharaCur[k] + 1 < NCharaCur[k]) {
               PATHWay[k] = CharaCur[k * n + (++ICharaCur[k])];
               DownStates(k + com.ns);
               break;
            }
            else { /* if (next equivocal node is not ancestor) update node k */
               for (i = j - 1; i >= 0; i--) if (Equivoc[(int)visit[i]]) break;
               if (i >= 0) {
                  for (it = k + com.ns, i = visit[i] + com.ns; ; it = nodes[it].father)
                     if (it == tree.root || nodes[it].father == i) break;
                  if (it == tree.root)
                     DownStatesOneNode(k + com.ns, nodes[k + com.ns].father);
               }
            }
         }
         if (j < 0) break;
      }
      fprintf(fout, " |%4d (%d)", npath, nchange);
   }   /* for (h) */
   free(PATHWay);
   return (0);
}

#endif



#if(BASEML || CODEML)


int BootstrapSeq(char* seqf)
{
   /* This is called from within ReadSeq(), right after the sequences are read
    and before the data are coded.
    jackknife if(lsb<com.ls && com.ngene==1).
    gmark[start+19] marks the position of the 19th site in that gene.
    */
   int iboot, nboot = com.bootstrap, h, is, ig, lg[NGENE] = { 0 }, j, start;
   int lsb = com.ls, n31 = 1, gap = 10, gpos[NGENE];
   int* sites = (int*)malloc(com.ls * sizeof(int)), * gmark = NULL;
   FILE* fseq = (FILE*)gfopen(seqf, "w");
   enum { PAML = 0, PAUP };
   char* datatype = (com.seqtype == AAseq ? "protein" : "dna");
   char* paupstart = "paupstart", * paupblock = "paupblock", * paupend = "paupend";
   int  format = 0;  /* 0: paml-phylip; 1:paup-nexus */

   if (com.readpattern) error2("work on bootstrapping pattern data.");

   printf("\nGenerating bootstrap samples in file %s\n", seqf);
   if (format == PAUP) {
      printf("%s, %s, & %s will be appended if existent.\n",
         paupstart, paupblock, paupend);
      appendfile(fseq, paupstart);
   }

   if (com.seqtype == CODONseq || com.seqtype == CODON2AAseq) { n31 = 3; gap = 1; }
   if (sites == NULL) error2("oom in BootstrapSeq");
   if (com.ngene > 1) {
      if (lsb < com.ls) error2("jackknife when #gene>1");
      if ((gmark = (int*)malloc(com.ls * sizeof(int))) == NULL)
         error2("oom in BootstrapSeq");

      for (ig = 0; ig < com.ngene; ig++)  com.lgene[ig] = gpos[ig] = 0;
      for (h = 0; h < com.ls; h++)  com.lgene[com.pose[h]]++;
      for (j = 0; j < com.ngene; j++) lg[j] = com.lgene[j];
      for (j = 1; j < com.ngene; j++) com.lgene[j] += com.lgene[j - 1];

      if (noisy && com.ngene > 1) {
         printf("Bootstrap uses stratefied sampling for %d partitions.", com.ngene);
         printf("\nnumber of sites in each partition: ");
         for (ig = 0; ig < com.ngene; ig++) printf(" %4d", lg[ig]);
         printf("\n");
      }

      for (h = 0; h < com.ls; h++) {     /* create gmark[] */
         ig = com.pose[h];
         start = (ig == 0 ? 0 : com.lgene[ig - 1]);
         gmark[start + gpos[ig]++] = h;
      }
   }

   for (iboot = 0; iboot < nboot; iboot++) {
      if (com.ngene <= 1)
         for (h = 0; h < lsb; h++) sites[h] = (int)(rndu() * com.ls);
      else {
         for (ig = 0; ig < com.ngene; ig++) {
            start = (ig == 0 ? 0 : com.lgene[ig - 1]);
            for (h = 0; h < lg[ig]; h++)
               sites[start + h] = gmark[start + (int)(rndu() * lg[ig])];
         }
      }

      /* print out the bootstrap sample */
      if (format == PAUP) {
         fprintf(fseq, "\n\n[Replicate # %d]\n", iboot + 1);
         fprintf(fseq, "\nbegin data;\n");
         fprintf(fseq, "   dimensions ntax=%d nchar=%d;\n", com.ns, lsb * n31);
         fprintf(fseq, "   format datatype=%s missing=? gap=-;\n   matrix\n", datatype);

         for (is = 0; is < com.ns; is++) {
            fprintf(fseq, "%-20s  ", com.spname[is]);
            for (h = 0; h < lsb; h++) {
               for (j = 0; j < n31; j++) fprintf(fseq, "%c", com.z[is][sites[h] * n31 + j]);
               if ((h + 1) % gap == 0) fprintf(fseq, " ");
            }
            fprintf(fseq, "\n");
         }

         fprintf(fseq, "   ;\nend;");
         /* site partitions */
         if (com.ngene > 1) {
            fprintf(fseq, "\n\nbegin paup;\n");
            for (ig = 0; ig < com.ngene; ig++)
               fprintf(fseq, "   charset partition%-2d = %-4d - %-4d;\n",
                  ig + 1, (ig == 0 ? 1 : com.lgene[ig - 1] + 1), com.lgene[ig]);
            fprintf(fseq, "end;\n");
         }
         appendfile(fseq, paupblock);
      }
      else {
         if (com.ngene == 1)
            fprintf(fseq, "%6d %6d\n", com.ns, lsb * n31);
         else {
            fprintf(fseq, "%6d %6d  G\nG %d  ", com.ns, lsb * n31, com.ngene);
            for (ig = 0; ig < com.ngene; ig++)
               fprintf(fseq, " %4d", lg[ig]);
            fprintf(fseq, "\n\n");
         }
         for (is = 0; is < com.ns; is++) {
            fprintf(fseq, "%-20s  ", com.spname[is]);
            for (h = 0; h < lsb; h++) {
               for (j = 0; j < n31; j++)
                  fprintf(fseq, "%c", com.z[is][sites[h] * n31 + j]);
               if ((h + 1) % gap == 0) fprintf(fseq, " ");
            }
            fprintf(fseq, "\n");
         }
      }

      if (noisy && (iboot + 1) % 10 == 0) printf("\rdid sample #%d", iboot + 1);
      fprintf(fseq, "\n");
   }  /* for(iboot) */
   free(sites);  if (com.ngene > 1) free(gmark);
   fclose(fseq);
   return(0);
}



int rell(FILE* flnf, FILE* fout, int ntree)
{
   /* This implements three methods for tree topology comparison.  The first
    tests the log likelihood difference using a normal approximation
    (Kishino and Hasegawa 1989).  The second does approximate bootstrap sampling
    (the RELL method, Kishino and Hasegawa 1989, 1993).  The third is a
    modification of the K-H test with a correction for multiple comparison
    (Shimodaira and Hasegawa 1999) .
    The routine reads input from the file lnf.

    fpattB[npatt] stores the counts of site patterns in the bootstrap sample,
    with sitelist[ls] listing sites by gene, for stratefied sampling.

    com.space[ntree*(npatt+nr+5)]:
    lnf[ntree*npatt] lnL0[ntree] lnL[ntree*nr] pRELL[ntree] pSH[ntree] vdl[ntree]
    btrees[ntree]
    */
   char* line, timestr[64];
   int nr = (com.ls < 100000 ? 10000 : (com.ls < 10000 ? 5000 : 50));
   //   int nr = (com.ls < 100000 ? 10000 : (com.ls < 10000 ? 5000 : 500));
   int lline = 16000, ntree0, ns0 = com.ns, ls0, npatt0;
   int itree, h, ir, j, k, ig, mltree, nbtree, * btrees, status = 0;
   int* sitelist, * fpattB, * lgeneB, * psitelist;
   double* lnf, * lnL0, * lnL, * pRELL, * lnLmSH, * pSH, * vdl, y, mdl, smallv = 1e-5;
   size_t s;

   fflush(fout);
   puts("\nTree comparisons (Kishino & Hasegawa 1989; Shimodaira & Hasegawa 1999)");
   fputs("\nTree comparisons (Kishino & Hasegawa 1989; Shimodaira & Hasegawa 1999)\n", fout);
   fprintf(fout, "Number of replicates: %d\n", nr);

   fscanf(flnf, "%d%d%d", &ntree0, &ls0, &npatt0);
   if (ntree0 != -1 && ntree0 != ntree)  error2("rell: input data file strange.  Check.");
   if (ls0 != com.ls || npatt0 != com.npatt)
      error2("rell: input data file incorrect.");
   s = ntree * (com.npatt + nr + 5) * sizeof(double);
   if (com.sspace < s) {
      if (noisy) printf("resetting space to %zu bytes in rell.\n", s);
      com.sspace = s;
      if ((com.space = (double*)realloc(com.space, com.sspace)) == NULL)
         error2("oom space");
   }
   lnf = com.space; lnL0 = lnf + ntree * com.npatt; lnL = lnL0 + ntree; pRELL = lnL + ntree * nr;
   pSH = pRELL + ntree; vdl = pSH + ntree; btrees = (int*)(vdl + ntree);
   fpattB = (int*)malloc((com.npatt + com.ls + com.ngene) * sizeof(int));
   if (fpattB == NULL) error2("oom fpattB in rell.");
   sitelist = fpattB + com.npatt;  lgeneB = sitelist + com.ls;

   lline = (com.seqtype == 1 ? ns0 * 8 : ns0) + 100;
   lline = max2(16000, lline);
   if ((line = (char*)malloc((lline + 1) * sizeof(char))) == NULL) error2("oom rell");

   /* read lnf from file flnf, calculates lnL0[] & find ML tree */
   for (itree = 0, mltree = 0; itree < ntree; itree++) {
      printf("\r\tReading lnf for tree # %d", itree + 1);
      fscanf(flnf, "%d", &j);
      if (j != itree + 1)
      {
         printf("\nerr: lnf, reading tree %d.", itree + 1); return(-1);
      }
      for (h = 0, lnL0[itree] = 0; h < com.npatt; h++) {
         fscanf(flnf, "%d%d%lf", &j, &k, &y);
         if (j != h + 1)
         {
            printf("\nlnf, patt %d.", h + 1); return(-1);
         }
         fgets(line, lline, flnf);
         lnL0[itree] += com.fpatt[h] * (lnf[itree * com.npatt + h] = y);
      }
      if (itree && lnL0[itree] > lnL0[mltree]) mltree = itree;
   }
   printf(", done.\n");
   free(line);

   /* calculates SEs (vdl) by sitewise comparison */

   printtime(timestr);
   printf("\r\tCalculating SEs by sitewise comparison");
   for (itree = 0; itree < ntree; itree++) {
      if (itree == mltree) { vdl[itree] = 0; continue; }
      mdl = (lnL0[itree] - lnL0[mltree]) / com.ls;
      for (h = 0, vdl[itree] = 0; h < com.npatt; h++) {
         y = lnf[itree * com.npatt + h] - lnf[mltree * com.npatt + h];
         vdl[itree] += com.fpatt[h] * (y - mdl) * (y - mdl);
      }
      vdl[itree] = sqrt(vdl[itree]);
   }
   printf(", %s\n", printtime(timestr));

   /* bootstrap resampling */
   for (ig = 0; ig < com.ngene; ig++)
      lgeneB[ig] = (ig ? com.lgene[ig] - com.lgene[ig - 1] : com.lgene[ig]);
   for (h = 0, k = 0; h < com.npatt; h++)
      for (j = 0; j < (int)com.fpatt[h]; j++)
         sitelist[k++] = h;

   zero(pRELL, ntree); zero(pSH, ntree); zero(lnL, ntree * nr);
   for (ir = 0; ir < nr; ir++) {
      for (h = 0; h < com.npatt; h++) fpattB[h] = 0;
      for (ig = 0, psitelist = sitelist; ig < com.ngene; psitelist += lgeneB[ig++]) {
         for (k = 0; k < lgeneB[ig]; k++) {
            j = (int)(lgeneB[ig] * rndu());
            h = psitelist[j];
            fpattB[h]++;
         }
      }
      for (h = 0; h < com.npatt; h++) {
         if (fpattB[h])
            for (itree = 0; itree < ntree; itree++)
               lnL[itree * nr + ir] += fpattB[h] * lnf[itree * com.npatt + h];
      }

      /* y is the lnL for the best tree from replicate ir. */
      for (j = 1, nbtree = 1, btrees[0] = 0, y = lnL[ir]; j < ntree; j++) {
         if (fabs(lnL[j * nr + ir] - y) < smallv)
            btrees[nbtree++] = j;
         else if (lnL[j * nr + ir] > y)
         {
            nbtree = 1; btrees[0] = j; y = lnL[j * nr + ir];
         }
      }

      for (j = 0; j < nbtree; j++)
         pRELL[btrees[j]] += 1. / (nr * nbtree);
      if (nr > 100 && (ir + 1) % (nr / 100) == 0)
         printf("\r\tRELL Bootstrapping.. replicate: %6d / %d %s", ir + 1, nr, printtime(timestr));

   }
   free(fpattB);

   if (fabs(1 - sum(pRELL, ntree)) > 1e-6) error2("sum pRELL != 1.");

   /* Shimodaira & Hasegawa correction (1999), working on lnL[ntree*nr] */
   printf("\nnow doing S-H test");
   if ((lnLmSH = (double*)malloc(nr * sizeof(double))) == NULL) error2("oom in rell");
   for (j = 0; j < ntree; j++)  /* step 3: centering */
      for (ir = 0, y = sum(lnL + j * nr, nr) / nr; ir < nr; ir++) lnL[j * nr + ir] -= y;
   for (ir = 0; ir < nr; ir++) {
      for (j = 1, lnLmSH[ir] = lnL[ir]; j < ntree; j++)
         if (lnL[j * nr + ir] > lnLmSH[ir]) lnLmSH[ir] = lnL[j * nr + ir];
   }
   for (itree = 0; itree < ntree; itree++) {  /* steps 4 & 5 */
      for (ir = 0; ir < nr; ir++)
         if (lnLmSH[ir] - lnL[itree * nr + ir] > lnL0[mltree] - lnL0[itree])
            pSH[itree] += 1. / nr;
   }

   fprintf(fout, "\n%6s %12s %9s %9s%8s%10s%9s\n\n",
      "tree", "li", "Dli", " +- SE", "pKH", "pSH", "pRELL");
   for (j = 0; j < ntree; j++) {
      mdl = lnL0[j] - lnL0[mltree];
      if (j == mltree || fabs(vdl[j]) < 1e-6) { y = -1; pSH[j] = -1; status = -1; }
      else y = 1 - CDFNormal(-mdl / vdl[j]);
      fprintf(fout, "%6d%c%12.3f %9.3f %9.3f%8.3f%10.3f%9.3f\n",
         j + 1, (j == mltree ? '*' : ' '), lnL0[j], mdl, vdl[j], y, pSH[j], pRELL[j]);
   }

   fprintf(frst1, "%3d %12.6f", mltree + 1, lnL0[mltree]);
   for (j = 0; j < ntree; j++) fprintf(frst1, " %5.3f", pRELL[j]);
   /*
    for(j=0;j<ntree;j++) if(j!=mltree) fprintf(frst1,"%9.6f",pSH[j]);
    */

   fputs("\npKH: P value for KH normal test (Kishino & Hasegawa 1989)\n", fout);
   fputs("pRELL: RELL bootstrap proportions (Kishino & Hasegawa 1989)\n", fout);
   fputs("pSH: P value with multiple-comparison correction (MC in table 1 of Shimodaira & Hasegawa 1999)\n", fout);
   if (status) fputs("(-1 for P values means N/A)\n", fout);

   printf("\n");
   free(lnLmSH);
   return(0);
}

#endif




#ifdef LFUNCTIONS
#ifdef RECONSTRUCTION


void ListAncestSeq(FILE* fout, char* zanc);

void ListAncestSeq(FILE* fout, char* zanc)
{
   /* zanc[nintern*com.npatt] holds ancestral sequences.
    Extant sequences are coded if cleandata.
    */
   int wname = 15, j, h, n31 = (com.seqtype == CODONseq || com.seqtype == CODON2AAseq ? 3 : 1);
   int lst = (com.readpattern ? com.npatt : com.ls);

   fputs("\n\n\nList of extant and reconstructed sequences\n\n", fout);
   if (!com.readpattern) fprintf(fout, "%6d %6d\n\n", tree.nnode, lst * n31);
   else                 fprintf(fout, "%6d %6d  P\n\n", tree.nnode, lst * n31);
   for (j = 0; j < com.ns; j++) {
      fprintf(fout, "%-*s   ", wname, com.spname[j]);
      print1seq(fout, com.z[j], lst, com.pose);
      fprintf(fout, "\n");
   }
   for (j = 0; j < tree.nnode - com.ns; j++) {
      fprintf(fout, "node #%-*d  ", wname - 5, com.ns + j + 1);
      print1seq(fout, zanc + j * com.npatt, lst, com.pose);
      fprintf(fout, "\n");
   }
   if (com.readpattern) {
      fprintf(fout, "\n");
      for (h = 0; h < com.npatt; h++) {
         fprintf(fout, " %4.0f", com.fpatt[h]);
         if ((h + 1) % 15 == 0) fprintf(fout, "\n");
      }
      fprintf(fout, "\n\n");
   }
}

int ProbSitePattern(double x[], double* lnL, double fhsiteAnc[], double ScaleC[]);
int AncestralMarginal(FILE* fout, double x[], double fhsiteAnc[], double Sir[]);
int AncestralJointPPSG2000(FILE* fout, double x[]);


int ProbSitePattern(double x[], double* lnL, double fhsiteAnc[], double ScaleC[])
{
   /* This calculates probabilities for observing site patterns fhsite[].
    The following notes are for ncatG>1 and method = 0.
    The routine calculates the scale factor common to all site classes (ir),
    that is, the greatest of the scale factors among the ir classes.
    The common scale factors will be used in scaling nodes[].conP for all site
    classes for all nodes in PostProbNode().  Small conP for some site classes
    will be essentially set to 0, which is fine.

    fhsite[npatt]
    ScaleSite[npatt]

    Ziheng Yang, 7 Sept, 2001
*/
   int ig, i, k, h, ir;
   double fh, S, y = 1;

   if (com.ncatG > 1 && com.method == 1) error2("don't need this?");
   if (SetParameters(x)) puts("par err.");
   for (h = 0; h < com.npatt; h++)
      fhsiteAnc[h] = 0;
   if (com.ncatG <= 1) {
      for (ig = 0, *lnL = 0; ig < com.ngene; ig++) {
         if (com.Mgene > 1) SetPGene(ig, 1, 1, 0, x);
         ConditionalPNode(tree.root, ig, x);
         for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
            for (i = 0; i < com.ncode; i++)
               fhsiteAnc[h] += com.pi[i] * nodes[tree.root].conP[h * com.ncode + i];
            *lnL -= log(fhsiteAnc[h]) * com.fpatt[h];
            if (com.NnodeScale)
               for (k = 0; k < com.NnodeScale; k++)
                  *lnL -= com.nodeScaleF[k * com.npatt + h] * com.fpatt[h];
         }
      }
   }
   else {
      for (ig = 0; ig < com.ngene; ig++) {
         if (com.Mgene > 1 || com.nalpha > 1)
            SetPGene(ig, com.Mgene > 1, com.Mgene > 1, com.nalpha > 1, x);
         for (ir = 0; ir < com.ncatG; ir++) {
#ifdef CODEML
            if (com.seqtype == 1 && com.NSsites /* && com.model */) IClass = ir;
#endif
            SetPSiteClass(ir, x);
            ConditionalPNode(tree.root, ig, x);

            for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
               for (i = 0, fh = 0; i < com.ncode; i++)
                  fh += com.pi[i] * nodes[tree.root].conP[h * com.ncode + i];

               if (com.NnodeScale) {
                  for (k = 0, S = 0; k < com.NnodeScale; k++)  S += com.nodeScaleF[k * com.npatt + h];
                  y = 1;
                  if (ir == 0)               ScaleC[h] = S;
                  else if (S <= ScaleC[h])   y = exp(S - ScaleC[h]);
                  else {    /* change of scale factor */
                     fhsiteAnc[h] *= exp(ScaleC[h] - S);  ScaleC[h] = S;
                  }
               }
               fhsiteAnc[h] += com.freqK[ir] * fh * y;
            }
         }
      }
      for (h = 0, *lnL = 0; h < com.npatt; h++)
         *lnL -= log(fhsiteAnc[h]) * com.fpatt[h];
      if (com.NnodeScale)
         for (h = 0; h < com.npatt; h++)
            *lnL -= ScaleC[h] * com.fpatt[h];
   }
   if (noisy) printf("\nlnL = %12.6f from ProbSitePattern.\n", -*lnL);

   return (0);
}


int updateconP(double x[], int inode);

int PostProbNode(int inode, double x[], double fhsiteAnc[], double ScaleC[],
   double* lnL, double pChar1node[], char za[], double pnode[])
{
   /* This calculates the full posterior distribution for node inode at each site.
    Below are special comments on gamma models and method = 0.

    Marginal reconstruction under gamma models, with complications arising from
    scaling on large trees (com.NnodeScale) and the use of two iteration algorithms
    (method).
    Z. Yang Sept 2001

    The algorithm is different depending on method, which makes the code clumsy.

    gamma method=0 or 2 (simultaneous updating):
    nodes[].conP overlap and get destroyed for different site classes (ir)
    The same for scale factors com.nodeScaleF.
    fhsite[npatt] and common scale factors ScaleC[npatt] are calculated for all
    nodes before this routine is called.  The common scale factors are then
    used to adjust nodes[].conP before they are summed across ir classes.

    gamma method=1 (one branch at a time):
    nodes[].conP (and com.nodeScaleF if node scaling is on) are separately
    allocated for different site classes (ir), so that all info needed is
    available.  Use of updateconP() saves computation on large trees.
    Scale factor Sir[] is of size ncatG and reused for each h.
    */
   int n = com.ncode, i, k, h, ir, it = -1, best, ig;
   double fh, y, pbest, * Sir = ScaleC, S;

   *lnL = 0;
   zero(pChar1node, com.npatt * n);

   /* nodes[].conP are reused for different ir, with or without node scaling */
   if (com.ncatG > 1 && com.method != 1) {
      ReRootTree(inode);
      for (ig = 0; ig < com.ngene; ig++) {
         if (com.Mgene > 1 || com.nalpha > 1)
            SetPGene(ig, com.Mgene > 1, com.Mgene > 1, com.nalpha > 1, x);
         for (ir = 0; ir < com.ncatG; ir++) {
#ifdef CODEML
            if (com.seqtype == 1 && com.NSsites)  IClass = ir;
#endif
            SetPSiteClass(ir, x);
            ConditionalPNode(tree.root, ig, x);

            for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
               if (!com.NnodeScale) S = 1;
               else {
                  for (k = 0, S = 0; k < com.NnodeScale; k++)
                     S += com.nodeScaleF[k * com.npatt + h];
                  S = exp(S - ScaleC[h]);
               }
               for (i = 0, fh = 0; i < n; i++) {
                  y = com.freqK[ir] * com.pi[i] * nodes[tree.root].conP[h * n + i] * S;
                  fh += y;
                  pChar1node[h * n + i] += y;
               }
            }
         }
      }
      for (h = 0; h < com.npatt; h++) {
         for (i = 0, y = 0; i < n; i++) y += (pChar1node[h * n + i] /= fhsiteAnc[h]);
         if (fabs(1 - y) > 1e-5)
            error2("PostProbNode: sum!=1");
         for (i = 0, best = -1, pbest = -1; i < n; i++)
            if (pChar1node[h * n + i] > pbest) {
               best = i;
               pbest = pChar1node[h * n + i];
            }
         za[(inode - com.ns) * com.npatt + h] = (char)best;
         pnode[(inode - com.ns) * com.npatt + h] = pbest;
         *lnL -= log(fhsiteAnc[h]) * com.fpatt[h];
         if (com.NnodeScale) *lnL -= ScaleC[h] * com.fpatt[h];
      }
   }
   else {  /* all other cases: (alpha==0 || method==1) */
      for (i = 0; i < tree.nnode; i++) com.oldconP[i] = 1;
      ReRootTree(inode);
      updateconP(x, inode);
      if (com.alpha == 0 && com.ncatG <= 1) { /* (alpha==0) (ngene>1 OK) */
         for (ig = 0; ig < com.ngene; ig++) {
            if (com.Mgene == 2 || com.Mgene == 4)
               xtoy(com.piG[ig], com.pi, n);
            for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
               for (i = 0, fh = 0, pbest = 0, best = -1; i < n; i++) {
                  y = com.pi[i] * nodes[tree.root].conP[h * n + i];
                  fh += y;
                  if (y > pbest)
                  {
                     pbest = y; best = i;
                  }
                  pChar1node[h * n + i] = y;
               }
               za[(inode - com.ns) * com.npatt + h] = (char)best;
               pnode[(inode - com.ns) * com.npatt + h] = (pbest /= fh);
               for (i = 0; i < n; i++)
                  pChar1node[h * n + i] /= fh;
               *lnL -= log(fh) * (double)com.fpatt[h];
               for (i = 0; i < com.NnodeScale; i++)
                  *lnL -= com.nodeScaleF[i * com.npatt + h] * com.fpatt[h];
            }
         }
      }
      else {  /* (ncatG>1 && method = 1)  This should work for NSsites? */
         for (ig = 0; ig < com.ngene; ig++) {
            if (com.Mgene == 2 || com.Mgene == 4)
               xtoy(com.piG[ig], com.pi, n);
            for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
               if (com.NnodeScale)
                  for (ir = 0, it = 0; ir < com.ncatG; ir++) {  /* Sir[it] is the biggest */
                     for (k = 0, Sir[ir] = 0; k < com.NnodeScale; k++)
                        Sir[ir] += com.nodeScaleF[ir * com.NnodeScale * com.npatt + k * com.npatt + h];
                     if (Sir[ir] > Sir[it]) it = ir;
                  }
               for (i = 0, fh = 0; i < n; i++) {
                  for (ir = 0; ir < com.ncatG; ir++) {
                     if (com.method == 1)
                        y = nodes[tree.root].conP[ir * (tree.nnode - com.ns) * com.npatt * n + h * n + i];
                     else
                        y = nodes[tree.root].conP[h * n + i]; /* wrong right now */
                     y *= com.pi[i] * com.freqK[ir];
                     if (com.NnodeScale) y *= exp(Sir[ir] - Sir[it]);

                     pChar1node[h * n + i] += y;
                     fh += y;
                  }
               }
               for (i = 0, best = 0; i < n; i++) {
                  pChar1node[h * n + i] /= fh;
                  if (i && pChar1node[h * n + best] < pChar1node[h * n + i])
                     best = i;
               }
               za[(inode - com.ns) * com.npatt + h] = (char)best;
               pnode[(inode - com.ns) * com.npatt + h] = pChar1node[h * n + best];
               *lnL -= log(fh) * com.fpatt[h];
               if (com.NnodeScale) *lnL -= Sir[it] * com.fpatt[h];
            }
         }
      }
   }
   return(0);
}


void getCodonNode1Site(char codon[], char zanc[], int inode, int site);

int AncestralMarginal(FILE* fout, double x[], double fhsiteAnc[], double Sir[])
{
   /* Ancestral reconstruction for each interior node.  This works under both
    the one rate and gamma rates models.
    pnode[npatt*nid] stores the prob for the best chara at a node and site.
    The best character is kept in za[], coded as 0,...,n-1.
    The data may be coded (com.cleandata==1) or not (com.cleandata==0).
    Call ProbSitePatt() before running this routine.
    pMAPnode[NS-1], pMAPnodeA[] stores the MAP probabilities (accuracy)
    for a site and for the entire sequence, respectively.

    The routine PostProbNode calculates pChar1node[npatt*ncode], which stores
    prob for each char at each pattern at each given node inode.  The rest of
    the routine is to output the results in different ways.

    Deals with node scaling to avoid underflows.  See above
    (Z. Yang, 2 Sept 2001)
    */
   char* pch = (com.seqtype == 0 ? BASEs : (com.seqtype == 2 ? AAs : BINs));
   char* zanc, str[4] = "", codon[2][4] = { "   ","   " }, aa[4] = "";
   char* sitepatt = (com.readpattern ? "pattern" : "site");
   int n = com.ncode, inode, ic = 0, b[3], i, j, k1 = -1, k2 = -1, k3, lsc = com.ls;
   int lst = (com.readpattern ? com.npatt : com.ls);
   int h, hp, ig, best, oldroot = tree.root;
   int nid = tree.nnode - com.ns, nchange;
   double lnL = 0, fh, y, pbest, * pChar1node, * pnode, p1 = -1, p2 = -1;
   double pMAPnode[NS - 1], pMAPnodeA[NS - 1], smallp = 0.001;

   char coding = 0, * bestAA = NULL;
   double pAA[21], * pbestAA = NULL, ns, na, nst, nat, S, N;
   /* bestAA[nid*npatt], pbestAA[nid*npatt]:
    To reconstruct aa seqs using codon or nucleotide seqs, universal code */

   if (noisy) puts("Marginal reconstruction.");

   fprintf(fout, "\n(1) Marginal reconstruction of ancestral sequences\n");
   fprintf(fout, "(eqn. 4 in Yang et al. 1995 Genetics 141:1641-1650).\n");
   pChar1node = (double*)malloc(com.npatt * n * sizeof(double));
   pnode = (double*)malloc((nid * com.npatt + 1) * (sizeof(double) + sizeof(char)));
   if (pnode == NULL || pChar1node == NULL)
      error2("oom pnode");
   zanc = (char*)(pnode + nid * com.npatt);

#ifdef BASEML
   if (com.seqtype == 0 && com.ls % 3 == 0 && com.coding) { coding = 1; lsc = com.ls / 3; }
#endif
   if (com.seqtype == 1) { coding = 1; lsc = com.npatt; }
   if (coding == 1) {
      if ((pbestAA = (double*)malloc(nid * lsc * 2 * sizeof(double))) == NULL)
         error2("oom pbestAA");
      bestAA = (char*)(pbestAA + nid * lsc);
   }

   if (SetParameters(x)) puts("par err.");

   if (com.verbose > 1)
      fprintf(fout, "\nProb distribs at nodes, those with p < %.3f not listed\n", smallp);

   /* This loop reroots the tree at inode & reconstructs sequence at inode */
   for (inode = com.ns; inode < tree.nnode; inode++) {
      PostProbNode(inode, x, fhsiteAnc, Sir, &lnL, pChar1node, zanc, pnode);
      if (noisy) printf("\tNode %3d: lnL = %12.6f\n", inode + 1, -lnL);

      /* print Prob distribution at inode if com.verbose>1 */
      if (com.verbose > 1) {
         fprintf(fout, "\nProb distribution at node %d, by %s\n", inode + 1, sitepatt);
         fprintf(fout, "\n%7s  Freq   Data\n\n", sitepatt);
         for (h = 0; h < lst; h++) {
            hp = (!com.readpattern ? com.pose[h] : h);
            fprintf(fout, "%7d%7.0f   ", h + 1, com.fpatt[hp]);
            print1site(fout, hp);
            fputs(": ", fout);
            for (j = 0; j < n; j++) {
               if (com.seqtype != CODONseq) {
                  str[0] = pch[j];
                  str[1] = 0;
               }
               else
                  strcpy(str, CODONs[j]);
               fprintf(fout, "%s(%5.3f) ", str, pChar1node[hp * n + j]);
            }
            fprintf(fout, "\n");
         }
      }     /* if (verbose) */


      /* find the best amino acid for coding seqs */
#ifdef CODEML
      if (com.seqtype == CODONseq)
         for (h = 0; h < com.npatt; h++) {
            for (j = 0; j < 20; j++) pAA[j] = 0;
            for (j = 0; j < n; j++) {
               i = GeneticCode[com.icode][FROM61[j]];
               pAA[i] += pChar1node[h * n + j];
            }
            /* matout(F0,pAA,1,20); */
            for (j = 0, best = 0, pbest = 0; j < 20; j++)
               if (pAA[j] > pbest) { pbest = pAA[j]; best = j; }
            bestAA[(inode - com.ns) * com.npatt + h] = (char)best;
            pbestAA[(inode - com.ns) * com.npatt + h] = pbest;
         }
#endif
      if (com.seqtype == 0 && coding) { /* coding seqs analyzed by baseml */
         for (h = 0; h < lsc; h++) {  /* h-th codon */
             /* sums up probs for the 20 AAs for each node. Stop codons are
              ignored, and so those probs are approxiamte. */
            for (j = 0, y = 0; j < 20; j++) pAA[j] = 0;
            for (k1 = 0; k1 < 4; k1++) for (k2 = 0; k2 < 4; k2++) for (k3 = 0; k3 < 4; k3++) {
               ic = k1 * 16 + k2 * 4 + k3;
               b[0] = com.pose[h * 3 + 0] * n + k1;
               b[1] = com.pose[h * 3 + 1] * n + k2;
               b[2] = com.pose[h * 3 + 2] * n + k3;
               fh = pChar1node[b[0]] * pChar1node[b[1]] * pChar1node[b[2]];
               if ((ic = GeneticCode[com.icode][ic]) == -1)
                  y += fh;
               else
                  pAA[ic] += fh;
            }
            if (fabs(1 - y - sum(pAA, 20)) > 1e-6)
               error2("AncestralMarginal strange?");

            for (j = 0, best = 0, pbest = 0; j < 20; j++)
               if (pAA[j] > pbest) { pbest = pAA[j]; best = j; }

            bestAA[(inode - com.ns) * com.ls / 3 + h] = (char)best;
            pbestAA[(inode - com.ns) * com.ls / 3 + h] = pbest / (1 - y);
         }
      }
   }        /* for (inode), This closes the big loop */

   for (i = 0; i < tree.nnode; i++)
      com.oldconP[i] = 0;
   ReRootTree(oldroot);

   if (com.seqtype == 0 && coding && !com.readpattern) { /* coding seqs analyzed by baseml */
      fputs("\nBest amino acids reconstructed from nucleotide model.\n", fout);
      fputs("Prob at each node listed by amino acid (codon) site\n", fout);
      fputs("(Please ignore if not relevant)\n\n", fout);
      for (h = 0; h < com.ls / 3; h++) {
         fprintf(fout, "%4d ", h + 1);
         for (j = 0; j < com.ns; j++) {
            getCodonNode1Site(codon[0], NULL, j, h);
            Codon2AA(codon[0], aa, com.icode, &i);
            fprintf(fout, " %s(%c)", codon[0], AAs[i]);
         }
         fprintf(fout, ": ");
         for (j = 0; j < tree.nnode - com.ns; j++) {
            fprintf(fout, " %1c (%5.3f)", AAs[(int)bestAA[j * com.ls / 3 + h]], pbestAA[j * com.ls / 3 + h]);
         }
         fprintf(fout, "\n");
      }
   }

   /* calculate accuracy measures */
   zero(pMAPnode, nid);  fillxc(pMAPnodeA, 1., nid);
   for (inode = 0; inode < tree.nnode - com.ns; inode++) {
      for (h = 0; h < com.npatt; h++) {
         pMAPnode[inode] += com.fpatt[h] * pnode[inode * com.npatt + h] / com.ls;
         pMAPnodeA[inode] *= pow(pnode[inode * com.npatt + h], com.fpatt[h]);
      }
   }

   fprintf(fout, "\nProb of best state at each node, listed by %s", sitepatt);
   if (com.ngene > 1) fprintf(fout, "\n\n%7s (g) Freq  Data: \n", sitepatt);
   else               fprintf(fout, "\n\n%7s   Freq   Data: \n", sitepatt);

   for (h = 0; h < lst; h++) {
      hp = (!com.readpattern ? com.pose[h] : h);
      fprintf(fout, "\n%4d ", h + 1);
      if (com.ngene > 1) {  /* which gene the site is from */
         for (ig = 1; ig < com.ngene; ig++)
            if (hp < com.posG[ig]) break;
         fprintf(fout, "(%d)", ig);
      }
      fprintf(fout, " %5.0f   ", com.fpatt[hp]);
      print1site(fout, hp);
      fprintf(fout, ": ");

      for (j = 0; j < nid; j++) {
         if (com.seqtype != CODONseq)
            fprintf(fout, "%c(%5.3f) ", pch[(int)zanc[j * com.npatt + hp]], pnode[j * com.npatt + hp]);
#ifdef CODEML
         else {
            ic = zanc[j * com.npatt + hp];
            Codon2AA(CODONs[ic], aa, com.icode, &i);
            fprintf(fout, " %s %1c %5.3f (%1c %5.3f)",
               CODONs[ic], AAs[i], pnode[j * com.npatt + hp], AAs[(int)bestAA[j * com.npatt + hp]], pbestAA[j * com.npatt + hp]);
         }
#endif
      }
      if (noisy && (h + 1) % 100000 == 0) printf("\r\tprinting, %d sites done", h + 1);
   }
   if (noisy && h >= 100000) printf("\n");

   /* Map changes onto branches
    k1 & k2 are the two characters; p1 and p2 are the two probs. */

   if (!com.readpattern) {
      fputs("\n\nSummary of changes along branches.\n", fout);
      fputs("Check root of tree for directions of change.\n", fout);
      if (!com.cleandata && com.seqtype == 1)
         fputs("Counts of n & s are incorrect along tip branches with ambiguity data.\n", fout);
      for (j = 0; j < tree.nbranch; j++) {
         inode = tree.branches[j][1];
         nchange = 0;
         fprintf(fout, "\nBranch %d:%5d..%-2d", j + 1, tree.branches[j][0] + 1, inode + 1);
         if (inode < com.ns) fprintf(fout, " (%s) ", com.spname[inode]);

         if (coding) {
            lsc = (com.seqtype == 1 ? com.ls : com.ls / 3);
            for (h = 0, nst = nat = 0; h < lsc; h++) {
               getCodonNode1Site(codon[0], zanc, inode, h);
               getCodonNode1Site(codon[1], zanc, tree.branches[j][0], h);
               if (strcmp(codon[0], codon[1])) {
                  difcodonNG(codon[0], codon[1], &S, &N, &ns, &na, 0, com.icode);
                  nst += ns;
                  nat += na;
               }
            }
            fprintf(fout, " (n=%4.1f s=%4.1f)", nat, nst);
         }
         fprintf(fout, "\n\n");
         for (h = 0; h < lst; h++) {
            hp = (!com.readpattern ? com.pose[h] : h);
            if (com.seqtype != CODONseq) {
               if (inode < com.ns)
                  k2 = pch[(int)com.z[inode][hp]];
               else {
                  k2 = pch[(int)zanc[(inode - com.ns) * com.npatt + hp]];
                  p2 = pnode[(inode - com.ns) * com.npatt + hp];
               }
               k1 = pch[(int)zanc[(tree.branches[j][0] - com.ns) * com.npatt + hp]];
               p1 = pnode[(tree.branches[j][0] - com.ns) * com.npatt + hp];
            }
#ifdef CODEML
            else {
               if (inode < com.ns) {
                  strcpy(codon[1], CODONs[(int)com.z[inode][hp]]);
                  k2 = GetAASiteSpecies(inode, hp);
               }
               else {
                  strcpy(codon[1], CODONs[(int)zanc[(inode - com.ns) * com.npatt + hp]]);
                  k2 = AAs[(int)bestAA[(inode - com.ns) * com.npatt + hp]];
                  p2 = pbestAA[(inode - com.ns) * com.npatt + hp];
               }
               strcpy(codon[0], CODONs[(int)zanc[(tree.branches[j][0] - com.ns) * com.npatt + hp]]);
               k1 = AAs[(int)bestAA[(tree.branches[j][0] - com.ns) * com.npatt + hp]];
               p1 = pbestAA[(tree.branches[j][0] - com.ns) * com.npatt + hp];

               /* best codon may not match best AA */
               if (strcmp(codon[0], codon[1])) {
                  if (inode < com.ns)
                     fprintf(fout, "\t%4d %s (%c %.3f) -> %s (%c)\n", h + 1, codon[0], k1, p1, codon[1], k2);
                  else
                     fprintf(fout, "\t%4d %s (%c %.3f) -> %s (%c) %.3f\n", h + 1, codon[0], k1, p1, codon[1], k2, p2);
               }
               k1 = k2 = 0;
            }
#endif
            if (k1 == k2) continue;
            fprintf(fout, "\t%4d ", h + 1);

#ifdef SITELABELS
            if (sitelabels) fprintf(fout, " %5s   ", sitelabels[h]);
#endif
            if (inode < com.ns) fprintf(fout, "%c %.3f -> %1c\n", k1, p1, k2);
            else             fprintf(fout, "%c %.3f -> %1c %.3f\n", k1, p1, k2, p2);
            nchange++;
         }
         fprintf(fout, "\n");
      }
   }

   ListAncestSeq(fout, zanc);
   fprintf(fout, "\n\nOverall accuracy of the %d ancestral sequences:", nid);
   matout2(fout, pMAPnode, 1, nid, 9, 5);  fputs("for a site.\n", fout);
   matout2(fout, pMAPnodeA, 1, nid, 9, 5); fputs("for the sequence.\n", fout);

   /* best amino acid sequences from codonml */
#ifdef CODEML
   if (com.seqtype == 1) {
      fputs("\n\nAmino acid sequences inferred by codonml.\n", fout);
      if (!com.cleandata)
         fputs("Results unreliable for sites with alignment gaps.\n", fout);
      for (inode = 0; inode < nid; inode++) {
         fprintf(fout, "\nNode #%-10d  ", com.ns + inode + 1);
         for (h = 0; h < lst; h++) {
            hp = (!com.readpattern ? com.pose[h] : h);
            fprintf(fout, "%c", AAs[(int)bestAA[inode * com.npatt + hp]]);
            if ((h + 1) % 10 == 0) fputc(' ', fout);
         }
      }
      fprintf(fout, "\n");
   }
#endif
   ChangesSites(fout, coding, zanc);

   free(pnode);
   free(pChar1node);
   if (coding) free(pbestAA);
   return (0);
}


void getCodonNode1Site(char codon[], char zanc[], int inode, int site)
{
   /* this is used to retrive the codon from a codon sequence for codonml
    or coding sequence in baseml, used in ancestral reconstruction
    zanc has ancestral sequences
    site is codon site
    */
   int i;

   for (i = 0; i < 3; i++)  /* to force crashes */
      codon[i] = -1;
   if (com.seqtype == CODONseq) {
#ifdef CODEML
      int hp = (!com.readpattern ? com.pose[site] : site);
      if (inode >= com.ns)
         strcpy(codon, CODONs[(int)zanc[(inode - com.ns) * com.npatt + hp]]);
      else
         strcpy(codon, CODONs[(int)com.z[inode][hp]]);
#endif
   }
   else {      /* baseml coding reconstruction */
      if (inode >= com.ns)
         for (i = 0; i < 3; i++)
            codon[i] = BASEs[(int)zanc[(inode - com.ns) * com.npatt + com.pose[site * 3 + i]]];
      else
         for (i = 0; i < 3; i++)
            codon[i] = BASEs[(int)com.z[inode][com.pose[site * 3 + i]]];
   }

}

int ChangesSites(FILE* frst, int coding, char* zanc)
{
   /* this lists and counts changes at sites from reconstructed ancestral sequences
    com.z[] has the data, and zanc[] has the ancestors
    For codon sequences (codonml or baseml with com.coding), synonymous and
    nonsynonymous changes are counted separately.
    Added in Nov 2000.
    */
   char* pch = (com.seqtype == 0 ? BASEs : (com.seqtype == 2 ? AAs : BINs));
   char codon[2][4] = { "   ","   " };
   int  h, hp, inode, k1, k2, d, ls1 = (com.readpattern ? com.npatt : com.ls);
   double S, N, Sd, Nd, S1, N1, Sd1, Nd1, b, btotal = 0, p, C;

   if (com.seqtype == 0 && coding) ls1 /= 3;
   if (coding) {
      fprintf(frst, "\n\nCounts of changes at sites, listed by %s\n\n",
         (com.readpattern ? "pattern" : "site"));
      fprintf(frst1, "\nList of sites with changes according to ancestral reconstruction\n");
      fprintf(frst1, "Suzuki-Gojobori (1999) style test\n");
      if (!com.cleandata)
         fprintf(frst, "(Counts of n & s are incorrect at sites with ambiguity data)\n\n");

      for (inode = 0; inode < tree.nnode; inode++)
         if (inode != tree.root) btotal += nodes[inode].branch;
      for (h = 0; h < ls1; h++) {
         fprintf(frst, "%4d ", h + 1);
         for (inode = 0, S = N = Sd = Nd = 0; inode < tree.nnode; inode++) {
            if (inode == tree.root) continue;
            b = nodes[inode].branch;
            getCodonNode1Site(codon[0], zanc, nodes[inode].father, h);
            getCodonNode1Site(codon[1], zanc, inode, h);

            difcodonNG(codon[0], codon[1], &S1, &N1, &Sd1, &Nd1, 0, com.icode);
            S += S1 * b / btotal;
            N += N1 * b / btotal;
            if (Sd1 || Nd1) {
               Sd += Sd1;
               Nd += Nd1;
               fprintf(frst, " %3s.%3s ", codon[0], codon[1]);
            }
         }
         b = S + N; S /= b;  N /= b;
         fprintf(frst, "(S N: %7.3f%7.3f Sd Nd: %6.1f %5.1f)\n", S * 3, N * 3, Sd, Nd);
         fprintf(frst1, "%4d S N: %7.3f%7.3f Sd Nd: %6.1f %5.1f ", h + 1, S * 3, N * 3, Sd, Nd);
         if (Sd + Nd) {
            if (Nd / (Sd + Nd) < N) {
               for (d = 0, p = 0, C = 1; d <= Nd; d++) {
                  p += C * pow(N, d) * pow(1 - N, Sd + Nd - d);
                  C *= (Sd + Nd - d) / (d + 1);
               }
               fprintf(frst1, " - p =%6.3f %s", p, (p < .01 ? "**" : (p < .05 ? "*" : "")));
            }
            else {
               for (d = 0, p = 0, C = 1; d <= Sd; d++) {
                  p += C * pow(S, d) * pow(1 - S, Sd + Nd - d);
                  C *= (Sd + Nd - d) / (d + 1);
               }
               fprintf(frst1, " + p =%6.3f %s", p, (p < .01 ? "**" : (p < .05 ? "*" : "")));
            }
         }
         fprintf(frst1, "\n");
      }
   }
   else {  /* noncoding nucleotide or aa sequences */
      fprintf(frst, "\n\nCounts of changes at sites%s\n\n",
         (com.readpattern ? ", listed by pattern" : ""));
      for (h = 0; h < ls1; h++) {
         hp = (!com.readpattern ? com.pose[h] : h);
         fprintf(frst, "%4d ", h + 1);
         for (inode = 0, d = 0; inode < tree.nnode; inode++) {
            if (inode == tree.root) continue;
            k1 = pch[(int)zanc[(nodes[inode].father - com.ns) * com.npatt + hp]];
            if (inode < com.ns)
               k2 = pch[(int)com.z[inode][hp]];
            else
               k2 = pch[(int)zanc[(inode - com.ns) * com.npatt + hp]];
            if (k1 != k2) {
               d++;
               fprintf(frst, " %c%c", k1, k2);
            }
         }
         fprintf(frst, " (%d)\n", d);
      }
   }
   return(0);
}


/* The code here for joint ancestral reconstruction keep not just the best state
   but the best 4 states (when NBESTANC = 4) at each stage of the dynamic programming
   algorithm (Pupko et al. 2000).
*/

#define  NBESTANC  4  /* use 1 2 3 or 4 */
int  parsimony = 0, * nBestScore, * icharNode[NBESTANC], * combIndex;
double* fhsiteAnc, * lnPanc[NBESTANC], * PMatTips, * combScore;
char* charNode[NBESTANC], * ancSeq, * ancState1site;
FILE* fanc;
int largeReconstruction;

void DownPassPPSG2000OneSite(int h, int inode, int inodestate, int ipath);
void PrintAncState1site(char ancState1site[], double prob);


double P0[16] = { 0, 1, 1.5, 1.5,
    1, 0, 1.5, 1.5,
    1.5, 1.5, 0, 1,
    1.5, 1.5, 1, 0 };

double piroot[NCODE] = { 0 };

/* combIndex[] uses two bits for each son to record the path that is taken by
 each reconstruction; for 32-bit integers, the maximum number of sons for
 each node is 16.

 lnPanc[3][(tree.nnode-com.ns)*npatt*n] uses the space of com.conP.
 It holds the ln(Pr) for the best reconstructions at the subtree down inode
 given the state of the father node.
 charNode[0,1,2] holds the corresponding state at inode.

 int nBestScore[maxnson];
 int   combIndex[2*n*ncomb];
 double *combScore[n*ncomb];
 char ancSeq[nintern*npatt], ancState1site[nintern];
 int  icharNode[NBESTANC][nintern*npatt*n];
 char  charNode[NBESTANC][nintern*npatt*n];
 */

void UpPassPPSG2000(int inode, int igene, double x[])
{
   /* The algorithm of PPSG2000, modified.  This routine is based on ConditionalPNode().
    lnPanc[h*n+i] is the best lnP, given that inode has state i.
    charNode[] stores the characters that achieved the best lnP.
    */
   int debug = 0;
   int n = com.ncode, it, ibest, i, j, k, h, ison, nson = nodes[inode].nson;
   int pos0 = com.posG[igene], pos1 = com.posG[igene + 1], ichar, jchar;
   int ncomb = 1, icomb, ipath;
   double t, y, psum1site = -1;

   if (com.ncode != 4) debug = 0;

   for (i = 0; i < nson; i++)
      if (nodes[nodes[inode].sons[i]].nson > 0)
         UpPassPPSG2000(nodes[inode].sons[i], igene, x);
   for (i = 0, ncomb = 1; i < nson; i++)
      ncomb *= (nBestScore[i] = (nodes[nodes[inode].sons[i]].nson > 0 ? NBESTANC : 1));
   if (debug) {
      printf("\n\nNode %2d has sons ", inode + 1);
      for (i = 0; i < nson; i++) printf(" %2d", nodes[inode].sons[i] + 1);
      printf("  ncomb=%2d: ", ncomb);
      for (i = 0; i < nson; i++) printf(" %2d", nBestScore[i]);
      printf("\n");
   }

   if (inode != tree.root) {    /* calculate log{P(t)} from father to inode */
      t = nodes[inode].branch * _rateSite;
      if (com.clock < 5) {
         if (com.clock)  t *= GetBranchRate(igene, (int)nodes[inode].label, x, NULL);
         else           t *= com.rgene[igene];
      }
      GetPMatBranch(PMat, x, t, inode);
      for (j = 0; j < n * n; j++)
         PMat[j] = (PMat[j] < 1e-300 ? 300 : -log(PMat[j]));
   }

   for (h = pos0; h < pos1; h++) {  /* loop through site patterns */
      if (h) debug = 0;
      /* The last round for inode==tree.root, shares some code with other nodes,
       and is thus embedded in the same loop.  Alternatively this round can be
       taken out of the loop with some code duplicated.
       */
      for (ichar = 0; ichar < (inode != tree.root ? n : 1); ichar++) { /* ichar for father */
          /* given ichar for the father, what are the best reconstructions at
           inode?  Look at n*ncomb possibilities, given father state ichar.
           */
         if (debug) {
            if (inode == tree.root) printf("\n\nfather is root\n");
            else  printf("\n\nichar = %2d  %c for father\n", ichar + 1, BASEs[ichar]);
         }

         for (icomb = 0; icomb < n * ncomb; icomb++) {
            jchar = icomb / ncomb;      /* jchar is for inode */
            if (inode == tree.root)
               combScore[icomb] = -log(com.pi[jchar] + 1e-300);
            else
               combScore[icomb] = PMat[ichar * n + jchar];

            if (inode == tree.root && parsimony) combScore[icomb] = 0;

            if (debug) printf("comb %2d %c", icomb + 1, BASEs[jchar]);

            for (i = 0, it = icomb % ncomb; i < nson; i++) { /* The ibest-th state in ison. */
               ison = nodes[inode].sons[i];
               ibest = it % nBestScore[i];
               it /= nBestScore[i];

               if (nodes[ison].nson)    /* internal node */
                  y = lnPanc[ibest][(ison - com.ns) * com.npatt * n + h * n + jchar];
               else if (com.cleandata)  /* tip clean: PMatTips[] has log{P(t)}. */
                  y = PMatTips[ison * n * n + jchar * n + com.z[ison][h]];
               else {                   /* tip unclean: PMatTips[] has P(t). */
                  for (k = 0, y = 0; k < nChara[(int)com.z[ison][h]]; k++)
                     y += PMatTips[ison * n * n + jchar * n + CharaMap[(int)com.z[ison][h]][k]];
                  y = -log(y);
               }

               combScore[icomb] += y;
               if (debug) printf("%*s son %2d #%2d %7.1f\n", (i ? 10 : 1), "", ison + 1, ibest + 1, y);
            }
         }  /* for(icomb) */

         if (debug) {
            printf("score ");
            for (i = 0; i < n * ncomb; i++) printf(" %4.1f", combScore[i]);
            printf("\n");
         }
         indexing(combScore, n * ncomb, combIndex, 0, combIndex + n * ncomb);
         if (debug) {
            printf("index ");
            for (i = 0; i < n * ncomb; i++) printf(" %4d", combIndex[i]);
            printf("\n");
         }

         /* print out reconstructions at the site if inode is root. */
         if (inode == tree.root) {
            fprintf(fanc, "%4d ", h + 1);
            if (com.ngene > 1) fprintf(fanc, "(%d) ", igene + 1);
            fprintf(fanc, " %6.0f  ", com.fpatt[h]);
            print1site(fanc, h);
            fprintf(fanc, ": ");
         }
         psum1site = 0;  /* used if inode is root */

         for (j = 0; j < (inode != tree.root ? NBESTANC : n * ncomb); j++) {
            jchar = (it = combIndex[j]) / ncomb; it %= ncomb;
            if (j < NBESTANC) {
               lnPanc[j][(inode - com.ns) * com.npatt * n + h * n + ichar] = combScore[combIndex[j]];
               charNode[j][(inode - com.ns) * com.npatt * n + h * n + ichar] = jchar;
            }
            if (debug) printf("\t#%d: %6.1f %c ", j + 1, combScore[combIndex[j]], BASEs[jchar]);

            for (i = 0, ipath = 0; i < nson; i++) {
               ison = nodes[inode].sons[i];
               ibest = it % nBestScore[i];
               it /= nBestScore[i];
               ipath |= ibest << (2 * i);
               if (debug) printf("%2d", ibest + 1);
            }
            if (j < NBESTANC)
               icharNode[j][(inode - com.ns) * com.npatt * n + h * n + ichar] = ipath;

            if (debug) printf(" (%o)", ipath);

            /* print if inode is root. */
            if (inode == tree.root) {
               ancState1site[inode - com.ns] = jchar;
               if (parsimony) y = combScore[combIndex[j]];
               else          psum1site += y = exp(-combScore[combIndex[j]] - fhsiteAnc[h]);

               for (i = 0; i < nson; i++) {
                  if (nodes[ison = nodes[inode].sons[i]].nson)
                     DownPassPPSG2000OneSite(h, tree.root, jchar, ipath);
               }
               PrintAncState1site(ancState1site, y);
               if (j > NBESTANC && y < .001) break;
            }
         }  /* for(j) */
      }     /* for(ichar) */
      if (inode == tree.root) fprintf(fanc, " (total %6.3f)\n", psum1site);

      if (largeReconstruction && (h + 1) % 2000 == 0)
         printf("\r\tUp pass for gene %d node %d sitepatt %d.", igene + 1, inode + 1, h + 1);

   }        /* for(h) */
   if (largeReconstruction)
      printf("\r\tUp pass for gene %d node %d.", igene + 1, inode + 1);
}


void DownPassPPSG2000OneSite(int h, int inode, int inodestate, int ipath)
{
   /* this puts the state in ancState1site[nintern], using
    int icharNode[NBESTANC][nintern*npatt*n],
    char charNode[NBESTANC][nintern*npatt*n].
    jchar is the state at inode, and ipath is the ipath code for inode.
    */
   int n = com.ncode, i, ison, ibest, sonstate;

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      if (nodes[ison].nson > 1) {
         ibest = (ipath & (3 << (2 * i))) >> (2 * i);

         ancState1site[ison - com.ns] = sonstate =
            charNode[ibest][(ison - com.ns) * com.npatt * n + h * n + inodestate];
         DownPassPPSG2000OneSite(h, ison, sonstate,
            icharNode[ibest][(ison - com.ns) * com.npatt * n + h * n + inodestate]);
      }
   }
}


void PrintAncState1site(char ancState1site[], double prob)
{
   int i;
   char* pch = (com.seqtype == 0 ? BASEs : (com.seqtype == 2 ? AAs : BINs));

   for (i = 0; i < tree.nnode - com.ns; i++) {
      if (com.seqtype == 1) {
#ifdef CODEML
         char codon[4] = "";
         fprintf(fanc, "%s ", getcodon(codon, FROM61[(int)ancState1site[i]]));
#endif
      }
      else
         fprintf(fanc, "%c", pch[(int)ancState1site[i]]);
   }
   fprintf(fanc, " (%5.3f) ", prob);
}

void DownPassPPSG2000(int inode)
{
   /* this reads out the best chara for inode from charNode[] into ancSeq[].
    */
   int i, ison, h;
   char c0 = 0;

   for (h = 0; h < com.npatt; h++) {
      if (inode != tree.root)
         c0 = ancSeq[(nodes[inode].father - com.ns) * com.npatt + h];
      ancSeq[(inode - com.ns) * com.npatt + h]
         = charNode[0][(inode - com.ns) * com.npatt * com.ncode + h * com.ncode + c0];
   }
   for (i = 0; i < nodes[inode].nson; i++)
      if (nodes[ison = nodes[inode].sons[i]].nson > 1)
         DownPassPPSG2000(ison);
}



int AncestralJointPPSG2000(FILE* fout, double x[])
{
   /* Ziheng Yang, 8 June 2000, rewritten on 8 June 2005.
    Joint ancestral reconstruction, taking character states for all nodes at a
    site as one entity, based on the algorithm of Pupko et al. (2000
    Mol. Biol. Evol. 17:890-896).

    fhsiteAns[]: fh[] for each site pattern
    nodes[].conP[] are destroyed and restored at the end of the routine.
    ancSeq[] stores the ancestral seqs, the best reconstruction.

    This outputs results by pattern.  I tried to print results by site (rather than by pattern),
    but gave up as some variables use the same memory (e.g., combIndex) for different site patterns.
    */
   int n = com.ncode, nintern = tree.nnode - com.ns, i, j, k, igene;
   int maxnson, maxncomb;
   double t;
   size_t sconPold = com.sconP, s;

   largeReconstruction = (noisy && (com.ns > 300 || com.ls > 1000000));

   if (noisy) puts("Joint reconstruction.");
   for (i = com.ns, maxnson = 0; i < tree.nnode; i++) {
      for (j = 0, k = 0; j < nodes[i].nson; j++)  /* k is the number of non-tip sons */
         if (nodes[nodes[i].sons[j]].nson > 0) k++;
      if (k > maxnson) maxnson = k;
   }
   if (maxnson > 16 || NBESTANC > 4) /* for int at least 32 bits */
      error2("NBESTANC too large or too many sons.");
   for (i = 0, maxncomb = 1; i < maxnson; i++) maxncomb *= NBESTANC;
   if ((PMatTips = (double*)malloc(com.ns * n * n * sizeof(double))) == NULL)
      error2("oom PMatTips");
   s = NBESTANC * nintern * (size_t)com.npatt * n * sizeof(double);
   if (s > sconPold) {
      com.sconP = s;
      printf("\n%9zu bytes for conP, adjusted\n", com.sconP);
      if ((com.conP = (double*)realloc(com.conP, com.sconP)) == NULL)
         error2("oom conP");
   }
   s = NBESTANC * nintern * com.npatt * n;
   s = ((s * sizeof(int) + (s + nintern) * sizeof(char) + 16) / sizeof(double)) * sizeof(double);
   if (s > com.sspace) {
      com.sspace = s;
      printf("\n%9zu bytes for space, adjusted\n", com.sspace);
      if ((com.space = (double*)realloc(com.space, com.sspace)) == NULL) error2("oom space");
   }
   for (i = 0; i < NBESTANC; i++) {
      lnPanc[i] = com.conP + (size_t)i * nintern * com.npatt * n;
      icharNode[i] = (int*)com.space + i * nintern * com.npatt * n;
      charNode[i] = (char*)((int*)com.space + NBESTANC * nintern * com.npatt * n) + i * nintern * com.npatt * n;
      ancState1site = charNode[0] + NBESTANC * nintern * com.npatt * n;
   }
   if ((ancSeq = (char*)malloc(nintern * com.npatt * n * sizeof(char))) == NULL)
      error2("oom charNode");

   if ((combScore = (double*)malloc((3 * n * maxncomb + com.ns) * sizeof(double))) == NULL)
      error2("oom combScore");
   nBestScore = (int*)(combScore + n * maxncomb);
   combIndex = nBestScore + com.ns;  /* combIndex[2*n*ncomb] contains work space */

   fanc = fout;
   fprintf(fout, "\n\n(2) Joint reconstruction of ancestral sequences\n");
   fprintf(fout, "(eqn. 2 in Yang et al. 1995 Genetics 141:1641-1650), using ");
   fprintf(fout, "the algorithm of Pupko et al. (2000 Mol Biol Evol 17:890-896),\n");
   fprintf(fout, "modified to generate sub-optimal reconstructions.\n");
   fprintf(fout, "\nReconstruction (prob.), listed by pattern (use the observed data to find the right site).\n");
   fprintf(fout, "\nPattern Freq   Data:\n\n");

   for (igene = 0; igene < com.ngene; igene++) {
      if (com.Mgene > 1) SetPGene(igene, 1, 1, 0, x);
      for (i = 0; i < com.ns; i++) {
         t = nodes[i].branch * _rateSite;
         if (com.clock < 5) {
            if (com.clock)  t *= GetBranchRate(igene, (int)nodes[i].label, x, NULL);
            else           t *= com.rgene[igene];
         }
         GetPMatBranch(PMatTips + i * n * n, x, t, i);
      }

      if (com.cleandata) {
         for (i = 0; i < com.ns * n * n; i++)
            PMatTips[i] = (PMatTips[i] < 1e-20 ? 300 : -log(PMatTips[i]));
      }
      if (parsimony)
         for (i = 0; i < com.ns; i++)
            xtoy(P0, PMatTips + i * n * n, n * n);

      UpPassPPSG2000(tree.root, igene, x); /* this prints into frst as well */
   }

   if (largeReconstruction) puts("\n\tDown pass.");
   DownPassPPSG2000(tree.root);

   ListAncestSeq(fout, ancSeq);

   free(ancSeq);
   free(PMatTips);
   free(combScore);
   com.sconP = sconPold;
   if ((com.conP = (double*)realloc(com.conP, com.sconP)) == NULL)
      error2("conP");
   PointconPnodes();
   return (0);
}



int AncestralSeqs(FILE* fout, double x[])
{
   /* Ancestral sequence reconstruction using likelihood (Yang et al. 1995).
    Marginal works with constant rate and variable rates among sites.
    Joint works only with constant rate among sites (ncatG=1).
    */
   int h, k;
   double lnL, * ScaleC = NULL;  /* collected scale factors */

   if (com.Mgene == 1)
      error2("When Mgene=1, use RateAncestor = 0.");
   if (tree.nnode == com.ns) {
      puts("\nNo ancestral nodes to reconstruct..\n");  return(0);
   }
   if (noisy) printf("\nReconstructed ancestral states go into file rst.\n");
   fprintf(fout, "\nAncestral reconstruction by %sML.\n",
      (com.seqtype == 0 ? "BASE" : (com.seqtype == 1 ? "CODON" : "AA")));
   fprintf(fout, "\n");
   OutTreeN(fout, 1, 1);  fprintf(fout, "\n\n");
   OutTreeN(fout, 0, 0);  fprintf(fout, "\n\n");
   OutTreeB(fout);        fprintf(fout, "\n");

   fputs("\ntree with node labels for Rod Page's TreeView\n", fout);
   OutTreeN(fout, 1, PrNodeNum);  fprintf(fout, "\n");

   fprintf(fout, "\nNodes %d to %d are ancestral\n", com.ns + 1, tree.nnode);
   if ((fhsiteAnc = (double*)malloc(com.npatt * sizeof(double))) == NULL)
      error2("oom fhsiteAnc");
   if (com.NnodeScale && com.ncatG > 1)
      if ((ScaleC = (double*)malloc(max2(com.npatt, com.ncatG) * sizeof(double))) == NULL)
         error2("oom ScaleC in AncestralSeqs");

   if (com.alpha)
      puts("Rates are variable among sites, marginal reconstructions only.");
   if (!com.cleandata) fputs("Unreliable at sites with alignment gaps\n", fout);

   if (com.ncatG <= 1 || com.method != 1)
      ProbSitePattern(x, &lnL, fhsiteAnc, ScaleC);

#ifdef BASEML
   if (com.nhomo <= 1)
#endif
      AncestralMarginal(fout, x, fhsiteAnc, ScaleC);

   fflush(fout);
   /* fhsiteAnc[] is modified by both Marginal and Joint. */
   if (com.ncatG <= 1 && tree.nnode > com.ns + 1) {
      ProbSitePattern(x, &lnL, fhsiteAnc, ScaleC);
      for (h = 0; h < com.npatt; h++) {
         fhsiteAnc[h] = log(fhsiteAnc[h]);
         for (k = 0; k < com.NnodeScale; k++)
            fhsiteAnc[h] += com.nodeScaleF[k * com.npatt + h];
      }
      /* AncestralJointPPSG2000 corrupts com.conP[] and fhsiteAnc[].
       */
      AncestralJointPPSG2000(fout, x);
   }
   fprintf(fout, "\n");
   free(fhsiteAnc);
   if (com.NnodeScale && com.ncatG > 1) free(ScaleC);

   return (0);
}


#endif


int SetNodeScale(int inode);
int NodeScale(int inode, int pos0, int pos1);
void InitializeNodeScale(void)
{
   /* This allocates memory to hold scale factors for nodes and also decide on the
    nodes for scaling by calling SetNodeScale().
    The scaling node is chosen before the iteration by counting the number of
    nodes visited in the post-order tree travesal algorithm (see the routine
    SetNodeScale).
    See Yang (2000 JME 51:423-432) for details.
    The memory required is  com.NnodeScale*com.npatt*sizeof(double).
    */
   int i, nS;

   if (com.clock >= 5) return;

   com.NnodeScale = 0;
   com.nodeScale = (char*)realloc(com.nodeScale, tree.nnode * sizeof(char));
   if (com.nodeScale == NULL) error2("oom");
   for (i = 0; i < tree.nnode; i++) com.nodeScale[i] = (char)0;
   SetNodeScale(tree.root);
   nS = com.NnodeScale * com.npatt;
   if (com.conPSiteClass) nS *= com.ncatG;
   if (com.NnodeScale) {
      if ((com.nodeScaleF = (double*)realloc(com.nodeScaleF, nS * sizeof(double))) == NULL)
         error2("oom nscale");
      memset(com.nodeScaleF, 0, nS * sizeof(double));

      if (noisy) {
         printf("\n%d node(s) used for scaling (Yang 2000 J Mol Evol 51:423-432):\n", com.NnodeScale);
         for (i = 0; i < tree.nnode; i++)
            if (com.nodeScale[i]) printf(" %2d", i + 1);
         printf("\n");
      }
   }
}


int SetNodeScale(int inode)
{
   /* This marks nodes for applying scaling factors when calculating f[h].
    */
   int i, ison, d = 0, every;

   if (com.seqtype == 0)       every = 100;   /* baseml */
   else if (com.seqtype == 1)  every = 15;    /* codonml */
   else                     every = 50;    /* aaml */

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      d += (nodes[ison].nson ? SetNodeScale(ison) : 1);
   }
   if (inode != tree.root && d > every) {
      com.nodeScale[inode] = 1;
      d = 1;
      com.NnodeScale++;
   }
   return(d);
}


int NodeScale(int inode, int pos0, int pos1)
{
   /* scale to avoid underflow
    */
   int h, k, j, n = com.ncode;
   double t;

   for (j = 0, k = 0; j < tree.nnode; j++)   /* k-th node for scaling */
      if (j == inode)
         break;
      else if (com.nodeScale[j])
         k++;

   for (h = pos0; h < pos1; h++) {
      for (j = 0, t = 0; j < n; j++)
         if (nodes[inode].conP[h * n + j] > t)
            t = nodes[inode].conP[h * n + j];

      if (t < 1e-300) {
         for (j = 0; j < n; j++)
            nodes[inode].conP[h * n + j] = 1;  /* both 0 and 1 fine */
         com.nodeScaleF[k * com.npatt + h] = -800;  /* this is problematic? */
      }
      else {
         for (j = 0; j < n; j++)
            nodes[inode].conP[h * n + j] /= t;
         com.nodeScaleF[k * com.npatt + h] = log(t);
      }
   }
   return(0);
}



static double* dfsites;

int fx_r(double x[], int np);


#if (BASEML || CODEML)

int HessianSKT2004(double xmle[], double lnLm, double g[], double H[])
{
   /* this calculates the hessian matrix of branch lengths using the approximation
    of Seo et al. (2004), especially useful for approximate likelihood calcualtion
    in divergence time estimation.
    df[0][i*com.npatt+h] has   d log(f_h)/d b_i.
    method = 0 uses difference approximation to first derivatives.
    method = 1 uses analytical calculation of first derivatives (Yang 2000).
    I am under the impression that method = 1 may be useful for very large datasets
    with >10M sites, but I have not implemented this method because the analytical
    calculation of first derivatives is possible for branch lengths only, and not
    available for other parameters.  Right now with method = 0, H and the SEs are
    calculated for all parameters although the H matrix in rst2 is a subset for
    branch lengths only.  More thought about what to do.  Ziheng's note on 8 March 2010.
    */
   int method = 0, backforth, h, i, j, lastround0 = LASTROUND, nzero = 0;
   double* x, * lnL[2], * df[2], eh0 = Small_Diff * 2, eh;

   if (com.np != tree.nbranch && method == 1)
      error2("I think HessianSKT2004 works for branch lengths only");
   df[0] = (double*)malloc((com.npatt * 2 + 1) * com.np * sizeof(double));
   if (df[0] == NULL) error2("oom space in HessianSKT2004");
   df[1] = df[0] + com.npatt * com.np;
   x = df[1] + com.npatt * com.np;
   lnL[0] = (double*)malloc(com.np * 2 * sizeof(double));
   lnL[1] = lnL[0] + com.np;

   LASTROUND = 2;

   for (backforth = 0; backforth < 2; backforth++) {
      for (i = 0; i < com.np; i++) {
         xtoy(xmle, x, com.np);
         eh = eh0 * (fabs(xmle[i]) + 1);
         if (backforth == 0) x[i] = xmle[i] - eh;
         else             x[i] = xmle[i] + eh;
         if (x[i] <= 4e-6)
            nzero++;
         dfsites = df[backforth] + i * com.npatt;
         lnL[backforth][i] = -com.plfun(x, com.np);
      }
   }

   for (i = 0; i < com.np; i++) {
      eh = eh0 * (fabs(xmle[i]) + 1);
      g[i] = (lnL[1][i] - lnL[0][i]) / (eh * 2);
   }
   /*
    printf("\nx gL g H");
    matout(F0, xmle, 1, com.np);
    matout(F0, g, 1, com.np);
    */
   zero(H, com.np * com.np);
   for (i = 0; i < com.np; i++) {
      eh = eh0 * (fabs(xmle[i]) + 1);
      for (h = 0; h < com.npatt; h++)
         df[0][i * com.npatt + h] = (df[1][i * com.npatt + h] - df[0][i * com.npatt + h]) / (eh * 2);
   }

   for (i = 0; i < com.np; i++) {
      for (j = 0; j < com.np; j++)
         for (h = 0; h < com.npatt; h++)
            H[i * com.np + j] -= df[0][i * com.npatt + h] * df[0][j * com.npatt + h] * com.fpatt[h];
   }

   if (nzero) printf("\nWarning: Hessian matrix may be unreliable for zero branch lengths\n");
   LASTROUND = lastround0;
   free(df[0]);
   free(lnL[0]);
   return(0);
}



int lfunRates(FILE* fout, double x[], int np)
{
   /* for dG, AdG or similar non-parametric models
    This distroys com.fhK[], and in return,
    fhK[<npatt] stores rates for conditional mean (re), and
    fhK[<2*npatt] stores the most probable rate category number.
    fhsite[npatt] stores fh=log(fh).
    */
   int ir, il, it, h, hp, j, nscale = 1, direction = -1;
   int lst = (com.readpattern ? com.npatt : com.ls);
   double lnL = 0, fh, fh1, t, re, mre, vre, b1[NCATG], b2[NCATG], * fhsite;

   if (noisy) printf("\nEstimated rates for sites go into file %s\n", ratef);
   if (SetParameters(x)) puts("par err. lfunRates");

   fprintf(fout, "\nEstimated rates for sites from %sML.\n",
      (com.seqtype == 0 ? "BASE" : (com.seqtype == 1 ? "CODON" : "AA")));
   OutTreeN(fout, 1, 1); fprintf(fout, "\n");
   fprintf(fout, "\nFrequencies and rates for categories (K=%d)", com.ncatG);
   fprintf(fout, "\nrate:");  for (j = 0; j < com.ncatG; j++) fprintf(fout, " %8.5f", com.rK[j]);
   fprintf(fout, "\nfreq:");  for (j = 0; j < com.ncatG; j++) fprintf(fout, " %8.5f", com.freqK[j]);
   fprintf(fout, "\n");

   if (com.rho) {
      fprintf(fout, "\nTransition prob matrix over sites");
      matout2(fout, com.MK, com.ncatG, com.ncatG, 8, 4);
   }

   if ((fhsite = (double*)malloc(com.npatt * sizeof(double))) == NULL) error2("oom fhsite");
   fx_r(x, np);
   if (com.NnodeScale) {
      for (h = 0; h < com.npatt; h++) {
         for (ir = 1, it = 0; ir < com.ncatG; ir++)
            if (com.fhK[ir * com.npatt + h] > com.fhK[it * com.npatt + h])
               it = ir;
         t = com.fhK[it * com.npatt + h];
         lnL -= com.fpatt[h] * t;
         for (ir = 0; ir < com.ncatG; ir++)
            com.fhK[ir * com.npatt + h] = exp(com.fhK[ir * com.npatt + h] - t);
      }
   }
   for (h = 0; h < com.npatt; h++) {
      for (ir = 0, fhsite[h] = 0; ir < com.ncatG; ir++)
         fhsite[h] += com.freqK[ir] * com.fhK[ir * com.npatt + h];
   }

   if (com.rho == 0) {     /* dG model */
      if (com.verbose > 1) {
         fprintf(fout, "\nPosterior probabilities for site classes, by %s\n\n",
            (com.readpattern ? "pattern" : "site"));
         for (h = 0; h < lst; h++) {
            fprintf(fout, " %5d  ", h + 1);
            hp = (!com.readpattern ? com.pose[h] : h);
            for (ir = 0; ir < com.ncatG; ir++)
               fprintf(fout, " %9.4f", com.freqK[ir] * com.fhK[ir * com.npatt + hp] / fhsite[hp]);
            fprintf(fout, "\n");
         }
      }

      fprintf(fout, "\n%7s  Freq   Data    Rate (posterior mean & category)\n\n",
         (com.readpattern ? "Pattern" : "Site"));
      for (h = 0, mre = vre = 0; h < com.npatt; h++) {
         for (ir = 0, it = 0, t = re = 0; ir < com.ncatG; ir++) {
            fh1 = com.freqK[ir] * com.fhK[ir * com.npatt + h];
            if (fh1 > t) { t = fh1; it = ir; }
            re += fh1 * com.rK[ir];
         }
         lnL -= com.fpatt[h] * log(fhsite[h]);

         re /= fhsite[h];
         mre += com.fpatt[h] * re / com.ls;
         vre += com.fpatt[h] * re * re / com.ls;
         com.fhK[h] = re;
         com.fhK[com.npatt + h] = it + 1.;
      }
      vre -= mre * mre;
      for (h = 0; h < lst; h++) {
         hp = (!com.readpattern ? com.pose[h] : h);
         fprintf(fout, "%7d %5.0f  ", h + 1, com.fpatt[hp]);
         print1site(fout, hp);
         fprintf(fout, " %8.3f%6.0f\n", com.fhK[hp], com.fhK[com.npatt + hp]);
      }
   }
   else {      /* Auto-dGamma model */
      fputs("\nSite Freq  Data  Rates\n\n", fout);
      h = (direction == 1 ? com.ls - 1 : 0);
      for (il = 0, mre = vre = 0; il < lst; h -= direction, il++) {
         hp = (!com.readpattern ? com.pose[h] : h);
         if (il == 0)
            for (ir = 0; ir < com.ncatG; ir++)
               b1[ir] = com.fhK[ir * com.npatt + hp];
         else {
            for (ir = 0; ir < com.ncatG; ir++) {
               for (j = 0, fh = 0; j < com.ncatG; j++)
                  fh += com.MK[ir * com.ncatG + j] * b1[j];
               b2[ir] = fh * com.fhK[ir * com.npatt + hp];
            }
            xtoy(b2, b1, com.ncatG);
         }
         if ((il + 1) % nscale == 0)
         {
            fh = sum(b1, com.ncatG); abyx(1 / fh, b1, com.ncatG); lnL -= log(fh);
         }

         for (ir = 0, it = -1, re = fh1 = t = 0; ir < com.ncatG; ir++) {
            re += com.freqK[ir] * b1[ir] * com.rK[ir];
            fh1 += com.freqK[ir] * b1[ir];
            if (b1[ir] > t) { it = ir; t = b1[ir]; }
         }
         re /= fh1;
         mre += re / com.ls;
         vre += re * re / com.ls;

         fprintf(fout, "%4d %5.0f  ", h + 1, com.fpatt[hp]);
         print1site(fout, hp);
         fprintf(fout, " %8.3f%6.0f\n", re, it + 1.);
      }  /* for(il) */
      vre -= mre * mre;
      for (ir = 0, fh = 0; ir < com.ncatG; ir++)  fh += com.freqK[ir] * b1[ir];
      lnL -= log(fh);
   }
   if (noisy) printf("lnL =%14.6f\n", -lnL);
   fprintf(fout, "\nlnL =%14.6f\n", -lnL);
   if (com.ngene == 1) {
      fprintf(fout, "\nmean(r^)=%9.4f  var(r^)=%9.4f", mre, vre);
      fprintf(fout, "\nAccuracy of rate prediction: corr(r^,r) =%9.4f\n",
         sqrt(com.alpha * vre));
   }
   free(fhsite);
   return (0);
}


double lfunAdG(double x[], int np)
{
   /* Auto-Discrete-Gamma rates for sites
    See notes in lfundG().
    */
   int  nscale = 1, h, il, ir, j, FPE = 0;
   int  direction = -1;  /* 1: n->1;  -1: 1->n */
   double lnL = 0, b1[NCATG], b2[NCATG], fh;

   NFunCall++;
   fx_r(x, np);
   if (com.NnodeScale)
      for (h = 0; h < com.npatt; h++) {
         fh = com.fhK[0 * com.npatt + h];
         lnL -= fh * com.fpatt[h];
         for (ir = 1, com.fhK[h] = 1; ir < com.ncatG; ir++)
            com.fhK[ir * com.npatt + h] = exp(com.fhK[ir * com.npatt + h] - fh);
      }
   h = (direction == 1 ? com.ls - 1 : 0);
   for (il = 0; il < com.ls; h -= direction, il++) {
      if (il == 0)
         for (ir = 0; ir < com.ncatG; ir++)
            b1[ir] = com.fhK[ir * com.npatt + com.pose[h]];
      else {
         for (ir = 0; ir < com.ncatG; ir++) {
            for (j = 0, fh = 0; j < com.ncatG; j++)
               fh += com.MK[ir * com.ncatG + j] * b1[j];
            b2[ir] = fh * com.fhK[ir * com.npatt + com.pose[h]];
         }
         xtoy(b2, b1, com.ncatG);
      }
      if ((il + 1) % nscale == 0) {
         fh = sum(b1, com.ncatG);
         if (fh < 1e-90) {
            if (!FPE) {
               FPE = 1; printf("h,fh%6d %12.4e\n", h + 1, fh);
               print1site(F0, h);
               printf("\n");
            }
            fh = 1e-300;
         }
         abyx(1 / fh, b1, com.ncatG); lnL -= log(fh);
      }
   }
   for (ir = 0, fh = 0; ir < com.ncatG; ir++)  fh += com.freqK[ir] * b1[ir];
   lnL -= log(fh);
   return (lnL);
}

#endif




#if (defined(BASEML))

int GetPMatBranch(double Pt[], double x[], double t, int inode)
{
   /* P(t) for branch leading to inode, called by routines ConditionalPNode()
    and AncestralSeq() in baseml and codeml.  x[] is not used by baseml.
    */
   int n = com.ncode, i;
   double space[NCODE * NCODE * 3] = { 0 };

   if (com.model <= K80)
      PMatK80(Pt, t, (com.nhomo == 2 ? *nodes[inode].pkappa : com.kappa));
   else {
      if (com.nhomo == 2)
         eigenTN93(com.model, *nodes[inode].pkappa, -1, com.pi, &nR, Root, Cijk);
      else if (com.nhomo > 2 && com.model <= TN93)
         eigenTN93(com.model, *nodes[inode].pkappa, *(nodes[inode].pkappa + 1), nodes[inode].pi, &nR, Root, Cijk);
      else if (com.nhomo > 2 && com.model == REV)
         eigenQREVbase(NULL, Pt, nodes[inode].pkappa, nodes[inode].pi, &nR, Root, Cijk);

      if (com.model <= REV || com.model == REVu)
         PMatCijk(Pt, t);
      else {
         QUNREST(NULL, Pt, x + com.ntime + com.nrgene, com.pi);
         for (i = 0; i < n * n; i++) Pt[i] *= t;
         matexp(Pt, n, 7, 5, space);
      }
   }
   return(0);
}

#elif (defined(CODEML))

int GetPMatBranch(double Pt[], double x[], double t, int inode)
{
   /* P(t) for branch leading to inode, called by routines ConditionalPNode()
    and AncestralSeq() in baseml and codeml.

    Qfactor in branch & site models (model = 2 or 3 and NSsites = 2 or 3):
    Qfactor scaling is applied here and not inside eigenQcodon().
    */
   int iUVR = 0, nUVR = NBTYPE + 2, ib = (int)nodes[inode].label;
   double* pkappa, w, mr = 1, Qfactor = 1;

   pkappa = (com.hkyREV || com.codonf == FMutSel ? x + com.ntime + com.nrgene : &com.kappa);

   if (com.seqtype == CODONseq && com.NSsites && com.model >= NSbranchB && com.model <= NSbranch3) {
      /* branch&site models (both NSsites & model):
       Usual likelihood calculation, no need to re-calculate UVRoot.
       Only need to point to the right place.
       */
      iUVR = Set_UVR_BranchSite(IClass, ib);
      Qfactor = Qfactor_NS_branch[ib];
   }
   else if (com.seqtype == CODONseq && BayesEB == 2 && com.model >= NSbranch2 && com.model <= NSbranch3) { /* BEB for A&C */
       /* branch&site models (both NSsites & model) BEB calculation:
        Need to calculate UVRoot, as w is different.  com.pomega points to wbranches[]
        in get_grid_para_like_M2M8() or get_grid_para_like_AC().

        Qfactor_NS_branch[] is fixed at the MLE:
        "we fix the branch lengths at the synonymous sites (i.e., the expected
        number of synonymous substitutions per codon) at their MLEs."
        */
      w = com.pomega[ib];
      eigenQcodon(1, -1, NULL, NULL, NULL, Root, U, V, &mr, pkappa, w, Pt);
      Qfactor = Qfactor_NS_branch[ib];
   }
   else if (com.seqtype == CODONseq && com.model >= NSbranchB && com.model <= NSbranch2 && com.nbtype <= nUVR) {
      /* branch model, also for AAClasses */
      iUVR = (int)nodes[inode].label;
      U = _UU[iUVR]; V = _VV[iUVR]; Root = _Root[iUVR];
   }
   else if (com.seqtype == CODONseq && com.model >= NSbranchB && com.model <= NSbranch3) {
      mr = 0;
      if (com.aaDist == AAClasses) { /* AAClass model */
         com.pomega = PointOmega(x + com.ntime, -1, inode, -1);
         eigenQcodon(1, -1, NULL, NULL, NULL, Root, U, V, &mr, pkappa, -1, Pt);
      }
      else if (com.nbtype > nUVR) {  /* branch models, with more than 8 omega */
         eigenQcodon(1, -1, NULL, NULL, NULL, Root, U, V, &mr, pkappa, nodes[inode].omega, Pt);
      }
   }

   if (com.seqtype == AAseq && com.model == Poisson)
      PMatJC69like(Pt, t, com.ncode);
   else {
      t *= Qfactor;
      PMatUVRoot(Pt, t, com.ncode, U, V, Root);
   }

   return(0);
}

#endif



void print_lnf_site(int h, double logfh)
{
#if(defined BASEML || defined CODEML)
   fprintf(flnf, "\n%6d %6.0f %16.10f %16.12f %12.4f  ",
      h + 1, com.fpatt[h], logfh, exp(logfh), com.ls * exp(logfh));
   print1site(flnf, h);

#endif
}

double lfundG(double x[], int np)
{
   /* likelihood function for site-class models.
      This deals with scaling for nodes to avoid underflow if(com.NnodeScale).
      The routine calls fx_r() to calculate com.fhK[], which holds log{f(x|r)}
      when scaling or f(x|r) when not.  Scaling factors are set and used for each
      site class (ir) to calculate log(f(x|r).  When scaling is used, the routine
      converts com.fhK[] into f(x|r), by collecting scaling factors into lnL.
      The rest of the calculation then becomes the same and relies on f(x|r).
      Check notes in fx_r.
      This is also used for NSsites models in codonml.
      Note that scaling is used between fx_r() and ConditionalPNode()
      When this routine is used under the multiple-gene site-class model, note
      that right now it assumes one set of com.freqK[] for the different genes,
      which may be an issue.
   */
   int h, ir, it, FPE = 0;
   double lnL = 0, fh = 0, t;

   NFunCall++;
   fx_r(x, np);

   for (h = 0; h < com.npatt; h++) {
      if (com.fpatt[h] <= 0 && com.print >= 0) continue;
      if (com.NnodeScale) { /* com.fhK[] has log{f(x|r}.  Note the scaling for nodes */
         for (ir = 1, it = 0; ir < com.ncatG; ir++) /* select term for scaling */
            if (com.fhK[ir * com.npatt + h] > com.fhK[it * com.npatt + h]) it = ir;
         t = com.fhK[it * com.npatt + h];
         for (ir = 0, fh = 0; ir < com.ncatG; ir++)
            fh += com.freqK[ir] * exp(com.fhK[ir * com.npatt + h] - t);
         fh = t + log(fh);
      }
      else {
         for (ir = 0, fh = 0; ir < com.ncatG; ir++)
            fh += com.freqK[ir] * com.fhK[ir * com.npatt + h];
         if (fh <= 0) {
            if (!FPE) {
               FPE = 1;  matout(F0, x, 1, np);
               printf("\nlfundG: h=%4d  fhK=%9.6e\ndata: ", h + 1, fh);
               print1site(F0, h);
               printf("\n");
            }
            fh = 1e-300;
         }
         fh = log(fh);
      }
      lnL -= fh * com.fpatt[h];
      if (LASTROUND == 2) dfsites[h] = fh;
      if (com.print < 0) print_lnf_site(h, fh);
   }

   return(lnL);
}


int SetPSiteClass(int iclass, double x[])
{
   /* This sets parameters for the iclass-th site class
    This is used by ConditionalPNode() and also updateconP in both algorithms
    For method=0 and 1.
    */
   _rateSite = com.rK[iclass];
#if CODEML
   double* pkappa = NULL, * xcom = x + com.ntime, mr;
   int k = com.nrgene + !com.fix_kappa;

   IClass = iclass;
   mr = 1 / Qfactor_NS;
   pkappa = (com.hkyREV || com.codonf == FMutSel ? xcom + com.nrgene : &com.kappa);
   if (com.seqtype == CODONseq && com.NSsites) {
      _rateSite = 1;
      if (com.model == 0) {
         if (com.aaDist) {
            if (com.aaDist < 10)       com.pomega = xcom + k + com.ncatG - 1 + 2 * iclass;
            else if (com.aaDist == 11) com.pomega = xcom + k + com.ncatG - 1 + 4 * iclass;
            else if (com.aaDist == 12) com.pomega = xcom + k + com.ncatG - 1 + 5 * iclass;
         }
         eigenQcodon(1, -1, NULL, NULL, NULL, Root, U, V, &mr, pkappa, com.rK[iclass], PMat);
      }
   }
#endif
   return (0);
}


extern int prt, Locus, Ir;


int fx_r(double x[], int np)
{
   /* This calculates f(x|r) if(com.NnodeScale==0) or log{f(x|r)}
    if(com.NnodeScale>0), that is, the (log) probability of observing data x
    at a site, given the rate r or dN/dS ratio for the site.  This is used by
    the discrete-gamma models in baseml and codeml as well as the NSsites models
    in codeml.
    The results are stored in com.fhK[com.ncatG*com.npatt].
    This deals with underflows with large trees using global variables
    com.nodeScale and com.nodeScaleF[com.NnodeScale*com.npatt].
    */
   int  h, ir, i, k, ig;
   double fh; /* for testing site class with w=0 */

   if (!BayesEB)
      if (SetParameters(x)) puts("\npar err..");

   for (ig = 0; ig < com.ngene; ig++) { /* alpha may differ over ig */
      if (com.Mgene > 1 || com.nalpha > 1)
         SetPGene(ig, com.Mgene > 1, com.Mgene > 1, com.nalpha > 1, x);
      for (ir = 0; ir < com.ncatG; ir++) {
         if (ir && com.conPSiteClass) {  /* shift com.nodeScaleF & conP */
            if (com.NnodeScale)
               com.nodeScaleF += (size_t)com.npatt * com.NnodeScale;
            for (i = com.ns; i < tree.nnode; i++)
               nodes[i].conP += (size_t)(tree.nnode - com.ns) * com.ncode * com.npatt;
         }
         SetPSiteClass(ir, x);
         ConditionalPNode(tree.root, ig, x);

         for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
            if (com.fpatt[h] <= 0 && com.print >= 0) continue;
            for (i = 0, fh = 0; i < com.ncode; i++)
               fh += com.pi[i] * nodes[tree.root].conP[h * com.ncode + i];

            if (fh <= 0) {
               if (fh < -1e-10) {
                  matout(F0, x, 1, np);
                  printf("\nfx_r: h = %d  r = %d fhK = %.5e ", h + 1, ir + 1, fh);
                  if (com.seqtype == 0 || com.seqtype == 2) {
                     printf("Data: ");
                     print1site(F0, h);
                     printf("\n");
                  }
               }
               fh = 1e-300;
            }
            if (!com.NnodeScale)
               com.fhK[ir * com.npatt + h] = fh;
            else {
               com.fhK[ir * com.npatt + h] = log(fh);
               for (k = 0; k < com.NnodeScale; k++)
                  com.fhK[ir * com.npatt + h] += com.nodeScaleF[k * com.npatt + h];
            }
         }  /* for (h) */
      }     /* for (ir) */

      if (com.conPSiteClass) {  /* shift pointers conP back */
         if (com.NnodeScale)
            com.nodeScaleF -= (com.ncatG - 1) * com.NnodeScale * (size_t)com.npatt;
         for (i = com.ns; i < tree.nnode; i++)
            nodes[i].conP -= (size_t)(com.ncatG - 1) * (tree.nnode - com.ns) * com.ncode * com.npatt;
      }
   }  /* for(ig) */
   return(0);
}


double lfun(double x[], int np)
{
   /* likelihood function for models of one rate for all sites including Mgene models.
   */
   int  h, i, k, ig, FPE = 0;
   double lnL = 0, fh;

   NFunCall++;
   if (SetParameters(x)) puts("\npar err..");
   for (ig = 0; ig < com.ngene; ig++) {
      if (com.Mgene > 1)
         SetPGene(ig, 1, 1, 0, x);
      ConditionalPNode(tree.root, ig, x);

      for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
         if (com.fpatt[h] <= 0 && com.print >= 0) continue;
         for (i = 0, fh = 0; i < com.ncode; i++)
            fh += com.pi[i] * nodes[tree.root].conP[h * com.ncode + i];
         if (fh <= 0) {
            if (fh < -1e-5 && noisy) {
               printf("\nfh = %.6f negative\n", fh);
               exit(-1);
            }
            if (!FPE) {
               FPE = 1;  matout(F0, x, 1, np);
               printf("lfun: h=%4d  fh=%9.6e\nData: ", h + 1, fh);
               print1site(F0, h);
               printf("\n");
               matout(F0, x, 1, com.np);
            }
            fh = 1e-80;
         }
         fh = log(fh);
         for (k = 0; k < com.NnodeScale; k++)
            fh += com.nodeScaleF[k * com.npatt + h];

         lnL -= fh * com.fpatt[h];
         if (LASTROUND == 2) dfsites[h] = fh;
         if (com.print < 0)
            print_lnf_site(h, fh);
      }
   }
   return (lnL);
}




#if(defined MINIMIZATION)

/* November, 1999, Minimization branch by branch */
int noisy_minbranches;
double* space_minbranches, * g_minbranches, * varb_minbranches, e_minbranches;

double minbranches(double xcom[], int np);
int lfunt(double t, int a, int b, double x[], double* l, double space[]);
int lfuntdd(double t, int a, int b, double x[], double* l, double* dl, double* ddl,
   double space[]);
int lfunt_SiteClass(double t, int a, int b, double x[], double* l, double space[]);
int lfuntdd_SiteClass(double t, int a, int b, double x[],
   double* l, double* dl, double* ddl, double space[]);

int minB(FILE* fout, double* lnL, double x[], double xb[][2], double e0, double space[])
{
   /* This calculates lnL for given values of common parameters by optimizing
    branch lengths, cycling through them.
    Z. Yang, November 1999
    This calls minbranches to optimize branch lengths and ming2 to
    estimate other paramters.
    At the end of the routine, there is a call to lfun to restore nodes[].conP.
    Returns variances of branch lengths in space[].
    space[] is com.space[].  com.space may be reallocated here, which may be unsafe
    as the pointers in the calling routine may not be pointing to the right places.

    return value: 0 convergent;  -1: not convergent.
    */
   int ntime0 = com.ntime, fix_blength0 = com.fix_blength;
   int status = 0, i, npcom = com.np - com.ntime;
   size_t s;
   double* xcom = x + com.ntime, lnL0 = *lnL, dl, e = 1e-5;
   double(*xbcom)[2] = xb + ntime0;
   int small_times = 0, max_small_times = 100, ir, maxr = (npcom ? 200 : 1);
   double small_improvement = 0.001;
   char timestr[64];

   if (com.conPSiteClass) {
      s = (2 * com.ncode * com.ncode + com.ncode * (size_t)com.npatt) * sizeof(double);
      if (com.sspace < s) {  /* this assumes that space is com.space */
         printf("\n%zu bytes in space, %zu bytes needed\n", com.sspace, s);
         printf("minB: reallocating memory for working space.\n");
         com.space = (double*)realloc(com.space, s);
         if (com.space == NULL) error2("oom space");
         com.sspace = s;
      }
   }
   g_minbranches = com.space;
   varb_minbranches = com.space + com.np;
   s = (3 * com.ncode * com.ncode + (com.conPSiteClass) * 4 * (size_t)com.npatt) * sizeof(double);
   if ((space_minbranches = (double*)malloc(s)) == NULL)
      error2("oom minB");
   if (com.ntime == 0) error2("minB: should not come here");

   if (*lnL <= 0)  *lnL = com.plfun(x, com.np);
   e = e_minbranches = (npcom ? 5.0 : e0);
   com.ntime = 0; com.fix_blength = 2;
#if(CODEML)
   if (com.NSsites == 0) com.pomega = xcom + com.nrgene + !com.fix_kappa;
#endif

   for (ir = 0; (npcom == 0 || com.method) && ir < maxr; ir++) {
      if (npcom) {
         if (noisy > 2) printf("\n\nRound %da: Paras (%d) (e=%.6g)", ir + 1, npcom, e);
         ming2(NULL, lnL, com.plfun, NULL, xcom, xbcom, com.space, e, npcom);
         if (noisy > 2) {
            printf("\n");
            for (i = 0; i < npcom; i++) printf(" %11.6f", xcom[i]);
            printf("%8s%s\n", "", printtime(timestr));
         }
      }
      noisy_minbranches = noisy;
      if (noisy > 2)
         printf("\nRound %db: Blengths (%d, e=%.6g)\n", ir + 1, tree.nbranch, e_minbranches);

      *lnL = minbranches(xcom, -1);
      for (i = 0; i < tree.nnode; i++)
         if (i != tree.root)
            x[nodes[i].ibranch] = nodes[i].branch;
      if (noisy > 2) printf("\n%s\n", printtime(timestr));

      if ((dl = fabs(*lnL - lnL0)) < e0 && e <= 0.02) break;
      if (dl < small_improvement) small_times++;
      else                        small_times = 0;
      if ((small_times > max_small_times && ntime0 < 200) || (com.method == 2 && ir == 1)) {
         if (noisy && com.method != 2) puts("\nToo slow, switching algorithm.");
         status = 2;
         break;
      }
      if (noisy && small_times > 5)
         printf("\n%d rounds of small improvement.", small_times);

      e /= 2;  if (dl < 1) e /= 2;
      if (dl < 0.5)     e = min2(e, 1e-3);
      else if (dl > 10) e = max2(e, 0.1);
      e_minbranches = max2(e, 1e-6);
      e = max2(e, 1e-6);

      lnL0 = *lnL;
      if (fout) {
         fprintf(fout, "%4d %12.5f x ", ir + 1, *lnL);
         for (i = 0; i < com.np; i++)
            fprintf(fout, " %8.5f", x[i]);
         fprintf(fout, "\n");
         fflush(fout);
      }
   }
   if (npcom && ir == maxr) status = -1;

   if (npcom && status == 2) {
      noisy_minbranches = 0;
      com.ntime = ntime0;
      com.fix_blength = fix_blength0;
      ming2(NULL, lnL, com.plfun, NULL, x, xb, com.space, e0, com.np);
      for (i = 0; i < tree.nnode; i++) space[i] = -1;
   }

   for (i = 0; i < tree.nnode; i++)
      if (i != tree.root) x[nodes[i].ibranch] = nodes[i].branch;

   if (noisy > 2) printf("\nlnL  = %12.6f\n", -*lnL);

   com.ntime = ntime0;
   com.fix_blength = fix_blength0;
   *lnL = com.plfun(x, com.np); /* restore things, for e.g. AncestralSeqs */
   if (fabs(*lnL - lnL0) > 1e-5)
      printf("%.6f != %.6f lnL error.  Something is wrong in minB\n", *lnL, lnL0);
   free(space_minbranches);

   return (status == -1 ? -1 : 0);
}


/*********************  START: Testing iteration algorithm ******************/

int minB2(FILE* fout, double* lnL, double x[], double xb[][2], double e0, double space[])
{
   /*
    */
   int ntime0 = com.ntime, fix_blength0 = com.fix_blength;
   int status = 0, i, npcom = com.np - com.ntime;
   size_t s;
   double* xcom = x + com.ntime;
   double(*xbcom)[2] = xb + ntime0;

   s = (3 * com.ncode * com.ncode + (com.conPSiteClass) * 4 * (size_t)com.npatt) * sizeof(double);
   if ((space_minbranches = (double*)malloc(s)) == NULL)  error2("oom minB2");
   if (com.ntime == 0 || npcom == 0) error2("minB2: should not come here");

   noisy_minbranches = 0;
   /* if(*lnL<=0)  *lnL=com.plfun(x,com.np); */
   com.ntime = 0; com.fix_blength = 2;
#if(CODEML)
   if (com.NSsites == 0) com.pomega = xcom + com.nrgene + !com.fix_kappa;
#endif

   ming2(NULL, lnL, minbranches, NULL, xcom, xbcom, space, e0, npcom);

   com.ntime = ntime0;  com.fix_blength = fix_blength0;
   for (i = 0; i < tree.nnode; i++)
      if (i != tree.root) x[nodes[i].ibranch] = nodes[i].branch;
   *lnL = com.plfun(x, com.np); /* restore things, for e.g. AncestralSeqs */
   free(space_minbranches);

   return (status == -1 ? -1 : 0);
}

/*********************  END: Testing iteration algorithm ******************/


int updateconP(double x[], int inode)
{
   /* update conP for inode.

      Confusing decision about x[] follows.  Think about redesign.

      (1) Called by PostProbNode for ancestral reconstruction, with com.clock = 0,
      1, 2: x[] is passed over and com.ntime is used to get xcom in
      SetPSiteClass()
      (2) Called from minbranches(), with com.clock = 0.  xcom[] is passed
      over by minbranches and com.ntime=0 is set.  So SetPSiteClass()
      can still get the correct substitution parameters.
      Also look at ConditionalPNode().

      Note that com.nodeScaleF and nodes[].conP are shifted if(com.conPSiteClass).
   */
   int ig, i, ir;

   if (com.conPSiteClass == 0)
      for (ig = 0; ig < com.ngene; ig++) {
         if (com.Mgene > 1 || com.nalpha > 1)
            SetPGene(ig, com.Mgene > 1, com.Mgene > 1, com.nalpha > 1, x);
         /* x[] needed by local clock models and if(com.aaDist==AAClasses).
          This is called from PostProbNode
          */

         ConditionalPNode(inode, ig, x);
      }
   else {  /* site-class models */
      for (ir = 0; ir < com.ncatG; ir++) {
#ifdef CODEML
         IClass = ir;
#endif
         if (ir) {
            if (com.NnodeScale)
               com.nodeScaleF += com.NnodeScale * (size_t)com.npatt;
            for (i = com.ns; i < tree.nnode; i++)
               nodes[i].conP += (tree.nnode - com.ns) * com.ncode * (size_t)com.npatt;
         }
         SetPSiteClass(ir, x);
         for (ig = 0; ig < com.ngene; ig++) {
            if (com.Mgene > 1 || com.nalpha > 1)
               SetPGene(ig, com.Mgene > 1, com.Mgene > 1, com.nalpha > 1, x);
            if (com.nalpha > 1) SetPSiteClass(ir, x);
            ConditionalPNode(inode, ig, x);
         }
      }

      /* shift positions */
      com.nodeScaleF -= (com.ncatG - 1) * com.NnodeScale * com.npatt;
      for (i = com.ns; i < tree.nnode; i++)
         nodes[i].conP -= (com.ncatG - 1) * (tree.nnode - com.ns) * com.ncode * (size_t)com.npatt;
   }
   return(0);
}


double minbranches(double x[], int np)
{
   /* Ziheng, November 1999.
    optimizing one branch at a time

    for each branch a..b, reroot the tree at b, and
    then calculate conditional probability for node a.
    For each branch, this routine determines the Newton search direction
    p = -dl/dll.  It then halves the steplength to make sure -lnL is decreased.
    When the Newton solution is correct, this strategy will waste one
    extra call to lfunt.  It does not seem possible to remove calculation of
    l (lnL) in lfuntddl().
    lfun or lfundG and thus SetParameters are called once beforehand to set up
    globals like com.pomega.
    This works with NSsites and NSbranch models.

    com.oldconP[] marks nodes that need to be updated when the tree is rerooted.
    The array is declared in baseml and codeml and used in the following
    routines: ReRootTree, minbranches, and ConditionalPNode.

    Note: At the end of the routine, nodes[].conP are not updated.
    */
   int ib, oldroot = tree.root, a, b;
   int icycle, maxcycle = 500, icycleb, ncycleb = 10, i;
   double lnL, lnL0 = 0, l0, l, dl, ddl = -1, t, t0, t00, p, step = 1, smallv = 1e-20, y;
   double tb[2] = { 1e-8,50 }, e = e_minbranches, * space = space_minbranches;
   double* xcom = x + com.ntime;  /* this is incorrect as com.ntime=0 */

   if (com.ntime) error2("ntime should be 0 in minbranches");
   lnL0 = l0 = l = lnL = com.plfun(xcom, -1);

   if (noisy_minbranches > 2) printf("\tlnL0 =    %14.6f\n", -lnL0);

   for (icycle = 0; icycle < maxcycle; icycle++) {
      for (ib = 0; ib < tree.nbranch; ib++) {
         t = t0 = t00 = nodes[tree.branches[ib][1]].branch;
         l0 = l;
         a = tree.branches[ib][0];
         b = tree.branches[ib][1];
         /* if a is the root, why do we want to reroot the tree at b?  Just switch a with b? */

         for (i = 0; i < tree.nnode; i++)
            com.oldconP[i] = 1;
         ReRootTree(b);
         updateconP(x, a);

         for (icycleb = 0; icycleb < ncycleb; icycleb++) {  /* iterating a branch */
            if (!com.conPSiteClass)
               lfuntdd(t, a, b, xcom, &y, &dl, &ddl, space);
            else
               lfuntdd_SiteClass(t, a, b, xcom, &y, &dl, &ddl, space);

            p = -dl / fabs(ddl);
            /* p = -dl/ddl; newton direction */
            if (fabs(p) < smallv) step = 0;
            else if (p < 0)       step = min2(1, (tb[0] - t0) / p);
            else                  step = min2(1, (tb[1] - t0) / p);

            if (icycle == 0 && step != 1 && step != 0)
               step *= 0.99; /* avoid border */
            for (i = 0; step > smallv; i++, step /= 4) {
               t = t0 + step * p;
               if (!com.conPSiteClass) lfunt(t, a, b, xcom, &l, space);
               else                   lfunt_SiteClass(t, a, b, xcom, &l, space);
               if (l < l0) break;
            }
            if (step <= smallv) { t = t0; l = l0; break; }
            if (fabs(t - t0) < e * fabs(1 + t) && fabs(l - l0) < e) break;
            t0 = t; l0 = l;
         }
         nodes[a].branch = t;

         g_minbranches[ib] = -dl;
         varb_minbranches[ib] = -ddl;
      }   /* for (ib) */
      lnL = l;
      if (noisy_minbranches > 2) printf("\tCycle %2d: %14.6f\n", icycle + 1, -l);
      if (fabs(lnL - lnL0) < e) break;
      lnL0 = lnL;
   }  /* for (icycle) */
   ReRootTree(oldroot);  /* did not update conP */
   for (i = 0; i < tree.nnode; i++)
      com.oldconP[i] = 0;
   return(lnL);
}



int lfunt(double t, int a, int b, double xcom[], double* l, double space[])
{
   /* See notes for lfunt_dd and minbranches
    */
   int i, j, k, h, ig, n = com.ncode, nroot = n;
   int n1 = (com.cleandata && b < com.ns ? 1 : n), xb;
   double expt, multiply;
   double* P = space, piqi, pqj, fh;

#if (CODEML)
   int nUVR = NBTYPE + 2;
   double* pkappa, mr = 0;

   pkappa = (com.hkyREV || com.codonf == FMutSel ? xcom + com.nrgene : &com.kappa);
   if (com.seqtype == CODONseq && com.model) {
      if ((com.model == NSbranchB || com.model == NSbranch2) && com.NSsites == 0 && com.nbtype <= nUVR) {
         U = _UU[(int)nodes[a].label];
         V = _VV[(int)nodes[a].label];
         Root = _Root[(int)nodes[a].label];
      }
      else {
         eigenQcodon(1, -1, NULL, NULL, NULL, Root, U, V, &mr, pkappa, nodes[a].omega, PMat);
      }
   }
#endif

#if (BASEML)
   if (com.nhomo == 2)
      eigenTN93(com.model, *nodes[a].pkappa, 1, com.pi, &nR, Root, Cijk);
   nroot = nR;
#endif

   * l = 0;
   for (ig = 0; ig < com.ngene; ig++) {
      if (com.Mgene > 1) SetPGene(ig, 1, 1, 0, xcom); /* com.ntime=0 */
      for (i = 0; i < n * n; i++) P[i] = 0;

      for (k = 0, expt = 1; k < nroot; k++) {
         multiply = com.rgene[ig] * Root[k];
         if (k) expt = exp(t * multiply);

#if (CODEML)  /* uses U & V */
         double uexpt = 0;
         for (i = 0; i < n; i++)
            for (j = 0, uexpt = U[i * n + k] * expt; j < n; j++)
               P[i * n + j] += uexpt * V[k * n + j];
#elif (BASEML) /* uses Cijk */
         for (i = 0; i < n; i++) for (j = 0; j < n; j++)
            P[i * n + j] += Cijk[i * n * nroot + j * nroot + k] * expt;
#endif
      }

      for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
         n1 = (b < com.ns ? nChara[(int)com.z[b][h]] : n);
         for (i = 0, fh = 0; i < n1; i++) {
            xb = i;
            if (b < com.ns)  piqi = com.pi[xb = CharaMap[(int)com.z[b][h]][i]];
            else             piqi = com.pi[i] * nodes[b].conP[h * n + i];

            for (j = 0, pqj = 0; j < n; j++)
               pqj += P[xb * n + j] * nodes[a].conP[h * n + j];
            fh += piqi * pqj;
         }
         if (noisy && fh < 1e-250)
            printf("a bit too small: fh[%d] = %10.6e\n", h, fh);
         if (fh < 0) fh = -500;
         else     fh = log(fh);

         *l -= fh * com.fpatt[h];
         for (i = 0; i < com.NnodeScale; i++)
            *l -= com.nodeScaleF[i * com.npatt + h] * com.fpatt[h];
      }
   }
   return(0);
}


int lfuntdd(double t, int a, int b, double xcom[], double* l, double* dl, double* ddl, double space[])
{
   /* Calculates lnL for branch length t for branch b->a.
    See notes in minbranches().
    Conditional probability updated correctly already.

    i for b, j for a?
    */
   int i, j, k, h, ig, n = com.ncode, nroot = n;
   int n1 = (com.cleandata && b < com.ns ? 1 : n), xb;
   double expt, multiply;
   double* P = space, * dP = P + n * n, * ddP = dP + n * n, piqi, pqj, dpqj, ddpqj, fh, dfh, ddfh;

#if(CODEML)
   int nUVR = NBTYPE + 2;
   double* pkappa, mr = 0, uexpt = 0;

   pkappa = (com.hkyREV || com.codonf == FMutSel ? xcom + com.nrgene : &com.kappa);
   if (com.seqtype == CODONseq && com.model) {
      if ((com.model == NSbranchB || com.model == NSbranch2) && com.NSsites == 0 && com.nbtype <= nUVR) {
         U = _UU[(int)nodes[a].label];
         V = _VV[(int)nodes[a].label];
         Root = _Root[(int)nodes[a].label];
      }
      else {
         eigenQcodon(1, -1, NULL, NULL, NULL, Root, U, V, &mr, pkappa, nodes[a].omega, PMat);
      }
   }
#endif

#if(BASEML)
   if (com.nhomo == 2)
      eigenTN93(com.model, *nodes[a].pkappa, 1, com.pi, &nR, Root, Cijk);
   nroot = nR;
#endif
   * l = *dl = *ddl = 0;
   for (ig = 0; ig < com.ngene; ig++) {
      if (com.Mgene > 1) SetPGene(ig, 1, 1, 0, xcom);  /* com.ntime=0 */
      for (i = 0; i < n * n; i++) P[i] = dP[i] = ddP[i] = 0;

      for (k = 0, expt = 1; k < nroot; k++) {
         multiply = com.rgene[ig] * Root[k];
         if (k) expt = exp(t * multiply);

#if (CODEML)  /* uses U & V */
         for (i = 0; i < n; i++)
            for (j = 0, uexpt = U[i * n + k] * expt; j < n; j++) {
               P[i * n + j] += uexpt * V[k * n + j];
               if (k) {
                  dP[i * n + j] += uexpt * V[k * n + j] * multiply;
                  ddP[i * n + j] += uexpt * V[k * n + j] * multiply * multiply;
               }
            }
#elif (BASEML) /* uses Cijk */
         for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
            P[i * n + j] += Cijk[i * n * nroot + j * nroot + k] * expt;
            if (k) {
               dP[i * n + j] += Cijk[i * n * nroot + j * nroot + k] * expt * multiply;
               ddP[i * n + j] += Cijk[i * n * nroot + j * nroot + k] * expt * multiply * multiply;
            }
         }
#endif
      }

      for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
         n1 = (b < com.ns ? nChara[(int)com.z[b][h]] : n);
         for (i = 0, fh = dfh = ddfh = 0; i < n1; i++) {
            xb = i;
            if (b < com.ns) piqi = com.pi[xb = CharaMap[(int)com.z[b][h]][i]];
            else            piqi = com.pi[i] * nodes[b].conP[h * n + i];
            for (j = 0, pqj = dpqj = ddpqj = 0; j < n; j++) {
               pqj += P[xb * n + j] * nodes[a].conP[h * n + j];
               dpqj += dP[xb * n + j] * nodes[a].conP[h * n + j];
               ddpqj += ddP[xb * n + j] * nodes[a].conP[h * n + j];
            }
            fh += piqi * pqj;
            dfh += piqi * dpqj;
            ddfh += piqi * ddpqj;
         }
         if (noisy && fh < 1e-250) {
            printf("too small: fh[%d] = %10.6e\n", h, fh);
            OutTreeN(F0, 0, 1);
         }
         *l -= log(fh) * com.fpatt[h];
         for (i = 0; i < com.NnodeScale; i++)
            *l -= com.nodeScaleF[i * com.npatt + h] * com.fpatt[h];
         *dl -= dfh / fh * com.fpatt[h];
         *ddl -= (fh * ddfh - dfh * dfh) / (fh * fh) * com.fpatt[h];
      }
   }  /* for(ig) */
   return(0);
}


int lfunt_SiteClass(double t, int a, int b, double xcom[], double* l, double space[])
{
   /* see notes in lfuntdd_SiteClass
    For branch&site models, look at the notes in GetPMatBranch()
    */
   int i, j, k, h, ig, ir, it, n = com.ncode, nroot = n;
   int n1 = (com.cleandata && b < com.ns ? 1 : n), xb;
   double y, expt, multiply, piqi, pqj;
   double* P = space, * fh = P + n * n;
   double* Sh = fh + com.npatt;  /* scale factor for each site pattern*/
   double* pK = com.fhK;  /* proportion for each site class after scaling */

#if (BASEML)
   if (com.nhomo == 2)
      eigenTN93(com.model, *nodes[a].pkappa, 1, com.pi, &nR, Root, Cijk);
   nroot = nR;
#endif

   if (com.NnodeScale == 0)
      for (ir = 0; ir < com.ncatG; ir++)
         for (h = 0; h < com.npatt; h++)
            pK[ir * com.npatt + h] = com.freqK[ir];
   else {
      for (h = 0; h < com.npatt; h++) {
         for (ir = 0, it = 0; ir < com.ncatG; ir++) {
            for (k = 0, y = 0; k < com.NnodeScale; k++)
               y += com.nodeScaleF[ir * com.NnodeScale * com.npatt + k * com.npatt + h];
            if ((pK[ir * com.npatt + h] = y) > pK[it * com.npatt + h])
               it = ir;
         }
         Sh[h] = pK[it * com.npatt + h];
         for (ir = 0; ir < com.ncatG; ir++)
            pK[ir * com.npatt + h] = com.freqK[ir] * exp(pK[ir * com.npatt + h] - Sh[h]);
      }
   }

   for (h = 0; h < com.npatt; h++) fh[h] = 0;
   for (ir = 0; ir < com.ncatG; ir++) {
      SetPSiteClass(ir, xcom);  /* com.ntime=0 */

#if CODEML  /* branch b->a */
        /* branch&site models */
      if (com.seqtype == CODONseq && com.NSsites && com.model)
         Set_UVR_BranchSite(ir, (int)nodes[a].label);
#endif

      if (ir) {
         for (i = com.ns; i < tree.nnode; i++)
            nodes[i].conP += (tree.nnode - com.ns) * n * (size_t)com.npatt;
      }
      for (ig = 0; ig < com.ngene; ig++) {
         if (com.Mgene > 1 || com.nalpha > 1)
            SetPGene(ig, com.Mgene > 1, com.Mgene > 1, com.nalpha > 1, xcom);  /* com.ntime=0 */
         if (com.nalpha > 1) SetPSiteClass(ir, xcom);    /* com.ntime=0 */

         for (i = 0; i < n * n; i++) P[i] = 0;
         for (k = 0, expt = 1; k < nroot; k++) {
            multiply = com.rgene[ig] * Root[k] * _rateSite;
#if (CODEML)
            if (com.seqtype == 1 && com.model >= 2)
               multiply *= Qfactor_NS_branch[(int)nodes[a].label];
#endif
            if (k) expt = exp(t * multiply);

#if (CODEML)  /* uses U & V */
            double uexpt = 0;
            for (i = 0; i < n; i++)
               for (j = 0, uexpt = U[i * n + k] * expt; j < n; j++)
                  P[i * n + j] += uexpt * V[k * n + j];
#elif (BASEML) /* uses Cijk */
            for (i = 0; i < n; i++)
               for (j = 0; j < n; j++)
                  P[i * n + j] += Cijk[i * n * nroot + j * nroot + k] * expt;
#endif
         }  /* for (k), look through eigenroots */
         for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
            n1 = (b < com.ns ? nChara[(int)com.z[b][h]] : n);
            for (i = 0; i < n1; i++) {
               xb = i;
               if (b < com.ns)
                  piqi = pK[ir * com.npatt + h] * com.pi[xb = CharaMap[(int)com.z[b][h]][i]];
               else
                  piqi = pK[ir * com.npatt + h] * com.pi[i] * nodes[b].conP[h * n + i];

               for (j = 0, pqj = 0; j < n; j++)
                  pqj += P[xb * n + j] * nodes[a].conP[h * n + j];
               fh[h] += piqi * pqj;
            }
         }  /* for (h) */
      }     /* for (ig) */
   }        /* for(ir) */

   for (i = com.ns; i < tree.nnode; i++)  /* shift position */
      nodes[i].conP -= (com.ncatG - 1) * (tree.nnode - com.ns) * n * (size_t)com.npatt;
   for (h = 0, *l = 0; h < com.npatt; h++) {
      if (fh[h] < 1e-250)
         printf("small (lfunt_SiteClass): fh[%d] = %10.6e\n", h, fh[h]);

      *l -= log(fh[h]) * com.fpatt[h];
      if (com.NnodeScale) *l -= Sh[h] * com.fpatt[h];
   }
   return(0);
}


int lfuntdd_SiteClass(double t, int a, int b, double xcom[],
   double* l, double* dl, double* ddl, double space[])
{
   /* dt and ddt for site-class models, modified from lfuntdd()
    nodes[].conP (and com.nodeScaleF if scaling is used) is shifted for ir,
    and moved back to the rootal place at the end of the routine.

    At the start of this routine, nodes[].conP has the conditional probabilties
    for each node, each site pattern, for each site class (ir).
    Scaling: When scaling is used, scale factors
    com.nodeScaleF[ir*com.NnodeScale*com.npatt + k*com.npatt+h] for all nodes
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
   int i, j, k, h, ig, ir, it, n = com.ncode, nroot = n;
   int n1 = (com.cleandata && b < com.ns ? 1 : n), xb;
   double y, expt, multiply, piqi, pqj, dpqj, ddpqj;
   double* P = PMat, * dP = space, * ddP = dP + n * n;
   double* fh = ddP + n * n, * dfh = fh + com.npatt, * ddfh = dfh + com.npatt;
   double* Sh = ddfh + com.npatt;  /* scale factor for each site pattern */
   double* pK = com.fhK;  /* proportion for each site class after scaling */

#if (BASEML)
   if (com.nhomo == 2)
      eigenTN93(com.model, *nodes[a].pkappa, 1, com.pi, &nR, Root, Cijk);
   nroot = nR;
#endif
   if (com.NnodeScale == 0)
      for (ir = 0; ir < com.ncatG; ir++)
         for (h = 0; h < com.npatt; h++)
            pK[ir * com.npatt + h] = com.freqK[ir];
   else {
      for (h = 0; h < com.npatt; h++) {
         for (ir = 0, it = 0; ir < com.ncatG; ir++) {
            for (k = 0, y = 0; k < com.NnodeScale; k++)
               y += com.nodeScaleF[ir * com.NnodeScale * com.npatt + k * com.npatt + h];
            if ((pK[ir * com.npatt + h] = y) > pK[it * com.npatt + h])
               it = ir;
         }
         Sh[h] = pK[it * com.npatt + h];
         for (ir = 0; ir < com.ncatG; ir++)
            pK[ir * com.npatt + h] = com.freqK[ir] * exp(pK[ir * com.npatt + h] - Sh[h]);
      }
   }

   for (h = 0; h < com.npatt; h++)
      fh[h] = dfh[h] = ddfh[h] = 0;
   for (ir = 0; ir < com.ncatG; ir++) {
      SetPSiteClass(ir, xcom);   /* com.ntime=0 */

#if CODEML  /* branch b->a */
        /* branch&site models */
      if (com.seqtype == CODONseq && com.NSsites && com.model)
         Set_UVR_BranchSite(ir, (int)nodes[a].label);
#endif

      if (ir) {
         for (i = com.ns; i < tree.nnode; i++)
            nodes[i].conP += (tree.nnode - com.ns) * n * (size_t)com.npatt;
      }
      for (ig = 0; ig < com.ngene; ig++) {
         if (com.Mgene > 1 || com.nalpha > 1)
            SetPGene(ig, com.Mgene > 1, com.Mgene > 1, com.nalpha > 1, xcom);   /* com.ntime=0 */
         if (com.nalpha > 1) SetPSiteClass(ir, xcom);   /* com.ntime=0 */

         for (i = 0; i < n * n; i++)
            P[i] = dP[i] = ddP[i] = 0;
         for (k = 0, expt = 1; k < nroot; k++) {   /* k loops through eigenroots */
            multiply = com.rgene[ig] * Root[k] * _rateSite;
#if (CODEML)
            if (com.seqtype == 1 && com.model >= 2)
               multiply *= Qfactor_NS_branch[(int)nodes[a].label];
#endif
            if (k) expt = exp(t * multiply);

#if (CODEML)  /* uses U & V */
            double uexpt = 0;
            for (i = 0; i < n; i++)
               for (j = 0, uexpt = U[i * n + k] * expt; j < n; j++) {
                  P[i * n + j] += uexpt * V[k * n + j];
                  if (k) {
                     dP[i * n + j] += uexpt * V[k * n + j] * multiply;
                     ddP[i * n + j] += uexpt * V[k * n + j] * multiply * multiply;
                  }
               }
#elif (BASEML) /* uses Cijk */
            for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
               P[i * n + j] += Cijk[i * n * nroot + j * nroot + k] * expt;
               if (k) {
                  dP[i * n + j] += Cijk[i * n * nroot + j * nroot + k] * expt * multiply;
                  ddP[i * n + j] += Cijk[i * n * nroot + j * nroot + k] * expt * multiply * multiply;
               }
            }
#endif
         }

         for (h = com.posG[ig]; h < com.posG[ig + 1]; h++) {
            n1 = (b < com.ns ? nChara[(int)com.z[b][h]] : n);
            for (i = 0; i < n1; i++) {
               xb = i;
               if (b < com.ns)
                  piqi = pK[ir * com.npatt + h] * com.pi[xb = CharaMap[(int)com.z[b][h]][i]];
               else
                  piqi = pK[ir * com.npatt + h] * com.pi[i] * nodes[b].conP[h * n + i];

               for (j = 0, pqj = dpqj = ddpqj = 0; j < n; j++) {
                  pqj += P[xb * n + j] * nodes[a].conP[h * n + j];
                  dpqj += dP[xb * n + j] * nodes[a].conP[h * n + j];
                  ddpqj += ddP[xb * n + j] * nodes[a].conP[h * n + j];
               }
               fh[h] += piqi * pqj;
               dfh[h] += piqi * dpqj;
               ddfh[h] += piqi * ddpqj;
            }
         }  /* for (h) */
      }     /* for (ig) */
   }        /* for(ir) */

   for (i = com.ns; i < tree.nnode; i++)
      nodes[i].conP -= (com.ncatG - 1) * (tree.nnode - com.ns) * n * (size_t)com.npatt;
   for (h = 0, *l = *dl = *ddl = 0; h < com.npatt; h++) {
      if (fh[h] < 1e-250)
         printf("small fh[%d] = %10.6e\n", h, fh[h]);

      *l -= log(fh[h]) * com.fpatt[h];
      if (com.NnodeScale) *l -= Sh[h] * com.fpatt[h];
      *dl -= dfh[h] / fh[h] * com.fpatt[h];
      *ddl -= (fh[h] * ddfh[h] - dfh[h] * dfh[h]) / (fh[h] * fh[h]) * com.fpatt[h];
   }

   return(0);
}

#endif


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
   int i, j, it, imin, fixt0 = 1;
   double la = birth, mu = death, rho = sample, tmin, r, t[NS - 1];
   double phi, eml, y;

   if (sample == 0)  /* coalescent model.  Check this!!!  */
      for (i = com.ns, y = 0; i > 1; i--)
         nodes[com.ns * 2 - i].age = y += -log(rndu()) / (i * (i - 1.) / 2.) * mut / 2;
   else {         /* BD with sampling */
      if (fixt0) t[com.ns - 2] = 1;
      if (fabs(la - mu) > 1e-6) {
         eml = exp(mu - la);
         phi = (rho * la * (eml - 1) + (mu - la) * eml) / (eml - 1);
         for (i = 0; i < com.ns - 1 - (fixt0); i++) {
            r = rndu();
            t[i] = log((phi - r * rho * la) / (phi - r * rho * la + r * (la - mu))) / (mu - la);
         }
      }
      else
         for (i = 0; i < com.ns - 1 - (fixt0); i++) {
            r = rndu();
            t[i] = r / (1 + la * rho * (1 - r));
         }
      /* bubble sort */
      for (i = 0; i < com.ns - 1 - 1; i++) {
         for (j = i + 1, tmin = t[i], imin = i; j < com.ns - 1; j++)
            if (tmin > t[j]) { tmin = t[j]; imin = j; }
         t[imin] = t[i];
         t[i] = tmin;
      }
      for (i = com.ns; i > 1; i--)
         nodes[com.ns * 2 - i].age = t[com.ns - i] * mut;
   }
   for (i = 0; i < com.ns; i++) nodes[i].age = 0;
   for (i = 0; i < tree.nnode; i++)
      if (i != tree.root)
         nodes[i].branch = nodes[nodes[i].father].age - nodes[i].age;
   if (!rooted) {
      it = nodes[tree.root].sons[2];
      nodes[it].branch = 2 * nodes[2 * com.ns - 2].age - nodes[tree.root].age - nodes[it].age;
   }
}

#endif


#ifdef NODESTRUCTURE
#ifdef EVOLVER

int RandomLHistory(int rooted, double space[])
{
   /* random coalescence tree, with each labeled history having equal probability.
    interior nodes are numbered ns, ns+1, ..., 2*ns-1-!rooted
    */
   int ns = com.ns, i, j, it = 0, * nodea = (int*)space;

   for (i = 0; i < 2 * ns - 1 - !rooted; i++)
      ClearNode(i);

   for (i = 0; i < ns; i++) nodea[i] = i;
   for (i = ns; i > (1 + !rooted); i--) {
      nodes[it = 2 * ns - i].nson = 2;
      j = (int)(i * rndu());
      nodes[nodea[j]].father = it;
      nodes[it].sons[0] = nodea[j];
      nodea[j] = nodea[i - 1];
      j = (int)((i - 1) * rndu());
      nodes[nodea[j]].father = it;
      nodes[it].sons[1] = nodea[j];
      nodea[j] = it;
      if (!rooted && i == 3) {
         nodes[it].nson++;
         nodes[nodea[1 - j]].father = it;
         nodes[it].sons[2] = nodea[1 - j];
      }
   }
   tree.root = it;
   tree.nnode = ns * 2 - 1 - !rooted;
   NodeToBranch();
   return (0);
}

#endif

#endif  /* NODESTRUCTURE */



#if (defined MCMCTREE)

int ProcessNodeAnnotation(int* haslabel)
{
   /* This processes fossil calibration information that has been read into
      nodes[].annotation, and copy the info into stree.nodes[].
      It uses both stree and nodes[], before nodes[] is destroyed.
      This is called before sequence alignments at loci are read.

      Further nodes:
      Other distributions such as G, SN, ST must be specified using the format 'G{alpha, beta}',
      say, and are processed here.  Simple bounds can also be specified using the format
      'L{0.5}', 'U{1.0}', or 'B{0.5, 1.0}', in which case they are processed here.
      I kept this complexity, as the notation of < and > is easy to understand, and because
      ReadTreeN can read node labels such as #, $ anyway.

      stree.duplication:
      stree.nodes[i].label = 0:  usual node, which may and may not have calibrations.
                            -1:  node i is the driving node, whose age is shared by other nodes.
      stree.nodes[j].label = i:  mirrored node j shares the age of node i, and can't have calibration.
   */
   int s = stree.nspecies, i, j, k, nfossiltype = 7;
   int* nodelabel, nmirror, mismatch;
   char* pch, * pch2, * braces = "({";
   double tailL = 0.025, tailR = 0.025, p_LOWERBOUND = 0.1, c_LOWERBOUND = 1.0;

   nodelabel = (int*)malloc(s * sizeof(int));
   if (nodelabel == NULL) error2("oom");
   for (i = 0; i < s - 1; i++)  stree.nodes[s + i].label = nodelabel[i] = 0;

   for (i = s; i < s * 2 - 1; i++) {
      if (nodes[i].annotation == NULL)
         continue;
      if ((pch = strchr(nodes[i].annotation, '#'))) {  /* for models for duplicated genes */
         *haslabel = 1;
         sscanf(pch + 1, "%lf", &stree.nodes[i].label);
         nodelabel[i - s] = (int)stree.nodes[i].label;
      }

      if (strchr(nodes[i].annotation, '>') && strchr(nodes[i].annotation, '<')) {
         stree.nfossil++;
         stree.nodes[i].fossil = BOUND_F;
         stree.nodes[i].pfossil[2] = tailL;
         stree.nodes[i].pfossil[3] = tailR;
         if ((pch = strchr(nodes[i].annotation, '>')))
            sscanf(pch + 1, "%lf", &stree.nodes[i].pfossil[0]);
         if ((pch = strchr(nodes[i].annotation, '<')))
            sscanf(pch + 1, "%lf", &stree.nodes[i].pfossil[1]);
         if (stree.nodes[i].pfossil[0] > stree.nodes[i].pfossil[1]) {
            printf("fossil bounds (%.4f, %.4f)", stree.nodes[i].pfossil[0], stree.nodes[i].pfossil[1]);
            error2("fossil bounds in tree incorrect");
         }
      }
      else if ((pch = strchr(nodes[i].annotation, '>'))) {
         stree.nfossil++;
         stree.nodes[i].fossil = LOWER_F;
         sscanf(pch + 1, "%lf", &stree.nodes[i].pfossil[0]);
         stree.nodes[i].pfossil[1] = p_LOWERBOUND;
         stree.nodes[i].pfossil[2] = c_LOWERBOUND;
         stree.nodes[i].pfossil[3] = tailL;
      }
      else if ((pch = strchr(nodes[i].annotation, '<'))) {
         stree.nfossil++;
         stree.nodes[i].fossil = UPPER_F;
         sscanf(pch + 1, "%lf", &stree.nodes[i].pfossil[1]);
         stree.nodes[i].pfossil[2] = tailR;
      }
      else {
         for (j = 1; j < nfossiltype + 1; j++)
            if ((pch = strstr(nodes[i].annotation, fossils[j]))) break;
         if (j == nfossiltype + 1) continue;

         stree.nfossil++;
         stree.nodes[i].fossil = j;
         pch2 = strchr(pch, braces[0]);
         if (pch2 == NULL) pch2 = strchr(pch, braces[1]);
         if (pch2 == NULL)
            error2("did not find left brace");
         pch = pch2 + 1;

         switch (j) {
         case (LOWER_F):
            /* truncated Cauchy default prior L(tL, p, c) */
            stree.nodes[i].pfossil[1] = p_LOWERBOUND;
            stree.nodes[i].pfossil[2] = c_LOWERBOUND;
            stree.nodes[i].pfossil[3] = tailL;
            sscanf(pch, "%lf,%lf,%lf,%lf", &stree.nodes[i].pfossil[0], &stree.nodes[i].pfossil[1],
               &stree.nodes[i].pfossil[2], &stree.nodes[i].pfossil[3]);
            break;
         case (UPPER_F):
            stree.nodes[i].pfossil[2] = tailR;
            sscanf(pch, "%lf,%lf", &stree.nodes[i].pfossil[1], &stree.nodes[i].pfossil[2]);
            break;
         case (BOUND_F):
            stree.nodes[i].pfossil[2] = tailL;
            stree.nodes[i].pfossil[3] = tailR;
            sscanf(pch, "%lf,%lf,%lf,%lf", &stree.nodes[i].pfossil[0], &stree.nodes[i].pfossil[1],
               &stree.nodes[i].pfossil[2], &stree.nodes[i].pfossil[3]);
            if (stree.nodes[i].pfossil[0] > stree.nodes[i].pfossil[1]) {
               printf("fossil bounds (%.4f, %.4f)", stree.nodes[i].pfossil[0], stree.nodes[i].pfossil[1]);
               error2("fossil bounds in tree incorrect");
            }
            break;
         case (GAMMA_F):
            sscanf(pch, "%lf,%lf", &stree.nodes[i].pfossil[0], &stree.nodes[i].pfossil[1]);
            break;
         case (SKEWN_F):
            sscanf(pch, "%lf,%lf,%lf", &stree.nodes[i].pfossil[0], &stree.nodes[i].pfossil[1], &stree.nodes[i].pfossil[2]);
            break;
         case (SKEWT_F):
            sscanf(pch, "%lf,%lf,%lf,%lf", &stree.nodes[i].pfossil[0], &stree.nodes[i].pfossil[1], &stree.nodes[i].pfossil[2], &stree.nodes[i].pfossil[3]);
            break;
         case (S2N_F):
            sscanf(pch, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &stree.nodes[i].pfossil[0], &stree.nodes[i].pfossil[1],
               &stree.nodes[i].pfossil[2], &stree.nodes[i].pfossil[3], &stree.nodes[i].pfossil[4],
               &stree.nodes[i].pfossil[5], &stree.nodes[i].pfossil[6]);
            break;
         }
      }
      stree.nodes[i].usefossil = 1;
      nodes[i].branch = 0;
      free(nodes[i].annotation);
      nodes[i].annotation = NULL;
   }

   /* With duplication dating, fossil calibrations are copied to the main node, while mirrored nodes
      have no calibrations.  If i is the main node, mirrored node j has stree.nodes[j].label = i.
   */
   if (stree.duplication) {
      if (*haslabel == 0) error2("need node labels for duplication dating...");
      for (i = s; i < s * 2 - 1; i++) printf("%4d", nodelabel[i - s]);
      printf("\n");

      for (i = s; i < s * 2 - 1; i++) {
         if (nodelabel[i - s] == 0)  continue;
         nmirror = 0;
         for (j = i + 1; j < s * 2 - 1; j++) {
            if (nodelabel[i - s] != nodelabel[j - s]) continue;
            nmirror++;
            stree.duplication++;
            /* node j mirros i. */
            if (stree.nodes[i].fossil && stree.nodes[j].fossil) {
               /* if both & j have calibrations, they must be the same */
               mismatch = 0;
               if (stree.nodes[i].fossil != stree.nodes[j].fossil) mismatch = 1;
               for (k = 0; k < npfossils[(int)stree.nodes[j].fossil]; k++)
                  if (fabs(stree.nodes[i].pfossil[k] - stree.nodes[j].pfossil[k]) > 1e-6) mismatch = 1;
               if (mismatch) {
                  printf("\nnodes %2d %2d share age but have different calibrations? ", i + 1, j + 1);
                  error2("err");
               }
               stree.nfossil--;
            }
            if (stree.nodes[i].fossil == 0 && stree.nodes[j].fossil) { /* copy fossil from j to i. */
               stree.nodes[i].usefossil = 1;
               stree.nodes[i].fossil = stree.nodes[j].fossil;
               for (k = 0; k < npfossils[(int)stree.nodes[j].fossil]; k++) {
                  stree.nodes[i].pfossil[k] = stree.nodes[j].pfossil[k];
                  stree.nodes[j].pfossil[k] = -1;
               }
            }
            nodelabel[j - s] = 0;
            stree.nodes[j].usefossil = 0;
            stree.nodes[j].fossil = 0;
            stree.nodes[j].label = i;
            stree.nodes[i].label = -1;

            if (noisy) printf("node %2d is mirrored by node %d.\n", i, j);
         }
         if (nmirror == 0) {
            printf("\nnode %2d has label %d, but is not mirrored by any other node? ", i, nodelabel[i]);
            error2("error");
         }
      }

      for (i = s; i < s * 2 - 1; i++) printf("%4d", nodelabel[i - s]);
      printf("\n");
      for (i = s; i < s * 2 - 1; i++) printf("%4.0f", stree.nodes[i].label);
      printf("\n");

      /* stree.duplication has the number of mirrored nodes, so that s - 1 - stree.duplication is the
      number of distinct ages on the tree */
      stree.duplication--;
      if (noisy)
         printf("\nduplication dating: %2d node ages are mirrored\n", stree.duplication);
   }
   free(nodelabel);
   return(0);
}

#endif

#if (defined BASEML || defined CODEML)

int ProcessNodeAnnotation(int* haslabel)
{
   int cladeLabels = 0, i, j;
   char* pch;

   *haslabel = 0;
   /* check and convert clade labels $ into # */
#if(defined(BASEML))
   if (com.seqtype == 0 && com.nhomo == 5) cladeLabels = 1;
#endif
   if (com.clock > 1 || (com.seqtype == 1 && com.model >= 2 && com.model <= 3))
      cladeLabels = 1;
   if (cladeLabels) {
      for (i = 0, j = 0; i < tree.nnode; i++) {
         if (nodes[i].label2 != -1) j++;
      }
      if (j) {  /* j is number of clade labels */
         DownTreeCladeLabel(tree.root, 0);
      }
   }
   for (i = 0; i < tree.nnode; i++) {
      /* default label is 0. */
      if (nodes[i].label == -1)
         nodes[i].label = 0;
      else {
         *haslabel = 1;
         if (nodes[i].annotation == NULL) continue;
         if ((pch = strchr(nodes[i].annotation, '#')))
            sscanf(pch + 1, "%lf", &nodes[i].label);
         free(nodes[i].annotation);
         nodes[i].annotation = NULL;
      }
   }
   for (i = 0, com.nbtype = 0; i < tree.nnode; i++) {
      if (i == tree.root) continue;
      j = (int)nodes[i].label;
      if (j + 1 > com.nbtype)
         com.nbtype = j + 1;
      if (j<0 || j>tree.nbranch - 1)
         error2("branch label in the tree (note labels start from 0 and are consecutive)");
   }
   if (com.nbtype > 1)
      printf("\n%d branch types are in tree. Stop if wrong.", com.nbtype);

   return(0);
}
#endif


/* routines for dating analysis of heterogeneous data */
#if (defined BASEML || defined CODEML || defined MCMCTREE)

void copy_from_Sptree(struct SPECIESTREE* stree, struct TREEN* nodes);
void copy_to_Sptree(struct SPECIESTREE* stree, struct TREEN* nodes);
int ReadMainTree(FILE *ftree, struct SPECIESTREE *stree);
int GenerateGtree_locus(int locus, int ns, int allocate_gnodes);

int ReadMainTree( FILE*ftree, struct SPECIESTREE *stree)
{
/* This reads main species tree into stree, and names into stree->nodes[].name.
   */
   int haslength, j, ntree;

   if (noisy) puts("\nReading main tree.");
   j = fscanf(ftree, "%d%d", &stree->nspecies, &ntree);
   if (ntree!=1) error2("maintree file can have only one tree");
   com.ns = stree->nspecies;
   if (com.ns > NS) error2("raise NS?");
   else if (com.ns < 2) error2("number of species <2?");
   /* to read species names into stree->nodes[].name */
   for (j = 0; j < stree->nspecies; j++)
      com.spname[j] = stree->nodes[j].name;
   nodes = nodes_t;
   ReadTreeN(ftree, &haslength, 1, 1);
   OutTreeN(F0, 1, 0); printf("\n");

   if (com.clock == 5 || com.clock == 6)
      for (j = 0; j < tree.nnode; j++)
         nodes[j].branch = nodes[j].label = 0;
   for (j = 0; j < tree.nnode; j++)
      if (nodes[j].label < 0) nodes[j].label = 0;  /* change -1 into 0 */

   /* copy main tree into stree */
#if(MCMCTREE)
   if (tree.nnode != 2 * com.ns - 1)
      error2("check and think about multifurcating trees.");
#endif
   copy_to_Sptree(stree, nodes);

   return(0);
}


int ReadTreeSeqs(FILE* fout)
{
   /* This reads the combined species tree, the fossil calibration information, and sequence
      data at each locus.  stree.nodes[].pfossil[] has tL, tU for bounds or alpha and beta
      for the gamma prior.

      This also constructs the gene tree at each locus, by pruning the main species tree.
   */
   FILE* fseq, * ftree;
   int i, j, locus, clean0 = com.cleandata;

   ftree = gfopen(com.treef, "r");
   ReadMainTree(ftree, &stree);
   /* process fossil calibration info 
     (((A1 , B1) '#1 B(0.1,0.3)', (A2 , B2) #1), C ) 'B(0.9, 1.1)';  
   */
#if (defined MCMCTREE)
   if (!tipdate.flag)
      ProcessNodeAnnotation(&i);
#else
   ProcessNodeAnnotation(&i);
#endif

   /* read sequences at each locus, construct gene tree by pruning stree */
   data.ngene = com.ndata;
   com.ndata = 1;
   fseq = gfopen(com.seqf, "r");
   if ((gnodes = (struct TREEN**)malloc(sizeof(struct TREEN*) * data.ngene)) == NULL)
      error2("oom");

   printf("\nReading sequence data..  %d loci\n", data.ngene);
   for (locus = 0; locus < data.ngene; locus++) {
      fprintf(fout, "\n\n*** Locus %d ***\n", locus + 1);
      printf("\n\n*** Locus %d ***\n", locus + 1);

      com.cleandata = (char)clean0;
      for (j = 0; j < stree.nspecies; j++)
         com.spname[j] = NULL;  /* points to nowhere */
#if (defined CODEML)
      if (com.seqtype == 1) {
         com.icode = data.icode[locus];
         setmark_61_64();
      }
#endif
      ReadSeq(fout, fseq, clean0, locus, 0);               /* allocates com.spname[] */
#if (defined CODEML)
      if (com.seqtype == 1) {
         if (com.sspace < max2(com.ngene + 1, com.ns) * (64 + 12 + 4) * sizeof(double)) {
            com.sspace = max2(com.ngene + 1, com.ns) * (64 + 12 + 4) * sizeof(double);
            if ((com.space = (double*)realloc(com.space, com.sspace)) == NULL)
               error2("oom space for #c");
         }
         InitializeCodon(fout, com.space);
      }
#endif

      data.ns[locus] = com.ns;
      data.ls[locus] = com.ls;
#if(MCMCTREE)
      if (data.datatype[locus] == MORPHC)
         ;
      else
#endif
      {
         if (com.seqtype == 0 || com.seqtype == 2)
            InitializeBaseAA(fout);
         fflush(fout);

#if(!defined MCMCTREE)
         if ((com.seqtype == 0 || com.seqtype == 2) && com.model == 0) {
            PatternWeightJC69like();
            if (fout) {
               fprintf(fout, "\n\nPrinting out site pattern counts\n");
               printPatterns(fout);
            }
         }
#endif

         xtoy(com.pi, data.pi[locus], com.ncode);
         data.cleandata[locus] = (char)com.cleandata;
         data.npatt[locus] = com.npatt;
         data.fpatt[locus] = com.fpatt; com.fpatt = NULL;
         for (i = 0; i < com.ns; i++) {
            data.z[locus][i] = com.z[i];
            com.z[i] = NULL;
         }
         printf("%3d patterns, %s\n", com.npatt, (com.cleandata ? "clean" : "messy"));
      }

      GenerateGtree_locus(locus, data.ns[locus], 1);    /* free com.spname[] */
   }  /* for(locus) */
   if (noisy) printf("\n");
   for (i = 0, com.cleandata = 1; i < data.ngene; i++)
      if (data.cleandata[i] == 0)
         com.cleandata = 0;

   fclose(ftree); fclose(fseq);
   SetMapAmbiguity(com.seqtype, 0);

#if(defined (MCMCTREE))
   if (tipdate.flag)
      GetTipDate();
#endif

   return(0);
}


int GenerateGtrees(FILE* fout, FILE* fseq, char *treesfilename)
{
   FILE* ftrees = gfopen(treesfilename, "w");
   int locus, s=stree.nspecies, i;

   for (locus = 0; locus < com.ndata; locus++) {
      for (i = 0; i < s; i++) com.spname[i] = NULL;
      ReadSeq(fout, fseq, 0, locus, 1);   /* allocates com.spname[] */
      GenerateGtree_locus(locus, com.ns, 0);
      printf("\nlocus %2d: ", locus + 1);
      OutTreeN(stdout, 1, PrLabel);  printf("\n");
      fprintf(ftrees, "\n%d  1\n", com.ns);
      OutTreeN(ftrees, 1, PrLabel);  fprintf(ftrees, "\n");
   }
   fclose(ftrees);
   return (0);
}

int GenerateGtree_locus(int locus, int ns, int allocate_gnodes)
{
/* construct the gene tree at locus in nodes_t.
   It prunes off tips in the main species tree stree.  
   com.spname[] have names of species at the current locus (read from the sequence alignment at the locus).  
   They are used by the routine to compare with stree.nodes[].name to decide which species to keep for 
   the locus.
   See GetSubTreeN() for more details.
   */
   int nnode0=stree.nnode, i, j, keep[NS], * newnodeNO;

   if (ns > stree.nspecies) 
      error2("ns in seqfile is > ns in main tree");
   newnodeNO = (int*)malloc(nnode0 * sizeof(int));
   if (newnodeNO == NULL) error2("oom GenerateGtree_locus");

   for (j = 0; j < stree.nspecies; j++) keep[j] = 0;
   for (i = 0; i < ns; i++) {
      for (j = 0; j < stree.nspecies; j++)
         if (!strcmp(com.spname[i], stree.nodes[j].name)) break;
      if (j == stree.nspecies) {
         printf("species %s not found in main tree\n", com.spname[i]);
         exit(-1);
      }
      if (keep[j]) {
         printf("\nspecies %s occurs twice in locus %d", com.spname[i], locus + 1);
         error2("\ngiving up...");
      }
      keep[j] = i + 1;
      free(com.spname[i]);
   }

   /* copy main species tree and then prune it. */
   nodes = nodes_t;
   copy_from_Sptree(&stree, nodes);
   GetSubTreeN(keep, newnodeNO);
   com.ns = ns;
   NodeToBranch();

#if(MCMCTREE)
   for (i = 0; i < stree.nnode; i++)
      if (newnodeNO[i] != -1) nodes[newnodeNO[i]].ipop = i;
   /* printGtree(0);  */
   if (allocate_gnodes) {
      gnodes[locus] = (struct TREEN*)malloc((ns * 2 - 1) * sizeof(struct TREEN));
      if (gnodes[locus] == NULL) error2("oom gtree");
      memcpy(gnodes[locus], nodes, (ns * 2 - 1) * sizeof(struct TREEN));
      data.root[locus] = tree.root;
   }
#endif
   free(newnodeNO);
   return(0);
}

int printGtree(int printBlength)
{
   int i, j;

   for (i = 0; i < com.ns; i++)
      com.spname[i] = stree.nodes[nodes[i].ipop].name;
   for (i = 0; i < tree.nnode; i++)
      if (i != tree.root)
         nodes[i].branch = nodes[nodes[i].father].age - nodes[i].age;
   printf("\nns = %d  nnode = %d", com.ns, tree.nnode);
   printf("\n%7s%7s %8s %7s%7s", "father", "node", "(ipop)", "nson:", "sons");
   for (i = 0; i < tree.nnode; i++) {
      printf("\n%7d%7d   (%2d) %7d  ",
         nodes[i].father + 1, i + 1, nodes[i].ipop + 1, nodes[i].nson);
      for (j = 0; j < nodes[i].nson; j++) printf(" %2d", nodes[i].sons[j] + 1);
   }
   printf("\n");
   OutTreeN(F0, 0, 0); printf("\n");
   OutTreeN(F0, 1, 0); printf("\n");
   if (printBlength) {
      OutTreeN(F0, 1, 1); printf("\n");
   }
   return(0);
}


void copy_from_Sptree(struct SPECIESTREE *stree, struct TREEN *nodes)
{
/* This copies stree into nodes, tree, & com.spname, for printing or editing
   */
   int i, j;

   com.ns = stree->nspecies;   tree.root = stree->root;
   tree.nnode = stree->nnode;  tree.nbranch = stree->nbranch;
   for (i = 0; i < stree->nnode; i++) {
      /* this is used by mcmctree */
      if (i < com.ns)
         com.spname[i] = stree->nodes[i].name; 
      /* The following may be needed by bpp.  Check carefully. */
      /*
       strcpy(com.spname[i], stree->nodes[i].name);

       com.spname[i] = nodes[i].name = stree->nodes[i].name; 
       */
      nodes[i].father = stree->nodes[i].father;
      nodes[i].nson = stree->nodes[i].nson;
      for (j = 0; j < nodes[i].nson; j++)
         nodes[i].sons[j] = stree->nodes[i].sons[j];
      nodes[i].fossil = stree->nodes[i].fossil;
      nodes[i].age = stree->nodes[i].age;
      nodes[i].label = stree->nodes[i].label;
      if (i != tree.root)
         nodes[i].branch = stree->nodes[nodes[i].father].age - stree->nodes[i].age;
   }
}


void copy_to_Sptree(struct SPECIESTREE* stree, struct TREEN* nodes)
{
   /* This copies nodes, tree, & com.spname into stree (main stree)  */
   int i, j;

   stree->nspecies = com.ns;
   for (j = 0; j < stree->nspecies; j++)
      strcpy(stree->nodes[j].name, com.spname[j]);
   stree->nnode = tree.nnode;  stree->nbranch = tree.nbranch;
   stree->root = tree.root;    stree->nfossil = 0;
   for (i = 0; i < stree->nnode; i++) {
      stree->nodes[i].father = nodes[i].father;
      stree->nodes[i].nson = nodes[i].nson;
      stree->nodes[i].label = nodes[i].label;
#if(MCMCTREE)
      if (nodes[i].nson != 0 && nodes[i].nson != 2)
         error2("main tree is not binary.");
#endif
      for (j = 0; j < stree->nodes[i].nson; j++)
         stree->nodes[i].sons[j] = nodes[i].sons[j];
   }
   for (j = 0; j < stree->nspecies; j++)
      com.spname[j] = NULL;
}


void printStree(void)
{
   int i, k, maxlname = 10;

   for (i = 0; i < stree.nspecies; i++)
      if ((k = (int)strlen(stree.nodes[i].name)) > maxlname)
         maxlname = k;

   printf("\n************\nSpecies tree\nns = %d  nnode = %d", stree.nspecies, stree.nnode);
   printf("\n%7s%7s  %-*s %9s %8s%16s\n", "father", "node", maxlname + 2, "name", "time", "sons", "fossil");
   for (i = 0; i < stree.nnode; i++) {
      printf("%7d%7d  %-*s %8.5f",
         stree.nodes[i].father + 1, i + 1, maxlname + 6, stree.nodes[i].name, stree.nodes[i].age);
      if (stree.nodes[i].nson)
         printf("  (%2d %2d)", stree.nodes[i].sons[0] + 1, stree.nodes[i].sons[1] + 1);

#ifdef MCMCTREE
      int j;

      if ((k = stree.nodes[i].fossil)) {
         printf("  %s ( ", fossils[k]);
         for (j = 0; j < npfossils[k]; j++) {
            printf("%6.4f", stree.nodes[i].pfossil[j + (k == UPPER_F)]);
            printf("%s", (j == npfossils[k] - 1 ? " ) " : ", "));
         }
      }
#endif

      printf("\n");
   }
   nodes = nodes_t;
   copy_from_Sptree(&stree, nodes);
   printf("\n");
   OutTreeN(F0, 0, 0);  printf("\n");
   OutTreeN(F0, 1, 0);  printf("\n");
   OutTreeN(F0, 1, 1);  printf("\n");
}


#endif




#if (defined BASEML || defined CODEML)

#if (defined CODEML)

int GetMemPUVR(int nc, int nUVR)
{
   /* this gets mem for nUVR sets of matrices
    */
   int i;

   PMat = (double*)malloc((nc * nc + nUVR * nc * nc * 2 + nUVR * nc) * sizeof(double));
   if (PMat == NULL) error2("oom getting P&U&V&Root");
   U = _UU[0] = PMat + nc * nc;
   V = _VV[0] = _UU[0] + nc * nc;
   Root = _Root[0] = _VV[0] + nc * nc;
   for (i = 1; i < nUVR; i++) {
      _UU[i] = _UU[i - 1] + nc * nc * 2 + nc; _VV[i] = _VV[i - 1] + nc * nc * 2 + nc;
      _Root[i] = _Root[i - 1] + nc * nc * 2 + nc;
   }
   return(0);
}

void FreeMemPUVR(void)
{
   free(PMat);
}


int GetUVRoot_codeml(void)
{
   /* This uses data.daafile[] to set up the eigen matrices U, V, Root for
    combined clock analyses of multiple protein data sets (clock = 5 or 6).
    */
   int locus, nc = (com.seqtype == 1 ? 64 : 20), nUVR = data.ngene;
   double mr = 0;

   if (com.seqtype == 1 && (!com.fix_kappa || !com.fix_omega)) nUVR = 1;
   GetMemPUVR(nc, nUVR);

   if (nUVR > 6) error2("The maximum number of proteins is set to 6.");
   if (com.seqtype == 2) {
      for (locus = 0; locus < data.ngene; locus++) {
         if (data.ngene > 1)
            strcpy(com.daafile, data.daafile[locus]);
         GetDaa(NULL, com.daa);
         if (com.model == Empirical_F)
            xtoy(data.pi[locus], com.pi, nc);
         eigenQaa(NULL, _Root[locus], _UU[locus], _VV[locus], NULL);
      }
   }
   else if (com.seqtype == 1 && com.fix_kappa & com.fix_omega) {
      for (locus = 0; locus < data.ngene; locus++) {
         if (com.seqtype == 1) {
            com.icode = data.icode[locus];
            setmark_61_64();
         }
         com.kappa = data.kappa[locus];
         com.omega = data.omega[locus];
         xtoy(data.pi[locus], com.pi, com.ncode);
         eigenQcodon(1, -1, NULL, NULL, NULL, _Root[locus], _UU[locus], _VV[locus], &mr,
            &com.kappa, com.omega, PMat);
      }
   }
   return(0);
}


#endif


int UseLocus(int locus, int copycondP, int setmodel, int setSeqName)
{
   /* This point nodes to the gene tree at locus gnodes[locus] and set com.z[]
    etc. for likelihood calculation for the locus.
    */
   int i, nS;

   com.ns = data.ns[locus]; com.ls = data.ls[locus];
   tree.root = data.root[locus];
   tree.nnode = 2 * com.ns - 1;  /* assumes binary tree */
   tree.nbranch = tree.nnode - 1;

   nodes = gnodes[locus];

   com.cleandata = data.cleandata[locus];
   com.npatt = com.posG[1] = data.npatt[locus];  com.posG[0] = 0;
   com.fpatt = data.fpatt[locus];
   for (i = 0; i < com.ns; i++) com.z[i] = data.z[locus][i];

   /* The following is model-dependent */
   if (setmodel) {

      com.kappa = data.kappa[locus];
      com.omega = data.omega[locus];
      com.alpha = data.alpha[locus];

#if(defined CODEML)
      if (com.seqtype == 1) {
         com.icode = data.icode[locus];
         setmark_61_64();
      }
#endif

#if(defined BASEML)
      if (com.seqtype == 0 && com.model != 0 && com.model != 1)
         xtoy(data.pi[locus], com.pi, com.ncode);
      if (com.model <= TN93)
         eigenTN93(com.model, com.kappa, com.kappa, com.pi, &nR, Root, Cijk);
      else if (com.model == REV)
         eigenQREVbase(NULL, PMat, &com.kappa, com.pi, &nR, Root, Cijk);
#else
      double mr = 0;
      if ((com.seqtype == 1 && com.codonf) || (com.seqtype == 2 && com.model == 3))
         xtoy(data.pi[locus], com.pi, com.ncode);

      if ((com.seqtype == 2 && (com.model == 2 || com.model == 3))
         || (com.seqtype == 1 && com.fix_kappa && com.fix_omega)) {
         Root = _Root[locus]; U = _UU[locus];  V = _VV[locus];
      }
      else {
         eigenQcodon(1, -1, NULL, NULL, NULL, Root, U, V, &mr, &com.kappa, com.omega, PMat);
      }

#endif
      if (com.alpha)
         DiscreteGamma(com.freqK, com.rK, com.alpha, com.alpha, com.ncatG, DGammaUseMedian);

      com.NnodeScale = data.NnodeScale[locus];
      com.nodeScale = data.nodeScale[locus];
      nS = com.NnodeScale * com.npatt * (com.conPSiteClass ? com.ncatG : 1);
      for (i = 0; i < nS; i++) com.nodeScaleF[i] = 0;
   }
   if (setSeqName)
      for (i = 0; i < com.ns; i++)
         com.spname[i] = stree.nodes[nodes[i].ipop].name;
   return(0);
}


void GetMemBC(void)
{
   /* This gets memory for baseml and codeml under local clock models for analysis
      of combined data from multiple loci.
      com.conP[] is shared across loci.
      fhK[] uses shared space for loci.
   */
   int j, locus, nc = (com.seqtype == 1 ? 64 : com.ncode);
   size_t maxsizeScale = 0, nS, sfhK = 0, s1, sizenode;

   for (locus = 0, com.sconP = 0; locus < data.ngene; locus++) {
      sizenode = nc * data.npatt[locus];
      s1 = sizenode * (data.ns[locus] - 1) * sizeof(double);
      if (com.alpha) {     /* this is for step 1, using method = 1 */
         com.conPSiteClass = 1;
         s1 *= com.ncatG;
      }
      if (s1 > com.sconP) com.sconP = s1;
      if (com.alpha && (size_t)data.npatt[locus] > sfhK)
         sfhK = data.npatt[locus];
   }

   com.conP = (double*)malloc(com.sconP);
   printf("\n%5zu bytes for conP\n", com.sconP);
   if (com.conP == NULL)
      error2("oom conP");
   if (com.alpha) {
      sfhK *= com.ncatG * sizeof(double);
      if ((com.fhK = (double*)realloc(com.fhK, sfhK)) == NULL) error2("oom");
   }

   /* set gnodes[locus][].conP for internal nodes */
   for (locus = 0; locus < data.ngene; locus++) {
      sizenode = nc * data.npatt[locus];
      for (j = data.ns[locus]; j < data.ns[locus] * 2 - 1; j++)
         gnodes[locus][j].conP = com.conP + (j - data.ns[locus]) * sizenode;
   }
   for (locus = 0; locus < data.ngene; locus++) {
      if (!data.cleandata[locus]) {
         UseLocus(locus, -1, 0, 0);
      }
   }

   if (stree.nspecies > 20) {
      for (locus = 0; locus < data.ngene; locus++) {
         UseLocus(locus, -1, 0, 0);
         com.NnodeScale = 0;
         com.nodeScale = data.nodeScale[locus] = (char*)malloc(tree.nnode * sizeof(char));
         if (com.nodeScale == NULL)  error2("oom");
         for (j = 0; j < tree.nnode; j++) com.nodeScale[j] = 0;

         SetNodeScale(tree.root);

         data.NnodeScale[locus] = com.NnodeScale;
         nS = com.NnodeScale * com.npatt;
         if (com.conPSiteClass) nS *= com.ncatG;
         maxsizeScale = max2(maxsizeScale, nS);

         if (com.NnodeScale) {
            printf("\n%d node(s) used for scaling at locus %d: \n", com.NnodeScale, locus + 1);
            for (j = 0; j < tree.nnode; j++)
               if (com.nodeScale[j]) printf(" %2d", j + 1);
            printf("\n");
         }
      }
      if (maxsizeScale) {
         if ((com.nodeScaleF = (double*)malloc(maxsizeScale * sizeof(double))) == NULL)
            error2("oom nscale");
         for (j = 0; j < (int)maxsizeScale; j++) com.nodeScaleF[j] = 0;
      }
   }

}

void FreeMemBC(void)
{
   int locus, j;

   for (locus = 0; locus < data.ngene; locus++)
      free(gnodes[locus]);
   free(gnodes);
   free(com.conP);
   for (locus = 0; locus < data.ngene; locus++) {
      free(data.fpatt[locus]);
      for (j = 0; j < data.ns[locus]; j++)
         free(data.z[locus][j]);
   }
   if (com.alpha)
      free(com.fhK);

   if (stree.nspecies > 20) {
      for (locus = 0; locus < data.ngene; locus++)
         free(data.nodeScale[locus]);
      if (com.nodeScaleF) free(com.nodeScaleF);
   }
}




double nu_AHRS = 0.001, * varb_AHRS;


double funSS_AHRS(double x[], int np);


double lnLfunHeteroData(double x[], int np)
{
   /* This calculates the log likelihood, the log of the probability of the data
    given gtree[] for each locus.  This is for step 3 of Yang (2004. Acta
    Zoologica Sinica 50:645-656)
    x[0,1,...s-k] has node ages in the species tree, followed by branch rates
    for genes 1, 2, ..., then kappa for genes, then alpha for genes
    */
   int i, k, locus;
   double lnL = 0, lnLt, * pbrate;

   /* ??? need more work for codon sequences */
   for (locus = 0, k = com.ntime - 1; locus < data.ngene; locus++)
      k += data.nbrate[locus];
   if (!com.fix_kappa) for (i = 0; i < data.ngene; i++) data.kappa[i] = x[k++];
   if (!com.fix_omega) for (i = 0; i < data.ngene; i++) data.omega[i] = x[k++];
   if (!com.fix_alpha) for (i = 0; i < data.ngene; i++) data.alpha[i] = x[k++];

   /* update node ages in species tree */
   nodes = nodes_t;
   copy_from_Sptree(&stree, nodes);
   SetBranch(x);
   for (i = 0; i < tree.nnode; i++) stree.nodes[i].age = nodes[i].age;

   for (locus = 0, pbrate = x + com.ntime - 1; locus < data.ngene; locus++) {

      UseLocus(locus, -1, 1, 1);
      /* copy node ages to gene tree */
      for (i = 0; i < tree.nnode; i++)   nodes[i].age = stree.nodes[nodes[i].ipop].age;
      for (i = 0; i < tree.nnode; i++) {
         if (i != tree.root) {
            nodes[i].branch = (nodes[nodes[i].father].age - nodes[i].age)
               * pbrate[(int)nodes[i].label];
            if (nodes[i].branch < -1e-4)
               puts("b<0");
         }
      }
      lnL += lnLt = com.plfun(x, -1);
      pbrate += data.nbrate[locus];
   }
   return(lnL);
}


double funSS_AHRS(double x[], int np)
{
   /* Function to be minimized in the ad hoc rate smoothing procedure:
    lnLb + lnLr
    nodes[].label has node rate.
    lnLb is weighted sum of squares using approximate variances for branch lengths.

    lnLr is the log of the prior of rates under the geometric Brownian motion
    model of rate evolution. There is no need for recursion as the order at
    which stree.nodes are visited is unimportant.  The rates are stored in
    gnodes[].label.
    The root rate is fixed to be the weighted average rate of its two sons,
    inversely weighted by the divergence times.
    */
   int locus, j, k, root, pa, son0, son1;
   double lnLb, lnLr, lnLbi, lnLri;  /* lnLb & lnLr are sum of squares for b and r */
   double b, be, t, t0, t1, r, rA, w, y, smallv = 1e-20, smallage = AgeLow[stree.root] * smallv;
   double nu = nu_AHRS, * varb = varb_AHRS;

   /* set up node ages in species tree */
   nodes = nodes_t; 
   copy_from_Sptree(&stree, nodes);
   SetBranch(x);
   for (j = 0; j < tree.nnode; j++)
      stree.nodes[j].age = nodes[j].age;

   k = com.ntime - 1;
   for (locus = 0, lnLb = lnLr = 0; locus < data.ngene; varb += com.ns * 2 - 1, locus++) {
      UseLocus(locus, -1, 0, 0);
      if (data.fix_nu == 2)      nu = x[np - 1];
      else if (data.fix_nu == 3) nu = x[np - 1 - (data.ngene - 1 - locus)];

      root = tree.root;
      son0 = nodes[root].sons[0];
      son1 = nodes[root].sons[1];
      /* copy node ages and rates into gene tree nodes[]. */
      for (j = 0; j < tree.nnode; j++) { /* age and rates */
         nodes[j].age = stree.nodes[nodes[j].ipop].age;
         if (j != root)
            nodes[j].label = x[k++];
      }
      t0 = nodes[root].age - nodes[son0].age;
      t1 = nodes[root].age - nodes[son1].age;
      if (t0 + t1 < 1e-7)
         error2("small root branch.  Think about what to do.");
      nodes[root].label = (nodes[son0].label * t1 + nodes[son1].label * t0) / (t0 + t1);

      for (j = 0, lnLbi = 0; j < tree.nnode; j++) {
         if (j == son0 || j == son1) continue;
         pa = nodes[j].father;
         if (j == root) {
            b = nodes[son0].branch + nodes[son1].branch;
            be = (nodes[j].age - nodes[son0].age) * (nodes[root].label + nodes[son0].label) / 2
               + (nodes[j].age - nodes[son1].age) * (nodes[root].label + nodes[son1].label) / 2;
         }
         else {
            b = nodes[j].branch;
            be = (nodes[pa].age - nodes[j].age) * (nodes[pa].label + nodes[j].label) / 2;
         }
         w = varb[j];
         if (w < smallv)
            puts("small variance");
         lnLbi -= square(be - b) / (2 * w);
      }

      for (j = 0, lnLri = 0; j < tree.nnode; j++) {
         if (j == root) continue;
         pa = nodes[j].father;
         t = nodes[pa].age - nodes[j].age;
         t = max2(t, smallage);
         r = nodes[j].label;
         rA = nodes[pa].label;

         if (rA < smallv || t < smallv || r < smallv)  puts("small r, rA, or t");
         y = log(r / rA) + t * nu / 2;
         lnLri -= y * y / (2 * t * nu) - log(r) - log(2 * Pi * t * nu) / 2;
      }

      if (data.fix_nu > 1) lnLri += -nu / nu_AHRS - log(nu);  /* exponential prior */
      lnLb -= lnLbi;
      lnLr -= lnLri;
   }
   return (lnLb + lnLr);
}


void SetBranchRates(int inode)
{
   /* this uses node rates to set branch rates, and is used only after the ad hoc
    rate smoothing iteration is finished.
    */
   int i;
   if (inode < com.ns)
      nodes[inode].label = (nodes[inode].label + nodes[nodes[inode].father].label) / 2;
   else
      for (i = 0; i < nodes[inode].nson; i++)
         SetBranchRates(nodes[inode].sons[i]);
}


int GetInitialsClock6Step1(double x[], double xb[][2])
{
   /* This is for clock 6 step 1.
    */
   int i, k;
   double tb[] = { .0001, 999 };

   com.ntime = k = tree.nbranch;
   GetInitialsTimes(x);

   com.plfun = (com.alpha == 0 ? lfun : lfundG);
   com.conPSiteClass = (com.method && com.plfun == lfundG);

   InitializeNodeScale();

   if (com.seqtype == 0)  com.nrate = !com.fix_kappa;

   com.np = com.ntime + !com.fix_kappa + !com.fix_alpha;
   if (com.seqtype == 1 && !com.fix_omega) com.np++;

   if (!com.fix_kappa) x[k++] = com.kappa;
   if (!com.fix_omega) x[k++] = com.omega;
   if (!com.fix_alpha) x[k++] = com.alpha;
   NodeToBranch();

   for (i = 0; i < com.ntime; i++)
   {
      xb[i][0] = tb[0]; xb[i][1] = tb[1];
   }
   for (; i < com.np; i++)
   {
      xb[i][0] = .001; xb[i][1] = 999;
   }

   if (noisy > 3 && com.np < 200) {
      printf("\nInitials (np=%d)\n", com.np);
      for (i = 0; i < com.np; i++) printf(" %10.5f", x[i]);
      printf("\n");
      for (i = 0; i < com.np; i++) printf(" %10.5f", xb[i][0]);
      printf("\n");
      for (i = 0; i < com.np; i++) printf(" %10.5f", xb[i][1]);
      printf("\n");
   }
   return (0);
}



int GetInitialsClock56Step3(double x[])
{
   /* This is for clock 5 or clock 6 step 3
    */
   int i, j, k = 0;

   if (com.clock == 5)
      GetInitialsTimes(x);

   com.plfun = (com.alpha == 0 ? lfun : lfundG);
   com.conPSiteClass = (com.method && com.plfun == lfundG);

   InitializeNodeScale();

   com.np = com.ntime - 1 + (1 + !com.fix_kappa + !com.fix_omega + !com.fix_alpha) * data.ngene;
   if (com.clock == 5)
      for (i = com.ntime - 1; i < com.np; i++) x[i] = .2 + rndu();
   else if (com.clock == 6) {
      for (j = 0, k = com.ntime - 1; j < data.ngene; k += data.nbrate[j], j++)
         com.np += data.nbrate[j] - 1;
      if (!com.fix_kappa)
         for (j = 0; j < data.ngene; j++) x[k++] = data.kappa[j];
      if (!com.fix_omega)
         for (j = 0; j < data.ngene; j++) x[k++] = data.omega[j];
      if (!com.fix_alpha)
         for (j = 0; j < data.ngene; j++) x[k++] = data.alpha[j];
      for (i = k; i < com.np; i++) x[i] = (.5 + rndu()) / 2;
   }
   return (0);
}


double GetMeanRate(void)
{
   /* This gets the rough average rate for the locus
    */
   int inode, i, j, k, ipop, nleft, nright, marks[NS], sons[2], nfossil;
   double mr, md;

   mr = 0; nfossil = 0;
   for (inode = com.ns; inode < tree.nnode; inode++) {
      ipop = nodes[inode].ipop;
      if (stree.nodes[ipop].fossil == 0) continue;
      sons[0] = nodes[inode].sons[0];
      sons[1] = nodes[inode].sons[1];
      for (i = 0, nleft = nright = 0; i < com.ns; i++) {
         for (j = i, marks[i] = 0; j != tree.root; j = nodes[j].father) {
            if (j == sons[0]) { marks[i] = 1; nleft++;  break; }
            else if (j == sons[1]) { marks[i] = 2; nright++; break; }
         }
      }
      if (nleft == 0 || nright == 0) {
         puts("this calibration is not in gene tree.");
         continue;
      }
      nfossil++;

      for (i = 0, md = 0; i < com.ns; i++) {
         for (j = 0; j < com.ns; j++) {
            if (marks[i] == 1 && marks[j] == 2) {
               for (k = i; k != inode; k = nodes[k].father)
                  md += nodes[k].branch;
               for (k = j; k != inode; k = nodes[k].father)
                  md += nodes[k].branch;
            }
         }
      }
      md /= (nleft * nright);
      mr += md / (stree.nodes[ipop].age * 2);

      /*
       printf("node age & mr n%-4d %9.5f%9.5f  ", inode, stree.nodes[ipop].age, md);
       if(com.ns<100) for (i= 0; i< com.ns; i++) printf("%d", marks[i]);
       printf("\n");
       */
   }
   mr /= nfossil;
   if (nfossil == 0) {
      printf("need fossils for this locus\n");
      exit(-1);
   }

   return(mr);
}


int AdHocRateSmoothing(FILE* fout, double x[NS * 3], double xb[NS * 3][2], double space[])
{
   /* ad hoc rate smoothing for likelihood estimation of divergence times.
    Step 1: Use JC69 to estimate branch lengths under no-clock model.
    Step 2: ad hoc rate smoothing, estimating one set of divergence times
    and many sets of branch rates for loci.  Rate at root is set to
    weighted average of rate at the two sons.
    */
   int model0 = com.model;  /* is this useful? */
   int fix_kappa0 = com.fix_kappa, fix_omega0 = com.fix_omega, fix_alpha0 = com.fix_alpha;
   int ib, son0, son1;
   double kappa0 = com.kappa, omega0 = com.omega, alpha0 = com.alpha, t0, t1, * varb;
   double f, e = 1e-8, pb = 0.00001, rb[] = { 0.001,99 }, lnL, lnLsum = 0;
   double mbrate[20], Rj[20], r, minr, maxr, beta, * pnu = &nu_AHRS, nu, mr[NGENE];
   int i, j, k, k0, locus, nbrate[20], maxnbrate = 20;
   char timestr[32];
   FILE* fBV = gfopen("in.BV", "w");
   FILE* fdist = gfopen("RateDist.txt", "w");
   FILE* finStep1 = fopen("in.ClockStep1", "r"),
      * finStep2 = fopen("in.ClockStep2", "r");

   noisy = 4;
   for (locus = 0, k = 0; locus < data.ngene; locus++)
      k += 2 * data.ns[locus] - 1;
   if ((varb_AHRS = (double*)malloc(k * sizeof(double))) == NULL)
      error2("oom AHRS");
   for (i = 0; i < k; i++)  varb_AHRS[i] = -1;


   /* Step 1: Estimate branch lengths without clock.  */
   printf("\nStep 1: Estimate branch lengths under no clock.\n");
   fprintf(fout, "\n\nStep 1: Estimate branch lengths under no clock.\n");
   com.clock = 0; com.method = 1;
   /*
    com.model=0;  com.fix_kappa=1; com.kappa=1;
    com.fix_alpha=1; com.alpha=0;
    */
   for (locus = 0; locus < data.ngene; locus++) {
      if (!com.fix_kappa) data.kappa[locus] = com.kappa;
      if (!com.fix_omega) data.omega[locus] = com.omega;
      if (!com.fix_alpha) data.alpha[locus] = com.alpha;
   }
   for (locus = 0, varb = varb_AHRS; locus < data.ngene; varb += com.ns * 2 - 1, locus++) {
      UseLocus(locus, -1, 1, 1);

      fprintf(fout, "\nLocus %d (%d sequences)\n", locus + 1, com.ns);

      son0 = nodes[tree.root].sons[0];
      son1 = nodes[tree.root].sons[1];

      GetInitialsClock6Step1(x, xb);

      lnL = 0;
      if (com.ns > 30) fprintf(frub, "\n\nLocus %d\n", locus + 1);
      if (finStep1) {
         puts("read MLEs from step 1 from file");
         for (i = 0; i < com.np; i++)
            fscanf(finStep1, "%lf", &x[i]);
      }
      else {
         j = minB((com.ns > 30 ? frub : NULL), &lnL, x, xb, e, space);
         for (j = 0; j < com.ns * 2 - 1; j++) {
            ib = nodes[j].ibranch;
            if (j != tree.root) varb[j] = (x[ib] > 1e-8 ? -1 / varb_minbranches[ib] : 999);
         }
         /*
          matout(F0, x, 1, com.ntime);
          matout2(F0, varb, 1, tree.nnode, 10, 7);
          fout = stdout;
          exit(0);
          */
      }

      if (!com.fix_kappa) data.kappa[locus] = x[com.ntime];
      if (!com.fix_omega) data.omega[locus] = x[com.ntime + !com.fix_kappa];
      if (!com.fix_alpha) data.alpha[locus] = x[com.ntime + !com.fix_kappa + !com.fix_omega];

      lnLsum += lnL;

      t0 = nodes[son0].branch;
      t1 = nodes[son1].branch;
      varb[tree.root] = varb[t0 > t1 ? son0 : son1];
      nodes[son0].branch = nodes[son1].branch = (t0 + t1) / 2;  /* arbitrary */
      mr[locus] = GetMeanRate();

      printf("   Locus %d: %d sequences, %d blengths, lnL = %15.6f mr=%.5f%10s\n",
         locus + 1, com.ns, com.np - 1, -lnL, mr[locus], printtime(timestr));
      fprintf(fout, "\nlnL = %.6f\n\n", -lnL);
      OutTreeB(fout);  fprintf(fout, "\n");
      for (i = 0; i < com.np; i++) fprintf(fout, " %8.5f", x[i]);
      fprintf(fout, "\n");
      for (i = 0; i < tree.nbranch; i++) fprintf(fout, " %8.5f", sqrt(varb[tree.branches[i][1]]));
      fprintf(fout, "\n\n");
      OutTreeN(fout, 1, 1);
      fprintf(fout, "\n");
      fflush(fout);

      fprintf(fBV, "\n\nLocus %d: %d sequences, %d+1 branches\nlnL = %15.6f\n\n",
         locus + 1, com.ns, tree.nbranch - 1, -lnL);
      OutTreeB(fBV);  fprintf(fBV, "\n");
      for (i = 0; i < tree.nbranch; i++) fprintf(fBV, " %12.9f", x[i]);
      fprintf(fBV, "\n");
      for (i = 0; i < tree.nbranch; i++) fprintf(fBV, " %12.9f", sqrt(varb[tree.branches[i][1]]));
      fprintf(fBV, "\n\n");  OutTreeN(fBV, 1, 1);  fprintf(fBV, "\n");
      fflush(fBV);
   }
   fclose(fBV);
   if (data.ngene > 1) fprintf(fout, "\nSum of lnL over loci = %15.6f\n", -lnLsum);

   /* Step 2: ad hoc rate smoothing to estimate branch rates.  */
   printf("\nStep 2: Ad hoc rate smoothing to estimate branch rates.\n");
   fprintf(fout, "\n\nStep 2: Ad hoc rate smoothing to estimate branch rates.\n");
   /* s - 1 - NFossils node ages, (2*s_i - 2) rates for branches at each locus */
   com.clock = 1;
   nodes = nodes_t; 
   copy_from_Sptree(&stree, nodes);
   GetInitialsTimes(x);

   for (locus = 0, com.np = com.ntime - 1; locus < data.ngene; locus++)
      com.np += data.ns[locus] * 2 - 2;
   if (data.fix_nu == 2) com.np++;
   if (data.fix_nu == 3) com.np += data.ngene;

   if (com.np > NS * 6) error2("change NP for ad hoc rate smoothing.");
   for (i = 0; i < com.ntime - 1; i++)
   {
      xb[i][0] = pb;  xb[i][1] = 1 - pb;
   }
   if (!nodes[tree.root].fossil)
   {
      xb[0][0] = AgeLow[tree.root] * 1.0001; xb[0][1] = max2(AgeLow[tree.root] * 10, 50);
   }
   for (; i < com.np; i++) { /* for rates */
      xb[i][0] = rb[0]; xb[i][1] = rb[1];
   }
   for (locus = 0, i = com.ntime - 1; locus < data.ngene; locus++)
      for (j = 0; j < data.ns[locus] * 2 - 2; j++)
         x[i++] = mr[locus] * (.8 + .4 * rndu());
   for (; i < com.np; i++)   /* nu */
      x[i] = 0.001 + 0.1 * rndu();

   if (noisy > 3) {
      for (i = 0; i < com.np; i++) {
         printf(" %10.5f", x[i]);
         if (i == com.ntime - 2) printf("\n");
      }
      printf("\n");
      if (com.np < 200) {
         for (i = 0; i < com.np; i++)  printf(" %10.5f", xb[i][0]);
         printf("\n");
         for (i = 0; i < com.np; i++)  printf(" %10.5f", xb[i][1]);
         printf("\n");
      }
   }

   if (data.fix_nu > 1)
      pnu = x + com.np - (data.fix_nu == 2 ? 1 : data.ngene);
   printf("  %d times, %d rates, %d parameters, ", com.ntime - 1, k, com.np);

   noisy = 3;
   f = funSS_AHRS(x, com.np);
   if (noisy > 2) printf("\nf0 = %12.6f\n", f);

   if (finStep2) {
      puts("read MLEs from step 2 from file");
      for (i = 0; i < com.np; i++) fscanf(finStep2, "%lf", &x[i]);
      matout(F0, x, 1, com.np);
   }
   else {
      j = ming2(frub, &f, funSS_AHRS, NULL, x, xb, space, 1e-9, com.np);

      /* generate output to in.clockStep2
       matout(fout,x,1,com.np);
       */

      if (j == -1)
      {
         puts("\nad hoc rate smoothing iteration may not have converged.\nEnter to continue; Ctrl-C to break.");
         getchar();
      }
   }
   free(varb_AHRS);

   fputs("\nEstimated divergence times from ad hoc rate smoothing\n\n", fout);
   nodes = nodes_t; 
   copy_from_Sptree(&stree, nodes);
   for (i = 0; i < tree.nnode; i++) nodes[i].branch *= 100;
   for (i = com.ns; i < tree.nnode; i++)
      fprintf(fout, "Node %2d   Time %9.5f\n", i + 1, nodes[i].age * 100);
   fprintf(fout, "\n");  OutTreeN(fout, 1, 1);  fprintf(fout, "\n");

   fprintf(fout, "\nEstimated rates from ad hoc rate smoothing\n");
   for (locus = 0, k = k0 = com.ntime - 1; locus < data.ngene; k0 += data.nbrate[locus++]) {

      UseLocus(locus, -1, 0, 1);
      for (i = 0; i < tree.nnode; i++)
         if (i != tree.root)  nodes[i].label = x[k++];
      son0 = nodes[tree.root].sons[0]; son1 = nodes[tree.root].sons[1];
      t0 = nodes[tree.root].age - nodes[son0].age;
      t1 = nodes[tree.root].age - nodes[son1].age;
      nodes[tree.root].label = (nodes[son0].label * t1 + nodes[son1].label * t0) / (t0 + t1);
      SetBranchRates(tree.root);  /* node rates -> branch rates */

      nu = (data.fix_nu == 3 ? *(pnu + locus) : *pnu);
      fprintf(fout, "\nLocus %d (%d sequences)\n\n", locus + 1, com.ns);
      fprintf(fout, "nu = %.6g\n", nu);

      /* this block can be deleted? */
      fprintf(fout, "\nnode \tage \tlength \trate\n");
      for (i = 0; i < tree.nnode; i++) {
         fprintf(fout, "%02d\t%.3f", i + 1, nodes[i].age);
         if (i != tree.root)
            fprintf(fout, "\t%.5f\t%.5f", nodes[i].branch, nodes[i].label);
         fprintf(fout, "\n");
      }

      fprintf(fout, "\nRates as labels in tree:\n");
      OutTreeN(fout, 1, PrLabel);
      fprintf(fout, "\n");
      fflush(fout);

      if (data.nbrate[locus] > maxnbrate) error2("too many rate classes?  Change source.");
      for (i = 0, minr = 1e6, maxr = 0; i < tree.nnode; i++)
         if (i != tree.root) {
            r = nodes[i].label;
            if (r < 0 && noisy)
               puts("node label<0?");
            minr = min2(minr, r);
            maxr = max2(maxr, r);
         }

      fprintf(fdist, "\n%6d\n", tree.nnode - 1);
      for (i = 0; i < tree.nnode; i++) {
         if (i == tree.root) continue;
         fprintf(fdist, "R%-10.7f  ", nodes[i].label);
         for (j = 0; j < i; j++)
            if (j != tree.root)
               fprintf(fdist, " %9.6f", fabs(nodes[i].label - nodes[j].label));
         fprintf(fdist, "\n");
      }
      fflush(fdist);
      /*
       for(j=0; j<data.nbrate[locus]; j++)
       Rj[j]=minr+(j+1)*(maxr-minr)/data.nbrate[locus];
       */
      beta = pow(1 / (data.nbrate[locus] + 1.), 1 / (data.nbrate[locus] - 1.));
      beta = 0.25 + 0.25 * log((double)data.nbrate[locus]);
      if (beta > 1) beta = 0.99;
      for (j = 0; j < data.nbrate[locus]; j++)
         Rj[j] = minr + (maxr - minr) * pow(beta, data.nbrate[locus] - 1. - j);

      printf("\nLocus %d: nu = %.6f, rate range (%.6f, %.6f)\n", locus + 1, nu, minr, maxr);
      printf("Cutting points:\n");
      for (j = 0; j < data.nbrate[locus]; j++)
         printf(" < %.6f, ", Rj[j]);
      printf("\nThe number of rate groups (0 for no change)? ");
      /* scanf("%d", &j); */
      j = 0;
      if (j) {
         data.nbrate[locus] = j;
         printf("input %d cutting points? ", data.nbrate[locus] - 1);
         for (j = 0, Rj[data.nbrate[locus] - 1] = maxr; j < data.nbrate[locus] - 1; j++)
            scanf("%lf", &Rj[j]);
      }

      for (i = 0; i < data.nbrate[locus]; i++) { mbrate[i] = 0; nbrate[i] = 0; }
      for (i = 0; i < tree.nnode; i++) {
         if (i == tree.root) continue;
         r = nodes[i].label;
         for (j = 0; j < data.nbrate[locus] - 1; j++)
            if (r < Rj[j]) break;
         mbrate[j] += r;
         nbrate[j] ++;
         nodes[i].label = j;
      }
      nodes[tree.root].label = -1;
      for (i = 0; i < data.nbrate[locus]; i++)
         mbrate[i] = (nbrate[i] ? mbrate[i] / nbrate[i] : -1);

      fprintf(fout, "\nCollapsing rates into groups\nRate range: (%.6f, %.6f)\n", minr, maxr);
      /*      fprintf(fout,"\nCollapsing rates into groups\nbeta = %.6g  Rate range: (%.6f, %.6f)\n", beta, minr,maxr);
       */
      for (j = 0; j < data.nbrate[locus]; j++)
         fprintf(fout, "rate group %d  (%2d): <%9.6f, mean %9.6f\n",
            j, nbrate[j], Rj[j], mbrate[j]);

      fprintf(fout, "\n");
      OutTreeN(fout, 1, PrLabel);
      fprintf(fout, "\n\nRough rates for branch groups at locus %d\n", locus + 1);
      for (i = 0; i < data.nbrate[locus]; i++)
         x[k0 + i] = mbrate[i];
   }

   printf("\n\n%d times, %d timerates from AHRS:\n", com.ntime - 1, k0);
   fprintf(fout, "\n\n%d times, %d timerates from AHRS\n", com.ntime - 1, k0);
   for (i = 0; i < k0; i++) {
      printf("%12.6f", x[i]);
      if (i == com.ntime - 2)  printf("\n");
      fprintf(fout, "%12.6f", x[i]);
      if (i == com.ntime - 2)  fprintf(fout, "\n");
   }
   printf("\n");
   fprintf(fout, "\n");

   for (i = 0; i < k0; i++) x[i] *= 0.9 + 0.2 * rndu();

   com.model = model0;  com.clock = 6;


   com.fix_kappa = fix_kappa0; com.kappa = kappa0;
   com.fix_omega = fix_omega0; com.omega = omega0;
   com.fix_alpha = fix_alpha0; com.alpha = alpha0;

#if 0
   /* fix parameters: value > 0, precise value unimportant */
   if (!fix_kappa0) { com.fix_kappa = 1; com.kappa = 0.1; }
   if (!fix_omega0) { com.fix_omega = 1; com.omega = 0.1; }
   if (!fix_alpha0) { com.fix_alpha = 1; com.alpha = 0.1; }
#endif

   fclose(fdist);
   fflush(fout);
   printf(" %10s\n", printtime(timestr));

   if (finStep1) fclose(finStep1);
   if (finStep2) fclose(finStep2);

   return(0);
}


void DatingHeteroData(FILE* fout)
{
   /* This is for clock and local-clock dating using heterogeneous data from
    multiple loci.  Some species might be missing at some loci.  Thus
    gnodes[locus] stores the gene tree at locus.  Branch lengths in the gene
    tree are constructed using the divergence times in the main species tree,
    and the rates for genes and branches.

    com.clock = 5: global clock
    6: local clock
    */
   char timestr[64];
   int i, j, k, s, np, locus;
   double x[NS * 6], xb[NS * 6][2], lnL, e = 1e-7, * var = NULL;
   int nbrate = 4;
   size_t maxnpML, maxnpADRS;

   data.fix_nu = 3;
   /*
    if(com.clock==6) {
    printf("nu (1:fix; 2:estimate one for all genes; 3:estimate one for every gene)? ");
    scanf("%d", &data.fix_nu);
    if(data.fix_nu==1) scanf("%lf", &nu_AHRS);
    }
    */
   ReadTreeSeqs(fout);
   com.nbtype = 1;
   for (j = 0; j < stree.nnode; j++) {
      stree.nodes[j].pfossil[0] = stree.nodes[j].pfossil[1] = -1;
   }
   for (j = stree.nspecies, com.ntime = j - 1, stree.nfossil = 0; j < stree.nnode; j++) {
      if (stree.nodes[j].fossil) {
         com.ntime--;
         stree.nfossil++;
         printf("node %2d age fixed at %.3f\n", j, stree.nodes[j].age);
      }
   }
   GetMemBC();
   s = stree.nspecies;
   maxnpML = s - 1 + (5 + 2) * data.ngene;
   maxnpADRS = s - 1 + (2 * s - 1) * data.ngene + 2 * data.ngene;
   com.sspace = max2(com.sspace, spaceming2(maxnpADRS));
   com.sspace = max2(com.sspace, maxnpML * (maxnpML + 1) * sizeof(double));
   if ((com.space = (double*)realloc(com.space, com.sspace)) == NULL)
      error2("oom space");

#if (defined CODEML)
   GetUVRoot_codeml();
#endif
   if (com.clock == 6) {
      if (data.fix_nu <= 1) {
         printf("nu & nbrate? ");
         scanf("%lf%d? ", &nu_AHRS, &nbrate);
      }
      for (locus = 0; locus < data.ngene; locus++)
         data.nbrate[locus] = nbrate;
      AdHocRateSmoothing(fout, x, xb, com.space);

      printf("\nStep 3: ML estimation of times and rates.");
      fprintf(fout, "\n\nStep 3: ML estimation of times and rates.\n");
   }
   else {   /* clock = 5, global clock */
      for (locus = 0; locus < data.ngene; locus++)
         for (i = 0, data.nbrate[locus] = 1; i < data.ns[locus] * 2 - 1; i++)
            gnodes[locus][i].label = 0;
   }

   noisy = 3;

   nodes = nodes_t; 
   copy_from_Sptree(&stree, nodes);
   GetInitialsClock56Step3(x);
   np = com.np;

   SetxBound(com.np, xb);
   lnL = lnLfunHeteroData(x, np);

   if (noisy) {
      printf("\nntime & nrate & np:%6d%6d%6d\n", com.ntime - 1, com.nrate, com.np);
      matout(F0, x, 1, np);
      printf("\nlnL0 = %12.6f\n", -lnL);
   }

   j = ming2(noisy > 2 ? frub : NULL, &lnL, lnLfunHeteroData, NULL, x, xb, com.space, e, np);

   if (noisy) printf("Out...\nlnL  = %12.6f\n", -lnL);

   LASTROUND = 1;
   for (i = 0, j = !stree.nodes[stree.root].fossil; i < stree.nnode; i++)
      if (i != stree.root && stree.nodes[i].nson && !stree.nodes[i].fossil)
         x[j++] = stree.nodes[i].age;       /* copy node ages into x[] */

   if (com.getSE) {
      if (np > 100 || (com.seqtype && np > 20)) puts("Calculating SE's");
      var = com.space + np;
      Hessian(np, x, lnL, com.space, var, lnLfunHeteroData, var + np * np);
      matinv(var, np, np, var + np * np);
   }
   nodes = nodes_t; 
   copy_from_Sptree(&stree, nodes);
   SetBranch(x);
   fprintf(fout, "\n\nTree:  ");  OutTreeN(fout, 0, 0);
   fprintf(fout, "\nlnL(ntime:%3d  np:%3d):%14.6f\n", com.ntime - 1, np, -lnL);
   OutTreeB(fout);   fprintf(fout, "\n");
   for (i = 0; i < np; i++) fprintf(fout, " %9.5f", x[i]);
   fprintf(fout, "\n");
   fflush(fout);

   if (com.getSE) {
      fprintf(fout, "SEs for parameters:\n");
      for (i = 0; i < np; i++) fprintf(fout, " %9.5f", (var[i * np + i] > 0. ? sqrt(var[i * np + i]) : -1));
      fprintf(fout, "\n");
      if (com.getSE == 2) matout2(fout, var, np, np, 15, 10);
   }

   fprintf(fout, "\nTree with node ages for TreeView\n");
   for (i = 0; i < tree.nnode; i++)  nodes[i].branch *= 100;
   fprintf(fout, "\n"); OutTreeN(fout, 1, 1);  fprintf(fout, "\n");
   fprintf(fout, "\n"); OutTreeN(fout, 1, PrNodeNum);  fprintf(fout, "\n");
   fprintf(fout, "\n"); OutTreeN(fout, 1, PrLabel | PrAge);  fprintf(fout, "\n");
   fprintf(fout, "\n"); OutTreeN(fout, 1, 0);  fprintf(fout, "\n");
   OutputTimesRates(fout, x, var);

   fprintf(fout, "\nSubstititon rates for genes (per time unit)\n");
   for (j = 0, k = com.ntime - 1; j < data.ngene; j++) {
      fprintf(fout, "   Gene %2d: ", j + 1);
      for (i = 0; i < data.nbrate[j]; i++, k++) {
         fprintf(fout, "%10.5f", x[k]);
         if (com.getSE) fprintf(fout, " +- %.5f", sqrt(var[k * np + k]));
      }
      if (com.clock == 6) fprintf(fout, " ");
      fprintf(fout, "\n");
   }
   if (!com.fix_kappa) {
      fprintf(fout, "\nkappa for genes\n");
      for (j = 0; j < data.ngene; j++, k++) {
         fprintf(fout, "%10.5f", data.kappa[j]);
         if (com.getSE) fprintf(fout, " +- %.5f", sqrt(var[k * np + k]));
      }
   }
   if (!com.fix_omega) {
      fprintf(fout, "\nomega for genes\n");
      for (j = 0; j < data.ngene; j++, k++) {
         fprintf(fout, "%10.5f", data.omega[j]);
         if (com.getSE) fprintf(fout, " +- %.5f", sqrt(var[k * np + k]));
      }
   }
   if (!com.fix_alpha) {
      fprintf(fout, "\nalpha for genes\n");
      for (j = 0; j < data.ngene; j++, k++) {
         fprintf(fout, "%10.5f", data.alpha[j]);
         if (com.getSE) fprintf(fout, " +- %.5f", sqrt(var[k * np + k]));
      }
   }
   fprintf(fout, "\n");
   FreeMemBC();
   printf("\nTime used: %s\n", printtime(timestr));
   exit(0);
}

#endif
