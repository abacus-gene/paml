/* treespace.c
   collection of tree perturbation routines
*/

#include"tools.h"

int MakeTreeIb (int ns, int Ib[], int rooted)
{
/* construct tree from Ib[] using the algorithm of adding species
   Ib[k] marks the branch to which the (k+3)_th species (or the root) 
   is added.  Ib[k] should be in the range [0,k+3]
*/
   int i,j,k, is,it;

   tree.nbranch=3;
   for (i=0; i<3; i++)  { tree.branches[i][0]=3;  tree.branches[i][1]=i; }
   for (k=0; k<ns-3; k++) {
      is=k+3;       /* add species (is) to branch Ib[k] */

      for (i=0; i<tree.nbranch; i++)  for (j=0; j<2; j++)
         if (tree.branches[i][j]>=is) tree.branches[i][j]+=2;
      it=tree.branches[Ib[k]][1];
      tree.branches[Ib[k]][1]=is+1;
      tree.branches[tree.nbranch][0]=is+1;
      tree.branches[tree.nbranch++][1]=it;
      tree.branches[tree.nbranch][0]=is+1;
      tree.branches[tree.nbranch++][1]=is;
   }
   tree.origin=tree.branches[0][0];
   BranchToNode ();
   
   if (rooted) {
      it=tree.branches[Ib[k]][0];
      tree.branches[Ib[k]][0]=tree.branches[tree.nbranch][0]=2*com.ns-2;
      tree.branches[tree.nbranch][1]=it;
      for (; it!=tree.origin;  it=nodes[it].father) {
         tree.branches[nodes[it].ibranch][0]=it;
         tree.branches[nodes[it].ibranch][1]=nodes[it].father;
      }
      tree.origin=2*com.ns-2;  tree.nbranch++;
      BranchToNode ();
   }
   return (0);
}

int GetTreeI (int itree, int ns, int rooted)
{
/* get the i_th tree.  Trees are ordered according to the algorithm of 
   adding species.
   returns a random tree if itree==-1, in which case ns can be large
*/
   int i, M[NS-2], nM=ns-3+rooted, Ib[NS-2];

   for (i=0; i<nM-1; i++) M[i]=2*i+5;
   for (i=0,M[nM-1]=1; i<nM-2; i++) M[nM-1-i-2]*=M[nM-1-i-1];

   if (itree==-1)  for (i=0; i<nM; i++) Ib[i]=(int)((2*i+3)*rndu());
   else            for (i=0; i<nM; i++) {Ib[i]=itree/M[i]; itree%=M[i]; } 
/*
   if (noisy>3) {
      FOR (i, nM) printf ("%5d ", M[i]);   FPN (F0);
      FOR (i, nM) printf ("%5d ", Ib[i]);  FPN (F0);
   }
*/
   MakeTreeIb (ns, Ib, rooted);
   return (0);
}

int NumberofTrees (int ns, int rooted)
{
   int i, ntree=1;

   if (ns>15) error ("ns too large in NumberofTrees().");
   for (i=4; i<=ns; i++)  ntree*=2*i-5;
   if (rooted) ntree*=2*i-3;
   return (ntree);
}


int ListTrees (FILE* fout, int ns, int rooted)
{
/* list trees by adding species, works fine with large ns
*/
   int NTrees, NTreeRoot=3;
   int i, Ib[NS-2], ns1=ns+rooted, nM=ns1-3,finish;

   printf ("%20s%20s%20s\n", "Taxa", "Unrooted trees", "Rooted trees");
   for (i=4,NTrees=1; i<=com.ns; i++)  
      printf ("%20d%20d%20d\n", i, (NTrees*=2*i-5), (NTreeRoot*=2*i-3));
   fprintf (fout, "%10d %10d\n", com.ns, (!rooted?NTrees:NTreeRoot));

   for (i=0;i<nM;i++) Ib[i]=0;
   for (NTrees=0; ; ) {
      MakeTreeIb (ns, Ib, rooted);
      OutaTreeN (fout, 0, 0);
      fprintf (fout, " %7d\n", NTrees++);

      for (i=nM-1,Ib[nM-1]++,finish=0; i>=0; i--) {
         if (Ib[i]<2*i+3) break;
         if (i==0) { finish=1; break; }
         Ib[i]=0; Ib[i-1]++; 
      }
      if (finish) break;
   }
   return (0);
}

int GetIofTree (int rooted, int keeptree, double space[])
{
/* Get the index of tree.
   tree.branches[] are destroyed for reconstruction, 
   and some branches are reversed which may affect nodes[] also.
   Both are restored before return if keeptree=1.
   Works with binary trees only.
   bA[nbranch*(ns-2)]
*/
   int M[NS-2], nM=com.ns-3+rooted;
   int i,j,k,is,it, b=0,a1,a2,Ib[NS-1], nid=tree.nnode-com.ns;
   char ns2=com.ns-2, *bA=(char*)space;  /* bA[b*ns2+j]: ancestors on branch b */
   struct TREEB tree0=tree;

   if (tree.nnode-com.ns!=com.ns-1-!rooted) error ("GetIofTree");
   if (com.ns>15) error("ns too large in GetIofTree");

   /* find new root. 
      Ib[]: No. of times inner nodes are visited on paths 1-2, 1-3, 2-3 */
   for (i=0; i<nid; i++) Ib[i]=0;
   for (i=0; i<2; i++)  for (j=i+1; j<3; j++)  {
      for (a1=i; a1!=-1; a1=nodes[a1].father) {
         for (a2=j; a2!=-1; a2=nodes[a2].father)  if (a1==a2) break;
         if (a1==a2) break;
      }
      for (k=0,Ib[a1-com.ns]++; k<2; k++) {  /* a1 is ancestor of i and j */
         a2=k?nodes[i].father:nodes[j].father;
         for (; a2!=a1; a2=nodes[a2].father)  Ib[a2-com.ns]++; 
      }
   }
   /* reset root of tree at it */
   for (it=com.ns; it<com.ns+nid; it++) if (Ib[it-com.ns]==3) break;
   ReRootTree (it);
   for (i=0,tree.nbranch=3; i<3; i++)  {
      tree.branches[i][0]=tree.origin;  tree.branches[i][1]=i;
      for (it=nodes[i].father,k=0; it!=tree.origin; it=nodes[it].father)
         bA[i*ns2+k++]=it;
      bA[i*ns2+k]=0;
   }
   
   for (is=3; is<com.ns; is++) { /* add species (is) on branch b at node it */
      for (it=nodes[is].father; it!=tree.origin; it=nodes[it].father) {
    for (b=0; b<tree.nbranch; b++) if (strchr(bA+b*ns2,it)) break;
    if (b<tree.nbranch) break;
      }
      Ib[is-3]=b;
      tree.branches[tree.nbranch][0]=it;
      tree.branches[tree.nbranch++][1]=tree.branches[b][1];
      tree.branches[tree.nbranch][0]=it;
      tree.branches[tree.nbranch++][1]=is;
      tree.branches[b][1]=it;

      if (is==com.ns-1) break;
      for (i=0; i<3; i++) {  /* modify bA[] for the 3 affected branches */
         if (i) b=tree.nbranch-3+i;
         it=nodes[tree.branches[b][1]].father; 
         for (k=0; it!=tree.branches[b][0]; it=nodes[it].father) 
            bA[b*ns2+k++]=it;
         bA[b*ns2+k]=0;
      }
   }  /* for (is) */
   if (rooted) {
      a1=nodes[k=tree0.origin].sons[0];  a2=nodes[tree0.origin].sons[1];
      if (nodes[a1].father==k)      k=a1;
      else if (nodes[a2].father==k) k=a2;
      else error ("rooooot");
      for (b=0; b<tree.nbranch; b++) if (tree.branches[b][1]==k) break;
      Ib[nM-1]=b;
   }
   if (keeptree)  { tree=tree0;  BranchToNode (); }

   for (i=0; i<nM-1; i++) M[i]=2*i+5;
   for (i=0,M[nM-1]=1; i<nM-2; i++) M[nM-1-i-2]*=M[nM-1-i-1];
   for (i=0,it=0; i<nM; i++) it+=Ib[i]*M[i];
   return (it);
}

void ReRootTree (int newroot)
{
/* reroot tree at newroot.  oldroot forgotten
*/
   int a,b;

   if (newroot==tree.origin) return;
   for (b=newroot,a=nodes[b].father; b!=tree.origin; b=a,a=nodes[b].father) {
      tree.branches[nodes[b].ibranch][0]=b;
      tree.branches[nodes[b].ibranch][1]=a;
   }
   tree.origin=newroot;
   BranchToNode ();
}

int NeighborNNI (int i_tree)
{
/* get the i_tree'th neighboring tree of tree by the nearest neighbor 
   interchange (NNI), each tree has 2*(# internal branches) neighbors.
   works with rooted and unrooted binary trees.

   Gives the ip_th neighbor for interior branch ib. 
   Involved branches are a..b, a..c, b..d,
   with a..b to be the internal branch.
   swap c with d, with d to be the ip_th son of b
*/ 
   int i, a,b,c,d, ib=i_tree/2, ip=i_tree%2;

   if (tree.nbranch!=com.ns*2-2-(nodes[tree.origin].nson==3)) 
      error ("err NeighborNNI: multificating tree.");

   /* locate a,b,c,d */
   for (i=0,a=0; i<tree.nbranch; i++)
      if (tree.branches[i][1]>=com.ns && a++==ib) break;
   ib=i;
   a=tree.branches[ib][0];  b=tree.branches[ib][1];
   c=nodes[a].sons[0];      if(c==b) c=nodes[a].sons[1];
   d=nodes[b].sons[ip];

   /* swap nodes c and d */
   tree.branches[nodes[c].ibranch][1]=d;
   tree.branches[nodes[d].ibranch][1]=c;
   BranchToNode ();
   return (0);
}

int GetLHistoryI (int iLH)
{
/* Get the ILH_th labelled history.  This function is rather similar to 
   GetTreeI which returns the I_th rooted or unrooted tree topology.
   The labeled history is recorded in the numbering of nodes: 
   node # increases as the node gets older: 
   node d corresponds to time 2*ns-2-d; tree.origin=ns*2-2;
   t0=1 > t1 > t2 > ... > t[ns-2]
   k ranges from 0 to i(i-1)/2 and indexes the pair (s1 & s2, with s1<s2),
   out of i lineages, that coalesce.
*/
   int i,k, inode, nodea[NS], s1, s2, it;

   tree.nnode=com.ns*2-1;
   for (i=0; i<tree.nnode; i++)  {
      nodes[i].father=nodes[i].nson=-1;  FOR (k,com.ns) nodes[i].sons[k]=-1;
   }
   for (i=0,inode=com.ns; i<com.ns; i++) nodea[i]=i;
   for (i=com.ns,it=iLH; i>=2; i--)  {
      k=it%(i*(i-1)/2);  it/=(i*(i-1)/2); 
      s2=(sqrt(1.+8*k)-1)/2+1;  s1=k-s2*(s2-1)/2; /* s1<s2, k=s2*(s2-1)/2+s1 */

      if (s1>=s2 || s1<0)  printf("\nijk%6d%6d%6d", s1, s2, k);

      nodes[nodea[s1]].father=nodes[nodea[s2]].father=inode;
      nodes[inode].nson=2;
      nodes[inode].sons[0]=nodea[s1];
      nodes[inode].sons[1]=nodea[s2];
      nodea[s1]=inode;  nodea[s2]=nodea[i-1]; 
      inode++;
   }
   tree.origin=inode-1;
   NodeToBranch();
   return (0);
}

int GetIofLHistory (void)
{
/* Get the index of the labelled history (rooted tree with nodes ordered
   according to time).  
   Numbering of nodes: node # increases as the node gets older:
   node d corresponds to time 2*ns-2-d; tree.origin=ns*2-2;
   t0=1 > t1 > t2 > ... > t[ns-2]
*/
   int index, i,j,k[NS+1], inode,nnode, nodea[NS], s[2];

   if (nodes[tree.origin].nson!=2 || tree.nnode!=com.ns*2-1
     || tree.origin!=com.ns*2-2)  error("IofLH");
   for (i=0; i<com.ns; i++) nodea[i]=i;
   for (inode=nnode=com.ns,index=0; inode<com.ns*2-1; inode++,nnode--) {
      FOR (i,2) FOR (j,nnode)  
         if (nodes[inode].sons[i]==nodea[j]) { s[i]=j; break; }
      k[nnode]=max(s[0],s[1]); s[0]=min(s[0],s[1]); s[1]=k[nnode];
      k[nnode]=s[1]*(s[1]-1)/2+s[0];
      nodea[s[0]]=inode; nodea[s[1]]=nodea[nnode-1];
   }
   for (nnode=2,index=0; nnode<=com.ns; nnode++)
      index=nnode*(nnode-1)/2*index+k[nnode];
   return (index);
}

int CountLHistory(char LHistories[], double space[])
{
/* counts the number of labeled histories corresponding to the same
   rooted tree topology.
   AA[level][] records all possible node numbers at level level and 
   depends on the value of AA[level-1][iA[level-1]]. 
   nA[] is the number of possibilities at each level, and iA[] indexes 
   node numbers along a path.  When level=nlevel-1, a path is formed.
   ipn is the node selected in the previous level and its sons are added 
   to the possilities for current level.
   LHistories[][com.ns-1] stores the orderings for each labeled history.
*/
   int level,nlevel, i,j;
   int *AA=(int*)space, ns1=com.ns-1, nA[NS-1], iA[NS-1], nLH, ipn;
   
   nlevel=tree.nnode-com.ns;
   if (com.ns-1!=tree.nnode-com.ns)  error ("binary tree?");
   FOR (i,com.ns-1) iA[i]=0;
/*
printf ("nlevel %2d\troot %2d\n", nlevel, tree.origin+1);
*/
   nA[0]=1;  ipn=AA[0*ns1+0]=tree.origin;
   for (level=1,nLH=0; ; level++) {
      nA[level]=0;
       /* inherate possiblities of level-1 */
      for (i=0,ipn=AA[(level-1)*ns1+iA[level-1]]; i<nA[level-1]; i++)
         if (i!=iA[level-1])  AA[level*ns1+nA[level]++]=AA[(level-1)*ns1+i];
      for (i=0; i<nodes[ipn].nson; i++) /* add new possibilities to level */
         { j=nodes[ipn].sons[i];  if (j>com.ns-1) AA[level*ns1+nA[level]++]=j; }

      if (level==nlevel-1)  {           /* last level */
         FOR (iA[nlevel-1], nA[nlevel-1])  {
            FOR (j, nlevel)  LHistories[nLH*(com.ns-1)+j]=AA[j*ns1+iA[j]];
            nLH++;
/*
printf ("\n%30s %3d: ", "path", nLH);
FOR (j, nlevel)  printf("%4d",AA[j*ns1+iA[j]]+1);
*/
         }
         for (level=nlevel-2; level>=0; level--) {
            if (++iA[level]<nA[level]) break;
            iA[level]=0;
         }
         if (level==-1) break;
      }
/*
printf("\nlevel %2d (ipn=%2d)(%2d/%2d):",level+1,ipn+1,iA[level]+1,nA[level]);
for (i=0; i<nA[level]; i++) printf ("%3d", AA[level*ns1+i]+1);
*/
   }
   return (nLH);
}

int ReorderNodes (char LHistory[])
{
/* changes interior node numbers so that the topology represent a labeled
   history.  LHistory recordes the order of interior nodes with [0] to be 
   root, and [ns-2] the youngest node.   
   Changes node LHistory[k] to com.ns*2-3-k.  
   Uses only tree but not nodes[] but transforms both tree and nodes into 
   the labeled history.
*/
   int i, j, k;

   if (tree.origin!=com.ns*2-2 || LHistory[0]!=com.ns*2-2) {
      tree.origin=com.ns*2-2;
/*      printf("\nRoot changed to %d in ReorderNodes..\n", com.ns*2-2+1); */
   }
   FOR (i, tree.nbranch) FOR (j,2) 
      if (tree.branches[i][j]>=com.ns) 
         FOR (k, com.ns-1)
            if (tree.branches[i][j]==LHistory[k]) 
               { tree.branches[i][j]=com.ns*2-2-k;  break; }
   BranchToNode();

   return (0);
}
