/*
    Copyright (C) 2015 Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#if defined(_MSC_VER)
static __declspec(align(32)) double lmatrix[4*4];
static __declspec(align(32)) double rmatrix[4*4];
#else
static double lmatrix[4*4] __attribute__ ((aligned (32)));
static double rmatrix[4*4] __attribute__ ((aligned (32)));
#endif

void pmat_jc69(double * pmatrix, double t)
{
  if (t < -0.0001)
    printf ("\nt = %.5f in pijJC69", t);
  
  if (t < 1e-100)
  {
    pmatrix[0]  = 1;
    pmatrix[1]  = 0;
    pmatrix[2]  = 0;
    pmatrix[3]  = 0;

    pmatrix[4]  = 0;
    pmatrix[5]  = 1;
    pmatrix[6]  = 0;
    pmatrix[7]  = 0;

    pmatrix[8]  = 0;
    pmatrix[9]  = 0;
    pmatrix[10] = 1;
    pmatrix[11] = 0;

    pmatrix[12] = 0;
    pmatrix[13] = 0;
    pmatrix[14] = 0;
    pmatrix[15] = 1;
  }
  else
  {
    double a =  (1 + 3*exp(-4*t/3) ) / 4;
    double b = (1 - a) / 3;

    pmatrix[0]  = a;
    pmatrix[1]  = b;
    pmatrix[2]  = b;
    pmatrix[3]  = b;

    pmatrix[4]  = b;
    pmatrix[5]  = a;
    pmatrix[6]  = b;
    pmatrix[7]  = b;

    pmatrix[8]  = b;
    pmatrix[9]  = b;
    pmatrix[10] = a;
    pmatrix[11] = b;

    pmatrix[12] = b;
    pmatrix[13] = b;
    pmatrix[14] = b;
    pmatrix[15] = a;
  }
}

static void tip_tip(int lchild, int rchild, double * clv)
{
  int sites = com.npatt;
  int states = 4;
  int n,i;

  __m256d ymm0,ymm1,ymm2;

  pmat_jc69(lmatrix, nodes[lchild].branch);
  pmat_jc69(rmatrix, nodes[rchild].branch);

  for (n = 0; n < sites; ++n)
  {
    double * lmat = lmatrix + com.z[lchild][n] * states;
    double * rmat = rmatrix + com.z[rchild][n] * states;

    // for (i = 0; i < states; ++i)
    //   clv[i] = lmat[i] * rmat[i];
    
    ymm0 = _mm256_load_pd(lmat);
    ymm1 = _mm256_load_pd(rmat);
    ymm2 = _mm256_mul_pd(ymm0,ymm1);

    _mm256_store_pd(clv,ymm2);

    clv += states;
  }
}

static void tip_tip_amb(int lchild, int rchild, double * clv)
{
  int sites = com.npatt;
  int states = 4;
  int n,k;

  double lvec[4];
  double rvec[4];

  pmat_jc69(lmatrix, nodes[lchild].branch);
  pmat_jc69(rmatrix, nodes[rchild].branch);
  
  __m256d xmm0,xmm1,ymm0,ymm1;

  for (n = 0; n < sites; ++n)
  {
    double * lmat = lmatrix + (CharaMap[com.z[lchild][n]][0] << 2); // *states;
    double * rmat = rmatrix + (CharaMap[com.z[rchild][n]][0] << 2); // *states;
    
    xmm0 = _mm256_load_pd(lmat);
    for (k = 1; k < nChara[com.z[lchild][n]]; ++k)
    {
      lmat = lmatrix + (CharaMap[com.z[lchild][n]][k] << 2); // *states;
      xmm1 = _mm256_load_pd(lmat);
      xmm0 = _mm256_add_pd(xmm0,xmm1);
    }

    ymm0 = _mm256_load_pd(rmat);
    for (k = 1; k < nChara[com.z[rchild][n]]; ++k)
    {
      rmat = rmatrix + (CharaMap[com.z[rchild][n]][k] << 2); // *states;
      ymm1 = _mm256_load_pd(rmat);
      ymm0 = _mm256_add_pd(ymm0,ymm1);
    }

    xmm1 = _mm256_mul_pd(xmm0,ymm0);
    _mm256_store_pd(clv,xmm1);

    clv += states;
  }
}

static void tip_inner(int lchild, int rchild, double * clv)
{
  int sites = com.npatt;
  int states = 4;
  int h;

  pmat_jc69(lmatrix, nodes[lchild].branch);
  pmat_jc69(rmatrix, nodes[rchild].branch);

  double * rclv = nodes[rchild].conP;
  
  __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
  __m256d xmm0,xmm4;

  for (h = 0; h < sites; ++h)
  {
    double * lmat = lmatrix + com.z[lchild][h] * states;
    double * rmat = rmatrix;

    xmm4 = _mm256_load_pd(lmat);

    ymm4 = _mm256_load_pd(rmat);
    ymm5 = _mm256_load_pd(rclv);
    ymm0 = _mm256_mul_pd(ymm4,ymm5);

    rmat += states;

    ymm4 = _mm256_load_pd(rmat);
    ymm1 = _mm256_mul_pd(ymm4,ymm5);

    rmat += states;

    ymm4 = _mm256_load_pd(rmat);
    ymm2 = _mm256_mul_pd(ymm4,ymm5);

    rmat += states;

    ymm4 = _mm256_load_pd(rmat);
    ymm3 = _mm256_mul_pd(ymm4,ymm5);

    /* compute y */
    ymm4 = _mm256_unpackhi_pd(ymm0,ymm1);
    ymm5 = _mm256_unpacklo_pd(ymm0,ymm1);

    ymm6 = _mm256_unpackhi_pd(ymm2,ymm3);
    ymm7 = _mm256_unpacklo_pd(ymm2,ymm3);

    ymm0 = _mm256_add_pd(ymm4,ymm5);
    ymm1 = _mm256_add_pd(ymm6,ymm7);

    ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
    ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
    ymm4 = _mm256_add_pd(ymm2,ymm3);

    /* compute x*y */
    xmm0 = _mm256_mul_pd(xmm4,ymm4);

    _mm256_store_pd(clv, xmm0);

    clv += states;
    rclv += states;
  }
}

static void tip_inner_amb(int lchild, int rchild, double * clv)
{
  int sites = com.npatt;
  int states = 4;
  int h,j,k;
  double y;
  double x;

  pmat_jc69(lmatrix, nodes[lchild].branch);
  pmat_jc69(rmatrix, nodes[rchild].branch);

  __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
  __m256d xmm0,xmm1;

  double * rclv = nodes[rchild].conP;

  for (h = 0; h < sites; ++h)
  {
    double * lmat = lmatrix + (CharaMap[com.z[lchild][h]][0] << 2); // *states;
    double * rmat = rmatrix;

    /* compute term for left child (tip) */
    xmm0 = _mm256_load_pd(lmat);
    for (k = 1; k < nChara[com.z[lchild][h]]; ++k)
    {
      lmat = lmatrix + (CharaMap[com.z[lchild][h]][k] << 2); // *states;
      xmm1 = _mm256_load_pd(lmat);
      xmm0 = _mm256_add_pd(xmm0,xmm1);
    }

    /* compute term for right child (inner) */
    ymm4 = _mm256_load_pd(rmat);
    ymm5 = _mm256_load_pd(rclv);
    ymm0 = _mm256_mul_pd(ymm4,ymm5);

    rmat += states;

    ymm4 = _mm256_load_pd(rmat);
    ymm1 = _mm256_mul_pd(ymm4,ymm5);

    rmat += states;

    ymm4 = _mm256_load_pd(rmat);
    ymm2 = _mm256_mul_pd(ymm4,ymm5);

    rmat += states;

    ymm4 = _mm256_load_pd(rmat);
    ymm3 = _mm256_mul_pd(ymm4,ymm5);

    ymm4 = _mm256_unpackhi_pd(ymm0,ymm1);
    ymm5 = _mm256_unpacklo_pd(ymm0,ymm1);

    ymm6 = _mm256_unpackhi_pd(ymm2,ymm3);
    ymm7 = _mm256_unpacklo_pd(ymm2,ymm3);

    ymm0 = _mm256_add_pd(ymm4,ymm5);
    ymm1 = _mm256_add_pd(ymm6,ymm7);

    ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
    ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
    ymm4 = _mm256_add_pd(ymm2,ymm3);

    ymm0 = _mm256_mul_pd(xmm0,ymm4);

    _mm256_store_pd(clv,ymm0);

    clv += states;
    rclv += states;
  }
  
}

static void inner_inner(int lchild, int rchild, double * clv)
{
  int sites = com.npatt;
  int states = 4;
  int h;

  __m256d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7;
  __m256d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7;

  pmat_jc69(lmatrix, nodes[lchild].branch);
  pmat_jc69(rmatrix, nodes[rchild].branch);

  double * lclv = nodes[lchild].conP;
  double * rclv = nodes[rchild].conP;

  for (h = 0; h < sites; ++h)
  {
    double * lmat = lmatrix;
    double * rmat = rmatrix;

    /* compute vector of x */
    xmm4 = _mm256_load_pd(lmat);
    xmm5 = _mm256_load_pd(lclv);
    xmm0 = _mm256_mul_pd(xmm4,xmm5);

    ymm4 = _mm256_load_pd(rmat);
    ymm5 = _mm256_load_pd(rclv);
    ymm0 = _mm256_mul_pd(ymm4,ymm5);

    lmat += states;
    rmat += states;

    xmm4 = _mm256_load_pd(lmat);
    xmm1 = _mm256_mul_pd(xmm4,xmm5);

    ymm4 = _mm256_load_pd(rmat);
    ymm1 = _mm256_mul_pd(ymm4,ymm5);

    lmat += states;
    rmat += states;

    xmm4 = _mm256_load_pd(lmat);
    xmm2 = _mm256_mul_pd(xmm4,xmm5);

    ymm4 = _mm256_load_pd(rmat);
    ymm2 = _mm256_mul_pd(ymm4,ymm5);

    lmat += states;
    rmat += states;

    xmm4 = _mm256_load_pd(lmat);
    xmm3 = _mm256_mul_pd(xmm4,xmm5);

    ymm4 = _mm256_load_pd(rmat);
    ymm3 = _mm256_mul_pd(ymm4,ymm5);

    /* compute x */
    xmm4 = _mm256_unpackhi_pd(xmm0,xmm1);
    xmm5 = _mm256_unpacklo_pd(xmm0,xmm1);

    xmm6 = _mm256_unpackhi_pd(xmm2,xmm3);
    xmm7 = _mm256_unpacklo_pd(xmm2,xmm3);

    xmm0 = _mm256_add_pd(xmm4,xmm5);
    xmm1 = _mm256_add_pd(xmm6,xmm7);

    xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));
    xmm3 = _mm256_blend_pd(xmm0,xmm1,12);
    xmm4 = _mm256_add_pd(xmm2,xmm3);

    /* compute y */
    ymm4 = _mm256_unpackhi_pd(ymm0,ymm1);
    ymm5 = _mm256_unpacklo_pd(ymm0,ymm1);

    ymm6 = _mm256_unpackhi_pd(ymm2,ymm3);
    ymm7 = _mm256_unpacklo_pd(ymm2,ymm3);

    ymm0 = _mm256_add_pd(ymm4,ymm5);
    ymm1 = _mm256_add_pd(ymm6,ymm7);

    ymm2 = _mm256_permute2f128_pd(ymm0,ymm1, _MM_SHUFFLE(0,2,0,1));
    ymm3 = _mm256_blend_pd(ymm0,ymm1,12);
    ymm4 = _mm256_add_pd(ymm2,ymm3);

    /* compute x*y */
    xmm0 = _mm256_mul_pd(xmm4,ymm4);

    _mm256_store_pd(clv, xmm0);

    clv += states;
    lclv += states;
    rclv += states;
  }
}

int ConditonalPNode (int inode)
{
  int states = 4;
  int sites = com.npatt;
  int i;
  int child, lchild, rchild;

  /* recursive call the ConditionalPNode on all children of the current node */
  for (i = 0; i < nodes[inode].nson; ++i)
  {
    /* get id of i-th child of current node */
    child = nodes[inode].sons[i];

    /* if child is an inner node then recursively compute its conditional
       probabilities vector */
    if (nodes[child].nson > 0 && (!mcmc.saveconP || !com.oldconP[child]))
      ConditonalPNode(nodes[inode].sons[i]);
  }

  /* initialize CLV entries of current node to 1 */
  int n = sites * states;
  for (i = 0; i < n; ++i)
    nodes[inode].conP[i] = 1;

  if (nodes[inode].nson == 0) return (0);

  lchild = nodes[inode].sons[0];
  rchild = nodes[inode].sons[1];

  int ltip = (nodes[lchild].nson == 0);
  int rtip = (nodes[rchild].nson == 0);

  if (ltip && rtip)
  {
    if (com.cleandata)
      tip_tip(lchild,rchild,nodes[inode].conP);
    else
      tip_tip_amb(lchild,rchild,nodes[inode].conP);
  }
  else if (ltip)
  {
    if (com.cleandata)
      tip_inner(lchild,rchild,nodes[inode].conP);
    else
      tip_inner_amb(lchild,rchild,nodes[inode].conP);
  }
  else if (rtip)
  {
    if (com.cleandata)
      tip_inner(rchild,lchild,nodes[inode].conP);
    else
      tip_inner_amb(rchild,lchild,nodes[inode].conP);
  }
  else
    inner_inner(lchild,rchild,nodes[inode].conP);

  return (0);
}
