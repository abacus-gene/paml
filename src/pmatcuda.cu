/* CUDA backend for PMatUVRoot and ConditionalPNode. See pmatcuda.h. */

#include "pmatcuda.h"

#include <cuda_runtime.h>
#include <stdio.h>
#include <string.h>

/* Device buffers, grown as needed and reused across calls. */
static double *d_P = NULL, *d_exptm1 = NULL, *d_U = NULL, *d_V = NULL;
static size_t cap_P = 0, cap_exptm1 = 0, cap_U = 0, cap_V = 0;
static int *d_slot = NULL;          /* batch slot -> node */
static size_t cap_slot = 0;
static int checked = 0, available = 0;

/* Resident conditional-probability state. d_conP mirrors com.conP, so a node's
   block is found at the same offset the host uses, and codeml's per-site-class
   pointer shift needs no special handling here. */
static double *d_conP = NULL;
static size_t cap_conP = 0;
static char *d_z = NULL;            /* [ns][npatt] tip data */
static size_t cap_z = 0;
static char *d_nChara = NULL, *d_CharaMap = NULL;
static double *d_scaleMax = NULL;
static size_t cap_scaleMax = 0;
static ConPWork *d_work = NULL;     /* all (parent, son) pairs, grouped */
static size_t *d_dst = NULL;        /* interior node offsets */
static size_t *d_sdst = NULL;       /* scaled node offsets, grouped */
static int *d_slot2 = NULL;         /* scale slot per scaled node */
static size_t cap_work = 0, cap_dst = 0, cap_slot2 = 0, cap_sdst = 0;
static double *h_pin = NULL;        /* staging for fetches */
static size_t cap_pin = 0;

/* Page-locked staging for U and V. They arrive as PAML's own malloc'd globals,
   and a pageable H2D costs far more than its size suggests; copying into pinned
   memory first and sending that is cheaper than letting the driver stage it. */
static double *h_UV = NULL;
static size_t cap_hUV = 0;

static int ensure(double **p, size_t *cap, size_t need)
{
    if (*cap >= need) return 0;
    if (*p) cudaFree(*p);
    if (cudaMalloc(p, need) != cudaSuccess) { *p = NULL; *cap = 0; return -1; }
    *cap = need;
    return 0;
}

static int ensure_raw(void **p, size_t *cap, size_t need)
{
    if (*cap >= need) return 0;
    if (*p) cudaFree(*p);
    if (cudaMalloc(p, need) != cudaSuccess) { *p = NULL; *cap = 0; return -1; }
    *cap = need;
    return 0;
}

/* Grid is (cells/blockDim, batch): blockIdx.y selects the matrix, blockIdx.x a
   chunk of its n*n outputs. One block per matrix would launch only `batch`
   blocks, which leaves most of the device idle at codeml's batch of ~23 and
   makes each thread walk many cells serially; splitting the cells widens the
   grid by ceil(n*n/blockDim) so a small batch still fills the SMs.
   Output goes to slot[b]'s node block, so the result is already in the layout
   the traversal reads and no host-side scatter is needed.
   exptm1[k] is hoisted out of the i loop as PAML does. Accumulation runs over k
   in PAML's order so the rounding matches. */
__global__ void KPMatUVRoot(double *Pall, const double *Eall, const double *Uall,
                            const double *Vall, const int *slot,
                            int n, int ustride, int vstride)
{
    extern __shared__ double sh[];
    double *e = sh;

    int b = blockIdx.y;
    const double *U = Uall + (size_t)b * ustride;
    const double *V = Vall + (size_t)b * vstride;
    const double *E = Eall + (size_t)b * n;
    double *P = Pall + (size_t)(slot ? slot[b] : b) * n * n;

    for (int k = threadIdx.x; k < n; k += blockDim.x)
        e[k] = E[k];
    __syncthreads();

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n * n) return;

    int i = idx / n, j = idx % n;
    double acc = 0.0;
    for (int k = 0; k < n; k++)
        acc += (U[i * n + k] * e[k]) * V[k * n + j];
    if (i == j) acc += 1.0;          /* PAML: P[i*n+i]++ after the sum */
    P[idx] = acc < 0.0 ? 0.0 : acc;  /* PAML: if (P[i] < smallp) P[i] = 0 */
}

int PMatCudaAvailable(void)
{
    int count = 0;
    if (checked) return available;
    checked = 1;
    if (cudaGetDeviceCount(&count) == cudaSuccess && count > 0)
        available = 1;
    return available;
}

int PMatUVRootBatchCuda(double *P, const double *exptm1, const double *U,
                        const double *V, int n, int batch,
                        int ustride, int vstride)
{
    return PMatUVRootBatchCudaSlot(P, exptm1, U, V, NULL, n, batch,
                                   ustride, vstride, 0);
}

int PMatUVRootBatchCudaSlot(double *P, const double *exptm1, const double *U,
                            const double *V, const int *slot, int n, int batch,
                            int ustride, int vstride, int keep_on_device)
{
    size_t mat = (size_t)n * n;
    size_t nU, nV, nP;
    int maxnode = batch;

    if (!PMatCudaAvailable() || batch <= 0 || n <= 0) return -1;

    if (slot) {
        maxnode = 0;
        for (int i = 0; i < batch; i++)
            if (slot[i] + 1 > maxnode) maxnode = slot[i] + 1;
    }
    nP = mat * maxnode * sizeof(double);

    /* stride 0 means one matrix shared by the whole batch */
    nU = (ustride ? (size_t)ustride * batch : mat) * sizeof(double);
    nV = (vstride ? (size_t)vstride * batch : mat) * sizeof(double);

    if (ensure(&d_P, &cap_P, nP) ||
        ensure(&d_exptm1, &cap_exptm1, (size_t)n * batch * sizeof(double)) ||
        ensure(&d_U, &cap_U, nU) || ensure(&d_V, &cap_V, nV))
        return -1;
    if (slot) {
        if (ensure_raw((void **)&d_slot, &cap_slot, (size_t)batch * sizeof(int)))
            return -1;
        if (cudaMemcpy(d_slot, slot, (size_t)batch * sizeof(int),
                       cudaMemcpyHostToDevice) != cudaSuccess)
            return -1;
    }

    if (cudaMemcpy(d_exptm1, exptm1, (size_t)n * batch * sizeof(double),
                   cudaMemcpyHostToDevice) != cudaSuccess)
        return -1;

    /* Stage U and V through page-locked memory. */
    if (cap_hUV < nU + nV) {
        if (h_UV) cudaFreeHost(h_UV);
        if (cudaHostAlloc(&h_UV, nU + nV, cudaHostAllocDefault) != cudaSuccess) {
            h_UV = NULL; cap_hUV = 0;
        } else {
            cap_hUV = nU + nV;
        }
    }
    if (h_UV) {
        memcpy(h_UV, U, nU);
        memcpy((char *)h_UV + nU, V, nV);
        if (cudaMemcpy(d_U, h_UV, nU, cudaMemcpyHostToDevice) != cudaSuccess ||
            cudaMemcpy(d_V, (char *)h_UV + nU, nV, cudaMemcpyHostToDevice) != cudaSuccess)
            return -1;
    } else {
        if (cudaMemcpy(d_U, U, nU, cudaMemcpyHostToDevice) != cudaSuccess ||
            cudaMemcpy(d_V, V, nV, cudaMemcpyHostToDevice) != cudaSuccess)
            return -1;
    }

    {
        int thr = 256;
        dim3 grd((n * n + thr - 1) / thr, batch);
        KPMatUVRoot<<<grd, thr, n * sizeof(double)>>>(
            d_P, d_exptm1, d_U, d_V, slot ? d_slot : NULL, n,
            ustride ? ustride : 0, vstride ? vstride : 0);
    }

    /* A failed launch leaves d_P holding the previous call's results, which
       would be copied back and reported as a valid answer. */
    if (cudaGetLastError() != cudaSuccess) return -1;

    /* The traversal reads P on the device, so the copy back is skipped and P is
       never returned to the host at all. That copy is 29.8 kB per 227 kflop of
       work and dominates the end-to-end cost when it is made. */
    if (!keep_on_device) {
        /* The D2H copy is stream-ordered after the kernel and blocks until it
           has completed, so no separate device sync is needed. */
        if (cudaMemcpy(P, d_P, nP, cudaMemcpyDeviceToHost) != cudaSuccess)
            return -1;
    }
    return 0;
}

/* ---- resident ConditionalPNode ------------------------------------------ */

/* Nodes at the same height are independent, so one launch covers all of them:
   blockIdx.y picks the work item, blockIdx.x a chunk of that node's outputs.
   A launch per node per son would issue ~570 launches per traversal on a
   192-taxon tree, each only 48 blocks against 142 SMs, which leaves the device
   idle and makes launch overhead the limit rather than arithmetic. */

/* conP[h*n+j] = val over every node in dst. */
__global__ void KConPSetMany(double *conP, const size_t *dst, int n,
                             int pos0, int pos1, double val)
{
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = (size_t)(pos1 - pos0) * n;
    if (i >= total) return;
    conP[dst[blockIdx.y] + (size_t)pos0 * n + i] = val;
}

/* Folds one son into its parent, for every (parent, son) pair in w. Pairs in
   one launch have distinct parents, so the read-modify-write cannot race; the
   sons of a single parent go in separate launches, in PAML's order. */
__global__ void KConPFold(double *conP, const double *Pall, const char *z,
                          const char *nChara, const char *CharaMap,
                          const ConPWork *w, int n, int pos0, int pos1,
                          int npatt, int cleandata)
{
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = (size_t)(pos1 - pos0) * n;
    if (i >= total) return;

    ConPWork it = w[blockIdx.y];
    int h = pos0 + (int)(i / n), j = (int)(i % n);
    const double *PMat = Pall + (size_t)it.pslot * n * n;
    const double *pr = PMat + (size_t)j * n;
    double t;

    if (it.tip >= 0) {
        unsigned char c = (unsigned char)z[(size_t)it.tip * npatt + h];
        if (cleandata) {
            t = pr[c];                       /* PMat[j*n + z[h]] */
        } else {
            int nc = (int)(unsigned char)nChara[c];
            t = 0;
            for (int k = 0; k < nc; k++)
                t += pr[(unsigned char)CharaMap[(size_t)c * 64 + k]];
        }
    } else {
        const double *cs = conP + it.src + (size_t)h * n;
        t = 0;
        for (int k = 0; k < n; k++)          /* PAML's k order */
            t += pr[k] * cs[k];
    }
    conP[it.dst + (size_t)h * n + j] *= t;
}

/* One block per (pattern, node): max over states, then divide. max is exact and
   order-independent, so the reduction shape does not affect the result.
   log() stays on the host: it is npatt calls against npatt*n divisions, and the
   device log does not agree with the host libm in the last bit. */
__global__ void KNodeScaleMany(double *conP, const size_t *dst, const int *slot,
                               double *out, int n, int pos0, int pos1, int npatt)
{
    extern __shared__ double red[];
    int h = pos0 + blockIdx.x;
    if (h >= pos1) return;
    double *c = conP + dst[blockIdx.y] + (size_t)h * n;

    double m = 0;
    for (int j = threadIdx.x; j < n; j += blockDim.x)
        if (c[j] > m) m = c[j];
    red[threadIdx.x] = m;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s && red[threadIdx.x + s] > red[threadIdx.x])
            red[threadIdx.x] = red[threadIdx.x + s];
        __syncthreads();
    }
    m = red[0];

    if (m < 1e-300) {
        for (int j = threadIdx.x; j < n; j += blockDim.x) c[j] = 1.0;
    } else {
        for (int j = threadIdx.x; j < n; j += blockDim.x) c[j] /= m;
    }
    if (threadIdx.x == 0) out[(size_t)slot[blockIdx.y] * npatt + h] = m;
}

int ConPCudaResize(size_t sconP, int nnode, int ns, int npatt, int n,
                   int nscale)
{
    size_t nw = (size_t)nnode;
    if (!PMatCudaAvailable()) return -1;
    if (ensure(&d_conP, &cap_conP, sconP)) return -1;
    if (ensure(&d_P, &cap_P, (size_t)nnode * n * n * sizeof(double))) return -1;
    if (ensure(&d_scaleMax, &cap_scaleMax,
               (size_t)(nscale ? nscale : 1) * npatt * sizeof(double)))
        return -1;
    if (ensure_raw((void **)&d_z, &cap_z, (size_t)ns * npatt)) return -1;
    if (ensure_raw((void **)&d_work, &cap_work, 2 * nw * sizeof(ConPWork)) ||
        ensure_raw((void **)&d_dst, &cap_dst, nw * sizeof(size_t)) ||
        ensure_raw((void **)&d_sdst, &cap_sdst, nw * sizeof(size_t)) ||
        ensure_raw((void **)&d_slot2, &cap_slot2, nw * sizeof(int)))
        return -1;
    if (d_nChara == NULL &&
        (cudaMalloc(&d_nChara, 256) != cudaSuccess ||
         cudaMalloc(&d_CharaMap, (size_t)256 * 64) != cudaSuccess))
        return -1;
    if (cap_pin < (size_t)npatt * n * sizeof(double)) {
        if (h_pin) cudaFreeHost(h_pin);
        if (cudaHostAlloc(&h_pin, (size_t)npatt * n * sizeof(double),
                          cudaHostAllocDefault) != cudaSuccess) {
            h_pin = NULL; cap_pin = 0; return -1;
        }
        cap_pin = (size_t)npatt * n * sizeof(double);
    }
    return 0;
}

int ConPCudaUploadTips(const char *const *z, int ns, int npatt,
                       const char *nChara, const char *CharaMap)
{
    for (int i = 0; i < ns; i++)
        if (cudaMemcpy(d_z + (size_t)i * npatt, z[i], npatt,
                       cudaMemcpyHostToDevice) != cudaSuccess)
            return -1;
    if (cudaMemcpy(d_nChara, nChara, 256, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_CharaMap, CharaMap, (size_t)256 * 64,
                   cudaMemcpyHostToDevice) != cudaSuccess)
        return -1;
    return 0;
}

#define GRID(total) (int)(((total) + 255) / 256)

/* The work lists go up once per traversal, not once per launch. A staged copy
   out of pageable memory costs ~15 us regardless of size, so uploading per
   launch cost more on a shallow tree than the batching saved. */
int ConPCudaUploadWork(const ConPWork *w, int nw)
{
    if (nw <= 0) return 0;
    return cudaMemcpy(d_work, w, (size_t)nw * sizeof(ConPWork),
                      cudaMemcpyHostToDevice) == cudaSuccess ? 0 : -1;
}

int ConPCudaUploadDst(const size_t *dst, int nw)
{
    if (nw <= 0) return 0;
    return cudaMemcpy(d_dst, dst, (size_t)nw * sizeof(size_t),
                      cudaMemcpyHostToDevice) == cudaSuccess ? 0 : -1;
}

int ConPCudaUploadScale(const size_t *dst, const int *slot, int nw)
{
    if (nw <= 0) return 0;
    if (cudaMemcpy(d_sdst, dst, (size_t)nw * sizeof(size_t),
                   cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_slot2, slot, (size_t)nw * sizeof(int),
                   cudaMemcpyHostToDevice) != cudaSuccess)
        return -1;
    return 0;
}

/* start indexes the arrays uploaded above; nw is that group's length. */
void ConPCudaSetFrom(int start, int nw, int n, int pos0, int pos1, double val)
{
    size_t total = (size_t)(pos1 - pos0) * n;
    if (nw <= 0 || total == 0) return;
    { dim3 g(GRID(total), nw);
      KConPSetMany<<<g, 256>>>(d_conP, d_dst + start, n, pos0, pos1, val); }
}

void ConPCudaFoldFrom(int start, int nw, int n, int pos0, int pos1,
                      int npatt, int cleandata)
{
    size_t total = (size_t)(pos1 - pos0) * n;
    if (nw <= 0 || total == 0) return;
    { dim3 g(GRID(total), nw);
      KConPFold<<<g, 256>>>(d_conP, d_P, d_z, d_nChara, d_CharaMap,
                            d_work + start, n, pos0, pos1, npatt, cleandata); }
}

void ConPCudaScaleFrom(int start, int nw, int n, int pos0, int pos1, int npatt)
{
    int nb = pos1 - pos0, thr = 64;
    if (nw <= 0 || nb <= 0) return;
    while (thr < n && thr < 256) thr <<= 1;
    { dim3 g(nb, nw);
      KNodeScaleMany<<<g, thr, thr * sizeof(double)>>>(
          d_conP, d_sdst + start, d_slot2 + start, d_scaleMax, n, pos0, pos1,
          npatt); }
}

/* Scale maxima for every scaled node, fetched once per traversal rather than
   once per node: each fetch is a synchronization, and a 192-taxon tree scales
   at 10 nodes. */
int ConPCudaFetchScale(double *dst, int nscale, int npatt)
{
    if (nscale <= 0) return 0;
    return cudaMemcpy(dst, d_scaleMax, (size_t)nscale * npatt * sizeof(double),
                      cudaMemcpyDeviceToHost) == cudaSuccess ? 0 : -1;
}

int ConPCudaFetch(double *dst, size_t off, size_t nelem)
{
    if (cudaGetLastError() != cudaSuccess) return -1;
    if (h_pin && nelem * sizeof(double) <= cap_pin) {
        if (cudaMemcpy(h_pin, d_conP + off, nelem * sizeof(double),
                       cudaMemcpyDeviceToHost) != cudaSuccess)
            return -1;
        memcpy(dst, h_pin, nelem * sizeof(double));
        return 0;
    }
    return cudaMemcpy(dst, d_conP + off, nelem * sizeof(double),
                      cudaMemcpyDeviceToHost) == cudaSuccess ? 0 : -1;
}

int ConPCudaSyncToHost(double *conP, size_t sconP)
{
    if (d_conP == NULL || sconP > cap_conP) return -1;
    return cudaMemcpy(conP, d_conP, sconP, cudaMemcpyDeviceToHost)
           == cudaSuccess ? 0 : -1;
}

int ConPCudaCheck(void)
{
    return cudaGetLastError() == cudaSuccess ? 0 : -1;
}

void *PMatCudaHostAlloc(size_t bytes)
{
    void *p = NULL;
    if (!PMatCudaAvailable()) return NULL;
    if (cudaHostAlloc(&p, bytes, cudaHostAllocDefault) != cudaSuccess) return NULL;
    return p;
}

void PMatCudaHostFree(void *p)
{
    if (p) cudaFreeHost(p);
}

void PMatCudaCleanup(void)
{
    if (d_P) cudaFree(d_P);
    if (d_exptm1) cudaFree(d_exptm1);
    if (d_U) cudaFree(d_U);
    if (d_V) cudaFree(d_V);
    if (d_slot) cudaFree(d_slot);
    if (d_conP) cudaFree(d_conP);
    if (d_z) cudaFree(d_z);
    if (d_nChara) cudaFree(d_nChara);
    if (d_CharaMap) cudaFree(d_CharaMap);
    if (d_scaleMax) cudaFree(d_scaleMax);
    if (d_work) cudaFree(d_work);
    if (d_dst) cudaFree(d_dst);
    if (d_sdst) cudaFree(d_sdst);
    if (d_slot2) cudaFree(d_slot2);
    if (h_pin) cudaFreeHost(h_pin);
    if (h_UV) cudaFreeHost(h_UV);
    d_P = d_exptm1 = d_U = d_V = d_conP = d_scaleMax = NULL;
    h_pin = h_UV = NULL;
    d_slot = NULL; d_dst = NULL; d_sdst = NULL; d_slot2 = NULL; d_work = NULL;
    d_z = d_nChara = d_CharaMap = NULL;
    cap_P = cap_exptm1 = cap_U = cap_V = cap_slot = 0;
    cap_conP = cap_z = cap_scaleMax = cap_pin = cap_hUV = 0;
    cap_work = cap_dst = cap_slot2 = cap_sdst = 0;
}
