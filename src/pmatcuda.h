/* CUDA backend for codeml's likelihood inner loops.
 *
 * Two pieces, both bit-identical to the C they replace:
 *
 *   PMatUVRoot      P(t) = U * exp{Root*t} * V          (tools.c:516)
 *   ConditionalPNode  the post-order traversal          (codeml.c:3690)
 *
 * PMatUVRoot alone is 39k-65k calls of a 61x61x61 dense double product per
 * codeml run on examples/HIVNSsites, and it is the whole cost on small trees.
 * On large trees it is not: at 192 taxa the traversal's own GEMM
 * (npatt x n x n per internal node) dominates, and returning P to the host
 * costs more than computing it there. Keeping P and conP resident removes both
 * the copy and the traversal's arithmetic.
 *
 * Results are bit-identical because every accumulation runs in PAML's own
 * order, one thread per output cell. expm1 and log stay on the host: each is
 * called n or npatt times against n^3 or npatt*n^2 for the products, so moving
 * them saves nothing, and CUDA's expm1 and log do not agree with the host libm
 * in the last bit.
 */

#ifndef PMATCUDA_H
#define PMATCUDA_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Non-zero if a usable device was found. Checked once; the caller falls back
   to the C path when this returns 0. */
int PMatCudaAvailable(void);

/* Batched P(t). Computes batch matrices of size n*n into P.
 *
 *   P     [batch][n*n]   output, row-major, as PMatUVRoot writes it
 *   exptm1[batch][n]     expm1(t*Root[k]) computed by the caller on the host
 *   U     [batch][n*n]   right eigenvectors
 *   V     [batch][n*n]   inverse eigenvectors
 *
 * U and V may be shared across the batch; pass ustride/vstride 0 to reuse a
 * single matrix for every entry, which is the common case since eigenQcodon
 * runs far less often than P(t).
 *
 * Returns 0 on success, non-zero on a CUDA error, in which case P is
 * untouched and the caller should use the C path.
 */
int PMatUVRootBatchCuda(double *P, const double *exptm1, const double *U,
                        const double *V, int n, int batch,
                        int ustride, int vstride);

/* As above, with two additions used by the resident traversal:
 *
 *   slot[batch]      writes matrix b to P + slot[b]*n*n instead of P + b*n*n,
 *                    so results land indexed by node and need no host scatter.
 *                    NULL keeps the plain batch layout.
 *   keep_on_device   non-zero skips the copy back; P is left on the device for
 *                    ConPCuda* to read, and the host pointer is untouched.
 */
int PMatUVRootBatchCudaSlot(double *P, const double *exptm1, const double *U,
                            const double *V, const int *slot, int n, int batch,
                            int ustride, int vstride, int keep_on_device);

/* Resident ConditionalPNode.
 *
 * ConPCudaResize mirrors com.conP on the device: a node's block sits at the
 * same offset the host uses, so codeml's per-site-class pointer shift needs no
 * translation. All offsets below are element offsets into that mirror,
 * i.e. nodes[i].conP - com.conP.
 *
 * Work is issued a tree level at a time. Nodes at the same height are
 * independent, so one launch covers all of them and the grid is wide enough to
 * fill the device; per-node launches are too small to do either.
 *
 * The calls are launches only and do not synchronize; ConPCudaCheck reports any
 * error since the last check, and the fetches synchronize implicitly.
 */

/* One (parent, son) fold. Pairs handed to ConPCudaFold together must have
   distinct parents, since each writes its parent's block. */
typedef struct {
   size_t dst;    /* parent conP offset */
   size_t src;    /* son conP offset; unused when tip >= 0 */
   int pslot;     /* son node index, selects its P(t) */
   int tip;       /* tip index into z, or -1 if the son is internal */
} ConPWork;

int ConPCudaResize(size_t sconP, int nnode, int ns, int npatt, int n,
                   int nscale);

/* Tip data, uploaded once per alignment. CharaMap is the flat [256][64] map. */
int ConPCudaUploadTips(const char *const *z, int ns, int npatt,
                       const char *nChara, const char *CharaMap);

/* The work lists go up once per traversal; the launchers then index into them
   by group. Uploading per launch costs more than the batching saves, since a
   staged copy out of pageable memory is ~15 us regardless of size. Buffers
   passed here should be page-locked. */
int ConPCudaUploadWork(const ConPWork *w, int nw);
int ConPCudaUploadDst(const size_t *dst, int nw);
int ConPCudaUploadScale(const size_t *dst, const int *slot, int nw);

/* conP[h*n+j] = val for dst[start .. start+nw). */
void ConPCudaSetFrom(int start, int nw, int n, int pos0, int pos1, double val);

/* Folds w[start .. start+nw). Those pairs must have distinct parents; the sons
   of one parent go in separate calls, in PAML's order. cleandata selects the
   observed-state lookup over the ambiguity sum, matching ConditionalPNode's two
   tip branches. */
void ConPCudaFoldFrom(int start, int nw, int n, int pos0, int pos1,
                      int npatt, int cleandata);

/* Scales as NodeScale does, writing each node's per-pattern maxima to its own
   row of the maxima buffer. */
void ConPCudaScaleFrom(int start, int nw, int n, int pos0, int pos1, int npatt);

/* All scale maxima, nscale*npatt doubles, from which the caller takes log(). */
int ConPCudaFetchScale(double *dst, int nscale, int npatt);

/* nelem doubles from the mirror at off. Used for the root block, the only part
   of conP the likelihood consumes on the host. */
int ConPCudaFetch(double *dst, size_t off, size_t nelem);

/* Whole mirror back to com.conP. Used once when the device stops owning conP,
   so host paths that read interior nodes see current values. */
int ConPCudaSyncToHost(double *conP, size_t sconP);

/* 0 if no CUDA error is pending. */
int ConPCudaCheck(void);

/* Page-locked host memory. Transfers from pageable memory stage through a
   driver bounce buffer, which dominates the cost at these batch sizes; buffers
   handed to PMatUVRootBatchCuda should come from here. Returns NULL on
   failure, in which case ordinary malloc still works, only slower. */
void *PMatCudaHostAlloc(size_t bytes);
void PMatCudaHostFree(void *p);

/* Frees device buffers. Safe to call when nothing was allocated. */
void PMatCudaCleanup(void);

#ifdef __cplusplus
}
#endif

#endif
