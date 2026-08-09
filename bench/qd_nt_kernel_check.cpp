/* Raw-limb identity check for the threaded NT gemm kernel (review2 §8.4, §10.1).
 *
 * B5's kernel test predates the NT commit and covered the NN kernel only, so
 * it cannot be cited for Rgemm_NT_omp. This harness is the NT-specific test:
 * it compares Rgemm_NT_omp (linked from the fork's own object file, exactly
 * the code the solver runs) against the serial NT body it replaced -- the
 * netlib dgemm loop copied verbatim from mpack/Rgemm.cpp's pre-dispatch NT
 * branch (git show 5e1cce0^:mpack/Rgemm.cpp) -- and requires every qd_real
 * limb of C to be memcmp-identical, with the padding rows of C untouched.
 *
 * Case axes: m/n/k from 1 to 129 (crossing the production gate both ways),
 * leading-dimension padding 0/3/7, alpha including 0 and negative, beta in
 * {0, 1, 0.37, -1} (0 and 1 are special-cased in the scaling pass), and
 * structural zeros placed in B by an INDEPENDENT predicate to exercise the
 * zero-skip branch. Fill values come from a fixed LCG: reproducible, no I/O.
 *
 * SELF-TESTS RUN FIRST, and the harness aborts if any fails -- a comparison
 * that cannot fail proves nothing (this project has shipped three such):
 *   - the fill really spans negatives (the v1 generator documented [-1,1) but
 *     returned [0,2), so its zero predicate `lcg() < -0.6` never fired and the
 *     advertised structural zeros were never placed: 0 of 405 cases had any);
 *   - a single perturbed limb in C is detected;
 *   - a write into C's padding rows is detected;
 *   - at least one zero was actually placed in B across the run, counted.
 *
 * The gate column reports the EXACT production predicate from
 * mpack_omp_tuning.h -- m*n*k >= MPACK_OMP_MIN_GEMM_WORK AND
 * n >= MPACK_OMP_MIN_GEMM_WIDTH -- not the work term alone (v1 labelled 18
 * n=1 cases "gate open" although production keeps them serial).
 *
 * Build (from the recipe repo root; FORK = fixes/sdpa-qd-omp, QD = the qd
 * install the fork was built against):
 *   g++-16 -fopenmp -O2 -D___MPACK_BUILD_WITH_QD___ -I$FORK/mpack -I$QD/include \
 *     patches/qd_nt_kernel_check.cpp $FORK/sdpa_qd-Rgemm_NT_omp.o \
 *     -L$QD/lib -lqd -o qd_nt_kernel_check
 * Run:  OMP_NUM_THREADS=<t> ./qd_nt_kernel_check
 * Output: one TSV row per case on stdout; exits nonzero on any limb mismatch,
 * padding clobber, or failed self-test.
 */
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <mblas_qd.h>
#include "mpack_omp_tuning.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void Rgemm_NT_omp(mpackint m, mpackint n, mpackint k, qd_real alpha,
    qd_real * A, mpackint lda, qd_real * B, mpackint ldb, qd_real beta,
    qd_real * C, mpackint ldc);

/* Serial NT body, copied verbatim from the pre-dispatch Rgemm.cpp NT branch. */
static void nt_serial(mpackint m, mpackint n, mpackint k, qd_real alpha,
    qd_real * A, mpackint lda, qd_real * B, mpackint ldb, qd_real beta,
    qd_real * C, mpackint ldc)
{
    qd_real Zero = 0.0, One = 1.0, temp;
    for (mpackint j = 0; j < n; j++) {
	if (beta == Zero) {
	    for (mpackint i = 0; i < m; i++)
		C[i + j * ldc] = Zero;
	} else if (beta != One) {
	    for (mpackint i = 0; i < m; i++)
		C[i + j * ldc] = beta * C[i + j * ldc];
	}
	for (mpackint l = 0; l < k; l++) {
	    if (B[j + l * ldb] != Zero) {
		temp = alpha * B[j + l * ldb];
		for (mpackint i = 0; i < m; i++)
		    C[i + j * ldc] = C[i + j * ldc] + temp * A[i + l * lda];
	    }
	}
    }
}

static unsigned long lcg_state = 12345;
static unsigned long lcg_raw()
{
    lcg_state = lcg_state * 6364136223846793005UL + 1442695040888963407UL;
    return lcg_state >> 33;                      /* 31 usable bits */
}
static double lcg()                              /* fill value in [-1, 1) */
{
    return (double)(long)(lcg_raw() % 200000) / 100000.0 - 1.0;
}
/* Zero selection is an INDEPENDENT draw, not a threshold on the fill value:
   a threshold couples "is zero" to the fill range and silently stops firing
   if that range moves, which is exactly how v1 lost its zero coverage. */
static int pick_zero()                           /* ~20% of entries */
{
    return (lcg_raw() % 5) == 0;
}

static const double SENTINEL = -777.125;         /* exact in binary; pad marker */

static int self_tests()
{
    /* 1. the fill spans negatives and stays inside [-1, 1) */
    double lo = 1e9, hi = -1e9;
    for (int i = 0; i < 100000; i++) { double v = lcg(); if (v < lo) lo = v; if (v > hi) hi = v; }
    if (!(lo < -0.9 && hi < 1.0 && hi > 0.9)) {
	fprintf(stderr, "SELF-TEST FAIL: fill range [%.5f, %.5f], want ~[-1, 1)\n", lo, hi);
	return 1;
    }
    /* 2. the zero predicate fires at roughly the advertised rate */
    int z = 0; for (int i = 0; i < 100000; i++) z += pick_zero();
    if (z < 15000 || z > 25000) {
	fprintf(stderr, "SELF-TEST FAIL: zero predicate fired %d/100000, want ~20%%\n", z);
	return 1;
    }
    /* 3. a single perturbed limb is detected, and 4. a padding write is detected */
    const int m = 4, n = 3, k = 5, ld = m + 2;
    std::vector<qd_real> A((size_t)ld * k), B((size_t)(n + 2) * k), C0((size_t)ld * n);
    for (auto &x : A) x = qd_real(lcg());
    for (auto &x : B) x = qd_real(lcg());
    for (auto &x : C0) x = qd_real(lcg());
    std::vector<qd_real> Ca = C0, Cb = C0;
    nt_serial(m, n, k, qd_real(1.0), A.data(), ld, B.data(), n + 2, qd_real(0.5), Ca.data(), ld);
    nt_serial(m, n, k, qd_real(1.0), A.data(), ld, B.data(), n + 2, qd_real(0.5), Cb.data(), ld);
    if (memcmp(Ca.data(), Cb.data(), Ca.size() * sizeof(qd_real)) != 0) {
	fprintf(stderr, "SELF-TEST FAIL: the serial body is not deterministic\n"); return 1;
    }
    Cb[0] = Cb[0] + qd_real(1e-60);              /* far below printed precision */
    if (memcmp(Ca.data(), Cb.data(), Ca.size() * sizeof(qd_real)) == 0) {
	fprintf(stderr, "SELF-TEST FAIL: a perturbed limb is NOT detected -- the comparison is vacuous\n");
	return 1;
    }
    Cb = Ca; Cb[m + 0 * ld] = qd_real(SENTINEL + 1.0);   /* first padding row */
    int pad_ok = 1;
    for (int j = 0; j < n; j++)
	for (int i = m; i < ld; i++)
	    if (Cb[i + (size_t)j * ld] != Ca[i + (size_t)j * ld]) pad_ok = 0;
    if (pad_ok) {
	fprintf(stderr, "SELF-TEST FAIL: a padding write is NOT detected\n"); return 1;
    }
    return 0;
}

int main()
{
    if (self_tests()) return 2;
    const int dims[] = {1, 2, 3, 5, 17, 33, 64, 65, 129};
    const int pads[] = {0, 3, 7};
    const double alphas[] = {1.0, -0.5, 0.0};
    const double betas[] = {0.0, 1.0, 0.37, -1.0};
    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif
    printf("m\tn\tk\tlda\tldb\tldc\talpha\tbeta\tthreads\tgate_open\tzeros_in_B\tresult\n");
    long cases = 0, failures = 0, total_zeros = 0, gated_on = 0;
    for (unsigned di = 0; di < sizeof(dims) / sizeof(*dims); di++)
    for (unsigned dj = 0; dj < sizeof(dims) / sizeof(*dims); dj += 2)
    for (unsigned dk = 0; dk < sizeof(dims) / sizeof(*dims); dk += 3) {
	int m = dims[di], n = dims[dj], k = dims[dk];
	int pad = pads[(di + dj + dk) % 3];
	int lda = m + pad, ldb = n + pad, ldc = m + pad;
	double a = alphas[(di + dk) % 3], b = betas[(di + dj) % 4];
	std::vector<qd_real> A((size_t)lda * k), B((size_t)ldb * k),
	    C0((size_t)ldc * n), Cs, Cp;
	for (auto &x : A) x = qd_real(lcg());
	/* Count only the region the kernel actually READS: B is ldb x k but the
	   loop touches B[j + l*ldb] for j < n, so zeros in the ldb-n padding rows
	   are never seen by the zero-skip branch and must not be claimed as
	   coverage (review2 §13.2 measured 2049 of the previously reported 28960
	   as padding). Padding is still filled, so an out-of-range read would
	   change the result. */
	int zeros = 0;
	for (mpackint l = 0; l < k; l++)
	    for (mpackint i = 0; i < ldb; i++) {
		bool z = pick_zero();
		B[i + (size_t)l * ldb] = z ? qd_real(0.0) : qd_real(lcg());
		if (z && i < n) zeros++;
	    }
	for (auto &x : C0) x = qd_real(lcg());
	/* mark the pad rows so an out-of-range write is caught, not missed */
	for (int j = 0; j < n; j++)
	    for (int i = m; i < ldc; i++)
		C0[i + (size_t)j * ldc] = qd_real(SENTINEL);
	Cs = C0; Cp = C0;
	nt_serial(m, n, k, qd_real(a), A.data(), lda, B.data(), ldb,
		  qd_real(b), Cs.data(), ldc);
	Rgemm_NT_omp(m, n, k, qd_real(a), A.data(), lda, B.data(), ldb,
		     qd_real(b), Cp.data(), ldc);
	int ok = memcmp(Cs.data(), Cp.data(),
			Cp.size() * sizeof(qd_real)) == 0;
	int pad_ok = 1;
	for (int j = 0; j < n && pad_ok; j++)
	    for (int i = m; i < ldc; i++)
		if (Cp[i + (size_t)j * ldc] != qd_real(SENTINEL)) { pad_ok = 0; break; }
	/* the EXACT production predicate, both terms */
	int gate = ((double)m * (double)n * (double)k >= MPACK_OMP_MIN_GEMM_WORK
		    && n >= MPACK_OMP_MIN_GEMM_WIDTH);
	cases++; total_zeros += zeros; gated_on += gate;
	if (!(ok && pad_ok)) failures++;
	printf("%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%d\t%d\t%d\t%s\n",
	       m, n, k, lda, ldb, ldc, a, b, threads, gate, zeros,
	       (ok && pad_ok) ? "identical" : (ok ? "PAD-CLOBBERED" : "LIMB-MISMATCH"));
    }
    if (total_zeros == 0) {
	fprintf(stderr, "FAIL: no structural zeros were placed -- the zero-skip branch was never exercised\n");
	failures++;
    }
    fprintf(stderr, "%ld cases (%ld with the production gate open, %ld serial), "
	    "%ld zeros placed in B, %ld failures, %d threads\n",
	    cases, gated_on, cases - gated_on, total_zeros, failures, threads);
    return failures ? 1 : 0;
}
