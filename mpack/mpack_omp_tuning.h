#ifndef MPACK_OMP_TUNING_H
#define MPACK_OMP_TUNING_H

/* NEW FILE (2026-08-05), for the threaded Rgemm NN kernel in this fork.
 *
 * This is deliberately NOT a copy of the dd fork's mplapack/mplapack_omp_tuning.h.
 * That header is 230 lines and calibrates Rtrsm_omp / Rtrmm_omp / Rsyrk_omp /
 * Rpotrf2, none of which exist here -- this fork's mpack/ has no threaded kernel
 * other than the Rgemm NN one. Importing the rest would import ~200 lines of dd
 * hardware tables describing code that is not in this tree, plus a
 * "generator refuses a stale header" version contract that has no generator on
 * this side. Only the two macros the gemm gate actually reads are carried over,
 * and BOTH VALUES HAD TO BE RE-DERIVED -- see below.
 */

/* Bump when macros are added or removed.
   1: initial -- MIN_GEMM_WORK / MIN_GEMM_WIDTH for Rgemm_NN_omp. */
#define MPACK_OMP_TUNING_VERSION 1

/* ---------------------------------------------------------------------------
   MPACK_OMP_MIN_GEMM_WORK -- minimum m*n*k before an OpenMP fork/join pays for
   itself in Rgemm_NN_omp.

   THE dd VALUE DOES NOT TRANSFER. The dd fork uses 20000.0, and that number
   counts *dd* multiply-adds. Its own header records a dd multiply-add at
   ~25-50 ns (later corrected to 2.6-4.1 ns on pi). A *qd* multiply-add is one
   to two orders of magnitude more expensive: measured in this tree on thanos
   (EPYC 7232P, instrumented build, 1 thread, uncontended),

       control1: 1248 calls of 10x10x10  =  0.1995 s  ->  160 us/call, ~160 ns/MAD
       gpp100 :   410 calls of 100^3     = 50.264 s   ->  123 ms/call, ~123 ns/MAD

   so a 10x10x10 qd gemm is already a 160 us serial call -- 30-80x a fork/join --
   yet m*n*k = 1000 is rejected by 20000. Carried over unchanged, dd's constant
   rejects 100% of control1's Rgemm (41.9% of that problem's wall) and 100% of
   truss5's NN bucket (13266 calls, 9.2% of wall). See
   patches/b5_notes/08_per_problem_gemm_census.md.

   MEASURED, 2026-08-05, thanos (EPYC 7232P, 8 physical cores, SMT2, schedutil).
   Method, and why it is not the obvious one: comparing an -fopenmp build of this
   kernel against a NO-OpenMP build produced speedups of 3.1x on TWO threads at
   m=n=k=5, which is arithmetically impossible when the split is over 5 columns.
   The two compilations of this source simply do not run at the same speed at tiny
   sizes, and that difference was being read as threading. The measurement below
   removes the confound: ONE binary, ONE set of flags, the kernel compiled twice
   from this same file with only MIN_GEMM_WORK differing (1e18 = gate always shut,
   0 = gate always open), the two called INTERLEAVED per repetition so a burst of
   foreign load lands on both arms, minimum over repetitions, threads bound one per
   physical core, 0.4 s warm-up per size against the 1500 MHz idle floor. That is
   exactly the decision this constant makes at run time, and the shut arm pays the
   same team-of-one cost it really pays when the gate rejects.

   speedup (gate shut / gate open), m = n = k, cubic:

       m=n=k   work      2 threads       4 threads       8 threads
                       meas / max      meas / max      meas / max
         5      125    1.59 / 1.67     2.12 / 2.50     2.98 / 5.00
         8      512    1.97 / 2.00     3.76 / 4.00     3.98 / 8.00
        10     1000    1.98 / 2.00     3.24 / 3.33     4.67 / 5.00
        12     1728    1.99 / 2.00     3.92 / 4.00     5.73 / 6.00
        16     4096    2.00 / 2.00     3.97 / 4.00     7.74 / 8.00
        20     8000    2.00 / 2.00     3.80 / 4.00     6.58 / 6.67
        32    32768    2.00 / 2.00     4.00 / 4.00     7.96 / 8.00
        50   125000    2.00 / 2.00     3.84 / 3.85     6.19 / 7.14
       100  1000000    1.51 / 2.00     2.67 / 4.00     3.83 / 7.69
       161  4173281    1.26 / 1.99     2.24 / 3.93     4.41 / 7.67

   "max" is the best a static split of n columns over min(T,n) threads can give.
   NOTHING IN THE MEASURED RANGE LOSES -- the crossover is below m=n=k=5, i.e.
   below work = 125 -- and the small sizes run at 95-100% of the achievable
   speedup. The 2-thread column reproduced to within 0.5% across repeated runs.

   THE BOX WAS SHARED while these were taken (another workflow held 4-6 of the 8
   cores; the harness records busy-core counts of 5-7.5 per row). That biases
   THREADED arms slow, so every number above is a LOWER bound and the conclusion
   "no tested shape loses" cannot be an artefact of it. It does mean the large-size
   rows (>= 64) are contention-limited and must NOT be quoted as this kernel's
   achievable scaling -- 161^3 at 8 threads reads 4.41x here and would read higher
   on an idle machine. Re-measure on an idle box before publishing any speedup.

   VALUE CHOSEN: 512, the smallest size measured to reach >= 98% of the achievable
   2-thread speedup. It is deliberately not set at the very edge of the measured
   range (125): thresholds should not be extrapolated from their last data point.
   It admits every NN shape this fork's corpus actually issues -- the smallest is
   10x10x10 (work 1000, truss5 and control1) which measures 1.98x/3.24x/4.67x.
   dd's 20000 would reject all of it.

   Note the (double) casts at the use site are load-bearing in THIS fork:
   mpack_config.h leaves USE64BITINT commented out, so mpackint is int32_t and a
   bare m*n*k overflows at m=n=k=1290. dd's mplapackint is int64_t and hides it. */
#ifndef MPACK_OMP_MIN_GEMM_WORK
#define MPACK_OMP_MIN_GEMM_WORK 512.0
#endif

/* Minimum width of the parallelised gemm loop. Rgemm_NN_omp splits the j loop,
   so the number of iterations available to share out is n. Work alone is not
   sufficient: a tall thin gemm (n = 1, large m*k) can clear MIN_GEMM_WORK and
   still offer exactly one iteration, paying team creation for no parallelism. */
#ifndef MPACK_OMP_MIN_GEMM_WIDTH
#define MPACK_OMP_MIN_GEMM_WIDTH 2
#endif

#endif
