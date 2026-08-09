#ifndef MPACK_OMP_TUNING_H
#define MPACK_OMP_TUNING_H

/* New file added by this fork, 2026-08-04. It lives in the bundled MPACK tree and is offered
   under the same GNU Lesser General Public License version 3 terms as the files around it.
   No upstream file corresponds to it; the notices below record its history. */

/* NEW FILE (LGPL-3 4a; this file is LGPL-3 like the MPACK tree it lives in), 2026-08-05,
   originally for the threaded Rgemm NN kernel; the NT
 * kernel (2026-08-08, mpack/Rgemm_NT_omp.cpp) shares the same two gemm gate
 * macros, and this fork ALSO carries the B3 triangulars (mpack/Rtrsm_omp.cpp,
 * mpack/Rtrmm_omp.cpp) whose gates live further down this file.
 *
 * This is deliberately NOT a copy of the dd fork's mplapack/mplapack_omp_tuning.h.
 * That header additionally calibrates Rsyrk_omp / Rpotrf2, which do not exist
 * here. Importing the rest would import ~200 lines of dd
 * hardware tables describing code that is not in this tree, plus a
 * "generator refuses a stale header" version contract that has no generator on
 * this side. Only the two macros the gemm gate actually reads are carried over,
 * and BOTH VALUES HAD TO BE RE-DERIVED -- see below.
 */

/* Bump when macros are added or removed.
   1: initial -- MIN_GEMM_WORK / MIN_GEMM_WIDTH for Rgemm_NN_omp.
   2: MIN_TRSM_WORK / MIN_TRMM_WORK / TRI_WORK / MIN_TRI_WIDTH for the B3 port,
      Rtrsm_omp (Left/Lower/NoTranspose) and Rtrmm_omp (Left/Lower/Transpose). */
#define MPACK_OMP_TUNING_VERSION 2

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
   the per-problem gemm census quoted in BENCHMARKS.md (recipe repo: results/qd_* tables).

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

/* ------------ Left-side triangular kernels: Rtrsm_omp / Rtrmm_omp ("B3") ------------

   These two are the Cholesky-inverse phase: Rtrsm against an IDENTITY right-hand side
   in Lal::getInvLowTriangularMatrix (sdpa_linear.cpp:446) and Rtrmm forming
   Z^-1 = L**T * L in Jal::getInvCholAndInv (sdpa_jordan.cpp:197).  Both are split over
   the n COLUMNS of B; see the header comment in either kernel for why columns and not
   rows, and why that distinction is silent if got wrong.

   The gate is on WORK, m*m*n (the multiply-add count to within a factor of two, and the
   same convention the dd fork uses so the two are comparable), not on one dimension, so
   a short wide B is judged by what it costs.  Keep the (double) casts at the use site:
   mpackint is int32_t in this fork and SDPLIB contains blocks of 1600, where a bare
   m*m*n is 4.1e9 and overflows.

   THE dd VALUE DOES NOT TRANSFER, AND THE ERROR WOULD HAVE BEEN INVISIBLE.  dd ships
   MIN_TRSM_WORK = MIN_TRMM_WORK = 15625 (n = 25 cubed).  Neither of dd's two
   justifications survives the move to this fork: its measured break-even was n = 10-13
   (thanos) / 12-18 (pi) for dd_real, and its "SDPLIB has a dead zone between 19 and 35"
   argument is about where dd's crossover happened to fall.  Note the DIRECTION of the
   error: for these Left kernels an over-large gate is not merely suboptimal, it is
   invisible -- the ported kernel becomes a no-op that still links and still passes a
   bit-identity ship gate.  Measured below: the qd crossover is at n = 4-5, not 25.

   MEASURED, 2026-08-05, pi (i9-13900K, 8 P-cores + 16 E-cores = 24 physical, idle:
   load 0.17 before, no foreign sdpa_ or regress process, busy_cores 1.8-3.6 per row out
   of 24).  Method is B5's, which is the one that worked in this codebase: ONE binary,
   ONE set of compiler flags, each kernel compiled TWICE from its own source with only
   the gate constant differing (1e18 = always delegate to the serial kernel, 0 = always
   thread), the two arms called INTERLEAVED per repetition so a burst of foreign load
   lands on both, minimum over repetitions, threads bound one per physical core
   (P-cores are logical 0,2,..,14 on this box; 24 threads uses 0-31), 0.4 s warm-up per
   size.  Do NOT instead compare a no-OpenMP build against an -fopenmp build: B5
   recorded that producing 3.1x on TWO threads at m=n=k=5, arithmetically impossible.

   OPERANDS ARE THE SOLVER'S.  This matters more here than for gemm because both
   kernels carry a zero-skip, so the zero pattern sets the cost profile: Rtrsm gets
   A = the lower Cholesky factor and B = the IDENTITY; Rtrmm gets A = B = L^-1, lower
   triangular; alpha = 1 for both (MONE is 1.0, sdpa_tool.h:68).  B is restored from a
   pristine copy before every call, outside the timed window.

   Two designs were run, because they disagree about trap 1 and it was not obvious which
   is right.  Design A is the interleaved one above: a GOMP thread pool IS alive when the
   gate rejects, which is what happens in the real solver.  Design B puts the serial arm
   in its own OMP_NUM_THREADS=1 process so no team is ever created.  They agree to within
   a few percent everywhere and agree on the sign at every size, so the choice does not
   matter here; A is the more pessimistic of the two at small n and is what the value
   below is set from.

   speedup (gate shut / gate open):

       Rtrsm: A = interleaved one process        B = serial arm alone at OMP_NUM_THREADS=1
       n    work        A 2t/4t/8t/24t              B 2t/4t/8t/24t
       3       27     0.79  0.72  0.73  0.73         0.95  0.80  0.85  0.66
       4       64     0.98  0.84  0.84  0.83         1.18  1.12  1.11  0.99
       5      125     1.14  1.23  1.18  1.21         1.30  1.57  1.32  1.37
       6      216     1.27  1.63  1.58  1.54         1.39  1.95  1.73  1.74
       8      512     1.53  2.24  2.27  2.23         1.66  2.44  2.66  2.68
      10     1000     1.69  2.72  3.05  2.76         1.74  2.96  3.52  3.26
      12     1728     1.80  3.05  3.90  3.44         1.85  3.18  4.23  3.78
      16     4096     1.90  3.44  5.45  4.74         1.90  3.56  5.43  4.02
      20     8000     1.95  3.65  6.33  6.06         1.96  3.73  6.34  5.18
      25    15625     1.96  3.81  6.77  6.20         1.99  3.85  6.84  5.96
      35    42875     2.00  3.92  7.51  8.60         2.01  3.94  6.81  8.36
      50   125000     2.03  3.99  7.51  8.83         2.03  4.02  7.15 11.00
      70   343000     1.97  3.75  7.20  8.08         2.04  4.04  7.29 13.66
     100  1000000     2.02  3.91  7.31 11.50         2.04  4.03  8.13 13.93
     161  4173281     2.05  4.05  7.95 14.34         2.05  4.03  7.08 14.10

       Rtrmm: A = interleaved one process        B = serial arm alone at OMP_NUM_THREADS=1
       n    work        A 2t/4t/8t/24t              B 2t/4t/8t/24t
       3       27     0.71  0.51  0.58  0.62         0.64  0.71  0.67  0.63
       4       64     0.82  0.74  0.76  0.74         1.01  0.93  0.84  0.86
       5      125     1.14  1.14  1.09  1.05         1.06  1.28  1.32  1.21
       6      216     1.18  1.42  1.37  1.40         1.34  1.74  1.50  1.65
       8      512     1.48  2.09  2.56  2.54         1.60  2.36  2.77  2.49
      10     1000     1.70  2.76  3.39  3.20         1.74  2.84  3.66  3.69
      12     1728     1.80  3.19  4.20  3.03         1.86  3.26  4.41  2.91
      16     4096     1.92  3.51  5.70  4.29         1.94  3.54  6.09  4.35
      20     8000     1.97  3.72  6.06  5.47         1.99  3.79  5.64  5.57
      25    15625     1.98  3.83  7.19  7.08         2.00  3.84  6.58  7.22
      35    42875     2.01  3.97  7.49 10.90         2.02  3.97  6.84 10.13
      50   125000     1.94  4.03  7.46  9.22         2.04  4.01  7.92 11.84
      70   343000     2.01  3.84  7.01  9.96         2.05  4.02  6.95 13.05
     100  1000000     2.04  4.00  7.66 12.31         2.05  4.03  7.38 13.53
     161  4173281     2.05  4.08  8.08 14.45         2.05  4.10  7.11 13.75

   Worst cell over all sixteen (kernel, thread-count, design) combinations, by size:

       n     3     4     5     6     8    10    12    16    20    25 .. 161
     worst 0.51  0.74  1.05  1.18  1.48  1.69  1.80  1.90  1.95  1.96 .. 2.05

   So the break-even is between n = 4 and n = 5, i.e. work between 64 and 125 -- a factor
   of 125-244 below dd's 15625, and in the same direction and of the same magnitude as
   the correction B5 already had to make to MIN_GEMM_WORK (dd's 20000 -> 512).

   REPRODUCIBILITY.  The whole grid was run twice on an idle pi.  Every row at n >= 8 --
   the rows the choice below actually rests on -- agreed to within 5.9%, and n=8 and n=10
   to within 1.6% and 2.2%.  The n <= 6 rows move by up to 19% run to run, but agree on
   the sign in both runs: n=3 and n=4 lose in both, n=5 and n=6 win in both.

   VALUE CHOSEN: 512, which is n = 8 cubed.  Three reasons.

     - Margin.  n=5 (work 125) is the first size that wins in all sixteen cells, but its
       worst cell is 1.05x -- inside the run-to-run spread at that size, so a gate cannot
       be sat on it.  n=8 is the first size with real headroom: worst cell 1.48x, best
       2.77x.  This is the same criterion dd applied when it rejected n=20 (worst 1.16x)
       in favour of n=25.

     - It costs nothing on the real corpus.  The dense SDP block sizes in this fork's
       10-problem set are, with multiplicity: 2x1, 150x3, 2x4, 1x5, 1x6, 34x10, 1x26,
       1x50, 3x100, 1x161 (censused directly from the .dat-s files, not inherited).
       Nothing lies between n=6 (work 216) and n=10 (work 1000), so EVERY threshold in
       (216, 1000] admits and rejects exactly the same real blocks; 512 sits inside that
       gap.  What it gives up against a threshold of 216 is hinf1's single n=6 block,
       measured 1.18-1.95x on a block whose work is 216 out of that problem's 344.

     - It is the value MIN_GEMM_WORK already carries, measured independently for the
       gemm kernel, so this fork carries one number rather than two that differ by less
       than their own uncertainty.

   THE TWO KERNELS GET THE SAME NUMBER.  Rtrsm crosses very slightly earlier than Rtrmm
   (at n=4 Rtrsm reads 0.83-1.18 against Rtrmm's 0.74-1.01), but the gap between them is
   smaller than the run-to-run spread at those sizes and no real block lands between the
   two crossovers, so splitting them would be false precision.

   WHAT dd's 15625 WOULD HAVE COST, quantified rather than asserted.  15625 falls in the
   corpus gap (1000, 17576], so it rejects every block with n <= 10 and admits only the
   six blocks at n = 26, 50, 100x3, 161.  By cube-weighted work that is 99.48% of the
   corpus still admitted, so the survey note that predicted "approximately nothing" from
   dd's constant is TOO STRONG and is corrected here.  The real cost is concentrated, not
   diffuse: it falls entirely on the problems whose blocks are all small.  truss5
   (1x1 + 33x10) and control1 (1x5 + 1x10) would get ZERO of their Cholesky-inverse phase
   threaded under 15625 and 99.997% / 88.9% of it under 512; truss6 (1x1 + 150x3) and
   hinf1 (2x4 + 1x6) get nothing under either, correctly, because n=3 and n=4 measure
   0.51-1.18x.  The other six problems are unaffected by the choice.  Measured speedups
   at the sizes 15625 throws away: 1.48-2.77x at n=8, 1.69-3.69x at n=10.

   MEASURED ON ONE MACHINE.  pi only.  thanos (EPYC 7232P, 8 physical) was carrying ~0.8
   of foreign load from other users for the duration and was not used for this grid; a
   second-machine confirmation is outstanding and this paragraph should not be deleted
   until it exists.  The value is nevertheless safe against a slower box: 512 sits three
   sizes above the measured crossover, and the quantity that would move it is the ratio
   of fork/join to a qd_real multiply-add, which on a slower core moves in the direction
   that makes threading pay SOONER, not later. */
/* Rtrsm_omp, Left/Lower/NoTranspose: minimum m*m*n to thread.  n = 8 cubed, MEASURED
   for qd_real on pi.  NOT dd's 15625, which is 30x too large for this fork. */
#ifndef MPACK_OMP_MIN_TRSM_WORK
#define MPACK_OMP_MIN_TRSM_WORK 512.0
#endif

/* Rtrmm_omp, Left/Lower/Transpose: same value, same measurement -- see above for why
   the two are not split. */
#ifndef MPACK_OMP_MIN_TRMM_WORK
#define MPACK_OMP_MIN_TRMM_WORK 512.0
#endif

/* Flat across the measured range (2..24 threads), which is what the table shows: no
   monotone function of the team size fits it -- at n=161 the 24-thread column is the
   loosest and at n=12 it is the tightest.  Grown linearly past 24 threads, where nothing
   has been measured on this fork, so an unvalidated team size errs towards the serial
   kernel rather than being silently extrapolated.  Form carried over from dd unchanged:
   the form is a guard, not a measurement. */
#ifndef MPACK_OMP_TRI_WORK
#define MPACK_OMP_TRI_WORK(w, nt) ((nt) <= 24 ? (double)(w) : (double)(w) * (double)(nt) / 24.0)
#endif

/* Minimum number of columns to share out.  Work alone is not sufficient: a tall thin
   solve with n=1 can clear the work gate and still offer the team a single iteration. */
#ifndef MPACK_OMP_MIN_TRI_WIDTH
#define MPACK_OMP_MIN_TRI_WIDTH 2
#endif

#endif
