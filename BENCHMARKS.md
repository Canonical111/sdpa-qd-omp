# Benchmarks — sdpa-qd-omp

Fork base: upstream `766eef3`. Measured at `443799f`.

**Read this first.** Upstream sdpa-qd's `main loop time` is process CPU time summed over
threads -- it *rises* with thread count even when the run gets faster, so any historical
threaded benchmark read from sdpa-qd's own output is suspect. This fork fixes the clock;
everything below is external wall time, taken outside the process.

## Summary

Two arms: the fork at `3376c79` (`base*`, the state the published tables were measured in)
and at `443799f` (`head*`). Both built from `git archive` of the named commit, on the machine
they run on, md5-distinct. Comparisons are **fork against fork at equal thread count**, over
the 4 of 5 published problems where the two arms walk the same path.

| | thanos, 8 cores | pi, 24 cores |
|---|---|---|
| port payoff at 1 thread (the neutrality control) | 1.031x | 1.021x |
| port payoff at 8 threads | **1.947x** | **2.329x** |
| port payoff at 24 threads | n/a | **2.660x** |
| thread scaling, pre-port arm (all 5 problems) | 1.384x (1->8) | 1.501x (1->24) |
| thread scaling, `443799f` (all 5 problems) | 2.862x (1->8) | 4.518x (1->24) |
| end to end, pre-port serial -> `443799f` threaded | 2.941x | 4.548x |

`gpp100` is excluded from those columns because its iteration count differs between the arms;
it is quoted per iteration in its own section, where it is also the source of this
regeneration's principal refutation.

## These tables replace the previously published qd tables. They do not refresh them.

Four things changed under the old tables. Each of them changes what a column *means*, which
is why the numbers were re-measured from scratch rather than re-run.

**1. `--enable-openmp` was inert, and is fixed in `f177644`.** At the published commit,
`AC_ARG_ENABLE(openmp)` declared the option with an empty action and the block that
consumed it expanded `$OPENMPFLAG`, which only `--with-openmpflag` ever sets; `AC_OPENMP`
was never called. So the switch added two preprocessor defines that nothing in the tree
reads, and no `-fopenmp`. The resulting binary still accepted `OMP_NUM_THREADS` and still
printed the same iteration table, so **it reported exactly 1.00x at every thread count and
was indistinguishable, in program output, from a threaded build.** What that means for the
old document is blunt: **its qd matrix could have been reproduced with a serial binary and
nothing in the table would have shown it.** The published binaries were in fact threaded,
but only because the build wrote `-fopenmp` into `CXXFLAGS` by hand; a reader who followed
the fork's own `--enable-openmp` instruction at that commit got the serial binary.

#### The six builds behind that claim

| arm | how OpenMP was requested | binary md5 | size (B) | OpenMP symbols in binary | compile lines with `-fopenmp` | TUs with OpenMP |
|---|---|---|---|---|---|---|
| baseline `3376c79` | CXXFLAGS=-O2 -funroll-all-loops -fopenmp | `71abe618bbf8e25b34923cdaec547cd3` | 6138496 | 169 | 61 | 1 of 60 |
| HEAD `443799f` | CXXFLAGS=-O2 -funroll-all-loops -fopenmp | `b8d52bfa70ae8b30632eb69ad37458a1` | 6226600 | 169 | 64 | 4 of 63 |
| baseline `3376c79` | --enable-openmp, no -fopenmp in CXXFLAGS | `28b4ced44fdcae1a7dfb55b0b25b4276` | 5527552 | 0 | 0 | — |
| baseline `3376c79` | no OpenMP flag at all | `28b4ced44fdcae1a7dfb55b0b25b4276` | 5527552 | 0 | 0 | — |
| HEAD `443799f` | --enable-openmp, no -fopenmp in CXXFLAGS | `b8d52bfa70ae8b30632eb69ad37458a1` | 6226600 | 169 | 64 | — |
| HEAD `443799f` | no OpenMP flag at all | `f3576332fb155e59db034fd21156a3f2` | 5611328 | 0 | 0 | — |

Read the two decisive pairs. At the baseline commit `--enable-openmp` yields
`28b4ced44fdcae1a7dfb55b0b25b4276` and so does passing no OpenMP flag at all — **byte-identical, zero
OpenMP symbols, `-fopenmp` on zero compile lines.** At `443799f`, `--enable-openmp` yields
`b8d52bfa70ae8b30632eb69ad37458a1`, which is byte-identical to the benchmarked `head` arm built with
`-fopenmp` in `CXXFLAGS`. The switch now delivers exactly what the manual flag delivers,
which is the strongest form the fix could take, and the no-flag build at HEAD
(`f3576332fb155e59db034fd21156a3f2`) is a different binary again, so the probe can tell the three
states apart.

**2. qd had no threaded BLAS at all; it now has a threaded `Rgemm` and threaded triangulars.**
A census of undefined `GOMP_*`/`omp_*` references across the translation units finds
**1 of 60** with any at the baseline commit — `sdpa_newton:10`,
the `k1` bMat loop in the Newton step, which is SDPA-level, not BLAS. `Rgemm`, `Rtrsm`,
`Rtrmm`, `Rsyrk`, `Rgemv`, `Rpotrf`, `Rsteqr`, `Rdot`, `Raxpy` and `Rcopy` had none.
At `443799f` the count is **4 of 63** —
`Rgemm_NN_omp:5,Rtrmm_omp:6,Rtrsm_omp:6,sdpa_newton:10` — from `f177644` (threaded `Rgemm` NN case) and `443799f` (threaded
X/Z inverse-Cholesky triangulars). Counted a second way, `-fopenmp` appears on
61 compile lines in the baseline build and 64 at HEAD: the same
three new translation units. The consequence is visible directly in the thread census below,
where the pre-port binary **never creates a second thread** on two of the five published
problems even at `OMP_NUM_THREADS=8` pinned to eight distinct physical cores.

**3. `63370f0` fixed an uninitialised block-size lookup in `iMlaenv`** that made iteration
counts depend on stack contents. The old table's `gpp100` row carried two different
iteration counts, 63 in the `pristine` column and 49 in the `hardened*` columns. What the fix did and did not do is set out under
`gpp100` below, and it is not what this document was originally expected to say.

**4. Termination phase was not recorded.** `gpp100` never reaches `pdOPT` in any build
measured here: it stops `pFEAS`, primal feasible with the dual not. The old tables therefore
presented a non-solve as a solve, and no ratio on that row can be rescued by dividing by
iteration count, because at least one of the two trajectories was not a solution at all.

## What the two arms are

| arm | commit | what it is |
|---|---|---|
| `base_*` / `base*` | `3376c79` | the fork **at the commit whose tree carries the published `BENCHMARKS.md` byte-for-byte** — pre-port, pre-`63370f0` |
| `head_*` / `head*` | `443799f` | measurement tip: `f177644` (threaded `Rgemm` NN, `--enable-openmp` fixed) plus `443799f` (threaded X/Z triangulars) |

> **Since these measurements** (`5e1cce0`, 2026-08-08): the `Rgemm` **NT** case is threaded as
> well. NT is `Rpotrf`'s blocked trailing update, i.e. the Schur-complement Cholesky, which
> every table on this page left serial. The NT work gate (shared with NN, `m*n*k >= 512`)
> **does admit** panels on this page's problems -- arch0's second panel alone is
> `46*64*64 = 188416` -- so the correct statement is not "gated off" but *measured unchanged*:
> at the current tip, `arch0`, `truss5` produce printed solutions IDENTICAL to the released
> binary's (same iterations, same full solution section; the raw record is retained in
> `results/qd_nt_solution_identity.tsv` in the recipe repository, as
> `bench/qd_nt_solution_identity.tsv` in the solver repository). The Schur phase is simply a
> small share of runtime at m <= 208, so the timing columns, measured at `443799f`, are left
> as they are rather than relabelled with a newer commit.
>
> Where NT matters is the inverse shape -- many constraints, small blocks. On a bootstrap
> problem with m = 2439 and 17 blocks all of order <= 30 -- two iterations, two repeats per
> cell, full binary/input/parameter hashes recorded in `results/qd_nt_m2439_probe.tsv`
> (recipe repository) / `bench/qd_nt_m2439_probe.tsv` (solver repository) -- the pre-NT build runs
> ~93 s wall at 1 AND 4 threads with CPU ~= wall: one core, however many are offered. The
> NT build runs 93 s at 1 thread and 34 s at 4. Old-versus-new at 4 threads is therefore
> **2.72x, measured** (92.85 / 34.13 s on medians), and the NT build's own 1-to-4-thread
> scaling is 2.71x, as it must be when the old build ignores threads. CPU-seconds
> are flat (92-95 s) across all 8 runs -- same arithmetic, spread across cores -- and
> the printed solution section is byte-identical across all 8 (one sha256). An earlier
> one-sample version of this table reported the old build's 4-thread run as a 120 s timeout
> and old-versus-new as a ">3.5x" bound; neither survived repeats, and the current table
> supersedes both claims.
>
> Two precision statements, kept distinct on purpose: the *kernel* is bit-identical to the
> serial body by construction and by its own raw-limb test -- 135 NT cases at each of
> 1/4/8 threads (74 with the production gate open, 61 deliberately
> below it), with padded leading dimensions, every beta special case and 28960 structural
> zeros placed in B, every one memcmp-identical (`bench/qd_nt_kernel_check.cpp`, also
> `patches/qd_nt_kernel_check.cpp` in the recipe repository; results in
> `results/qd_nt_kernel_identity.tsv` / `bench/qd_nt_kernel_identity.tsv`). B5's earlier
> kernel test covered the NN case only and predates the NT commit, so it is not cited for
> NT. The *solver-level* checks above compare the complete printed solution section, which
> is identity at the output's printed precision, not a raw-limb claim.

Note what this is *not*: the published tables' `pristine` column was **upstream `766eef3`**,
a different comparison entirely. Nothing in this document is upstream. Every ratio here is
fork-against-fork, which is why the numbers are not comparable to the published ones
column by column and are given a dedicated section below instead.

## Methodology

External wall-clock seconds, **median of 3 repeats**, spread reported where it exceeds
rounding. Runs are pinned to physical cores on Linux (`taskset` + `OMP_PROC_BIND=true
OMP_PLACES=cores`, CPU set recorded per row); macOS exposes no per-process affinity, so those
runs use OpenMP binding plus one unrecorded warmup. The parameter file is passed explicitly
and its SHA-256 recorded. Every repeat is a row in the raw TSVs published alongside this
document; a run is counted only if it exited cleanly and every field parsed. Iteration counts
and objectives are checked across repeats and across configurations -- a wall-time ratio
between builds with different iteration counts measures path length, not speed, and is
flagged.

Solver-internal timers are elapsed time in sdpa-dd/sdpa-gmp; upstream sdpa-qd's timer
reports process CPU time (summed over threads), which this fork fixes -- all qd numbers here
use external wall time and the corrected clock.

Repeat counts, stated because they are **not** uniform and the table captions say so: 3 per
cell for four of the five published problems on both machines; **10** for `control1` on
thanos, which was re-measured behind an explicit 20 s all-core pre-warm after its first
configuration turned out to be reading the schedutil governor's ramp off a 1500 MHz idle
floor rather than the solver (0.64 / 0.60 / 0.45 s falling monotonically, against every later
configuration flat to the millisecond); and 2 per cell on the extension, where a median is a
mean of two and carries no dispersion information.

Added for this campaign, and absent from the previously published qd tables:

* **Termination phase is recorded per run.** A run that stops primal-feasible-only is not a
  solve, and an iterations-and-objective table cannot show the difference.
* **Foreign CPU load is sampled throughout and joined to each repeat's own time window.**
  This project has twice published a threading ratio that was really a measurement of
  contention -- 0.68x where the truth was 2.07x, and 1.07x then 0.59x from identical
  binaries. On thanos the worst per-row sibling-CPU load over the 60 measured repeats was
  1.39 cores and the largest deviation of any row from its own configuration's median was
  +1.09%; on pi the worst per-row foreign load was 0.050 cores. No row was discarded.
* **Threads actually created are sampled from `/proc` inside each repeat's window.** This is
  the observable that distinguishes a threaded binary whose parallel region never fires from
  a binary that has no parallel region to fire, and the published tables had no way to tell
  those apart.
* **Both arms are built from `git archive` of the named commit**, in isolated trees, on the
  machine they are measured on, and their md5s are checked to differ before the campaign and
  again after. The fork's working tree was dirty at generation time, including in one of the
  two files `443799f` touches; building from the object store rather than the checkout is
  what keeps the arms honest.

## thanos — EPYC 7232P, 8 physical cores

## thanos-epyc7232p — qd — `wall_s`

Median of 3-10 (NOT UNIFORM) repeats, seconds (`wall_s`). Spread = (max-min)/median.

| problem | m | base_t1 | base_t8 | head_t1 | head_t8 |
|---|---|---|---|---|---|
| control1 | 21 | 0.360 | 0.370 ±5% | 0.350 | 0.280 ±7% |
| gpp100 ‡ | 101 | 62.860 | 62.850 | 60.820 | 20.880 |
| theta1 | 104 | 5.740 ±1% | 5.500 ±1% | 5.530 ±1% | 3.400 |
| truss5 | 208 | 20.810 | 16.350 | 20.710 | 14.990 |
| arch0 | 174 | 260.570 | 168.110 | 252.120 | 79.080 |
| **total** | | **350.3** | **253.2** | **339.5** | **118.6** |

‡ = `base_t8` and `head_t8` take **different numbers of iterations** on this problem (gpp100 49 vs 50). A ratio between those two cells measures path length as well as speed; the per-iteration figure below is the one to quote. Counts are listed under Integrity.

**head_t8 vs base_t8, per iteration: 2.14x** (4.931 s -> 2.307 s per iteration, summed over the 5 of 5 problems with a determinate iteration count on both sides).
Geometric mean of the per-problem per-iteration ratios, which weights a 0.03 s problem equally with a 120 s one: **1.72x**.

End-to-end wall time totals 253.2 s -> 118.6 s (2.13x). **That ratio is not a speed statement**: the two builds disagree on iteration count on 1 of 5 problems (‡ above), so it mixes speed with path length.
Restricted to the 4 problem(s) where both builds take the *same* number of iterations, wall time is like-for-like: 190.3 s -> 97.8 s = **1.95x** (control1, theta1, truss5, arch0).

### Integrity

- gpp100: ITERATION COUNT differs BETWEEN configs (49, 50) -- compare per-iteration cost, not wall time

### Peak RSS (MB, max over repeats)

| problem | base_t1 | base_t8 | head_t1 | head_t8 |
|---|---|---|---|---|
| control1 | 11.4 | 11.4 | 11.3 | 11.3 |
| gpp100 | 11.4 | 11.3 | 11.4 | 11.3 |
| theta1 | 11.3 | 11.4 | 11.4 | 11.3 |
| truss5 | 11.4 | 11.3 | 11.4 | 11.3 |
| arch0 | 14.2 | 14.1 | 14.2 | 14.2 |

#### thanos — every ratio the two arms support

| comparison | `control1` | `gpp100` | `theta1` | `truss5` | `arch0` | total, same-iteration subset | total, all |
|---|---|---|---|---|---|---|---|
| pre-port fork, 1 -> 8 threads | 0.973x | 1.000x | 1.044x | 1.273x | 1.550x | **1.384x** | 1.384x |
| `443799f`, 1 -> 8 threads | 1.250x | 2.913x | 1.626x | 1.382x | 3.188x | **2.862x** | 2.862x |
| port payoff at 1 thread | 1.029x | 1.034x ‡ (per-it 1.055x) | 1.038x | 1.005x | 1.034x | **1.031x** | 1.032x |
| port payoff at 8 threads | 1.321x | 3.010x ‡ (per-it 3.071x) | 1.618x | 1.091x | 2.126x | **1.947x** | 2.134x |
| end to end: pre-port serial -> `443799f` at 8 threads | 1.286x | 3.011x ‡ (per-it 3.072x) | 1.688x | 1.388x | 3.295x | **2.941x** | 2.953x |

‡ the two configurations take **different numbers of iterations** on `gpp100`, so the wall ratio there mixes speed with path length; the per-iteration figure beside it is the one to quote, and the bolded total excludes it.

#### thanos — iterations, phase, objective

| problem | m | `base_t1` iters | `base_t8` iters | `head_t1` iters | `head_t8` iters | phase | objective identical across configs |
|---|---|---|---|---|---|---|---|
| `control1` | 21 | 46 | 46 | 46 | 46 | pdOPT | yes |
| `gpp100` ‡ | 101 | 49 | 49 | 50 | 50 | **pFEAS** | yes |
| `theta1` | 104 | 43 | 43 | 43 | 43 | pdOPT | yes |
| `truss5` | 208 | 48 | 48 | 48 | 48 | pdOPT | yes |
| `arch0` | 174 | 53 | 53 | 53 | 53 | pdOPT | yes |

‡ = the iteration count is **not** the same in every configuration of this problem (`gpp100`), so any wall ratio across that row measures path length as well as speed.

A phase other than `pdOPT` means the run stopped **without an optimal certificate** and is not a solve, however long it took: `gpp100` (pFEAS). The previously published qd tables recorded no phase and so could not show this.

#### thanos — threads the process actually created

| problem | `base_t1` | `base_t8` | `head_t1` | `head_t8` |
|---|---|---|---|---|
| `control1` | 1 | 1 | 1 | *n/s* |
| `gpp100` | 1 | 1 | 1 | 8 |
| `theta1` | 1 | 8 | 1 | 8 |
| `truss5` | 1 | 8 | 1 | 8 |
| `arch0` | 1 | 8 | 1 | 8 |

*n/s* = no sample fell wholly inside that repeat's window, so the cell is reported as unmeasured rather than given a number the data does not support (`control1`).

Every affected cell is on `control1`, which solves in well under a second — shorter than the sampler's own window. Its cells rest on far fewer samples than the other problems' (which carry tens to hundreds each) and should be read as corroboration, not as the census's evidence.

Two results in that census. First, **every 1-thread configuration has exactly one thread**,
so no idle OpenMP team is alive during a serial measurement — the mechanism that manufactured
12-14x phantom speedups on this box in an earlier campaign is absent by construction, not by
hope. Second, **at `OMP_NUM_THREADS=8`, pinned to eight distinct physical cores, the pre-port
binary stays single-threaded on `gpp100`**, which is why its `gpp100` column is flat at
62.86 s and
62.85 s — 1.000x from 1 thread to 8. That is the
published table's `gpp100` row, explained, and it is not what the published document assumed:
the problem does not resist threading, the binary had nothing to thread on that path. The
`control1` row says the same thing on the weaker evidence its short runtime allows.

### Headline, thanos

Over the 4 problems where both arms take the same number of iterations
(`control1`, `theta1`, `truss5`, `arch0`), at equal thread count:

* **port payoff at 8 threads: 1.947x** (190.3 s -> 97.8 s)
* **port payoff at 1 thread: 1.031x** — the control. These are threading commits;
  serially they must be neutral, and they are.
* **end to end, pre-port serial -> `443799f` at 8 threads: 2.941x**
  (287.5 s -> 97.8 s)

Thread scaling within each arm, over all five problems: the pre-port fork buys
**1.384x** from 1 to 8 threads, `443799f` buys **2.862x**. That gap
is the whole content of the port.

## pi — i9-13900K, 24 physical cores

pi is the only 24-core machine, and 24 threads is part of the determinism guarantee, so it
carries the 24-thread column. Its 8-thread column is pinned to the 8 P-cores
(`0,2,4,6,8,10,12,14`); the 24-thread column uses `0-31`, which the build phase verified
lands 24 threads on 24 distinct physical cores, 8 P + 16 E, never two on one core.

## pi-i9-13900k — qd — `wall_s`

Median of 3 repeats, seconds (`wall_s`). Spread = (max-min)/median.

| problem | m | base1 | head1 | base8 | head8 | base24 | head24 |
|---|---|---|---|---|---|---|---|
| control1 | 21 | 0.120 | 0.120 | 0.120 | 0.080 | 0.120 | 0.080 |
| gpp100 ‡ | 101 | 23.480 | 25.490 | 23.480 | 6.310 | 23.480 | 5.350 |
| theta1 | 104 | 1.860 ±1% | 1.810 ±1% | 1.760 ±1% | 0.960 | 1.770 ±1% | 0.960 ±1% |
| truss5 | 208 | 7.060 | 7.060 | 5.090 | 4.600 | 5.000 | 4.550 |
| arch0 | 174 | 87.290 | 85.390 | 51.420 | 19.430 ±1% | 49.440 | 15.590 |
| **total** | | **119.8** | **119.9** | **81.9** | **31.4** | **79.8** | **26.5** |

‡ = `base24` and `head24` take **different numbers of iterations** on this problem (gpp100 57 vs 64). A ratio between those two cells measures path length as well as speed; the per-iteration figure below is the one to quote. Counts are listed under Integrity.

**head24 vs base24, per iteration: 3.01x** (1.493 s -> 0.497 s per iteration, summed over the 5 of 5 problems with a determinate iteration count on both sides).
Geometric mean of the per-problem per-iteration ratios, which weights a 0.03 s problem equally with a 120 s one: **2.16x**.

End-to-end wall time totals 79.8 s -> 26.5 s (3.01x). **That ratio is not a speed statement**: the two builds disagree on iteration count on 1 of 5 problems (‡ above), so it mixes speed with path length.
Restricted to the 4 problem(s) where both builds take the *same* number of iterations, wall time is like-for-like: 56.3 s -> 21.2 s = **2.66x** (control1, theta1, truss5, arch0).

### Integrity

- gpp100: ITERATION COUNT differs BETWEEN configs (57, 64) -- compare per-iteration cost, not wall time

### Peak RSS (MB, max over repeats)

| problem | base1 | head1 | base8 | head8 | base24 | head24 |
|---|---|---|---|---|---|---|
| control1 | 2.7 | 2.7 | 2.7 | 2.7 | 2.7 | 2.7 |
| gpp100 | 7.4 | 7.4 | 7.4 | 7.4 | 7.4 | 7.4 |
| theta1 | 4.2 | 4.2 | 4.2 | 4.4 | 4.2 | 4.4 |
| truss5 | 6.2 | 6.2 | 6.2 | 6.2 | 6.2 | 7.9 |
| arch0 | 14.2 | 14.2 | 14.2 | 14.2 | 14.2 | 14.2 |

#### pi — every ratio the two arms support

| comparison | `control1` | `gpp100` | `theta1` | `truss5` | `arch0` | total, same-iteration subset | total, all |
|---|---|---|---|---|---|---|---|
| pre-port fork, 1 -> 8 threads | 1.000x † | 1.000x | 1.057x | 1.387x | 1.698x | **1.463x** | 1.463x |
| pre-port fork, 1 -> 24 threads | 1.000x † | 1.000x | 1.051x | 1.412x | 1.766x | **1.501x** | 1.501x |
| `443799f`, 1 -> 8 threads | 1.500x † | 4.040x | 1.885x | 1.535x | 4.395x | **3.820x** | 3.820x |
| `443799f`, 1 -> 24 threads | 1.500x † | 4.764x | 1.885x | 1.552x | 5.477x | **4.518x** | 4.518x |
| port payoff at 1 thread | 1.000x † | 0.921x ‡ (per-it 1.034x) | 1.028x | 1.000x | 1.022x | **1.021x** | 0.999x |
| port payoff at 8 threads | 1.500x † | 3.721x ‡ (per-it 4.178x) | 1.833x | 1.107x | 2.646x | **2.329x** | 2.609x |
| port payoff at 24 threads | 1.500x † | 4.389x ‡ (per-it 4.928x) | 1.844x | 1.099x | 3.171x | **2.660x** | 3.008x |
| end to end: pre-port serial -> `443799f` at 24 threads | 1.500x † | 4.389x ‡ (per-it 4.928x) | 1.938x | 1.552x | 5.599x | **4.548x** | 4.516x |

‡ the two configurations take **different numbers of iterations** on `gpp100`, so the wall ratio there mixes speed with path length; the per-iteration figure beside it is the one to quote, and the bolded total excludes it.

† **resolution-limited, do not quote as a speedup.** `wall_s` carries two decimals, and on `control1` at least one side of the ratio is small enough that one 0.01 s quantum is more than 5% of it — the ratio rests on a handful of clock ticks. These cells still enter the totals, where they are a negligible share of the time.

#### pi — iterations, phase, objective

| problem | m | `base1` iters | `base8` iters | `base24` iters | `head1` iters | `head8` iters | `head24` iters | phase | objective identical across configs |
|---|---|---|---|---|---|---|---|---|---|
| `control1` | 21 | 46 | 46 | 46 | 46 | 46 | 46 | pdOPT | yes |
| `gpp100` ‡ | 101 | 57 | 57 | 57 | 64 | 64 | 64 | **pFEAS** | yes |
| `theta1` | 104 | 43 | 43 | 43 | 43 | 43 | 43 | pdOPT | yes |
| `truss5` | 208 | 48 | 48 | 48 | 48 | 48 | 48 | pdOPT | yes |
| `arch0` | 174 | 53 | 53 | 53 | 53 | 53 | 53 | pdOPT | yes |

‡ = the iteration count is **not** the same in every configuration of this problem (`gpp100`), so any wall ratio across that row measures path length as well as speed.

A phase other than `pdOPT` means the run stopped **without an optimal certificate** and is not a solve, however long it took: `gpp100` (pFEAS). The previously published qd tables recorded no phase and so could not show this.

The pre-port arm's `gpp100` row reads 23.48 s at 1
thread, 23.48 s at 8 and
23.48 s at 24 — **three identical medians, 1.000x
scaling**. That is the thanos thread census reproduced on the other machine, on a different
compiler, without needing a census to see it: a binary that creates no second thread cannot
get slower or faster when told to use more.

### Headline, pi

Over the four problems where both arms take the same number of iterations, at equal thread
count:

* **port payoff at 24 threads: 2.660x** (56.3 s -> 21.2 s)
* **port payoff at 8 P-cores: 2.329x**
* **port payoff at 1 thread: 1.021x** — again the control, again neutral
* **end to end, pre-port serial -> `443799f` at 24 threads: 4.548x**
  (96.3 s -> 21.2 s)

Thread scaling within each arm, all five problems: pre-port **1.501x** from 1
to 24 threads, `443799f` **4.518x**.

## Extension — problems larger than the published set (pi)

The published set is entirely `m <= 208`, which the audit noted makes threading gains hard
to see. This section adds larger problems. It is an **extension**, not part of the set the
headline figures above are computed over, and it is quoted separately for that reason.

Repeats here are **2 per cell**, not 3. A median of 2 is a
mean of 2 with no dispersion information, so the spread column is not meaningful on this
table and no figure from it is quoted as a headline.

This campaign was launched before pi dropped off the network and kept running through the
outage under `nohup`; its rows were collected afterwards, unchanged and un-rerun, with the
binary md5s re-checked at both ends of the run. Every one of its 24 rows is
joined to its own window of the 1 Hz per-PID load trace
([`bench/qd_ext_pi_load.tsv`](bench/qd_ext_pi_load.tsv)): all 24 have a real
window (these are long runs), and the worst per-row **foreign** load over the whole campaign
is 0.050 cores of 24, contributed by `gnome-shell`.

## pi-i9-13900k — qd — `wall_s`

Median of 2 repeats, seconds (`wall_s`). Spread = (max-min)/median.

| problem | m | base1 | head1 | base8 | head8 | base24 | head24 |
|---|---|---|---|---|---|---|---|
| control5 | 351 | 128.895 | 127.185 | 34.155 | 32.595 | 27.455 | 25.825 ±1% |
| gpp250-1 ‡ | 251 | 304.750 | 341.400 | 304.785 | 72.715 | 304.770 | 59.565 |
| **total** | | **433.6** | **468.6** | **338.9** | **105.3** | **332.2** | **85.4** |

‡ = `base24` and `head24` take **different numbers of iterations** on this problem (gpp250-1 50 vs 58). A ratio between those two cells measures path length as well as speed; the per-iteration figure below is the one to quote. Counts are listed under Integrity.

**head24 vs base24, per iteration: 4.33x** (6.634 s -> 1.533 s per iteration, summed over the 2 of 2 problems with a determinate iteration count on both sides).
Geometric mean of the per-problem per-iteration ratios, which weights a 0.03 s problem equally with a 120 s one: **2.51x**.

End-to-end wall time totals 332.2 s -> 85.4 s (3.89x). **That ratio is not a speed statement**: the two builds disagree on iteration count on 1 of 2 problems (‡ above), so it mixes speed with path length.
Restricted to the 1 problem(s) where both builds take the *same* number of iterations, wall time is like-for-like: 27.5 s -> 25.8 s = **1.06x** (control5).

### Integrity

- gpp250-1: ITERATION COUNT differs BETWEEN configs (50, 58) -- compare per-iteration cost, not wall time

### Peak RSS (MB, max over repeats)

| problem | base1 | head1 | base8 | head8 | base24 | head24 |
|---|---|---|---|---|---|---|
| control5 | 9.4 | 9.4 | 9.7 | 9.7 | 9.7 | 9.7 |
| gpp250-1 | 31.8 | 31.8 | 31.7 | 31.7 | 31.7 | 31.7 |

#### extension — every ratio the two arms support

| comparison | `control5` | `gpp250-1` | total, same-iteration subset | total, all |
|---|---|---|---|---|
| pre-port fork, 1 -> 24 threads | 4.695x | 1.000x | **1.305x** | 1.305x |
| `443799f`, 1 -> 24 threads | 4.925x | 5.732x | **5.488x** | 5.488x |
| port payoff at 8 threads | 1.048x | 4.192x ‡ (per-it 4.862x) | **1.048x** | 3.218x |
| port payoff at 24 threads | 1.063x | 5.117x ‡ (per-it 5.935x) | **1.063x** | 3.891x |
| end to end: pre-port serial -> `443799f` at 24 threads | 4.991x | 5.116x ‡ (per-it 5.935x) | **4.991x** | 5.078x |

‡ the two configurations take **different numbers of iterations** on `gpp250-1`, so the wall ratio there mixes speed with path length; the per-iteration figure beside it is the one to quote, and the bolded total excludes it.

#### extension — iterations, phase, objective

| problem | m | `base1` iters | `base8` iters | `base24` iters | `head1` iters | `head8` iters | `head24` iters | phase | objective identical across configs |
|---|---|---|---|---|---|---|---|---|---|
| `control5` | 351 | 51 | 51 | 51 | 51 | 51 | 51 | pdOPT | yes |
| `gpp250-1` ‡ | 251 | 50 | 50 | 50 | 58 | 58 | 58 | **pFEAS** | yes |

‡ = the iteration count is **not** the same in every configuration of this problem (`gpp250-1`), so any wall ratio across that row measures path length as well as speed.

A phase other than `pdOPT` means the run stopped **without an optimal certificate** and is not a solve, however long it took: `gpp250-1` (pFEAS). The previously published qd tables recorded no phase and so could not show this.

### What the larger problems actually show

| problem | m | pre-port, 1 thread | pre-port 1 -> 24 threads | port payoff at 24 threads, wall | per iteration |
|---|---|---|---|---|---|
| `control5` | 351 | 128.89 s | 4.695x | 1.063x | same path |
| `gpp250-1` | 251 | 304.75 s | 1.000x | 5.117x ‡ | **5.935x** (50 vs 58 iterations) |

‡ = iteration counts differ between the arms on this problem, so the wall figure is not a speed statement and the per-iteration column is the one to quote.

This is not the result adding larger problems was expected to produce. The audit's point was
that the published set (`m <= 208`) is too small for threading gains to show. Size is not
what separates these two. **`control5` (m=351) scales 4.695x from 1 to 24 threads in the pre-port
arm; `gpp250-1` (m=251) scales 1.000x — that is, not at all** — and the pre-port arm's only parallel
region is the `k1` bMat loop, whose work is the Schur complement. Where that loop dominates,
the pre-port fork already threads well and the port's marginal payoff is small; where it does
not, the pre-port fork does not thread at all and the port's payoff is the whole of the
speedup. The published set contains one of each (`arch0` at 1.766x pre-port scaling, `gpp100`
at 1.000x), so the effect was visible there too; the larger problems make it unmistakable and
show that the discriminant is problem *structure*, not problem size.

## `gpp100`

`gpp100` is the one problem in the published set whose iteration count is not identical
between the two arms, on either machine. Every wall ratio on that row therefore mixes speed
with path length and is quoted per iteration:

| | pre-port | `443799f` | wall | per iteration |
|---|---|---|---|---|
| thanos, 1 thread | 62.86 s / 49 it | 60.82 s / 50 it | 1.034x | **1.055x** |
| thanos, 8 threads | 62.85 s / 49 it | 20.88 s / 50 it | 3.010x | **3.071x** |
| pi, 1 thread | 23.48 s / 57 it | 25.49 s / 64 it | 0.921x | **1.034x** |
| pi, 24 threads | 23.48 s / 57 it | 5.35 s / 64 it | 4.389x | **4.928x** |

The pi 1-thread line is the reason this matters. `443799f` is **3.4%
faster per iteration** there and takes **12% more
iterations**, so its wall time looks 9% *worse*. Quoting
0.921x as a slowdown would repeat exactly the error the old 49/63 row made
in the other direction.

Both arms stop at `pFEAS` on both machines — **none of these trajectories is a
converged solve**: not 49 or 50 iterations on thanos, not
57 or 64 on pi. The objective nevertheless agrees to all 17
printed digits in every configuration on both machines.

### The iteration count is a property of the build, not of the source

#### Every iteration count this project has recorded for `gpp100`

| source file | machine | config | iterations | phase | objective |
|---|---|---|---|---|---|
| `qd_v2_thanos.tsv` — published, thanos | thanos-epyc7232p | `hardened1` | 49 | not recorded | -4.4943550775891147e+01 |
| `qd_v2_thanos.tsv` — published, thanos | thanos-epyc7232p | `hardened4` | 49 | not recorded | -4.4943550775891147e+01 |
| `qd_v2_thanos.tsv` — published, thanos | thanos-epyc7232p | `hardened8` | 49 | not recorded | -4.4943550775891147e+01 |
| `qd_v2_thanos.tsv` — published, thanos | thanos-epyc7232p | `pristine` | 63 | not recorded | -4.4943550775891147e+01 |
| `qd_v2_pi.tsv` — published, pi | pi-i9-13900k | `fork1` | 49 | not recorded | -4.4943550775891147e+01 |
| `qd_v2_pi.tsv` — published, pi | pi-i9-13900k | `fork24` | 49 | not recorded | -4.4943550775891147e+01 |
| `qd_v2_pi.tsv` — published, pi | pi-i9-13900k | `fork8P` | 49 | not recorded | -4.4943550775891147e+01 |
| `qd_v2_pi.tsv` — published, pi | pi-i9-13900k | `pristine` | 49 | not recorded | -4.4943550775891147e+01 |
| `qd_v2_mac.tsv` — published, Mac | mac-m1max | `fork1` | 56 | not recorded | -4.4943550775891147e+01 |
| `qd_v2_mac.tsv` — published, Mac | mac-m1max | `fork8` | 56 | not recorded | -4.4943550775891147e+01 |
| `qd_v2_mac.tsv` — published, Mac | mac-m1max | `pristine` | 56 | not recorded | -4.4943550775891147e+01 |
| `qd_v3_thanos_fixed.tsv` — post-`63370f0` check, thanos | thanos-epyc7232p | `fork1` | 52 | not recorded | -4.4943550775891147e+01 |
| `qd_v3_thanos_fixed.tsv` — post-`63370f0` check, thanos | thanos-epyc7232p | `fork8` | 52 | not recorded | -4.4943550775891147e+01 |
| `qd_v3_pi_fixed.tsv` — post-`63370f0` check, pi | pi-i9-13900k | `fork1` | 52 | not recorded | -4.4943550775891147e+01 |
| `qd_v3_pi_fixed.tsv` — post-`63370f0` check, pi | pi-i9-13900k | `fork24` | 52 | not recorded | -4.4943550775891147e+01 |
| `qd_v3_pi_fixed.tsv` — post-`63370f0` check, pi | pi-i9-13900k | `fork8P` | 52 | not recorded | -4.4943550775891147e+01 |
| `qd_v3_pi_bench.tsv` — post-`63370f0` bench, pi | pi-i9-13900k | `fixed1` | 52 | not recorded | -4.4943550775891147e+01 |
| `qd_v3_pi_bench.tsv` — post-`63370f0` bench, pi | pi-i9-13900k | `fixed24` | 52 | not recorded | -4.4943550775891147e+01 |
| `qd_v3_pi_bench.tsv` — post-`63370f0` bench, pi | pi-i9-13900k | `fixed8P` | 52 | not recorded | -4.4943550775891147e+01 |
| `qd_v3_pi_bench.tsv` — post-`63370f0` bench, pi | pi-i9-13900k | `published` | 60 | not recorded | -4.4943550775891147e+01 |
| `qd_regen_thanos_v2.tsv` — this campaign, thanos | thanos-epyc7232p | `base_t1` | 49 | pFEAS | -4.4943550775891147e+01 |
| `qd_regen_thanos_v2.tsv` — this campaign, thanos | thanos-epyc7232p | `base_t8` | 49 | pFEAS | -4.4943550775891147e+01 |
| `qd_regen_thanos_v2.tsv` — this campaign, thanos | thanos-epyc7232p | `head_t1` | 50 | pFEAS | -4.4943550775891147e+01 |
| `qd_regen_thanos_v2.tsv` — this campaign, thanos | thanos-epyc7232p | `head_t8` | 50 | pFEAS | -4.4943550775891147e+01 |
| `qd_v4_pi_bench.tsv` — this campaign, pi | pi-i9-13900k | `base1` | 57 | pFEAS | -4.4943550775891147e+01 |
| `qd_v4_pi_bench.tsv` — this campaign, pi | pi-i9-13900k | `base24` | 57 | pFEAS | -4.4943550775891147e+01 |
| `qd_v4_pi_bench.tsv` — this campaign, pi | pi-i9-13900k | `base8` | 57 | pFEAS | -4.4943550775891147e+01 |
| `qd_v4_pi_bench.tsv` — this campaign, pi | pi-i9-13900k | `head1` | 64 | pFEAS | -4.4943550775891147e+01 |
| `qd_v4_pi_bench.tsv` — this campaign, pi | pi-i9-13900k | `head24` | 64 | pFEAS | -4.4943550775891147e+01 |
| `qd_v4_pi_bench.tsv` — this campaign, pi | pi-i9-13900k | `head8` | 64 | pFEAS | -4.4943550775891147e+01 |

This is the refutation, and it is the main thing this regeneration found that was not
expected. The claim carried into the campaign was that `63370f0` had settled `gpp100` at
**a stable 52 on both machines**.

That claim was **true of the builds it was measured on**, and the survey above says so:
`52` on thanos and `52` on pi in the post-fix checks, at every
thread count, in both files. It is **not** true of `443799f`:

* at `443799f`, thanos gives 50 and pi gives
  64. The count **moved again** between the post-fix
  build and HEAD, and moved *differently on each machine* — so 52 was a property of one pair
  of builds, not a constant of the fixed source, and a document that printed "now a stable 52"
  would be wrong about the very binary it was describing;
* across every build this project has measured, `gpp100` has taken
  49, 50, 52, 56, 57, 60, 63, 64 iterations — 8 distinct values — all converging to the same
  objective, and stopping `pFEAS` in every run whose phase was recorded (the published
  tables recorded none, which is its own finding);
* the sharpest case is pi, where **the source did not change at all**. The published `fork1`
  column recorded 49 iterations at 0.4008 s each; this
  campaign's freshly built pre-port arm on the same machine, with the same parameter file,
  records 57 at 0.4119 s each — **2.8%
  apart in cost per iteration and 16%
  apart in path length.** The build phase found the mechanism: pi carries two entirely
  different QD library builds (`/usr/local/lib/libqd.a` at 408646 B against another at
  4002406 B), and three qd binaries already on that box give 49, 60 and 64 iterations on this
  problem with the same parameter file.

So the published `49` and `63` cannot be attributed to the uninitialised read on the evidence
available. Those two numbers came from **two different source trees** — `63` was upstream
`766eef3`, `49` was the fork — and every build measured since has produced its own value from
the same source. A trajectory that moves with the arithmetic library is not a trajectory that
identifies one defect.

What `63370f0` **did** do is real and is worth stating precisely: within a machine and within
a build, the count is deterministic — constant across repeats *and* across 1, 8 and 24
threads, in both arms — and the fix moved the count on both machines
(49 -> 50
on thanos, 57 ->
64 on pi), which is what fixing an uninitialised
read that feeds a block-size lookup should do. What it did not do is make the count a
constant of the source. A second observation to record rather than smooth over: the
pre-`63370f0` arm was **not** observed to be unstable here. Both arms are perfectly
deterministic across repeats and thread counts. That is what an uninitialised read does when
the value it picks up happens to be stable for a fixed binary and environment; "stable here"
is not "not a bug", and it means this observable cannot demonstrate the fix on these two
machines.

## Every published qd number, against its re-measurement

#### thanos — published against re-measured

`pristine` was **upstream** `766eef3`; nothing in the new campaign is upstream, so it has no counterpart and is listed with none rather than paired with something it is not. `hardened*` was the fork at `3376c79`, which is exactly what `base_*` is, so those lines are like-for-like.

| problem | published config | published s | new config | new s | change | published iters | new iters |
|---|---|---|---|---|---|---|---|
| `control1` | `pristine` | 0.560 | *(no counterpart)* | — | — | 46 | — |
| `gpp100` | `pristine` | 80.820 | *(no counterpart)* | — | — | 63 | — |
| `theta1` | `pristine` | 5.800 | *(no counterpart)* | — | — | 43 | — |
| `truss5` | `pristine` | 21.270 | *(no counterpart)* | — | — | 48 | — |
| `arch0` | `pristine` | 264.300 | *(no counterpart)* | — | — | 53 | — |
| `control1` | `hardened1` | 0.360 | `base_t1` | 0.360 | +0.0% | 46 | 46 |
| `gpp100` | `hardened1` | 63.290 | `base_t1` | 62.860 | -0.7% | 49 | 49 |
| `theta1` | `hardened1` | 5.810 | `base_t1` | 5.740 | -1.2% | 43 | 43 |
| `truss5` | `hardened1` | 20.990 | `base_t1` | 20.810 | -0.9% | 48 | 48 |
| `arch0` | `hardened1` | 262.840 | `base_t1` | 260.570 | -0.9% | 53 | 53 |
| `control1` | `hardened4` | 0.360 | *(no counterpart)* | — | — | 46 | — |
| `gpp100` | `hardened4` | 63.240 | *(no counterpart)* | — | — | 49 | — |
| `theta1` | `hardened4` | 5.570 | *(no counterpart)* | — | — | 43 | — |
| `truss5` | `hardened4` | 17.480 | *(no counterpart)* | — | — | 48 | — |
| `arch0` | `hardened4` | 188.470 | *(no counterpart)* | — | — | 53 | — |
| `control1` | `hardened8` | 0.370 | `base_t8` | 0.370 | +0.0% | 46 | 46 |
| `gpp100` | `hardened8` | 63.270 | `base_t8` | 62.850 | -0.7% | 49 | 49 |
| `theta1` | `hardened8` | 5.540 | `base_t8` | 5.500 | -0.7% | 43 | 43 |
| `truss5` | `hardened8` | 16.480 | `base_t8` | 16.350 | -0.8% | 48 | 48 |
| `arch0` | `hardened8` | 171.840 | `base_t8` | 168.110 | -2.2% | 53 | 53 |

#### pi — published against re-measured

The published `pristine` column on pi was run at **one thread** (`threads=1`, `cpuset=0`), so the published headline `fork24 vs pristine` was a 24-thread column over a 1-thread column. It has no counterpart here either.

| problem | published config | published s | new config | new s | change | published iters | new iters |
|---|---|---|---|---|---|---|---|
| `control1` | `pristine` | 0.110 | *(no counterpart)* | — | — | 46 | — |
| `gpp100` | `pristine` | 19.640 | *(no counterpart)* | — | — | 49 | — |
| `theta1` | `pristine` | 1.800 | *(no counterpart)* | — | — | 43 | — |
| `truss5` | `pristine` | 6.890 | *(no counterpart)* | — | — | 48 | — |
| `arch0` | `pristine` | 84.590 | *(no counterpart)* | — | — | 53 | — |
| `control1` | `fork1` | 0.110 | `base1` | 0.120 | +9.1% | 46 | 46 |
| `gpp100` | `fork1` | 19.640 | `base1` | 23.480 | +19.6% | 49 | 57 ‡ |
| `theta1` | `fork1` | 1.800 | `base1` | 1.860 | +3.3% | 43 | 43 |
| `truss5` | `fork1` | 6.820 | `base1` | 7.060 | +3.5% | 48 | 48 |
| `arch0` | `fork1` | 84.530 | `base1` | 87.290 | +3.3% | 53 | 53 |
| `control1` | `fork8P` | 0.110 | `base8` | 0.120 | +9.1% | 46 | 46 |
| `gpp100` | `fork8P` | 19.640 | `base8` | 23.480 | +19.6% | 49 | 57 ‡ |
| `theta1` | `fork8P` | 1.710 | `base8` | 1.760 | +2.9% | 43 | 43 |
| `truss5` | `fork8P` | 4.930 | `base8` | 5.090 | +3.2% | 48 | 48 |
| `arch0` | `fork8P` | 49.790 | `base8` | 51.420 | +3.3% | 53 | 53 |
| `control1` | `fork24` | 0.110 | `base24` | 0.120 | +9.1% | 46 | 46 |
| `gpp100` | `fork24` | 19.640 | `base24` | 23.480 | +19.6% | 49 | 57 ‡ |
| `theta1` | `fork24` | 1.710 | `base24` | 1.770 | +3.5% | 43 | 43 |
| `truss5` | `fork24` | 4.840 | `base24` | 5.000 | +3.3% | 48 | 48 |
| `arch0` | `fork24` | 47.800 | `base24` | 49.440 | +3.4% | 53 | 53 |

‡ = the iteration count moved, so the wall-time change on that line is not a speed change.

#### thanos — peak RSS, published against re-measured

| problem | `hardened1` -> `base_t1` | `hardened8` -> `base_t8` |
|---|---|---|
| `control1` | 11.7 -> 11.4 MB (-2.6%) | 11.7 -> 11.4 MB (-2.6%) |
| `gpp100` | 11.7 -> 11.4 MB (-2.6%) | 11.7 -> 11.3 MB (-3.4%) |
| `theta1` | 11.7 -> 11.3 MB (-3.4%) | 11.7 -> 11.4 MB (-2.6%) |
| `truss5` | 11.8 -> 11.4 MB (-3.4%) | 11.7 -> 11.3 MB (-3.4%) |
| `arch0` | 14.2 -> 14.2 MB (+0.0%) | 14.1 -> 14.1 MB (+0.0%) |

#### pi — peak RSS, published against re-measured

| problem | `fork1` -> `base1` | `fork8P` -> `base8` | `fork24` -> `base24` |
|---|---|---|---|
| `control1` | 2.7 -> 2.7 MB (+0.0%) | 2.7 -> 2.7 MB (+0.0%) | 2.7 -> 2.7 MB (+0.0%) |
| `gpp100` | 7.4 -> 7.4 MB (+0.0%) | 7.4 -> 7.4 MB (+0.0%) | 7.4 -> 7.4 MB (+0.0%) |
| `theta1` | 4.2 -> 4.2 MB (+0.0%) | 4.4 -> 4.2 MB (-4.5%) | 4.2 -> 4.2 MB (+0.0%) |
| `truss5` | 6.2 -> 6.2 MB (+0.0%) | 6.2 -> 6.2 MB (+0.0%) | 6.2 -> 6.2 MB (+0.0%) |
| `arch0` | 14.2 -> 14.2 MB (+0.0%) | 14.2 -> 14.2 MB (+0.0%) | 14.2 -> 14.2 MB (+0.0%) |

Peak RSS moved by at most 4.5% on any cell of either
machine, in both directions, on problems whose footprint is 2-14 MB. Nothing here is a
memory result; it is recorded because "every number that changed" includes the ones that
barely did.

### The published headline figures

| published claim | recomputed from its own TSV | status now |
|---|---|---|
| thanos `hardened8 vs pristine: 1.45x` | 372.8 s -> 257.5 s = 1.448x | superseded — `pristine` is upstream `766eef3`, which this campaign has no arm for |
| thanos "the honest headline is 1.50x" (gpp100 excluded) | 1.503x | same caveat; the fork-against-fork figure that replaces it is 1.947x at 8 threads |
| pi `fork24 vs pristine: 1.53x` | 113.0 s -> 74.1 s = 1.525x | **mis-framed.** That `pristine` column was run at **one thread**. `pristine` against `fork1`, both serial, is 1.001x — essentially none of the 1.53x was the fork. The fork's own 1 -> 24 thread scaling accounts for it: 1.524x |

The pi line is the one to read twice. The published document put a 24-thread column against
a 1-thread column and labelled the ratio as what the fork bought. Re-measured as an
intra-arm number on the same commit, the pre-port fork scales 1.501x from 1 to
24 threads on this set, and that is the whole of it.

### What the two porting commits claimed for themselves

Both figures below are quoted from the commit messages, and this document refuses to print
them unless they are still in those messages. Neither is the same measurement as the
regenerated one, and the differences are recorded as differences rather than reconciled.

| claim | where | this campaign |
|---|---|---|
| `time-weighted @8    1.80x` on the published set, thanos | `f177644` (`Rgemm` NN port alone) | 1.947x same-iteration subset / 2.134x all five, for `f177644`+`443799f` together |
| `gpp100 2.48x (best case)` at 8 threads, thanos | `f177644` alone | 3.010x wall, 3.071x per iteration |
| B3 payoff 1.220x at 8 threads and 1.287x at 24, pi | `results/qd_b3_payoff_pi.tsv` | **not separable here** — no arm was built at `f177644`, so this campaign cannot isolate `443799f` |

The `head` arm carries both commits, so a like-for-like check is only available as a
consistency argument: 1.80x for `f177644` alone against 1.220x claimed for `443799f` alone
multiplies to 2.20x, and the combined measurement over all five problems on thanos is
2.134x. Those are different aggregations over different
runs and their agreement is a sanity check, not a confirmation.

## Carried over, NOT re-measured: Apple M1 Max

The table below is the **pre-port fork** and predates every commit this document is about.
`443799f` has never been benchmarked on the Mac; its `fork8` column here is the same
pre-port binary as `base*` above, under a different name, and must not be read as HEAD.
It is retained because deleting measured data to tidy a document is its own distortion.

## mac-m1max — qd — `wall_s`

Median of 3 repeats, seconds (`wall_s`). Spread = (max-min)/median.

| problem | m | pristine | fork1 | fork8 |
|---|---|---|---|---|
| control1 | 21 | 0.230 ±4% | 0.230 | 0.230 |
| gpp100 | 101 | 22.520 | 22.520 | 22.520 |
| theta1 | 104 | 2.350 ±1% | 2.330 | 2.230 |
| truss5 | 208 | 9.430 | 9.230 | 7.130 |
| arch0 | 174 | 89.110 | 89.040 | 51.360 |
| **total** | | **123.6** | **123.4** | **83.5** |

**fork8 vs pristine: 1.48x**  (totals 123.6 s -> 83.5 s)

Iteration counts are identical between these two configurations on all 5 problems, so this is a like-for-like speed ratio.

### Integrity

- all repeats `ok`; iteration count and objective identical across repeats and across configs for every problem

### Peak RSS (MB, max over repeats)

| problem | pristine | fork1 | fork8 |
|---|---|---|---|
| control1 | 8.5 | 8.5 | 8.5 |
| gpp100 | 8.5 | 8.5 | 8.5 |
| theta1 | 8.5 | 8.5 | 8.5 |
| truss5 | 8.5 | 8.5 | 8.6 |
| arch0 | 14.8 | 15.0 | 15.2 |


On the Mac, `gpp100`'s iteration count was identical in every configuration measured (56) --
a fourth distinct value for the same problem, and further evidence for the build-dependence
argued above rather than against it.

## What is still missing

* **No upstream arm.** Every ratio here is fork-against-fork. The published `pristine`
  columns are upstream `766eef3` and are not comparable to anything measured here.
* **No Mac arm at `443799f`.**
* **`f177644` and `443799f` are not separated.** The `head` arm carries both, and the
  15 commits between `3376c79` and `443799f` besides. The parallel gain is attributable to
  those two commits — they are the only ones that add OpenMP translation units, and the
  1-thread aggregate is neutral to within 0.1% on pi — but isolating each one needs a third
  arm at `f177644`, which was not built.
* **The published set is small.** All five problems are `m <= 208`, which is why the
  extension section exists; the extension itself is two problems at 2 repeats and is not
  a substitute for measuring the larger corpus properly. It is also pi-only — the same two
  problems were not run on thanos, where a serial `gpp250-1` run alone costs ~5 minutes and
  the four configurations would have needed roughly an hour of exclusive machine time.
* **`gpp100`'s machine-to-machine iteration difference is uncontrolled by this campaign.**
  Both thanos arms link the same `libqd.a` (`40ebea2292b78afc344146ec9e410731`, no `QD_FMA`),
  so the 49 -> 50 shift there is attributable to the source diff; the thanos-against-pi
  difference is not, and separating it needs a matched-`libqd` experiment that was not run.

Raw data: [`bench/qd_regen_thanos_v2.tsv`](bench/qd_regen_thanos_v2.tsv) (thanos, with the
pre-warmed `control1` rows spliced in),
[`bench/qd_regen_thanos.tsv`](bench/qd_regen_thanos.tsv) (thanos as measured),
[`bench/qd_regen_thanos_load.tsv`](bench/qd_regen_thanos_load.tsv) (per-repeat foreign load
and thread census), [`bench/qd_v4_pi_bench.tsv`](bench/qd_v4_pi_bench.tsv) and
[`bench/qd_v4_pi_bench_load.tsv`](bench/qd_v4_pi_bench_load.tsv) (pi),
[`bench/qd_ext_pi.tsv`](bench/qd_ext_pi.tsv) and
[`bench/qd_ext_pi_load.tsv`](bench/qd_ext_pi_load.tsv) (the extension),
[`bench/qd_build_provenance_thanos.tsv`](bench/qd_build_provenance_thanos.tsv) (the six
builds behind the inertness table). Superseded but retained:
[`bench/qd_v2_thanos.tsv`](bench/qd_v2_thanos.tsv),
[`bench/qd_v2_pi.tsv`](bench/qd_v2_pi.tsv), [`bench/qd_v2_mac.tsv`](bench/qd_v2_mac.tsv).
