# Benchmarks — sdpa-qd-omp

Fork base: upstream `766eef3`.

**Read this first.** Upstream sdpa-qd's `main loop time` is process CPU time summed over
threads -- it *rises* with thread count even when the run gets faster, so any historical
threaded benchmark read from sdpa-qd's own output is suspect. This fork fixes the clock;
everything below is external wall time.

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

## thanos — EPYC 7232P, 8 physical cores

## thanos-epyc7232p — qd — `wall_s`

Median of 3 repeats, seconds (`wall_s`). Spread = (max-min)/median.

| problem | m | pristine | hardened1 | hardened4 | hardened8 |
|---|---|---|---|---|---|
| control1 | 21 | 0.560 ±45% | 0.360 ±3% | 0.360 | 0.370 ±3% |
| gpp100 | 101 | 80.820 | 63.290 | 63.240 | 63.270 |
| theta1 | 104 | 5.800 ±1% | 5.810 ±1% | 5.570 | 5.540 ±1% |
| truss5 | 208 | 21.270 | 20.990 | 17.480 | 16.480 |
| arch0 | 174 | 264.300 ±2% | 262.840 | 188.470 | 171.840 |
| **total** | | **372.8** | **353.3** | **275.1** | **257.5** |

**hardened8 vs pristine: 1.45x**  (totals 372.8 s -> 257.5 s)

Iteration counts differ on 1 problem(s) (gpp100), so the wall-time ratio above mixes speed with path length. **Per iteration: 1.37x** (6.860 s -> 5.014 s).

### Integrity

- gpp100: ITERATION COUNT differs BETWEEN configs (49, 63) -- compare per-iteration cost, not wall time

### Peak RSS (MB, max over repeats)

| problem | pristine | hardened1 | hardened4 | hardened8 |
|---|---|---|---|---|
| control1 | 11.8 | 11.7 | 11.8 | 11.7 |
| gpp100 | 11.7 | 11.7 | 11.8 | 11.7 |
| theta1 | 11.7 | 11.7 | 11.7 | 11.7 |
| truss5 | 11.7 | 11.8 | 11.8 | 11.7 |
| arch0 | 13.7 | 14.2 | 14.1 | 14.1 |


### The honest headline is 1.50x, not the table total

`gpp100` is the one problem where the patch is numerically non-neutral: upstream takes
63 iterations, the patched build 49, converging to the
same objective. Its wall-time gain is a shorter path, not faster execution -- per iteration
the two builds are within 1%. Excluding it, the four problems with identical trajectories
total 291.9 s -> 194.2 s = **1.50x**, concentrated in
the large F3-dominated `arch0` (1.54x). Why quad-double shifts one trajectory when
double-double shifts none (20/20 trajectory-neutral) is an open question, documented in the
companion repository rather than hidden.

## pi — i9-13900K, 24 physical cores (fork binary built per the README, fresh clone)

## pi-i9-13900k — qd — `wall_s`

Median of 3 repeats, seconds (`wall_s`). Spread = (max-min)/median.

| problem | m | pristine | fork1 | fork8P | fork24 |
|---|---|---|---|---|---|
| control1 | 21 | 0.110 | 0.110 | 0.110 | 0.110 |
| gpp100 | 101 | 19.640 | 19.640 | 19.640 | 19.640 |
| theta1 | 104 | 1.800 | 1.800 | 1.710 | 1.710 |
| truss5 | 208 | 6.890 | 6.820 | 4.930 | 4.840 |
| arch0 | 174 | 84.590 | 84.530 | 49.790 | 47.800 |
| **total** | | **113.0** | **112.9** | **76.2** | **74.1** |

**fork24 vs pristine: 1.53x**  (totals 113.0 s -> 74.1 s)

### Integrity

- all repeats `ok`; iteration count and objective identical across repeats and across configs for every problem

### Peak RSS (MB, max over repeats)

| problem | pristine | fork1 | fork8P | fork24 |
|---|---|---|---|---|
| control1 | 2.4 | 2.7 | 2.7 | 2.7 |
| gpp100 | 7.2 | 7.4 | 7.4 | 7.4 |
| theta1 | 3.9 | 4.2 | 4.4 | 4.2 |
| truss5 | 5.7 | 6.2 | 6.2 | 6.2 |
| arch0 | 13.9 | 14.2 | 14.2 | 14.2 |


## Mac — Apple M1 Max (fork binary built by this README's macOS instructions)

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


On the Mac, `gpp100`'s iteration count is identical in every configuration (56), so the
whole table is clean like-for-like there.

### The `gpp100` trajectory, measured on three platforms

| | EPYC 7232P | i9-13900K | M1 Max |
|---|---|---|---|
| pristine upstream, iterations | 63 | 49 | 56 |
| this fork, iterations | 49 | 49 | 56 |

*Pristine upstream by itself* takes three different paths on three microarchitectures — the
problem is trajectory-fragile to low-bit differences of any origin (platform, QD library
build, or this patch), and the one fork-induced shift (thanos) lands exactly on another
platform's pristine value. Same objective everywhere. Consequence: on this problem,
wall-time ratios measure path length as much as speed; per-iteration cost is the stable
metric, and by it the fork is neutral on `gpp100` on every platform measured.

Raw data: [`bench/qd_v2_thanos.tsv`](bench/qd_v2_thanos.tsv), [`bench/qd_v2_pi.tsv`](bench/qd_v2_pi.tsv), [`bench/qd_v2_mac.tsv`](bench/qd_v2_mac.tsv).
