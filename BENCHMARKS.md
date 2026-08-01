# Benchmarks — sdpa-qd-omp

Fork base: upstream `766eef3`. thanos: AMD EPYC 7232P, 8 physical cores, Ubuntu.

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

## Results

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

Raw data: [`bench/qd_v2_thanos.tsv`](bench/qd_v2_thanos.tsv).
