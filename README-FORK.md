# sdpa-qd-omp

OpenMP-threaded fork of sdpa-qd, with its benchmark clock fixed.

Fork of [nakatamaho/sdpa-qd](https://github.com/nakatamaho/sdpa-qd) at `766eef3`, carrying two
patches. Reported upstream; not adopted there.

1. **Timer fix.** `Time::rGetUseTime()` measured process CPU time summed over threads, so any
   parallel speedup *reported itself as a slowdown* (arch0: an 84.6 → 49.8 s run printed 90.9).
   Now `std::chrono::steady_clock` elapsed time. Anyone who ever benchmarked threaded sdpa-qd
   against its own output should re-check their conclusions.
2. **Schur-complement threading** — the first OpenMP in this codebase (upstream has none, and
   its `--enable-openmp` flags are vestigial: they do not add `-fopenmp`).

**Measured** (EPYC 7232P, external wall clock, median of 3 pinned repeats): **1.50×** over the
four problems with unchanged trajectories (arch0 264.3 → 171.8 s). `gpp100` excluded honestly:
the patch shortens its path (63 → 49 iterations, same objective), so per iteration it is 1.00×.

Every modified file carries an in-file, dated change notice (GPLv2 §2a): `sdpa_tool.cpp`,
`sdpa_newton.cpp`, `sdpa_parts.cpp`.

Build (QD library required, e.g. `libqd-dev`; `-fopenmp` must be passed explicitly):

```bash
./configure --with-qd-includedir=/usr/include --with-qd-libdir=/usr/lib/x86_64-linux-gnu \
            CXXFLAGS='-O2 -funroll-all-loops -fopenmp' CFLAGS='-O2 -funroll-all-loops -fopenmp' \
            LDFLAGS='-fopenmp'
make -j$(nproc)
nm sdpa_qd | grep -c GOMP     # the binary is static; ldd cannot verify OpenMP
```

License: GPL v2, unchanged (`COPYING`); original SDPA authors retain copyright.
