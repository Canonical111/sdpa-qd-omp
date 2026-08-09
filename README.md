# sdpa-qd-omp

OpenMP-threaded fork of sdpa-qd, with its benchmark clock fixed.

Fork of [nakatamaho/sdpa-qd](https://github.com/nakatamaho/sdpa-qd) (upstream README preserved as [README-UPSTREAM.md](README-UPSTREAM.md)) at `766eef3`, carrying four
patches. Reported upstream; not adopted there.

1. **Timer fix.** `Time::rGetUseTime()` measured process CPU time summed over threads, so any
   parallel speedup *reported itself as a slowdown* (arch0: an 84.6 → 49.8 s run printed 90.9).
   Now `std::chrono::steady_clock` elapsed time. Anyone who ever benchmarked threaded sdpa-qd
   against its own output should re-check their conclusions.
2. **Schur-complement threading** — the first OpenMP in this codebase (upstream has none).
3. **Threaded `Rgemm`, NN and NT cases** (`mpack/Rgemm_NN_omp.cpp`,
   `mpack/Rgemm_NT_omp.cpp`). Patch 2 threads one loop; every BLAS kernel underneath it
   stayed serial, and `Rgemm` is where the time is — measured on an EPYC 7232P at one
   thread, **79.8%** of `gpp100`'s wall and **36.9%** of `arch0`'s. On `gpp100` at
   `OMP_NUM_THREADS=8` the process did not create a second thread at all. NN is 100% of
   `gpp100`'s `Rgemm` time and 97.6% of `arch0`'s; NT is the blocked Cholesky's trailing
   update — negligible on those problems but the whole runtime for
   many-constraints/small-blocks shapes, where the `m x m` Schur Cholesky ran one-core flat
   before the NT kernel (measurements and raw-limb identity evidence in
   [BENCHMARKS.md](BENCHMARKS.md)). Both kernels are bit-identical to the serial code by
   construction and by a raw-limb kernel test, and independent of thread count — each
   column of `C` is written by exactly one thread and its accumulation order is unchanged.
   TN and TT still run the serial bodies; neither shows in any per-problem gemm census.
4. **`--enable-openmp` actually enables OpenMP.** It used to be inert: it added
   `-DENABLE_OPENMP -DNUM_OF_THREADS=` (neither of which anything reads) and no `-fopenmp`,
   so `./configure --enable-openmp` produced a binary byte-identical to a serial one that
   still accepted `OMP_NUM_THREADS` and still reported exactly 1.00× at every thread count.
   It now probes for a working OpenMP flag and **fails the configure** if it cannot find one.
   CI has a matrix row that builds with `--enable-openmp` and no `-fopenmp` in `CXXFLAGS` and
   asserts `nm` finds GOMP symbols — nothing in the program's output can tell those two
   binaries apart, so the check has to be at the object level.

**Measured** (external wall clock, median of 3 pinned repeats) — **these numbers predate
patch 3 and describe patches 1-2 only; they have not been re-taken with the threaded
`Rgemm`**: **1.50×** on an EPYC 7232P
over the four problems with unchanged trajectories (arch0 264.3 → 171.8 s); **1.53×** on an
i9-13900K and **1.48×** on an M1 Max over all five problems — trajectories are identical
everywhere on those two machines, `gpp100` included. On the EPYC, `gpp100` is excluded
honestly: the patch shortens its path there (63 → 49 iterations, same objective), so its
wall-time gain is path length, not speed — per iteration it is 1.00×. The i9 and M1 fork
binaries are the ones built by this README's own instructions from a fresh clone. Full
tables, methodology and raw per-repeat data: [BENCHMARKS.md](BENCHMARKS.md) and
[`bench/`](bench/).

Every modified source file carries an in-file, dated change notice — GPLv2 §2a for the
SDPA sources, LGPL-3 §4a for `mpack/Rgemm_NT_omp.cpp`, which is LGPL-3-only like the MPACK
file it was split from. The complete, always-current list is
`git diff --stat <upstream-base>..HEAD` — an enumeration here went stale twice and is
deliberately not repeated.

Build (QD library required, e.g. `libqd-dev`). Either pass `--enable-openmp`, which since
patch 4 really does add a working OpenMP flag and fails the configure if it cannot find
one, or spell the flags out as below — both are covered by CI matrix rows:

```bash
./configure --with-qd-includedir=/usr/include --with-qd-libdir=/usr/lib/x86_64-linux-gnu \
            CXXFLAGS='-O2 -funroll-all-loops -fopenmp' CFLAGS='-O2 -funroll-all-loops -fopenmp' \
            LDFLAGS='-fopenmp'
make -j$(nproc)
nm sdpa_qd | grep -c GOMP     # the binary is static; ldd cannot verify OpenMP
```

Prefer it automated? This repository packages an agent skill —
[`.claude/skills/install-sdpa-omp/`](.claude/skills/install-sdpa-omp/) — that performs the
whole installation (compiler detection, QD-from-source on macOS, the SPOOLES rescue, OpenMP
verification, smoke test) on Linux and macOS, verified on x86-64 and Apple Silicon. Claude
Code discovers it automatically in a clone;
`bash .claude/skills/install-sdpa-omp/scripts/install.sh` also works standalone.

### macOS (Apple Silicon) — verified on an M1 Max

Historically this platform had three traps; two are now fixed in-tree and only one remains.
**Homebrew's `qd` bottle cannot be used**: it is built with Apple clang against libc++, and
linking with g++ fails on `operator<<(std::ostream&, qd_real const&)` — QD must be built
from source with the same GCC. (The other two are gone: the tracked
`config.guess`/`config.sub` are current automake copies that know arm64, and the SPOOLES
compiler/flags fix is applied automatically by `spooles/Makefile` from
`spooles/patches/patch-Make.inc`, with a build-log assertion that it took effect — no
manual `Make.inc` surgery, no re-untar clobbering.)

```bash
brew install gcc autoconf automake
GCC=$(ls $(brew --prefix gcc)/bin/gcc-[0-9]* | head -1)   # resolve the current version
GXX=$(ls $(brew --prefix gcc)/bin/g++-[0-9]* | head -1)   # (brew install may have just upgraded it)
# QD from source, with the SAME compiler as the solver
curl -LO https://www.davidhbailey.com/dhbsoftware/qd-2.3.24.tar.gz
tar xzf qd-2.3.24.tar.gz
( cd qd-2.3.24 && ./configure CC="$GCC" CXX="$GXX" \
    --prefix=$HOME/qd-gcc --enable-fortran=no && make -j8 && make install )

./configure CC="$GCC" CXX="$GXX" \
            --with-qd-includedir=$HOME/qd-gcc/include --with-qd-libdir=$HOME/qd-gcc/lib \
            CXXFLAGS='-O2 -funroll-all-loops -fopenmp' \
            CFLAGS='-O2 -funroll-all-loops -fopenmp' LDFLAGS='-fopenmp'
make -j8
nm sdpa_qd | grep -c GOMP       # non-zero when OpenMP is really in
```

License: GPL v2, unchanged (`COPYING`); original SDPA authors retain copyright.

## Exit status

Scripts that loop over problems can rely on the exit code:

| outcome | exit |
|---|---|
| solver ran to any stopping condition (`pdOPT`, `pdFEAS`, `pFEAS`, `dFEAS`, `pdINF`, `pINF_dFEAS`, `pFEAS_dINF`, `pUNBD`, `dUNBD`, `noINFO`) | **0** |
| iteration limit reached | **0** |
| infeasibility / unboundedness detected | **0** |
| malformed input, unreadable file, invalid parameter | **1** with a diagnostic (line-numbered for data files) |
| numerical failure with nothing valid to print -- no iteration completed, or the updated `X`/`Z` left the cone **and the rollback could not be refactored** | **2**, `solveStatus = FAILURE` in the result file, no solution section |
| recoverable late failure after `k` good iterations -- a Schur factorisation failure, or an `X`/`Z` update that left the cone and was **rolled back**; the last valid iterate is printed and labelled | **3**, `solveStatus = PARTIAL`, `failureIteration = k` in the result file |

Infeasibility and the iteration limit are valid mathematical results, not errors. Upstream
exited 0 on *every* path -- including fatal errors -- so a crashed run was indistinguishable
from a solved one in any harness that checks exit codes.
