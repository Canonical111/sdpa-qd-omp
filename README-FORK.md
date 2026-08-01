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

### macOS (Apple Silicon) — verified on an M1 Max

Three traps, all with one-line causes. (1) **Homebrew's `qd` bottle cannot be used**: it is
built with Apple clang against libc++, and linking with g++ fails on
`operator<<(std::ostream&, qd_real const&)` — QD must be built from source with the same
GCC. (2) The 2009 `config.guess` predates arm64. (3) A failed `make` re-untars SPOOLES and
**clobbers any fix you made** — edit only after the first failure, and satisfy the
`libspooles.a` target so it is not re-extracted.

```bash
brew install gcc autoconf automake
GCC=$(ls $(brew --prefix gcc)/bin/gcc-[0-9]* | head -1)   # resolve the current version
GXX=$(ls $(brew --prefix gcc)/bin/g++-[0-9]* | head -1)   # (brew install may have just upgraded it)
# QD from source, with the SAME compiler as the solver
curl -LO https://www.davidhbailey.com/dhbsoftware/qd-2.3.24.tar.gz
tar xzf qd-2.3.24.tar.gz
( cd qd-2.3.24 && ./configure CC="$GCC" CXX="$GXX" \
    --prefix=$HOME/qd-gcc --enable-fortran=no && make -j8 && make install )

# arm64-aware config.guess/config.sub
cp "$(ls -d /opt/homebrew/share/automake-* | head -1)"/config.{guess,sub} .
chmod u+w config.guess config.sub

./configure CC="$GCC" CXX="$GXX" \
            --with-qd-includedir=$HOME/qd-gcc/include --with-qd-libdir=$HOME/qd-gcc/lib \
            CXXFLAGS='-O2 -funroll-all-loops -fopenmp' \
            CFLAGS='-O2 -funroll-all-loops -fopenmp' LDFLAGS='-fopenmp'
make -j8 || true   # first pass stops inside SPOOLES -- expected; the untar creates spooles/build
M=spooles/build/Make.inc
sed -i '' 's|^# CC = gcc|  CC = '$GCC'|' $M
sed -i '' 's|^  CFLAGS += -O2 -funroll-all-loops|  CFLAGS += -O2 -funroll-all-loops -Wno-error=int-conversion -Wno-error=implicit-function-declaration -Wno-error=incompatible-pointer-types|' $M
( cd spooles/build && find . -name '*.o' -delete && rm -f spooles.a && make global -f makefile )
cp spooles/build/spooles.a spooles/build/libspooles.a   # satisfies the target; prevents re-untar
make -j8
nm sdpa_qd | grep -c GOMP       # non-zero when OpenMP is really in
```

License: GPL v2, unchanged (`COPYING`); original SDPA authors retain copyright.
