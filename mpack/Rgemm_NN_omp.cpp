/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 *
 * Copyright 2008 by Nakata, Maho
 *
 * $Id: Rgemm.cpp,v 1.5 2009/09/25 04:00:39 nakatamaho Exp $
 *
 * MPACK - multiple precision arithmetic library
 *
 * This file is part of MPACK.
 *
 * MPACK is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License version 3
 * only, as published by the Free Software Foundation.
 *
 * MPACK is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License version 3 for more details
 * (a copy is included in the LICENSE file that accompanied this code).
 *
 * You should have received a copy of the GNU Lesser General Public License
 * version 3 along with MPACK.  If not, see
 * <http://www.gnu.org/licenses/lgpl.html>
 * for a copy of the LGPLv3 License.
 *
 ************************************************************************/
/*
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.

- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* NEW FILE (2026-08-05), split out of mpack/Rgemm.cpp and threaded. See git log.
 *
 * WHAT THIS IS.  The body below is the "Form C := alpha*A*B + beta*C" branch of
 * mpack/Rgemm.cpp -- netlib dgemm's NN case, including its zero-skip on B --
 * lifted into its own translation unit and parallelised over j, structured the
 * way the dd fork's mplapack/Rgemm_NN_omp.cpp is. Before this file existed, the
 * whole of this fork's mpack/ was serial: of the 60 translation units in an
 * -fopenmp build, 59 contained zero GOMP_ / omp_ references, and on gpp100 at
 * OMP_NUM_THREADS=8 the process never created a second thread at all.
 *
 * TWO TRANSFORMATIONS, BOTH BIT-NEUTRAL, AND WHY:
 *
 *  1. Loop fission. mpack/Rgemm.cpp scales column j by beta and then immediately
 *     accumulates into that same column. Here the beta scaling is a complete
 *     serial pass over all j, and the accumulation is a second pass. This is
 *     what makes the accumulation loop a bare canonical OpenMP for -- an
 *     `if (beta == Zero) ... else if (beta != One) ...` prologue inside the
 *     parallel body would be re-tested per thread for no reason, and the two
 *     phases have different cost profiles.
 *
 *  2. Parallelising the accumulation over j.
 *
 * For a fixed element C[i + j*ldc] the sequence of qd_real operations is
 * unchanged by either transformation: exactly one `= Zero` or `= beta * C`,
 * then one `= C + temp * A` per l in ascending order, with the same
 * `B[l + j*ldb] != Zero` skip. Fission only reorders operations on *different*
 * elements, which never interact. The j-split gives every column to exactly one
 * thread, so no element is written by two threads and no accumulation order
 * changes. The result is therefore bit-identical to the serial code AND
 * independent of the thread count -- which is what makes this shippable at all,
 * since the ship gate for this fork is bit-identity, not speed.
 *
 * The one assumption both transformations rest on is the standard BLAS one: C
 * does not alias A or B. Under aliasing the fission alone would already change
 * the answer (it zeroes/scales every column of C before any column is read,
 * where the serial code interleaves the two per column), and the j-split would
 * additionally race. Every call site in this fork passes a distinct destination
 * matrix -- sdpa_linear.cpp:896/:989/:1081 all take retMat, aMat, bMat as three
 * separate DenseMatrix objects -- so the assumption holds here; it is recorded
 * because it is the only way the bit-identity argument above can fail.
 *
 * ACCUMULATION FORM.  Written `C = C + temp * A`, not dd's `C += temp * A`.
 * For qd_real the two happen to be the same thing today, but only because
 * qd_inline.h:432 defines operator+= as literally `*this = *this + a`. That is
 * luck of the header, not design: a future qd release that gave += a fused or
 * reordered implementation would silently change the arithmetic trajectory of
 * every gemm in the solver. The explicit form is correct by construction.
 *
 * PRIVACY OF `temp`.  mpack/Rgemm.cpp declares `qd_real temp` at function scope
 * and shares it across all four transpose cases. That declaration must NOT be
 * inherited here without the `private(...)` clause below -- a function-scope
 * `temp` under an `omp parallel for` without `private` is a silent data race
 * that would corrupt results non-deterministically. The clause is the guard;
 * do not delete it, and do not move `temp` out of the private list while it is
 * declared at function scope.
 *
 * NESTING.  `!omp_in_parallel()` in the if-clause is carried over from dd. It
 * currently rejects nothing measurable in this fork -- instrumentation on arch0,
 * truss5 and gpp100 recorded zero Rgemm calls made from inside a live team, and
 * this fork's one other parallel region (the k1 bMat loop in sdpa_newton.cpp)
 * does not route its sparse work through Rgemm. But the call path
 * run_k1 -> Lal::let -> Lal::multiply -> sdpa_linear.cpp:896 IS an NN Rgemm, so
 * the guard is live insurance for any input where that gate does fire, and
 * without it a nested call would still pay the num_threads evaluation and the
 * fork/join per constraint even though GCC's default OMP_NESTED=false would run
 * it on one thread anyway. Keep it.
 */

#include <mblas_qd.h>
#include "mpack_omp_tuning.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void
Rgemm_NN_omp(mpackint m, mpackint n, mpackint k, qd_real alpha, qd_real * A,
    mpackint lda, qd_real * B, mpackint ldb, qd_real beta, qd_real * C,
    mpackint ldc)
{
    mpackint i, j, l;
    qd_real temp;
    qd_real Zero = 0.0, One = 1.0;

//Form C := alpha*A*B + beta*C.
//beta scaling, serial: one complete pass over the columns (see header note 1).
    for (j = 0; j < n; j++) {
	if (beta == Zero) {
	    for (i = 0; i < m; i++) {
		C[i + j * ldc] = Zero;
	    }
	} else if (beta != One) {
	    for (i = 0; i < m; i++) {
		C[i + j * ldc] = beta * C[i + j * ldc];
	    }
	}
    }
//main loop.  (double) casts are load-bearing: mpackint is int32_t in this fork
//(mpack_config.h leaves USE64BITINT commented out), so a bare m*n*k overflows at
//m=n=k=1290.  Do not remove them because dd, whose mplapackint is int64_t, can.
#ifdef _OPENMP
#pragma omp parallel for private(i, j, l, temp) if ((double)m * (double)n * (double)k >= MPACK_OMP_MIN_GEMM_WORK && n >= MPACK_OMP_MIN_GEMM_WIDTH && !omp_in_parallel()) num_threads(omp_get_max_threads() < (int)n ? omp_get_max_threads() : (int)n)
#endif
    for (j = 0; j < n; j++) {
	for (l = 0; l < k; l++) {
	    if (B[l + j * ldb] != Zero) {
		temp = alpha * B[l + j * ldb];
		for (i = 0; i < m; i++) {
		    C[i + j * ldc] =
			C[i + j * ldc] + temp * A[i + l * lda];
		}
	    }
	}
    }
    return;
}
