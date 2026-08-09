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

/* MODIFICATION NOTICE (LGPL-3 section 4a; this file is LGPL-3-only, NOT the
   GPLv2 the SDPA sources around it carry -- the earlier notice here said
   "GPLv2 2a" by copy-paste from the solver files and named the wrong licence),
   2026-08-08: new file. Threaded NT case of Rgemm, split out of
   mpack/Rgemm.cpp the same way Rgemm_NN_omp was. See git log.

   The notices ABOVE this line are a verbatim copy of mpack/Rgemm.cpp's, and
   are byte-identical to it including trailing whitespace. That whitespace is
   deliberate fidelity to the copied copyright and licence text, not an
   oversight: do not "clean" it, because editing a copyright notice to tidy
   whitespace is exactly the kind of silent alteration the notice exists to
   prevent. (Rgemm_NN_omp.cpp's copy was tidied before this was understood;
   that is the inconsistent one.) */

/* WHY THIS FILE EXISTS.  Rpotrf's blocked trailing update on the "Lower" path
 * -- the only path SDPA ever calls -- is Rgemm("No transpose","Transpose"),
 * and until this file that case ran the serial body.  For solver runs in the
 * published qd benchmark set that did not matter: NN is 100% of gpp100's Rgemm
 * time and 97.6% of arch0's (the per-problem gemm census quoted in BENCHMARKS.md (recipe repo: results/qd_* tables)).
 * It matters enormously for the inverse shape -- many constraints, small
 * blocks -- where the m x m Schur complement Cholesky IS the runtime and its
 * trailing update is this kernel.  First observed on a bootstrap problem with
 * m = 2439 and 17 blocks all of order <= 30: two iterations cost 94.8 s wall
 * and 94.3 s CPU at any thread count, i.e. one core, because this case was
 * serial.
 *
 * AXIS.  C := alpha*A*B' + beta*C.  Column j of C is written from B row j and
 * A only: every write inside the j-iteration is C[i + j*ldc], A is read-only,
 * and no other column of C is read.  Splitting over j therefore leaves each
 * output element's accumulation order (l outer, i inner, exactly the serial
 * body's order) untouched, which is what makes the result bit-identical to
 * the serial kernel at any thread count.  In Rpotrf's trailing update the
 * column count is the panel width jb (= 64 for every blocked call in this
 * fork, via the fixed iMlaenv), so the split hands 64 columns to the team;
 * num_threads is clamped to n so a narrower call never creates idle members.
 *
 * BETA SCALING is fissioned into a serial pass first, exactly as in
 * Rgemm_NN_omp: per element the scale still happens before every
 * accumulation, in the same order, so the fission is bit-neutral; keeping it
 * out of the parallel loop keeps the parallel body pure accumulation.
 *
 * ZERO-SKIP.  The `if (B[j + l*ldb] != Zero)` guard is the serial body's own
 * (netlib dgemm's) skip, kept verbatim.  In the trailing update B is the
 * current panel, which after the panel factorisation has no structural zeros
 * to speak of, so the guard rarely fires there -- it is kept because the
 * kernel must remain interchangeable with the serial body on EVERY input, not
 * because it wins time.
 *
 * ACCUMULATION FORM.  `C = C + temp * A`, not `+=`, for the reason recorded
 * in Rgemm_NN_omp.cpp: qd_inline.h happens to define += as `*this = *this+a`
 * today, and relying on that would let a future qd release silently change
 * every gemm trajectory.
 *
 * GATE.  Shares MPACK_OMP_MIN_GEMM_WORK / MPACK_OMP_MIN_GEMM_WIDTH with the
 * NN kernel: the per-madd cost is the same qd_real multiply-add (123-160 ns
 * measured on thanos), the fork/join cost is the same, and the split axis is
 * the same, so there is no basis for a different crossover.  The (double)
 * casts are load-bearing -- mpackint is int32_t in this fork and m*n*k
 * overflows at m=n=k=1290.
 *
 * PRIVACY.  i, j, l, temp are declared at function scope to mirror the serial
 * body; the private(...) clause below is what keeps them per-thread.  Do not
 * remove it.
 */

#include <mblas_qd.h>
#include "mpack_omp_tuning.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void
Rgemm_NT_omp(mpackint m, mpackint n, mpackint k, qd_real alpha, qd_real * A,
    mpackint lda, qd_real * B, mpackint ldb, qd_real beta, qd_real * C,
    mpackint ldc)
{
    mpackint i, j, l;
    qd_real temp;
    qd_real Zero = 0.0, One = 1.0;

//Form C := alpha*A*B' + beta*C.
//beta scaling, serial: one complete pass over the columns (see header).
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
//main loop.
#ifdef _OPENMP
#pragma omp parallel for private(i, j, l, temp) if ((double)m * (double)n * (double)k >= MPACK_OMP_MIN_GEMM_WORK && n >= MPACK_OMP_MIN_GEMM_WIDTH && !omp_in_parallel()) num_threads(omp_get_max_threads() < (int)n ? omp_get_max_threads() : (int)n)
#endif
    for (j = 0; j < n; j++) {
	for (l = 0; l < k; l++) {
	    if (B[j + l * ldb] != Zero) {
		temp = alpha * B[j + l * ldb];
		for (i = 0; i < m; i++) {
		    C[i + j * ldc] =
			C[i + j * ldc] + temp * A[i + l * lda];
		}
	    }
	}
    }
    return;
}
