/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rtrsm.cpp,v 1.4 2009/09/24 07:25:57 nakatamaho Exp $ 
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

/*
Based on http://www.netlib.org/blas/dtrsm.f
Rtrsm solves one of the matrix equations
 op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
where alpha is a scalar, X and B are m by n matrices, A is a unit, or
non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
   op( A ) = A   or   op( A ) = A'.
The matrix X is overwritten on B.
*/

/* NEW FILE (2026-08-05), ported from the dd fork's mplapack/Rtrsm_omp.cpp (bcb7801, "B3").
 * Column-parallel Rtrsm for the SINGLE case this fork's solver uses --
 * Left / Lower / NoTranspose, i.e. Lal::getInvLowTriangularMatrix at
 * sdpa_linear.cpp:446 -- and every other case, plus every sub-threshold call, is
 * handed to the serial Rtrsm in mpack/Rtrsm.cpp.
 *
 * WHY A SEPARATE ENTRY POINT RATHER THAN A PRAGMA INSIDE Rtrsm.cpp.
 * Two reasons, and the first is weaker in this fork than it is in dd, so it is
 * stated honestly.  (a) In dd, Rtrsm is also reached from Rpotrf2's recursion and
 * Rpotrf's UPPER path at shapes for which nothing is measured.  Here the internal
 * callers are only mpack/Rpotrf.cpp:120 ("Left","Upper","Transpose") and
 * mpack/Rpotrf.cpp:143 ("Right","Lower","Transpose") -- this fork has no
 * Rpotrf2.cpp -- and NEITHER is Left/Lower/NoTranspose, so an in-place pragma
 * confined to that one branch would in fact not leak.  (b) What does still argue
 * for the separate file: it keeps the permitted-call-site contract in one auditable
 * block (see mpack/mblas_qd.h), and mpack/Rpotrf.cpp:143 is exactly the
 * Right/Lower/Transpose shape dd's B1 threaded, so this symbol is the seam a later
 * B1 port plugs into.  Do not add that case here by analogy with the Left one --
 * see the axis note below; it is a different axis.
 *
 * THE AXIS: COLUMNS.  For side == "Left" the outermost loop of every branch of
 * Rtrsm runs over the n COLUMNS of B, and every read and write inside it is
 * B[.. + j*ldb] for that one j; A is read-only and no operation mixes two values of
 * j.  So the columns of B are independent, the arithmetic within a column is
 * untouched by the split, and the result is bit-identical to the serial kernel at
 * any thread count.  This is NOT true for side == "Right", where column j is
 * updated from column k for every k < j (mpack/Rtrsm.cpp:246-303) -- the columns
 * there are a dependent chain and it is the m ROWS that are independent instead.
 * A column split of a Right-side call compiles, converges, and returns wrong
 * answers.  This file implements Left only and refuses everything else, so that
 * distinction cannot be got wrong by accident.
 *
 * FOUR THINGS THAT DIFFER FROM THE dd SOURCE AND WOULD BE WRONG IF COPIED:
 *  1. This fork's kernels are 0-BASED (`for (j = 0; j < n; j++)`, `B[i + j*ldb]`);
 *     dd's are a 1-based transcription with `b[(i-1) + (j-1)*ldb]`.  The body below
 *     is a transcription of mpack/Rtrsm.cpp:159-177, NOT of dd's file.
 *  2. mpackint is int32_t here (mpack_config.h leaves USE64BITINT commented out),
 *     so the (double) casts in the work expression are load-bearing: a bare
 *     m*m*n overflows at m=n=1290, and SDPLIB contains blocks of 1600
 *     (1600^3 = 4.1e9 > INT32_MAX).  dd's int64_t mplapackint hides this.
 *  3. Mlsame_qd returns int, not bool; it is used bare in the conditions below,
 *     the way mpack/Rtrsm.cpp:90-107 uses it.
 *  4. max is the bare `using std::max` from mutils_qd.h:31, reached only via
 *     <mblas_qd.h> (which includes <qd/qd_real.h> first).  Include the umbrella
 *     header, never mutils_qd.h directly.
 *
 * `Zero` and `One` are const here where mpack/Rtrsm.cpp declares them mutable.
 * They are read-only inside the parallel region either way, so this is not a fix
 * for a race -- it is so that a future edit cannot introduce one.  The loop
 * indices i, j, k are declared inside the region, so they are private by
 * construction rather than by a private() clause that could drift from them.
 */
#include <mblas_qd.h>
#include "mpack_omp_tuning.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void
Rtrsm_omp(const char *side, const char *uplo, const char *transa,
    const char *diag, mpackint m, mpackint n, qd_real alpha, qd_real * A,
    mpackint lda, qd_real * B, mpackint ldb)
{
    const qd_real Zero = 0.0, One = 1.0;

    //The only case this file implements: B := alpha*inv(A)*B with A lower
    //triangular, not transposed, applied from the left.  Each test is POSITIVE
    //(Mlsame against the wanted letter), never "not the other one", so a malformed
    //argument string falls through to the serial kernel instead of being silently
    //treated as this case.
    int handled = Mlsame_qd(side, "L") && Mlsame_qd(uplo, "L")
	&& Mlsame_qd(transa, "N")
	&& (Mlsame_qd(diag, "N") || Mlsame_qd(diag, "U"));

    //Argument errors, the m==0/n==0 quick return and the alpha==Zero path are left
    //to Rtrsm so that Mxerbla_qd reporting lives in exactly one place.  nrowa is m
    //because side is "Left"; do not reuse this lda test for a Right-side case.
    int valid = (m > 0) && (n > 0) && (lda >= max((mpackint) 1, m))
	&& (ldb >= max((mpackint) 1, m)) && (alpha != Zero);

    int parallel = handled && valid;
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    parallel = parallel && (nthreads > 1) && !omp_in_parallel()
	&& ((double) m * (double) m * (double) n >=
	MPACK_OMP_TRI_WORK(MPACK_OMP_MIN_TRSM_WORK, nthreads))
	&& (n >= MPACK_OMP_MIN_TRI_WIDTH);
#else
    //Without OpenMP there is nothing to gate; everything goes to the serial kernel.
    parallel = 0;
#endif
    if (!parallel) {
	Rtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb);
	return;
    }
#ifdef _OPENMP
    mpackint nounit = Mlsame_qd(diag, "N");
    //Form B := alpha*inv(A)*B.  (lower, no transpose)
    //
    //schedule(dynamic,1), not the default static split: the B this call is given is
    //the IDENTITY (sdpa_linear.cpp:445 calls retMat.setIdentity() immediately
    //before), so the `B[k + j*ldb] != Zero` test at the top of the k loop skips the
    //leading zeros and column j costs O((m-j)^2).  A contiguous block split would
    //hand the first thread roughly a third of the total work.
    //
    //num_threads() is clamped to n so a narrow B never creates idle team members.
    //Neither the schedule nor the team size can change the result: every column is
    //assigned whole to exactly one thread, so no element is touched twice and no
    //accumulation order changes.
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads < (int)n ? nthreads : (int)n)
    for (mpackint j = 0; j < n; j++) {
	if (alpha != One) {
	    for (mpackint i = 0; i < m; i++) {
		B[i + j * ldb] = alpha * B[i + j * ldb];
	    }
	}
	for (mpackint k = 0; k < m; k++) {
	    if (B[k + j * ldb] != Zero) {
		if (nounit)
		    B[k + j * ldb] = B[k + j * ldb] / A[k + k * lda];
		for (mpackint i = k + 1; i < m; i++) {
		    B[i + j * ldb] =
			B[i + j * ldb] - B[k + j * ldb] * A[i + k * lda];
		}
	    }
	}
    }
#endif
    return;
}
