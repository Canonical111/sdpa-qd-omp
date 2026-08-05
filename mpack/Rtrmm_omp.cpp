/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rtrmm.cpp,v 1.4 2009/09/24 07:25:57 nakatamaho Exp $ 
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
Based on http://www.netlib.org/blas/dtrmm.f
Rtrmm performs one of the matrix-matrix operations
 B := alpha*op(A)*B, or B := alpha*B*op(A),
where alpha is a scalar, B is an m by n matrix, A is a unit, or
non-unit, upper or lower triangular matrix and op(A) is one  of
 op(A) = A  or op(A) = A'.
*/

/* NEW FILE (2026-08-05), ported from the dd fork's mplapack/Rtrmm_omp.cpp (bcb7801, "B3").
 * Column-parallel Rtrmm for the SINGLE case this fork's solver uses --
 * Left / Lower / Transpose, i.e. forming Z^-1 = L**T * L in Jal::getInvCholAndInv
 * at sdpa_jordan.cpp:190 -- and every other case, plus every sub-threshold call,
 * is handed to the serial Rtrmm in mpack/Rtrmm.cpp.
 *
 * WHY A SEPARATE ENTRY POINT.  mpack/Rlarfb.cpp calls Rtrmm 24 times and EVERY one
 * of them is side == "Right", so, as with Rtrsm_omp, an in-place pragma on the
 * Left/Lower/Transpose branch would not actually leak into them today.  The
 * separate file is kept for two other reasons.  First, the permitted-call-site
 * contract then lives in one auditable block (mpack/mblas_qd.h).  Second, and
 * specific to this file: mpack/Rtrmm.cpp declares ONE function-scope `qd_real temp`
 * at line 85 and shares it across all EIGHT branches.  An `omp parallel for` in
 * that function would need a private(temp) clause guarding a variable declared
 * three branches away from the pragma -- the kind of coupling that survives one
 * edit and not two.  Here `temp` is declared INSIDE the parallel body, so it is
 * private by construction and cannot be un-privatised by a later edit.
 *
 * THE AXIS: COLUMNS.  For side == "Left" the outermost loop of every branch of
 * Rtrmm runs over the n COLUMNS of B, and every read and write inside it is
 * B[.. + j*ldb] for that one j; A is read-only and no operation mixes two values
 * of j.  Within a column, i ascends and the inner k loop reads rows k > i, which
 * have not been written yet -- so the split changes neither the operands nor their
 * order, and the result is bit-identical to the serial kernel at any thread count.
 * This is NOT true for side == "Right", where the B*A loops read B[.. + k*ldb] from
 * OTHER columns (mpack/Rtrmm.cpp:203-285) and it is the m ROWS that are independent
 * instead.  This file implements Left only and refuses everything else.
 *
 * THE ZERO-SKIP MUST STAY IN STEP WITH mpack/Rtrmm.cpp.  The `B[k + j*ldb] != Zero`
 * test in the inner loop below is NOT in netlib dtrmm; this fork added it to
 * Rtrmm.cpp on 2026-08-04 (see the MODIFIED note at Rtrmm.cpp:76) because the B the
 * solver passes here is lower triangular.  It is reproduced verbatim because the
 * ship gate for this port is bit-identity against THIS fork's serial kernel, not
 * against netlib.  If Rtrmm.cpp's branch ever changes, this body changes with it.
 *
 * FOUR DIFFERENCES FROM THE dd SOURCE, as in Rtrsm_omp.cpp: 0-based indexing
 * (this body transcribes mpack/Rtrmm.cpp:187-199, not dd's 1-based file);
 * mpackint is int32_t so the (double) casts in the work expression are
 * load-bearing; Mlsame_qd returns int; and `max` is the bare using-declaration
 * from mutils_qd.h reached only through <mblas_qd.h>.
 *
 * ACCUMULATION FORM.  `temp = temp + A*B`, matching mpack/Rtrmm.cpp:194 exactly,
 * not dd's `temp += A*B`.  For qd_real the two are the same call today only because
 * qd_inline.h:432 defines operator+= as literally `*this = *this + a`; the explicit
 * form is correct by construction rather than by luck of the header.
 */
#include <mblas_qd.h>
#include "mpack_omp_tuning.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void
Rtrmm_omp(const char *side, const char *uplo, const char *transa,
    const char *diag, mpackint m, mpackint n, qd_real alpha, qd_real * A,
    mpackint lda, qd_real * B, mpackint ldb)
{
    const qd_real Zero = 0.0;

    //The only case this file implements: B := alpha*A'*B with A lower triangular,
    //applied from the left.  Each test is POSITIVE (Mlsame against the wanted
    //letter), never "not the other one", so a malformed argument string falls
    //through to the serial kernel instead of being silently treated as this case.
    //"C" is accepted alongside "T" because Rtrmm itself treats them identically for
    //a real matrix (mpack/Rtrmm.cpp:136 tests only for "N").
    int handled = Mlsame_qd(side, "L") && Mlsame_qd(uplo, "L")
	&& (Mlsame_qd(transa, "T") || Mlsame_qd(transa, "C"))
	&& (Mlsame_qd(diag, "N") || Mlsame_qd(diag, "U"));

    //Argument errors, the m==0/n==0 quick return and the alpha==Zero path are left
    //to Rtrmm so that Mxerbla_qd reporting lives in exactly one place.  nrowa is m
    //because side is "Left".
    int valid = (m > 0) && (n > 0) && (lda >= max((mpackint) 1, m))
	&& (ldb >= max((mpackint) 1, m)) && (alpha != Zero);

    int parallel = handled && valid;
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    parallel = parallel && (nthreads > 1) && !omp_in_parallel()
	&& ((double) m * (double) m * (double) n >=
	MPACK_OMP_TRI_WORK(MPACK_OMP_MIN_TRMM_WORK, nthreads))
	&& (n >= MPACK_OMP_MIN_TRI_WIDTH);
#else
    //Without OpenMP there is nothing to gate; everything goes to the serial kernel.
    parallel = 0;
#endif
    if (!parallel) {
	Rtrmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb);
	return;
    }
#ifdef _OPENMP
    mpackint nounit = Mlsame_qd(diag, "N");
    //Form B := alpha*A'*B.  (lower)
    //
    //schedule(dynamic,1) for the same reason as Rtrsm_omp: the B this call is given
    //is inv(L), lower triangular (sdpa_jordan.cpp:189 copies it from invCholMat), so
    //the zero-skip makes column j cost O((m-j)^2) and a contiguous block split is
    //badly unbalanced.  num_threads() is clamped to n so a narrow B never creates
    //idle team members.  Neither the schedule nor the team size can change the
    //result: every column is assigned whole to exactly one thread.
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads < (int)n ? nthreads : (int)n)
    for (mpackint j = 0; j < n; j++) {
	qd_real temp;
	for (mpackint i = 0; i < m; i++) {
	    temp = B[i + j * ldb];
	    if (nounit)
		temp = temp * A[i + i * lda];
	    for (mpackint k = i + 1; k < m; k++) {
		if (B[k + j * ldb] != Zero) {
		    temp = temp + A[k + i * lda] * B[k + j * ldb];
		}
	    }
	    B[i + j * ldb] = alpha * temp;
	}
    }
#endif
    return;
}
