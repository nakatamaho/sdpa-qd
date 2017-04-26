/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: mlapack_qd.h,v 1.6 2009/09/22 20:27:18 nakatamaho Exp $ 
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

#ifndef _MLAPACK_QD_H_
#define _MLAPACK_QD_H_

/* this is a subset of mpack for SDPA-GMP only */
/* http://mplapack.sourceforge.net/ */

/* mlapack prototypes */
void Rsteqr(const char *compz, mpackint n, qd_real * d, qd_real * e,
    qd_real * Z, mpackint ldz, qd_real * work, mpackint *info);
void
    Rsyev(const char *jobz, const char *uplo, mpackint n, qd_real * A,
    mpackint lda, qd_real * w, qd_real * work, mpackint *lwork, mpackint *info);
void Rpotrf(const char *uplo, mpackint n, qd_real * A, mpackint lda, mpackint *info);
mpackint iMlaenv_qd(mpackint ispec, const char *name, const char *opts, mpackint n1, mpackint n2,
    mpackint n3, mpackint n4);
qd_real Rlamch_qd(const char *cmach);
qd_real Rlansy(const char *norm, const char *uplo, mpackint n, qd_real * A,
    mpackint lda, qd_real * work);
void Rlascl(const char *type, mpackint kl, mpackint ku, qd_real cfrom, qd_real cto,
    mpackint m, mpackint n, qd_real * A, mpackint lda, mpackint *info);
void Rsytrd(const char *uplo, mpackint n, qd_real * A, mpackint lda, qd_real * d,
    qd_real * e, qd_real * tau, qd_real * work, mpackint lwork, mpackint *info);
void Rsytd2(const char *uplo, mpackint n, qd_real * A, mpackint lda, qd_real * d,
    qd_real * e, qd_real * tau, mpackint *info);
qd_real Rlanst(const char *norm, mpackint n, qd_real * d, qd_real * e);
void Rlae2(qd_real a, qd_real b, qd_real c, qd_real * rt1,
    qd_real * rt2);
qd_real Rlapy2(qd_real x, qd_real y);
void Rlasrt(const char *id, mpackint n, qd_real * d, mpackint *info);
void Rorgql(mpackint m, mpackint n, mpackint k, qd_real * A, mpackint lda, qd_real * tau,
    qd_real * work, mpackint lwork, mpackint *info);
void Rorgqr(mpackint m, mpackint n, mpackint k, qd_real * A, mpackint lda, qd_real * tau,
    qd_real * work, mpackint lwork, mpackint *info);
void Rlarfg(mpackint N, qd_real * alpha, qd_real * x, mpackint incx,
    qd_real * tau);
void Rlassq(mpackint n, qd_real * x, mpackint incx, qd_real * scale,
    qd_real * sumsq);
void Rorg2l(mpackint m, mpackint n, mpackint k, qd_real * A, mpackint lda, qd_real * tau,
    qd_real * work, mpackint *info);
void Rlarft(const char *direct, const char *storev, mpackint n, mpackint k,
    qd_real * v, mpackint ldv, qd_real * tau, qd_real * t, mpackint ldt);
void Rlarfb(const char *side, const char *trans, const char *direct,
    const char *storev, mpackint m, mpackint n, mpackint k, qd_real * V, mpackint ldv,
    qd_real * T, mpackint ldt, qd_real * C, mpackint ldc, qd_real * work,
    mpackint ldwork);
void Rorg2r(mpackint m, mpackint n, mpackint k, qd_real * A, mpackint lda, qd_real * tau,
    qd_real * work, mpackint *info);
void Rlarf(const char *side, mpackint m, mpackint n, qd_real * v, mpackint incv,
    qd_real tau, qd_real * C, mpackint ldc, qd_real * work);
void Rpotf2(const char *uplo, mpackint n, qd_real * A, mpackint lda, mpackint *info);
void Rlaset(const char *uplo, mpackint m, mpackint n, qd_real alpha, qd_real beta,
    qd_real * A, mpackint lda);
void Rlaev2(qd_real a, qd_real b, qd_real c, qd_real * rt1,
    qd_real * rt2, qd_real * cs1, qd_real * sn1);
void Rlasr(const char *side, const char *pivot, const char *direct, mpackint m,
    mpackint n, qd_real * c, qd_real * s, qd_real * A, mpackint lda);
void Rlartg(qd_real f, qd_real g, qd_real * cs, qd_real * sn,
    qd_real * r);
void Rlatrd(const char *uplo, mpackint n, mpackint nb, qd_real * A, mpackint lda, qd_real * e, qd_real * tau, qd_real * w, mpackint ldw);
void Rsterf(mpackint n, qd_real * d, qd_real * e, mpackint *info);
void Rorgtr(const char *uplo, mpackint n, qd_real * a, mpackint lda, qd_real * tau,
    qd_real * work, mpackint lwork, mpackint *info);
#endif
