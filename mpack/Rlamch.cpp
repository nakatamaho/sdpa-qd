/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2008 by Nakata, Maho
 * 
 * $Id: Rlamch.cpp,v 1.7 2009/09/18 23:01:08 nakatamaho Exp $ 
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

#include <mblas_qd.h>
#include <mlapack_qd.h>

//"E" denots we always calculate relative machine precision (e).
//where 1+e > 1, minimum of e.
qd_real RlamchE_qd(void)
{
  //about 1e-(16+16+16+16)=1e-64
  return qd_real::_eps;
}

//"S" denots we always calculate `safe minimum, such that 1/sfmin does not overflow'.
//cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchS_qd(void)
{
  //about 1e-(308-16-16-16)=1e-260
    return qd_real::_min_normalized;
}
//"B" base  = base of the machine
//cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchB_qd(void)
{
    qd_real two;
    two = 2.0;
    return two;
}

//"P" prec = eps*base
//cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchP_qd(void)
{
    qd_real base, eps, prec;

    base = RlamchB_qd();
    eps = RlamchE_qd();
    prec = eps * base;
    return prec;
}

//"N" t = number of digits in mantissa
//cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchN_qd(void)
{
  return (qd_real)209.0; //52*4
}

//"R" rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
//cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchR_qd(void)
{
    qd_real mtmp;
    mtmp = 1.0;
    return mtmp;
}

//"M"
//cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchM_qd(void)
{
  return qd_real(-1022.0 + 3.0*53.0);
}

//"U"
//cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchU_qd(void)
{
  return qd_real::_min_normalized;
}

//"L"
//cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchL_qd(void)
{
  return (qd_real)1024.0;
}

//"O"
//cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchO_qd(void)
{
    return qd_real::_max; //approx 1.7976931348623157E+308 in float.h
}

//"Z" :dummy
//cf.http://www.netlib.org/blas/dlamch.f
qd_real RlamchZ_qd(void)
{
    qd_real mtemp = 0.0;
    return mtemp;
}

qd_real Rlamch_qd(const char *cmach)
{
    if (Mlsame_qd(cmach, "E"))
	return RlamchE_qd();
    if (Mlsame_qd(cmach, "S"))
	return RlamchS_qd();
    if (Mlsame_qd(cmach, "B"))
	return RlamchB_qd();
    if (Mlsame_qd(cmach, "P"))
	return RlamchP_qd();
    if (Mlsame_qd(cmach, "N"))
	return RlamchN_qd();
    if (Mlsame_qd(cmach, "R"))
	return RlamchR_qd();
    if (Mlsame_qd(cmach, "M"))
	return RlamchM_qd();
    if (Mlsame_qd(cmach, "U"))
	return RlamchU_qd();
    if (Mlsame_qd(cmach, "L"))
	return RlamchL_qd();
    if (Mlsame_qd(cmach, "O"))
	return RlamchO_qd();

    Mxerbla_qd("Rlamch", 1);
    return RlamchZ_qd();
}
