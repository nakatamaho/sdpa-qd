/*************************************************************************
 *
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 * 
 * Copyright 2009 by Nakata, Maho
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

#ifndef _MUTILS_QD_H_
#define _MUTILS_QD_H_

using std::max;
using std::min;

qd_real Msign(qd_real a, qd_real b);
double cast2double(qd_real a);
int M2int(qd_real a);
//void mpf_pow(mpf_t ans, mpf_t x, mpf_t y);
qd_real mpf_approx_log(qd_real x);
qd_real mpf_approx_log2(qd_real x);
qd_real mpf_approx_log10(qd_real x);
qd_real mpf_approx_pow(qd_real x, qd_real y);
qd_real mpf_approx_cos(qd_real x);
qd_real mpf_approx_sin(qd_real x);
qd_real mpf_approx_exp(qd_real x);
qd_real mpf_approx_pi();

//implementation of sign transfer function.
inline qd_real Msign(qd_real a, qd_real b)
{
  qd_real mtmp;
  mtmp=abs(a);
  if (b<0.0) {
    mtmp=-mtmp;
  }
  return mtmp;
}

inline double
cast2double(qd_real a)
{
    return a.x[0];
}

inline int
M2int(qd_real a)
{
    int i;
    qd_real tmp;
    a = a + 0.5;
    tmp = floor(a);
    i = (int)tmp.x[0];
    return i;
}

#endif
