# -*- coding: utf-8 -*-
from __future__ import division
import cmath
import numpy as np
cimport numpy as np

DTYPE = np.int
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.int_t DTYPE_t

def _thetaI_np(double k, double h_s,
        double complex alpha, double complex Lambda, complex eta,
        double complex E_s, double complex E_eff):

    cdef double complex sk, sn, ck, cn
    if k*h_s > 8:
        return (E_s / E_eff * (-Lambda +
            (1 + alpha - Lambda*(1 - 8*np.exp(-h_s*(eta + k))
                                 +4*Lambda*np.exp(-2*eta*h_s)))/
            (alpha - Lambda + 1)))
    else:
        sk = np.sinh(k * h_s)
        sn = np.sinh(eta * h_s)
        ck = np.cosh(k * h_s)
        cn = np.cosh(eta * h_s)
        return (E_s / E_eff * (-Lambda/ np.tanh(eta*h_s) +
                (sk*sn + alpha*ck*sn - Lambda * (ck*cn - 2 + Lambda*sk/sn)) /
                (ck * sn + alpha*sk*sn - Lambda*sk*cn)))

def _thetaI_math(double k, double h_s,
        double complex alpha, double complex Lambda, complex eta,
        double complex E_s, double complex E_eff):

    cdef double complex sk, sn, ck, cn
    if k*h_s > 8:
        return (E_s / E_eff * (-Lambda +
            (1 + alpha - Lambda*(1 - 8*cmath.exp(-h_s*(eta + k))
                                 +4*Lambda*cmath.exp(-2*eta*h_s)))/
            (alpha - Lambda + 1)))
    else:
        sk = cmath.sinh(k * h_s)
        sn = cmath.sinh(eta * h_s)
        ck = cmath.cosh(k * h_s)
        cn = cmath.cosh(eta * h_s)
        return (E_s / E_eff * (-Lambda/ cmath.tanh(eta*h_s) +
                (sk*sn + alpha*ck*sn - Lambda * (ck*cn - 2 + Lambda*sk/sn)) /
                (ck * sn + alpha*sk*sn - Lambda*sk*cn)))
