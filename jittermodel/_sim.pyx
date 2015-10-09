# -*- coding: utf-8 -*-ystem
from __future__ import division
cimport cython
import cmath
import numpy as np

# define the complex c functions, which are roughly a factor of 
# 4 faster than using python's cmath functions 
cdef extern from "complex.h":
    double complex I
    double complex csqrt(double complex x)
    double complex csinh(double complex x)
    double complex ccosh(double complex x)
    double complex ctanh(double complex x)
    double complex cexp(double complex x)
    double creal(double complex x)
    double cimag(double complex x)

def _thetaI_np(double k, double h_s,
        double complex alpha, double complex Lambda, double complex eta,
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

cpdef double complex _thetaI_math(double k, double h_s,
        double complex alpha, double complex Lambda, double complex eta,
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


cpdef double complex _thetaI_c(double k, double h_s,
        double complex alpha, double complex Lambda, double complex eta,
        double complex E_s, double complex E_eff):

    cdef double complex sk, sn, ck, cn
    # TODO document this approximation; link to the ipython notebook
    if k*h_s > 8:
        return (E_s / E_eff * (-Lambda +
            (1 + alpha - Lambda*(1 - 8*cexp(-h_s*(eta + k))
                                 +4*Lambda*cexp(-2*eta*h_s)))/
            (alpha - Lambda + 1)))
    else:
        sk = csinh(k * h_s)
        sn = csinh(eta * h_s)
        ck = ccosh(k * h_s)
        cn = ccosh(eta * h_s)
        return (E_s / E_eff * (-Lambda/ ctanh(eta*h_s) +
                (sk*sn + alpha*ck*sn - Lambda * (ck*cn - 2 + Lambda*sk/sn)) /
                (ck * sn + alpha*sk*sn - Lambda*sk*cn)))


cpdef double complex _thetaII_math(double k, double h, double complex E_s,
                    double complex E_d, double complex E_eff,
                    double complex Lambda):
    return (E_s / E_d * ((E_eff + (1 - Lambda) * E_d / cmath.tanh(k*h)) /
                         (E_eff / cmath.tanh(k*h) + (1 - Lambda) * E_d)))

cpdef double complex _thetaII_c(double k, double h, double complex E_s,
                    double complex E_d, double complex E_eff,
                    double complex Lambda):
    return (E_s / E_d * ((E_eff + (1 - Lambda) * E_d / ctanh(k*h)) /
                         (E_eff / ctanh(k*h) + (1 - Lambda) * E_d)))

@cython.cdivision(True)
cdef double complex _eta_c(double k, double complex kappa, double complex E_s,
                            double D, double omega):
    return csqrt(k**2 + kappa**2 / E_s + omega/D*1j)


@cython.cdivision(True)
cdef double complex _lambda_c(double k, double complex eta,
                             double complex E_eff, double complex E_s):
    """Helper function for calculating the correlation integrand.
    See Lekkala, et al., 2013, Eq. 19"""
    return k/eta*(1 - E_eff/E_s)

def _im_dielectric(k, h_diel, h_trans, E_s, E_i, mu, omega, rho, T, k_B, q, E_0,
                   model):
    sigma = mu * rho * q 
    E_eff = E_s - sigma / (E_0 * omega) * 1j
    kappa = (2 * rho * q ** 2 / (E_0 * k_B * T)) ** 0.5
    diff = mu * k_B * T / q
    if model == 1:
        E_d = E_i
        h = h_diel + h_trans
        alpha = E_eff / E_d
        eta = _eta_c(k, kappa, E_s, diff, omega)
        Lambda = _lambda_c(k, eta, E_eff, E_s)
        theta = _thetaI_c(k, h, alpha, Lambda, eta, E_s, E_eff)
    elif model == 2:
        E_d = E_s
        h = h_diel
        eta = _eta_c(k, kappa, E_s, diff, omega)
        Lambda = _lambda_c(k, eta, E_eff, E_s)
        theta = _thetaII_c(k, h, E_s, E_d, E_eff, Lambda)
    elif model == 20:
        E_d = E_i
        h = h_diel
        eta = _eta_c(k, kappa, E_s, diff, omega)
        Lambda = _lambda_c(k, eta, E_eff, E_s)
        theta = _thetaII_c(k, h, E_s, E_d, E_eff, Lambda)

    result = (E_s - theta) / (E_s + theta)
    return result.imag

@cython.cdivision(True)
cpdef double _im_dielectric_c(double k, double h_diel, double h_trans,
                            double complex E_s, double complex E_i,
                            double mu, double omega, double rho, double T,
                            double k_B, double q, double E_0,
                            int model):
    cdef double sigma = mu * rho * q
    cdef double complex E_eff = E_s - sigma / (E_0 * omega) * I
    cdef double kappa = (2 * rho * q ** 2 / (E_0 * k_B * T)) ** 0.5
    cdef double diff = mu * k_B * T / q

    cdef double h
    cdef double complex E_d, alpha, eta, Lambda, theta
    if model == 1:
        E_d = E_i
        h = h_diel + h_trans
        alpha = E_eff / E_d
        eta = _eta_c(k, kappa, E_s, diff, omega)
        Lambda = _lambda_c(k, eta, E_eff, E_s)
        theta = _thetaI_c(k, h, alpha, Lambda, eta, E_s, E_eff)
    elif model == 2:
        E_d = E_s
        h = h_diel
        eta = _eta_c(k, kappa, E_s, diff, omega)
        Lambda = _lambda_c(k, eta, E_eff, E_s)
        theta = _thetaII_c(k, h, E_s, E_d, E_eff, Lambda)

    cdef double complex result = (E_s - theta) / (E_s + theta)
    return cimag(result)
