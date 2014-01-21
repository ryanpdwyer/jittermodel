"""
Unit Simulation
================

The simulation code, now with units.

Simulation code."""

from __future__ import division
import numpy as np
from numpy import (pi, sinh, tanh, arccosh,
                   exp, log10, arctanh)
from scipy.integrate import quad
from scipy.special import jn
from numdifftools import Derivative
import math
from copy import copy
from jittermodel import u

E_0 = 8.854e-3
k_B = 1.38065e-2
q = 1.602e-1

inf = float('+infinity')


def coth(x):
    """The hyperpolic cotanget of x."""
    return (1 / tanh(x))


def int_sum_sinh(alpha, eps):
    """Determines the number of terms needed to be within eps
    of the sum, using the integral test."""
    if alpha < 5:
        n_min = np.ceil(2 / alpha * arctanh(exp(-eps * alpha / sinh(alpha))))
    else:
        n_min = 1
    return n_min


def sum_sinh(alpha, eps=1e-8):
    """This calculates the infinite sum, sinh(alpha) / sinh(alpha * n),
    for n = 1 to infinity.

    We manually require at least 4 terms so that the derivative is
    numerically stable. We use math.fsum to give a numerically stable sum."""
    summand = lambda n: sinh(alpha) / sinh(alpha * n)

    N_max = max(4, int_sum_sinh(alpha, eps))
    terms = summand(np.arange(1, N_max + 1))
    return math.fsum(terms)


class UnitSimulation(object):

    """This calculates experimental parameters such as capacitance
    and sample-induced friction for a given cantilever, sample,
    and experiment.
    """

    def __init__(self, cantilever, sample, experiment):
        """Initialize the simulation with the values from the given
        cantilever, sample and experiment. It also calculates
        parameters used in the simulation"""
        self.Cant = copy(cantilever)
        self.Samp = copy(sample)
        self.Expt = copy(experiment)

        self.func_dict = {'friction': self.calc_gamma_s,
                          'jitter': self.calc_jitter}

    def assign(self, attr, val):
        """Assign the attribute 'attr' to value 'val'. Behind the
        scenes, this function finds which property came from which
        class, and then does the appropriate assignment."""
        for item in (self.Cant, self.Samp, self.Expt):
            if attr in item._all_attributes:
                item.assign(attr, val)

    def lookup(self, attr):
        """Looks for an attribute in the cantilever, sample, or
        experiment namespaces, and returns it."""
        for item in (self.Cant, self.Samp, self.Expt):
            if attr in item._all_attributes:
                return item.lookup(attr)
        return '{attr} not found'.format(attr=attr)

    # Capacitance Calculations  #################
    def C_sphere(self, d=None):
        """Capacitance between a sphere and a thin sample, calculated
        according to a modified version of the equation from Brus et al.,
        http://dx.doi.org/10.1021/jp0265438. The capacitance is
        approximated by truncating the exact infinite sum.
        """
        if d is None:
            d = self.Expt.d
        if d < 0:
            # Prevent numerical derivatives from giving an error.
            d = 0

        cant = self.Cant
        samp = self.Samp

        alpha = arccosh(1 + d / cant.R_tip + samp.h /
                        (samp.E_s1 * cant.R_tip))

        return (4 * pi * E_0 * cant.R_tip * sum_sinh(alpha))

    def Cd_sphere(self, d=None):
        """Returns the numerical derivative of the sphere capacitance
        C_sphere at a distance d. Uses numdifftools
        (https://code.google.com/p/numdifftools/) to calculate the
        derivative to high precision automatically."""
        if d is None:
            d = self.Expt.d
        Cd = Derivative(self.C_sphere)
        return float(Cd(d))

    def Cd2_sphere(self, d=None):
        """Returns the second derivative of the sphere capacitance
        C_sphere at a distance d. Uses numdifftools
        (https://code.google.com/p/numdifftools/) to calculate the
        derivative to high precision automatically."""
        if d is None:
            d = self.Expt.d
        Cd2 = Derivative(self.C_sphere, n=2)
        return float(Cd2(d))

    #

    # def calc_rho(self, V_g = None):
        # """Calculates the charge density rho when a given gate voltage
        # is placed on the sample. This function may be obsolete with the
        # newly refactored code."""
    #     if V_g is None:
    #         V_g = self.Expt.V_g
    #
    #     rho = 1 * self.Samp.C_i * V_g / ( self.Samp.h_trans * q)
    #     return rho

    # Correlation Function related functions. #####

    def _prefactor(self, omega):
        """This is the prefactor for correlation function integrals.
        See Eq. 16 in Lekkala et al. 2013."""
        return (-k_B * self.Samp.T / (4 * pi * E_0 * omega))

    def _im_dielectric(self, k, omega=None, model=None):
        """Computes the complicated expression containing
        :math:`\epsilon_\mathrm{rel}(\omega)` and
        :math:`\theta(k,\omega)` in the correlation function
        from Lekkala et al. 2013. Currently uses model II, but
        support for model I can be easily added. See Eq. 16.

        TODO: Add unittest, even with just hard-coded comparisons
        to Mathematica."""
        if omega is None:
            omega = self.Cant.omega_c
        if model is None:
            model = 'II'
        elif model == 'I':
            raise NotImplementedError('This model has not beed coded yet.')
        else:
            raise ValueError("Model must be either 'I' or 'II'.")

        samp = self.Samp
        if model == 'II':
            kappa = samp.kappa
            E_s = samp.E_s
            E_eff = samp.E_eff(omega)
            E_d = samp.E_s
            h = samp.h_diel

        nu = ((k ** 2 + kappa ** 2 / E_s + omega / self.Samp.diff * 1j)
              ** (0.5 + 0j))
        llambda = (1 - E_eff / E_s) * k / nu
        thetaII = E_s / E_d * ((E_eff + (1 - llambda) * E_d * coth(k * h)) /
                              (E_eff * coth(k * h) + (1 - llambda) * E_d))

        result = (E_s - thetaII) / (E_s + thetaII)
        return result.imag

    def _corr_integrand(self, r1, r2, z1, z2, k, omega=None, n=0):
        """The integrand of the correlation function from
        Lekkala et al., 2013, at the cantilever resonance frequency."""
        if omega is None:
            omega = self.Cant.omega_c
        return ((-k) ** n * jn(0, k * abs(r1 - r2)) *
                exp(-1 * k * (z1 + z2)) * self._im_dielectric(k, omega))

    def _base_corr(self, d, omega, n_derivs):
        """This is the base function all of the other correlation
        functions which are derivatives in z are defined using."""
        if omega is None:
            omega = self.Cant.omega_c

        integral, error = quad(lambda x: self._corr_integrand(
            0, 0, d, d, x, omega, n=n_derivs), 0, inf)

        return self._prefactor(omega) * integral

    def corr_PP(self, d, omega=None):
        """Autocorrelation function of potential fluctuations,
        defined according to Lekkala et al., 2013. See _base_corr for
        more information."""
        return self._base_corr(d, omega, 0)

    def corr_PEz(self, d, omega=None):
        """"""
        return self._base_corr(d, omega, 1)

    def corr_PEzz(self, d, omega=None):
        """"""
        return self._base_corr(d, omega, 2)

    def corr_EzEz(self, d, omega=None):
        """"""
        return self._base_corr(d, omega, 2)

    def corr_EzEzz(self, d, omega=None):
        """"""
        return self._base_corr(d, omega, 3)

    def corr_EzzEzz(self, d, omega=None):
        """"""
        return self._base_corr(d, omega, 4)

    # Calculates the one corr_ExxExx related function we need.
    def corr_ExxExx(self, d, omega=None):
        """Calculate the autocorrelation function of x-direction
        electric field gradient fluctuations.

        This formula uses the result from taking the derivative of
        the potential autocorrelation fluctuations. This result is
        summerized below.
        :math:`\\left.\\frac{\\partial^4}{\\partial x_1^2 \\partial x_2^2}C(r_1,d,r_2,d;\\omega)\\right|_{r_1 = r_2 = 0}`

        :math:`\\left.\\frac{\\partial^4}{\\partial x_1^2 \\partial x_2^2}J_0(k|r_1 - r_2|)|_{r_1 = r_2 = 0} = \\frac{3 k^4}{8}`

        In short, we get a factor of 3k^4/8 from this derivative.
        """

        _integrand = lambda k: (3 * k ** 4 / 8 * exp(-2 * k * d) *
                                self._im_dielectric(k, omega))
        return (self.prefactor(omega) * quad(_integrand, 0, inf))

    def _power_spectrum_parallel(self, omega, d=None):
        """Calculates the power spectrum of frequency fluctuations
        at the angular frequency omega for a parallel geometry
        cantilever. See Eq. 6, Lekkala 2013."""
        if d is None:
            d = self.Expt.d
        _prefactor = (self.V_ts * self.C_sphere(d) / self.k_c) ** 2

        return (_prefactor * self.corr_ExxExx(d, omega))

    def _power_spectrum_perpendicular(self, omega, d=None):
        """Calculates the power spectrum of frequency fluctuations
        at angular frequency omega for a perpendicular geometry
        cantilever. See Eq. 9, Lekkala 2013."""
        if d is None:
            d = self.Expt.d
        power_spectrum_prefactor = (
            self.Cant.f_c * self.Expt.V_ts / self.Cant.k_c) ** 2
        c = self.C_sphere(d)
        cd = self.Cd_sphere(d)
        cdd = self.Cd2_sphere(d)

        terms = [
            c ** 2 * self.corr_EzzEzz(d, omega),
            -4 * c * cd * self.corr_EzEzz(d, omega),
            4 * cd ** 2 * self.corr_EzEz(d, omega),
            -2 * c * cdd * self.corr_PEzz(d, omega),
            -4 * cd * cdd * self.corr_PEz(d, omega),
            cdd ** 2 * self.corr_PP(d, omega)
        ]
        return (power_spectrum_prefactor * sum(terms))

    def _gamma_perpendicular(self, d=None):
        """Calculates the perpendicular geometry friction,
        according to the formula in Lekkala et al 2013. Which formula?"""
        if d is None:
            d = self.Expt.d
        omega_c = self.Cant.omega_c

        gamma_prefactor = self.Expt.V_ts ** 2 / (k_B * self.Samp.T)

        c = self.C_sphere(d)
        cd = self.Cd_sphere(d)

        terms = [
            c ** 2 * self.corr_EzEz(d, omega_c),
            cd ** 2 * self.corr_PP(d, omega_c),
            -2 * c * cd * self.corr_PEz(d, omega_c)
        ]

        return (gamma_prefactor * sum(terms))

    def _gamma_parallel(self, d=None):
        """Calculates the parellel geometry friction,
        according to the formula in Lekkala et al 2013. Which formula?"""
        if d is None:
            d = self.Expt.d

        _prefactor = - \
            self.C_sphere(d) ** 2 * self.Expt.V_ts ** 2 / \
            (8 * pi * E_0 * self.Cant.omega_c)
        _integrand = lambda k: (
            k ** 2 * exp(-2 * k * d) * self._im_dielectric(k))
        _integral, _error = quad(_integrand, 0, inf)

        return _prefactor * _integral

    def calc_gamma_s(self, d=None):
        """Calculates sample-induced friction according to the
        appropriate equation, depending on whether the cantilever
        moves in the perpendicular or parallel geometry."""
        if d is None:
            d = self.Expt.d
        cant = self.Cant

        if cant.geometry_c == 'perpendicular':
            Gamma_s = self._gamma_perpendicular(d)
        elif cant.geometry_c == 'parallel':
            Gamma_s = self._gamma_parallel(d)
        else:
            raise ValueError("""geometry_c must be either 'parallel'\
                                or 'perpendicular'""")

        self.Gamma_s = Gamma_s
        return Gamma_s

    def calc_Pdf(self, omega, d=None):
        """Calculates the frequency power spectrum according
        to the appropriate equation, depending on whether the
        cantilever moves in the perpendicular or parallel geometry.
        Uses Lekkala 2013, et al., Eq. 6 and 10."""
        if d is None:
            d = self.Expt.d

        if self.Cant.geometry_c == 'perpendicular':
            Pdf = self._power_spectrum_perpendicular(omega, d)
        elif self.Cant.geometry_c == 'parallel':
            Pdf = self._power_spectrum_parallel(omega, d)
        else:
            raise ValueError("geometry_c must be either 'parallel'\
or 'perpendicular'")

        return Pdf

    def calc_power_spectrum(self, f_i=1e-5, f_f=1e3, n_pts=100, d=None):
        """Returns lists of f, Pdf(f) for frequencies in the given
        range. These can be used in plotting."""
        if d is None:
            d = self.Expt.d

        f_vals = np.logspace(log10(f_i), log10(f_f), n_pts)
        Pdf_vals = [self.calc_Pdf(2 * pi * f, d) for f in f_vals]
        self.f_vals = f_vals
        self.Pdf_vals = Pdf_vals
        return f_vals, Pdf_vals

    def calc_jitter(self, f_i=None, f_f=None, d=None):
        """Integrates the power spectrum over some
        frequency range f_i to f_f."""
        if f_i is None:
            f_i = self.Expt.jitter_f_i
        if f_f is None:
            f_f = self.Expt.jitter_f_f
        if d is None:
            d = self.Expt.d

        if f_i > f_f:
            raise ValueError("'jitter_f_i' must be less than 'jitter_f_f'.")

        calc_Pdf = lambda f: self.calc_Pdf(2 * pi * f, d)
        integral, error = quad(calc_Pdf, f_i, f_f)
        return integral

# def sum_csch(alpha):
# """This calculates the infinite sum, csch(alpha * n), for n = 1 to
# infinity. Testing the speed of this code, it is 1000 times slower
# than the equivalent version implemented using numpy functions."""
#     Nmax = int( np.max([4, np.ceil(20 / alpha)+1]) )
#     f = lambda n: sinh(alpha) / sinh(alpha * n)
#     terms = [f(n) for n in xrange(1,Nmax)]
#     return fsum(terms)

# def sum_csch(alpha, eps = 1e-8):
#     """This calculates the infinite sum, sinh(alpha) / sinh(alpha * n), for n = 1 to infinity."""
#     N_max_min = 1 + int(np.ceil(2 / alpha * arctanh( exp(-eps * alpha / sinh(alpha)) ) ))
#     N_max = np.max([5, N_max_min])
#     f = lambda n: sinh(alpha) / sinh(alpha * n)
#     terms = f(np.arange(1,N_max))
#     return np.sum(terms)

# def to_time(f):
#     return [f(alpha) for alpha in np.logspace(-4,1,6)]


# if __name__ == '__main__':
#     from timeit import timeit
# print(timeit("to_time(sum_csch)", setup="from __main__ import to_time, sum_csch", number = 1000))
#     print(timeit("to_time(sum_csch2)", setup="from __main__ import to_time, sum_csch2", number = 10000))
# from jittermodel.base import Cantilever, Sample, Experiment
"""Test Case - Tip Sample Capacitance"""
# ct1 = Cantilever(f_c = 46, Q = 2500, k_c = 8.5e9) # Should specify every parameter.
# st1 = Sample(h = 63e-3) # Here too
# et1 = Experiment() # Here three!
# simt1 = Simulation(ct1, st1, et1)
# print simt1.Cd2_sphere(100e-3)
# print simt1.C_sphere(100e-3)
"""Should give 6.9e-2, 5.1e-3 for the two capacitance derivatives.
This matches what Nik and Swapna reported in their paper
(Hoepker2011oct, Fig. 4)."""


"""Test Case - Im Dielectric """
# et2 = Experiment(V_ts = 3e3, d = 300e-3)
# ct2 = Cantilever(f_c = 65, Q = 2500, k_c = 3.5e9,
                   # geometry_c = 'perpendicular')
# st2 = Sample(h = 72e-3, V_g = 40e3, E_s1 = 3.4, mobility = 2.7e-4,
    # E_s2 = -0.05, E_i1 = 4.65)
# simt2 = Simulation(ct2, st2, et2)
# print simt2._im_dielectric(1)
# print simt2.corr_PP(0.1)
# print simt2.Samp.rho, simt2.Samp.sigma, simt2.Samp.kappa
# print simt2.Samp.C_i
"""Should give -0.003635 and 0.0019477. See Mathematica notebook 20130815."""


# print simt2.gamma_s()
"""Test Case - Perpendicular Geometry Friction"""
# et3 = Experiment(V_ts = 3e3, d = 300e-3)
# ct3 = Cantilever(f_c = 65, Q = 2500, k_c = 3.5e9,
    # geometry_c = 'perpendicular')
# st3 = Sample(h = 72e-3, E_s1 = 3.4, mobility = 2.7e-4, E_s2 = -0.05,
    # E_i1 = 3.4, E_i2 = -0.05,  = 40e3)
# simt3 = Simulation(ct3, st3, et3)
# print simt3.calc_gamma_s()