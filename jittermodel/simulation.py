# -*- coding: utf-8 -*-
"""
Simulation
==========

The simulation code, but with units.

It is more difficult to add units to the simulation code than
the simple base classes. This is because there are complicated integrals,
derivatives, sums, and other pieces of code which are both numerically
sensitive and demanding of execution time. The most naive addition of
units to the simulation is likely to be far too slow, as well as numerically
unstable.

Therefore, the best tactic to ensure stable, fast execution of this code is
to use some sort of input type checking in order to convert unitted quantities
to the appropriate unitless variable.

In order to do this, I should start with a dictionary containing the unitless
version of each variable I will encounter. This dictionary should be able to
transform any unit into the appropriate unitless variable needed in the
simulation. From there, I can simply use this dictionary to find the correct
unit, and keep all the math in this code unitfree.

Cython
------

The innermost loop (``Simulation._im_dielectric``) is now implemented in cython;
this speeds up the model 2 (the model which was
already working well) inner loop by a factor of 20, and the model 1 inner loop
by a factor of 40, enough to make their speeds roughly equal
(model 2 is now about one third faster).

If additional performance gains are desired, (``Simulation._corr_integrand``
could also be converted to Cython code). There is still very significant
overhead in all of the critical functions being python functions rather than
Cython; overall, we achieved a factor of 2 to 4 speed-up moving to Cython.

The scipy quad routine could also be converted to Cython; this could be done
directly by ``cdef extern from "gsl/gsl_file.h"``, or using the ``CythonGSL``
library.
"""

from __future__ import division
import numpy as np
from numpy import (pi, sinh, cosh, tanh, arccosh,
                   exp, log10, arctanh)
from scipy.integrate import quad
from scipy.special import jn
from scipy.misc import derivative
import math
from copy import copy
import jittermodel._sim
from jittermodel import u, q2unitless

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
    # For alpha greater than 37, the sum is 1 to within double precision,
    # so short-circuit the calculation to avoid overflow errors.
    if alpha > 37:
        return 1
    else:
        summand = lambda n: sinh(alpha) / sinh(alpha * n)

        N_max = max(4, int_sum_sinh(alpha, eps))
        terms = summand(np.arange(1, N_max + 1))
        return math.fsum(terms)


def _alpha(d, R_tip, h, E_s1):
    return arccosh(1 + d / R_tip + h / (E_s1 * R_tip))


class SphereCapacitance(object):
    """Implement a simple model of AFM cantilever-sample capacitance
    contributed by a spherical tip."""

    def __init__(self, simulation):
        """Initialize with an owner class containing the other methods
        and data used by the capacitance model.

        See http://stackoverflow.com/q/6556744/2823213"""
        self.sim = simulation

    def C(self, d=None):
        """Capacitance between a sphere and a thin sample, calculated
        according to a modified version of the equation from Brus et al.,
        http://dx.doi.org/10.1021/jp0265438. The capacitance is
        approximated by truncating the exact infinite sum."""
        sim = self.sim
        if d is None:
            d = sim.Expt.d
        if d < 0:
            # Prevent numerical derivatives from giving an error.
            d = 0

        cant = sim.Cant
        samp = sim.Samp

        alpha = _alpha(d, cant.R_tip, samp.h, samp.E_s1)

        return (4 * pi * sim.E_0 * cant.R_tip * sum_sinh(alpha))

    def Cd(self, d=None):
        """Returns the numerical derivative of the sphere capacitance
        C at a distance d. Uses scipy to calculate the derivative."""
        if d is None:
            d = self.sim.Expt.d
        return derivative(self.C, d, dx=1e-5, n=1)

    def Cd2(self, d=None):
        """Returns the second derivative of the sphere capacitance
        C at a distance d. Uses scipy to calculate the derivative."""
        if d is None:
            d = self.sim.Expt.d
        return derivative(self.C, d, dx=1e-5, n=2)


def _eta(k, kappa, E_s, D, omega):
    return np.sqrt(k**2 + kappa**2 / E_s + omega/D*1j)


def _lambda(k, eta, E_eff, E_s):
    """Helper function for calculating the correlation integrand.
    See Lekkala, et al., 2013, Eq. 19"""
    return k/eta*(1 - E_eff/E_s)


def _thetaI(k, h_s, alpha, Lambda, eta, E_s, E_eff):
    """Large term in brackets at the bottom of page 17 of
    Lekkala, 2013, *et al.*."""
    if k*h_s > 8:
        return (E_s / E_eff * (-Lambda +
            (1 + alpha - Lambda*(1 - 8*exp(-h_s*(eta + k))
                                 +4*Lambda*exp(-2*eta*h_s)))/
            (alpha - Lambda + 1)))
    else:
        sk = sinh(k * h_s)
        sn = sinh(eta * h_s)
        ck = cosh(k * h_s)
        cn = cosh(eta * h_s)
        return (E_s / E_eff * (-Lambda/ tanh(eta*h_s) +
                (sk*sn + alpha*ck*sn - Lambda * (ck*cn - 2 + Lambda*sk/sn)) /
                (ck * sn + alpha*sk*sn - Lambda*sk*cn)))


def _thetaII(k, h, E_s, E_d, E_eff, Lambda):
    return (E_s / E_d * ((E_eff + (1 - Lambda) * E_d / tanh(k*h)) /
                         (E_eff / tanh(k*h) + (1 - Lambda) * E_d)))

# TODO: test this function!!!

def _im_dielectric(k, h_diel, h_trans, E_s, E_i, mu, omega, rho, T, k_B, q, E_0,
                   model):
    sigma = mu * rho * q
    E_eff = E_s - sigma / (E_0 * omega) * 1j
    kappa = (2 * rho * q ** 2 / (E_0 * k_B * T)) ** 0.5
    diff = mu * k_B * T / q
    if int(model) == 1:
        E_d = E_i
        h = h_diel + h_trans
        alpha = E_eff / E_d
        eta = _eta(k, kappa, E_s, diff, omega)
        Lambda = _lambda(k, eta, E_eff, E_s)
        theta = _thetaI(k, h, alpha, Lambda, eta, E_s, E_eff)
    elif int(model) == 2:
        E_d = E_s
        h = h_diel
        eta = _eta(k, kappa, E_s, diff, omega)
        Lambda = _lambda(k, eta, E_eff, E_s)
        theta = _thetaII(k, h, E_s, E_d, E_eff, Lambda)
    else:
        raise ValueError("Model must be 1 or 2, not {0}".format(model))

    result = (E_s - theta) / (E_s + theta)
    return result.imag


class Simulation(object):
    """This calculates experimental parameters such as capacitance
    and sample-induced friction for a given cantilever, sample,
    and experiment.
    """
    units = {"[mass]": u.pg, "[length]": u.um, "[time]": u.ms,
             "[current]": u.aC / u.ms, "[temperature]": u.K, "[angle]": u.rad}

    # Store fundamental constants in the units defined above.
    E_0 = q2unitless(jittermodel.E_0, units)
    q = q2unitless(jittermodel.q, units)
    k_B = q2unitless(jittermodel.k_B, units)

    def __init__(self, cantilever, sample, experiment, model=2):
        """Initialize the simulation with the values from the given
        cantilever, sample and experiment. It also calculates
        parameters used in the simulation"""
        self.UCant = copy(cantilever)
        self.USamp = copy(sample)
        self.UExpt = copy(experiment)
        self.UCant._unitless_units = self.units
        self.USamp._unitless_units = self.units
        self.UExpt._unitless_units = self.units

        self.Cant = self.UCant.to_unitless()
        self.Samp = self.USamp.to_unitless()
        self.Expt = self.UExpt.to_unitless()

        self.model = model

        self.Sphere = SphereCapacitance(self)

        self.func_dict = {'friction': self.calc_gamma_s,
                          'jitter': self.calc_jitter,
                          'power spectrum': self.calc_Pdf}

    def assign(self, attr, val):
        """Assign the attribute 'attr' to value 'val'. Behind the
        scenes, this function finds which property came from which
        class, and then does the appropriate assignment."""
        if attr == 'omega':
            # Unlike other attributes, omega is owned by the simulation.
            # TODO: Why is kHz hard-coded here!
            self.omega = val.to(u.kHz).magnitude
        if attr == 'model':
            self.model = val
        else:
            found = False
            for item in (self.UCant, self.USamp, self.UExpt):
                # TODO: searching in a hidden, private attribute is bad form;
                # This could be nicely implemented with __contains__, which
                # would make this check look far saner.
                if attr in item._all_attributes:
                    found = True
                    item.assign(attr, val)
                    # Reinitialize the unitless version of the item
                    item.to_unitless()
                    break

            if not found:
                raise AttributeError(
                    "The attribute {0} was not found.".format(attr))

    def lookup(self, attr):
        """Looks for an attribute in the cantilever, sample, or
        experiment namespaces, and returns it."""
        for item in (self.Cant, self.Samp, self.Expt):
            if attr in item._all_attributes:
                return item.lookup(attr)
        return '{attr} not found'.format(attr=attr)

    # Correlation Function related functions. #

    def _prefactor(self, omega):
        """This is the prefactor for correlation function integrals.
        See Eq. 16 in Lekkala et al. 2013."""
        return (-self.k_B * self.Samp.T / (4 * pi * self.E_0 * omega))

    def _im_dielectric(self, k, omega):
        """Computes the complicated expression containing
        :math:`\\epsilon_\\mathrm{rel}(\\omega)` and
        :math:`\\theta(k,\\omega)` in the correlation function
        from Lekkala et al. 2013. Currently uses model II, but
        support for model I can be easily added. See Eq. 16.

        TODO: Add unittest, even with just hard-coded comparisons
        to Mathematica."""

        if self.model not in (1, 2, 20):
            raise ValueError("Model must be either 1, 2 or 20.")
        elif self.model == 1 or self.model == 2:
            return jittermodel._sim._im_dielectric_c(k, self.Samp.h_diel, self.Samp.h_trans,
                            self.Samp.E_s, self.Samp.E_i,
                            self.Samp.mobility, omega, self.Samp.rho, self.Samp.T,
                            self.k_B, self.q, self.E_0,
                            self.model)
        else:
            return jittermodel._sim._im_dielectric_c(k, self.Samp.h_d,
                self.Samp.h_s, self.Samp.E_s, self.Samp.E_i,
                self.Samp.mobility, omega, self.Samp.rho, self.Samp.T,
                self.k_B, self.q, self.E_0,
                self.model)


    def _corr_integrand(self, r1, r2, z1, z2, k, omega, n=0):
        """The integrand of the correlation function from
        Lekkala et al., 2013, at the cantilever resonance frequency."""
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

    # Correlation functions defined in terms of the function _base_corr

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

    # Power Spectrum calculations
    def _power_spectrum_parallel(self, omega, d=None):
        """Calculates the power spectrum of frequency fluctuations
        at the angular frequency omega for a parallel geometry
        cantilever. See Eq. 6, Lekkala 2013."""
        if d is None:
            d = self.Expt.d
        _prefactor = (self.V_ts * self.Sphere.C(d) / self.k_c) ** 2

        return (_prefactor * self.corr_ExxExx(d, omega))

    def _power_spectrum_perpendicular(self, omega, d=None):
        """Calculates the power spectrum of frequency fluctuations
        at angular frequency omega for a perpendicular geometry
        cantilever. See Eq. 9, Lekkala 2013."""
        if d is None:
            d = self.Expt.d
        power_spectrum_prefactor = (
            self.Cant.f_c * self.Expt.V_ts / self.Cant.k_c) ** 2
        c = self.Sphere.C(d)
        cd = self.Sphere.Cd(d)
        cdd = self.Sphere.Cd2(d)

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

        gamma_prefactor = self.Expt.V_ts ** 2 / (self.k_B * self.Samp.T)

        c = self.Sphere.C(d)
        cd = self.Sphere.Cd(d)

        terms = np.array([
            c ** 2 * self.corr_EzEz(d, omega_c),
            cd ** 2 * self.corr_PP(d, omega_c),
            -2 * c * cd * self.corr_PEz(d, omega_c)
        ])

        return (gamma_prefactor * np.sum(terms))

    def _gamma_parallel(self, d=None):
        """Calculates the parellel geometry friction,
        according to the formula in Lekkala et al 2013. Which formula?"""
        if d is None:
            d = self.Expt.d

        _prefactor = (-self.Sphere.C(d) ** 2 * self.Expt.V_ts ** 2 /
                      (8 * pi * self.E_0 * self.Cant.omega_c))
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
            raise ValueError(
                "geometry_c must be either 'parallel' or 'perpendicular'")

        self.Gamma_s = Gamma_s
        return Gamma_s

    def calc_Pdf(self, omega=None, d=None):
        """Calculates the frequency power spectrum according
        to the appropriate equation, depending on whether the
        cantilever moves in the perpendicular or parallel geometry.
        Uses Lekkala 2013, et al., Eq. 6 and 10."""
        if d is None:
            d = self.Expt.d
        if omega is None:
            omega = self.omega

        if self.Cant.geometry_c == 'perpendicular':
            Pdf = self._power_spectrum_perpendicular(omega, d)
        elif self.Cant.geometry_c == 'parallel':
            Pdf = self._power_spectrum_parallel(omega, d)
        else:
            raise ValueError(
                "geometry_c must be either 'parallel' or 'perpendicular'")

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


"""Test Case - Tip Sample Capacitance"""
# ct1 = Cantilever(f_c = 46, Q = 2500, k_c = 8.5e9) # Should specify every parameter.
# st1 = Sample(h = 63e-3) # Here too
# et1 = Experiment() # Here three!
# simt1 = Simulation(ct1, st1, et1)
# print simt1.Sphere.Cd2(100e-3)
# print simt1.Sphere.C(100e-3)
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
