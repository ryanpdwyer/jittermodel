# -*- coding: utf-8 -*-
"""
Base
====

This file defines the individual classes used to build a simulation, including
a ``Cantilever``, ``Transistor`` and ``Experiment``

NOTE: To use units in your own file, import the
unitregistry (``u``) from ``jittermodel``!

Class structure
---------------

These classes should be conceptually much simpler; they should be simple classes,
that perhaps contain additional logic buried underneath.

For example, the cantilever should *only* expose simple attributes, 

"""

from __future__ import division
from numpy import pi
from jittermodel import (u, UnitAssigner, NoUnitAssigner,
                         get_defaults, E_0, k_B, q)


class SimpleCantilever(UnitAssigner):
    """Implement a simple unit cantilever, and require all inputs."""

    k_B = k_B

    def __init__(self, f_c, k_c, Q):
        """Initialize the cantilever."""
        self.f_c = f_c
        self.k_c = k_c
        self.Q = Q

        self._default_units = {'f_c': u.kHz, 'k_c': u.N/u.m,
                               'Q': u.dimensionless}
        self._check_dimensionality_units()
        self._check_number_inputs_positive()

    # Properties of the cantilver
    @property
    def omega_c(self):
        """Return the angular resonance frequency of the cantilever."""
        return self.f_c * 2 * pi

    @property
    def Gamma_i(self):
        """Return the cantilever's intrinsic dissipation."""
        return (self.k_c / (self.omega_c * self.Q))

    def F_min(self, T, bandwidth=1*u.Hz):
        """Return the thermally limited minimum detectable
        force (pN).

        The optional bandwidth parameter allows determining
        a miniumun force over a broader or narrower bandwidth
        than 1 Hz."""
        return ((4 * self.k_B * self.Gamma_i * T * bandwidth) ** 0.5)

    # Representations of the cantilever
    def __str__(self):
        """Write out the cantilever as its most important parameters:
        resonance frequency, spring constant, quality factor and
        intrinsic friction."""
        return "f_c = {self.f_c}, k_c = {self.k_c}, Q = {self.Q}\
Gamma_i = {self.Gamma_i}".format(self=self)

    def __repr__(self):
        """Return a representation of the cantilever. Rounds the cantilever
        to 9 digits, so eval(repr(cantilever)) is not necessarily equal to
        cantilever."""
        return "Cantilever(f_c = {self.f_c}, k_c = {self.k_c}, \
Q = {self.Q})".format(self=self)


class Cantilever(SimpleCantilever):
    """Implement a Cantilever class with support for units."""

    def __init__(self, f_c=50*u.kHz, k_c=3*u.N/u.m, Q=1000*u.dimensionless,
                 R_tip=40*u.nm, L_tip=15*u.um, theta_tip=16*u.degrees,
                 geometry_c='perpendicular'):
        """Initialize the cantilever."""
        self.UnitlessClass = UnitlessCantilever

        self.f_c = f_c
        self.k_c = k_c
        self.Q = Q
        self.R_tip = R_tip
        self.L_tip = L_tip
        self.theta_tip = theta_tip
        self.geometry_c = geometry_c

        self._default_units = {'f_c': u.kHz, 'k_c': u.N/u.m,
                               'Q': u.dimensionless, 'R_tip': u.nm,
                               'L_tip': u.um, 'theta_tip': u.degrees}
        self._check_dimensionality_units()
        self._check_number_inputs_positive()
        self._check_theta_less_than_90()
        self._check_geometry()

    # Functions to check the inputs to the cantilever.
    def _check_geometry(self):
        """Raises an error if the geometry of the cantilever geometry_c
        is not either 'perpendicular' or 'parallel'."""
        if self.geometry_c is not ('perpendicular' or 'parallel'):
            raise ValueError("""geometry_c must be either 'perpendicular'\
                or 'parallel'""")

    def _check_theta_less_than_90(self):
        """Return a ValueError if theta_tip >= 90 degrees
        since this is unphysical."""
        if self.theta_tip >= 90 * u.degrees:
            raise ValueError("'theta_tip' must be less than 90 degrees.")

    # Representations of the cantilever
    def __repr__(self):
        """Return a representation of the cantilever. Rounds the cantilever
        to 9 digits, so eval(repr(cantilever)) is not necessarily equal to
        cantilever."""
        return "Cantilever(f_c = {self.f_c}, k_c = {self.k_c}, Q = {self.Q},\
R_tip = {self.R_tip}, L_tip = {self.L_tip},\
theta_tip = {self.theta_tip},\
geometry_c = '{self.geometry_c}')".format(self=self)


class UnitlessCantilever(NoUnitAssigner, Cantilever):
    pass


class Experiment(UnitAssigner):
    """Stores parameters set by the experimenter. Now with units!"""

    def __init__(self, d=100 * u.nm, V_ts=5 * u.V,
                 jitter_f_i=0.2 * u.Hz, jitter_f_f=3 * u.Hz):
        self.UnitlessClass = UnitlessExperiment

        self.d = d
        self.V_ts = V_ts
        self.jitter_f_i = jitter_f_i
        self.jitter_f_f = jitter_f_f

        self._default_units = {'d': u.nm, 'V_ts': u.V, 'jitter_f_i': u.Hz,
                               'jitter_f_f': u.Hz}

        self._check_dimensionality_units()
        self._check_number_inputs_positive()
        # Check for errors in the experimental parameters
        if self.V_ts < 0 * u.V:
            raise ValueError("The voltages 'V_g' and 'V_ts' must be positive.")
        if self.jitter_f_i > self.jitter_f_f:
            raise ValueError("'jitter_f_i' must be less than 'jitter_f_f'.")

    def __str__(self):
        """A nice string representation of the experiment."""
        return """Tip-sample: {self.d}, {self.V_ts}""".format(self=self)


class UnitlessExperiment(NoUnitAssigner, Experiment):
    pass


class Transistor(UnitAssigner):
    """A transistor sample, now with units."""
    E_0 = E_0
    k_B = k_B
    q = q

    def __init__(self, semiconductor='TPD',
                 h=70 * u.nm, h_trans=1 * u.nm, h_i=300 * u.nm,
                 E_s1=3.5, E_s2=-0.0005,
                 E_i1=4.65, E_i2=0,
                 mobility=3e-6 * u.cm ** 2 / u.V / u.s, T=298 * u.K,
                 V_g=10 * u.V, rho=None):
        """Initialize the sample with all of the experimentally
        relevant sample parameters."""
        self.UnitlessClass = UnitlessTransistor

        self.semiconductor = semiconductor
        self.h = h
        self.h_trans = h_trans
        self.h_i = h_i
        self.E_s1 = E_s1
        self.E_s2 = E_s2
        self.E_i1 = E_i1
        self.E_i2 = E_i2
        self.mobility = mobility
        self.T = T

        self._check_V_g_rho_defined(V_g, rho)

        self._default_units = {'h': u.nm, 'h_trans': u.nm, 'h_i': u.nm,
                      'mobility': u.cm ** 2 / u.V / u.s, 'T': u.K,
                      'V_g': u.V, 'rho': u.cm ** -3}

        self._check_dimensionality_units()
        self._check_number_inputs_positive()

    def _check_V_g_rho_defined(self, V_g, rho):
        """Checks to determine whether one, both, or none of V_g and rho
        were given when the sample was initialized, and properly assigns
        V_g and rho or throws an error as appropriate."""

        default_dict = get_defaults(self.__init__)

        if rho is None:
            self.V_g = V_g
        elif V_g == default_dict['V_g']:
            self.rho = rho
        else:
            raise ValueError("The provided values of 'V_g' and 'rho' are \
incompatible. Only specify one of 'V_g' or 'rho' when defining a Sample.")

    @property
    def diff(self):
        """Diffusion constant defined according to the Einstein relation."""
        return self.mobility * self.k_B * self.T / self.q

    @property
    def C_i(self):
        """Capacitance per unit area between the
        transistor gate and sample."""
        return self.E_i1 * self.E_0 / self.h_i

    @property
    def h_diel(self):
        """Layer of the sample which is acting as a pure dielectric layer."""
        return self.h - self.h_trans

    @property
    def E_s(self):
        """Total dielectric constant, adjusted for conductivity."""
        return self.E_s1 + self.E_s2*1j

    @property
    def E_i(self):
        """Total dielectric constant of the insulator layer,
        assuming the same lossy part of the spectrum as the
        sample layer."""
        return self.E_i1 + self.E_i2*1j

    #---------------------------------------------------------
    # The gate voltage, charge density pair are defined below.
    @property
    def V_g(self):
        """Get the gate voltage."""
        return self._V_g

    @V_g.setter
    def V_g(self, value):
        """Set the gate voltage. Updates the semiconductor
        carrier density hidden variable _rho to match the
        new gate voltage."""
        self._V_g = value
        self._rho = self.C_i * self._V_g / (self.h_trans * self.q)

    @property
    def rho(self):
        """Get the semiconductor carrier density rho."""
        return self._rho

    @rho.setter
    def rho(self, value):
        """Set the semiconductory carrier density rho.
        Updates the gate voltage hidden variable _V_g to match the new
        carrier density."""
        self._rho = value
        self._V_g = self.q * self._rho * self.h_trans / self.C_i

    #---------------------------------------------------------
    # Relevant properties derived from gate voltage / charge density.
    @property
    def sigma(self):
        """The conductivity sigma of the sample."""
        return self.mobility * self.rho * self.q

    @property
    def kappa(self):
        """Define the inverse Debye length screening length kappa,
        used in the Lekkala Loring theory. See Lekkala et al.,
        p4, at http://dx.doi.org/10.1063/1.4754602."""
        return (2 * self.rho * self.q ** 2 /
                (self.E_0 * self.k_B * self.T)) ** 0.5

    def E_eff(self, omega):
        """Defines the effective dielectric constant,
        which corrects for the effect of Ohm's law (carrier drift),
        at a particular angular frequency. See Eq. 14 in
        Lekkala et al., 2013."""
        return self.E_s - self.sigma / (self.E_0 * omega) * 1j

    def __str__(self):
        """Write out the sample, using its semiconductor material,
        height and mobility."""
        return "{self.semiconductor:P}  {self.h_nm:P}\
mobility {self.mobility_cm2:P}".format(self=self)

    def __repr__(self):
        return "Transistor(semiconductor = '{self.semiconductor}', \
h = {self.h}, h_trans = {self.h_trans}, \
h_i = {self.h_i}, E_s1 = {self.E_s1}, \
E_s2 = {self.E_s2}, E_i1 = {self.E_i1}, \
E_i2 = {self.E_i2}, mobility = {self.mobility}, \
T = {self.T}, V_g = {self.V_g})".format(self=self)


class UnitlessTransistor(NoUnitAssigner, Transistor):
    pass


class BareTransistor(UnitAssigner):
    """A transistor sample, now with units."""
    E_0 = E_0
    k_B = k_B
    q = q

    def __init__(self, semiconductor='Si',
                 h_d=300 * u.nm, h_s=500 * u.um,
                 E_d1=4.65, E_d2=-0.005,
                 E_s1=11.65, E_s2=-0.05,
                 mobility=1 * u.cm ** 2 / u.V / u.s,
                 T=298 * u.K,
                 rho=None):
        """Initialize the sample with all of the experimentally
        relevant sample parameters."""
        self.UnitlessClass = UnitlessTransistor

        self.semiconductor = semiconductor
        self.h_d = h_d
        self.h_s = h_s
        self.E_s1 = E_s1
        self.E_s2 = E_s2
        self.E_d1 = E_d1
        self.E_d2 = E_d2
        self.mobility = mobility
        self.T = T
        self.rho = rho

        self._default_units = {'h_d': u.nm, 'h_s': u.nm,
                      'mobility': u.cm ** 2 / u.V / u.s, 'T': u.K,
                      'rho': u.cm ** -3}

        self._check_dimensionality_units()
        self._check_number_inputs_positive()

    def _check_V_g_rho_defined(self, V_g, rho):
        """Checks to determine whether one, both, or none of V_g and rho
        were given when the sample was initialized, and properly assigns
        V_g and rho or throws an error as appropriate."""

        default_dict = get_defaults(self.__init__)

        if rho is None:
            self.V_g = V_g
        elif V_g == default_dict['V_g']:
            self.rho = rho
        else:
            raise ValueError("The provided values of 'V_g' and 'rho' are \
incompatible. Only specify one of 'V_g' or 'rho' when defining a Sample.")

    @property
    def diff(self):
        """Diffusion constant defined according to the Einstein relation."""
        return self.mobility * self.k_B * self.T / self.q

    @property
    def E_s(self):
        """Total dielectric constant, adjusted for conductivity."""
        return self.E_s1 + self.E_s2*1j

    @property
    def E_d(self):
        """Total dielectric constant of the insulator layer,
        assuming the same lossy part of the spectrum as the
        sample layer."""
        return self.E_d1 + self.E_d2*1j

    #---------------------------------------------------------
    # Relevant properties derived from gate voltage / charge density.
    @property
    def sigma(self):
        """The conductivity sigma of the sample."""
        return self.mobility * self.rho * self.q

    @property
    def kappa(self):
        """Define the inverse Debye length screening length kappa,
        used in the Lekkala Loring theory. See Lekkala et al.,
        p4, at http://dx.doi.org/10.1063/1.4754602."""
        return (2 * self.rho * self.q ** 2 /
                (self.E_0 * self.k_B * self.T)) ** 0.5

    def E_eff(self, omega):
        """Defines the effective dielectric constant,
        which corrects for the effect of Ohm's law (carrier drift),
        at a particular angular frequency. See Eq. 14 in
        Lekkala et al., 2013."""
        return self.E_s - self.sigma / (self.E_0 * omega) * 1j

    def __str__(self):
        """Write out the sample, using its semiconductor material,
        height and mobility."""
        return "{self.semiconductor:P}  {self.h_nm:P}\
mobility {self.mobility_cm2:P}".format(self=self)

    def __repr__(self):
        return "Transistor(semiconductor = '{self.semiconductor}', \
h = {self.h}, h_trans = {self.h_trans}, \
h_i = {self.h_i}, E_s1 = {self.E_s1}, \
E_s2 = {self.E_s2}, E_i1 = {self.E_i1}, \
E_i2 = {self.E_i2}, mobility = {self.mobility}, \
T = {self.T}, V_g = {self.V_g})".format(self=self)


class UnitlessTransistor(NoUnitAssigner, Transistor):
    pass
