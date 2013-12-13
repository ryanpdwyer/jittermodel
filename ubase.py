"""
UnitCantilever
Cantilever Object with units
2013-12-13
Ryan Dwyer

NOTE: To use units in your own file, import the unitregistry (u)
from jittermodel!
"""

from __future__ import division
from numpy import pi
from autoassign import autoassign
from . import u, UnitAssigner


# Universal Constants
E_0 = 8.854e-12 * u.F / u.m
k_B = 1.38065e-23 * u.J / u.K
q = 1.602e-19 * u.C


class UnitCantilever(UnitAssigner):
    """Implement a Cantilever class with support for units."""
    @autoassign
    def __init__(self, f_c=50*u.kHz, k_c=3*u.N/u.m, Q=1000*u.dimensionless,
                 R_tip=40*u.nm, L_tip=15*u.um, theta_tip=16*u.degrees,
                 geometry_c='perpendicular'):
        """Initialize the cantilever."""
        self.units = {'f_c': u.kHz, 'k_c': u.N/u.m, 'Q': u.dimensionless,
                      'R_tip': u.nm, 'L_tip': u.um, 'theta_tip': u.degrees}

        self._check_dimensionality_units()
        self._check_number_inputs_positive()
        self._check_theta_less_than_90()
        self._check_geometry()

    # Properties of the cantilver
    @property
    def omega_c(self):
        """Return the angular resonance frequency of the cantilever."""
        return self.f_c * 2 * pi

    @property
    def Gamma_i(self):
        """Return the cantilever's intrinsic dissipation."""
        return (self.k_c / (self.omega_c * self.Q)).ito(u.pN * u.s / u.m)

    def F_min(self, T, bandwidth=1*u.Hz):
        """Return the thermally limited minimum detectable
        force (pN).

        The optional bandwidth parameter allows determining
        a miniumun force over a broader or  narrower bandwidth
        than 1 Hz."""
        return ((4 * k_B * self.Gamma_i * T * bandwidth) ** 0.5).ito(u.pN)

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
        return "Cantilever(f_c = {self.f_c}, k_c = {self.k_c}, Q = {self.Q},\
R_tip = {self.R_tip}, L_tip = {self.L_tip},\
theta_tip = {self.theta_tip},\
geometry_c = '{self.geometry_c}')".format(self=self)


class UnitExperiment(UnitAssigner):
    """Stores parameters set by the experimenter. Now with units!"""
    @autoassign
    def __init__(self, d=100 * u.nm, V_ts=5 * u.V,
                 jitter_f_i=0.2 * u.Hz, jitter_f_f=3 * u.Hz):
        self.units = {'d': u.nm, 'V_ts': u.V, 'jitter_f_i': u.Hz,
                      'jitter_f_f': u.Hz}

        self._check_dimensionality_units()
        self._check_number_inputs_positive()
        # Check for errors in the experimental parameters
        if self.V_ts < 0:
            raise ValueError("The voltages 'V_g' and 'V_ts' must be positive.")
        if self.jitter_f_i > self.jitter_f_f:
            raise ValueError("'jitter_f_i' must be less than 'jitter_f_f'.")

    def __str__(self):
        """A nice string representation of the experiment."""
        return """Tip-sample: {self.d}, {self.V_ts}""".format(
            self=self)
