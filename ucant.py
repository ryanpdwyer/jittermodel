"""
UnitCantilever
Cantilever Object with units
2013-12-13
Ryan Dwyer

NOTE: To use units in your own file, import the unitregistry from
jittermodel.ucant!
"""

from __future__ import division
from numpy import pi
from autoassign import autoassign

from pint import UnitRegistry

u = UnitRegistry()

# Universal Constants
E_0 = 8.854e-12 * u.F / u.m 
k_B = 1.38065e-23 * u.J / u.K
q = 1.602e-19 * u.C


class UnitCantilever(Assigner):
    """Implement a Cantilever class with support for units."""
    @autoassign
    def __init__(self, f_c=50*u.kHz, k_c=3*u.N/u.m, Q=1000*u.dimensionless,
                 R_tip=40*u.nm, L_tip=15*u.um, theta_tip=16*u.degrees,
                 geometry_c='perpendicular'):
        self._check_number_inputs_positive()
        self._check_theta_less_than_90()
        self._check_geometry()

    @property
    def omega_c(self):
        """Return the angular resonance frequency of the cantilever."""
        return self.f_c * 2 * pi

    @property
    def Gamma_i(self):
        """Return the cantilever's intrinsic dissipation."""
        return (self.k_c / (self.omega_c * self.Q)).ito(u.pN * u.s / u.m)

    def _check_geometry(self):
        """Raises an error if the geometry of the cantilever geometry_c
        is not either 'perpendicular' or 'parallel'."""
        if self.geometry_c is not ('perpendicular' or 'parallel'):
            raise ValueError("""geometry_c must be either 'perpendicular'\
                or 'parallel'""")

    def _check_number_inputs_positive(self):
        """Returns a ValueError if the number inputs are not positive."""
        greater_than_zero = ('f_c', 'k_c', 'Q',
                             'R_tip', 'L_tip', 'theta_tip')
        
        for attr in greater_than_zero:
            if self.lookup(attr).magnitude <= 0:
                raise ValueError("""The attribute '{attr}'\
must be positive.""".format(attr=attr))

    def _check_theta_less_than_90(self):
        """Return a ValueError if theta_tip >= 90 degrees
        since this is unphysical."""
        if self.theta_tip >= 90 * u.degrees:
            raise ValueError("'theta_tip' must be less than 90 degrees.")

    def F_min(self, T, bandwidth=1*u.Hz):
        """Return the thermally limited minimum detectable
        force (pN).

        The optional bandwidth parameter allows determining
        a miniumun force over a broader or  narrower bandwidth
        than 1 Hz."""
        return ((4 * k_B * self.Gamma_i * T * bandwidth) ** 0.5).ito(u.pN)

    def __str__(self):
        """Write out the cantilever as its
        three most important parameters: resonance frequency,
        spring constant, and quality factor."""
        return "f_c = {self.f_c:!p}, k_c = {self.k_c:!p}, Q = {self.Q:!p}\
Gamma_i = {self.Gamma_i:!p}".format(self=self)

    def __repr__(self):
        """Return a representation of the cantilever which can be used
        to generate an equivalent copy of the cantilever."""
        return "Cantilever(f_c = {self.f_c}, k_c = {self.k_c}, Q = {self.Q},\
R_tip = {self.R_tip}, L_tip = {self.L_tip},\
theta_tip = {self.theta_tip},\
geometry_c = '{self.geometry_c}')".format(self=self)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__