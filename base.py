

"""
jittermodel.base

Created by Ryan Dwyer on 2013-05-18.

This file collects all of the functions for running capacitance
calculations. The classes Cantilever, Experiment and Sample
define the various parameters of the Simulation. Simulation
collects parameters from Cantilever, Experiment and Sample, and
has methods to calculate sample-induced friction and jitter
according to the Lekkala, Loring theory (Lekkala *et al.*,
2013, *in press*).

GenerateFrictionPlot is the object that can be used to generate
various plots, varying different cantilever, experiment or sample
parameters.

These are the units used in the following calculations. They are
chosen so that the integrals and derivatives that must be
evaluated are of order 1.

Units:
frequency: kHz
time: ms
distance: um
charge: aC
mass: pg
capacitance: fF
capacitance / length: fF / um = 1e-9 F/m
capacitance / length^2: fF / um^2 = 1e-3 F/m^2

voltage: mV

Friction: pg / ms = 1 pN s / m
From swapna's paper, we had about 10 pN s/ m = 1e-11 N s/m = 1e-5 uN s / m

mobility: um^2 / mV ms = (1e-6 m^2 / Vs)
Spring constant: pg / ms^2 = 1e-9 N/m
energy: J = Nm =
            kg m^2 / s^2 = pg um^2 / ms^2 = -15 - 12 -(-6) = 1e-21 J = zJ?


Boltzmann: 1.38065e-2
q = 1.602e-1
E_0 = 8.854e-3
"""

from __future__ import division
from numpy import pi
from autoassign import autoassign
from . import Assigner

# Universal Constants
E_0 = 8.854e-3
k_B = 1.38065e-2
q = 1.602e-1


class Cantilever(Assigner):
    """This defines a cantilever. The user sets the parameters listed below,
    with default values in parentheses. All values are input in the units
    listed at the top of the page, and in the list below.

    f_c
        The resonance frequency of the cantilever (50 kHz).

    k_c
        The spring constant of the cantilever (3e9 nN/m = 3 N/m).

    Q
        The quality factor of the cantilever (1000).

    R_tip
        The radius of the cantilever tip (0.04 :math:`\mu`m).

    L_tip
        The length of the cantilever cone (15 :math:`\mu`m).

    ThetaDegrees_tip
        The half-angle of the cantilever tip cone (:math:`16^\circ`).

    geometry_c
        The direction that the cantilever moves with respect to the sample.
        Ultra-sensitive cantilevers vibrate 'parallel' and commercial
        cantilevers vibrate 'perpendicular'.

In addition, the following quantities are defined in the class
    """
    @autoassign
    def __init__(self, f_c=50, k_c=3e9, Q=1000,
                 R_tip=0.04, L_tip=15, ThetaDegrees_tip=16,
                 geometry_c='perpendicular'):
        """Initialize the cantilever class with values, check to
        make sure the geometry is valid, and derive other useful
        cantilever properties from the given values.
        See class help for more initalization information (units).
        @autoassign automatically assigns the input quantities to self.
        See http://code.activestate.com/recipes/551763/ for more
        information."""

        self._check_number_inputs_positive()
        self._check_theta_less_than_90()
        self._check_geometry()

    # Define other useful cantilever properties.
    @property
    def Theta_tip(self):
        """Return the tip angle in radians"""
        return self.ThetaDegrees_tip * pi / 180

    @property
    def omega_c(self):
        """Return the angular resonance frequency of the cantilever."""
        return self.f_c * 2 * pi

    @property
    def Gamma_i(self):
        """Return the cantilever's intrinsic dissipation."""
        return self.k_c / (self.omega_c * self.Q)

    def _check_geometry(self):
        """Raises an error if the geometry of the cantilever geometry_c
        is not either 'perpendicular' or 'parallel'."""
        if self.geometry_c is not ('perpendicular' or 'parallel'):
            raise ValueError("""geometry_c must be either 'perpendicular'\
                or 'parallel'""")

    def _check_number_inputs_positive(self):
        """Returns a ValueError if the number inputs are not positive."""
        greater_than_zero = ('f_c', 'k_c', 'Q',
                             'R_tip', 'L_tip', 'ThetaDegrees_tip')
        
        for attr in greater_than_zero:
            if self.lookup(attr) <= 0:
                raise ValueError("""The attribute '{attr}'\
                                 must be positive.""".format(attr=attr))

    def _check_theta_less_than_90(self):
        """Return a ValueError if ThetaDegrees_tip >= 90,
        since this is unphysical."""
        if self.ThetaDegrees_tip >= 90:
            raise ValueError("'ThetaDegrees_tip' must be less than 90 degrees.")

    def F_min(self, T, bandwidth=0.001):
        """Return the thermally limited minimum detectable
        force (pN), or pN / Hz^0.5.

        The optional bandwidth parameter allows determining
        a miniumun force over a broader or  narrower bandwidth
        than 1 Hz (= 0.001 kHz, in the units used here.)"""
        return (4 * k_B * self.Gamma_i * T * bandwidth) ** 0.5

    def __str__(self):
        """Write out the cantilever as its
        three most important parameters: resonance frequency,
        spring constant, and quality factor."""
        return """f_c = {self.f_c:.2g} kHz,\
                  k_c = {k_c:.2g} N/m,\
                  Q = {self.Q:.2g}""".format(self=self, k_c=self.k_c / 1e9)

    def __repr__(self):
        """Return a representation of the cantilever which can be used
        to generate an equivalent copy of the cantilever."""
        return "Cantilever(f_c = {self.f_c}, k_c = {self.k_c},\
Q = {self.Q},  R_tip = {self.R_tip},\
L_tip = {self.L_tip},\
ThetaDegrees_tip = {self.ThetaDegrees_tip},\
geometry_c = '{self.geometry_c}')".format(self=self)



class Sample(Assigner):
    """A transistor sample.

    This class defines a transistor sample according to the
    parameters listed below. The temperature is considered a
    sample property so that the diffusion constant of the
    sample can be calculated according to the Einstein
    relation. Defaults are given in parentheses and all values
    are input in the units given below (see the module help
    for more information).

    semiconductor
        The semiconductor material of the transistor ('TPD').

    E_s2
        The static relative dielectric constant of the semiconductor (3.5).

    E_s2
        The imaginary part of the dielectic constant of the semiconductor
        E_s1, adjusted for the contribution of conductivity
        (:math:`\frac{\sigma}{\epsilon_0 D}`) (-0.0005).

    E_i1
        The static relative dielectric constant of the insulator (3.5).

    E_i2
        The imaginary part of the dielectric constant of the insulator (0).

    h
        The thickness of the semiconductor (0.07 um = 70 nm).

    h_trans
        The thickness of the semiconductor in which charge motion is
        confined (3e-3 um = 3 nm).

    h_i
        The thickness of the insulator (0.315 um = 315 nm).

    mobility
        The mobility of the sample (3e-4 um^2/mV/ms = 3e-10 m^2/V/s =
        3e-6 cm^2/V/s).
    T
        The temperature of the sample (298 K).

    V_g
        The gate voltage of the transistor (10e3 mV = 10 V).

    """
    # Exclude the two variables handled as properties from the autoassign
    @autoassign(exclude=('V_g', 'rho'))
    def __init__(self, semiconductor='TPD',
                 h=70e-3, h_trans=1e-3, h_i=300e-3,
                 E_s1=3.5, E_s2=-0.0005,
                 E_i1=4.65, E_i2=0,
                 mobility=3e-4, T=293,
                 V_g=10e3, rho=None):
        """Initialize the sample with all of the experimentally
        relevant sample parameters.

        @autoassign automatically assigns the input quantities
        to self. See http://code.activestate.com/recipes/551763/
        for more information."""

        if rho is None:
            self.V_g = V_g
        elif V_g == 10e3:
            self.rho = rho
        else:
            raise ValueError("""\
                The provided values of 'V_g' and 'rho'are incompatible.
                Only specify one of 'V_g' or 'rho' when defining a Sample.""")

        if V_g <= 0:
            raise ValueError("The voltage 'V_g' must be positive.")

    # TODO: Add units to all of these properties!
    # All of the properties here are functions of the initial parameters.
    @property
    def diff(self):
        """Diffusion constant defined according to the Einstein relation."""
        return self.mobility * k_B * self.T / q

    @property
    def C_i(self):
        """Capacitance per unit area between the
        transistor gate and sample."""
        return self.E_i1 * E_0 / self.h_i

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
        self._rho = self.C_i * self._V_g / (self.h_trans * q)

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
        self._V_g = q * self._rho * self.h_trans / self.C_i

    #---------------------------------------------------------
    # Relevant properties derived from gate voltage / charge density.
    @property
    def sigma(self):
        """The conductivity sigma of the sample."""
        return self.mobility * self.rho * q

    @property
    def kappa(self):
        """Define the inverse Debye length screening length kappa,
        used in the Lekkala Loring theory. See Lekkala et al.,
        p4, at http://dx.doi.org/10.1063/1.4754602."""
        return (2 * self.rho * q ** 2 / (E_0 * k_B * self.T)) ** 0.5

    def E_eff(self, omega):
        """Defines the effective dielectric constant,
        which corrects for the effect of Ohm's law (carrier drift),
        at a particular angular frequency. See Eq. 14 in
        Lekkala et al., 2013."""
        return (self.E_s - self.sigma / (E_0 * omega) * 1j)

    def __str__(self):
        """Write out the sample, using its semiconductor material,
        height and mobility."""
        return """{self.semiconductor}  {h_nm:.2g} nm
                   mobility {mobility_cm2:.2g} cm^2/Vs""".format(
            self=self, h_nm=self.h * 1000,
            mobility_cm2=self.mobility * 1e-2)

    def __repr__(self):
        """Machine readable representation of the object,
        that can return an identical copy of the object with"""
        return """Sample(semiconductor = '{self.semiconductor}',\
            h = {self.h}, h_trans = {self.h_trans},
            h_i = {self.h_i}, E_s1 = {self.E_s1},\
            E_s2 = {self.E_s2}, E_i1 = {self.E_i1},\
            E_i2 = {self.E_i2}, mobility = {self.mobility},\
            T = {self.T}, V_g = {self.V_g})""".format(self=self)


class Experiment(Assigner):
    """Stores parameters set by the experimenter.

    This sets the tip-sample voltage, transistor gate voltage
    and the distance between the cantilever and the sample.
    Defaults are given in parentheses and all values are input
    in the units given below (see the module help for more
    information).

    d
        Distance between the cantilever and the sample . (0.1 um = 100 nm)

    V_ts
        The tip-sample voltage (5e3 mV = 5 V)

    jitter_f_i
        The lower frequency endpoint used to calcualte
        jitter (0.2e-3 kHz = 0.2 Hz)

    jitter_f_f
        The upper frequency endpoint used to calcualte
        jitter (3e-3 kHz = 3 Hz)
    """
    @autoassign
    def __init__(self, d=0.1, V_ts=5e3,
                 jitter_f_i=0.2e-3, jitter_f_f=3e-3):
        """Initialize the experimental parameters, and then check for
        obvious errors, raising a ValueError."""

        # Check for errors in the experimental parameters
        if self.V_ts < 0:
            raise ValueError("The voltages 'V_g' and 'V_ts' must be positive.")
        if self.jitter_f_i > self.jitter_f_f:
            raise ValueError("'jitter_f_i' must be less than 'jitter_f_f'.")

    def __str__(self):
        """A nice string representation of the experiment."""
        return """Tip-sample: {d_nm:.2g} nm, {V_ts_v:.2g}""".format(
            d_nm=self.d * 1000, V_ts_v=self.V_ts / 1e3)

    def __repr__(self):
        """Machine readable representation of the object."""
        return """Experiment(d = {self.d}, V_ts = {self.V_ts},\
            jitter_f_i = {self.jitter_f_i},\
            jitter_f_f = {self.jitter_f_f})""".format(self=self)
