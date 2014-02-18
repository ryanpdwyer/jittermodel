"""
Test UnitCantilever Class
2013-12-12

Ryan Dwyer

"""

from jittermodel.ubase import (u, SUCantilever, UnitCantilever,
                               UnitExperiment, UnitTransistor)
from nose.tools import assert_raises, assert_almost_equals
from pint import DimensionalityError
import unittest

# TO DO
# Pint Helper function!


# ----- SUCantilever Tests ----------------------------------------------

def test_SUnitCantilever_bad_init():
    """Test that SUnit Cantilever returns an error when
    too few arguements are given."""
    to_test = [{'f_c': -10 * u.Hz, 'Q': 39 * u.dimensionless},
               {'Q': 39 * u.dimensionless, 'k_c': 1 * u.N / u.m},
               {'f_c': 40 * u.kHz, 'k_c': 1 * u.N / u.m}]
    for kwargs in to_test:
        assert_raises(TypeError, SUCantilever, **kwargs)

# ---- Unit Cantilever Tests -------------------------------------------- 


def test_UnitCantilever_input():
    """Make sure that defining a UnitCantilever with an incorrect geometry, or
    negative number raises a ValueError."""
    to_test = [{'f_c': -10 * u.Hz},
               {'Q': -39 * u.dimensionless},
               {'f_c': 40 * u.kHz, 'R_tip': -0.023 * u.nm},
               {'theta_tip': -10 * u.degrees},
               {'theta_tip': 100 * u.degrees},
               {'geometry_c': 'not perpendicular or parallel'}]

    for kwargs in to_test:
        assert_raises(ValueError, UnitCantilever, **kwargs)


def test_UnitCantilever_input_units():
    to_test = [{'f_c': -10 * u.s},
               {'Q': -39},
               {'f_c': 40 * u.kHz, 'R_tip': -0.023 * u.C},
               {'theta_tip': 12 * u.degK},
               {'theta_tip': 100}]

    for kwargs in to_test:
        assert_raises(DimensionalityError, UnitCantilever, **kwargs)


def test_UnitCantilever_init():
    """Make sure the unit cantilever initializes properly."""
    c1 = UnitCantilever(f_c=50 * u.kHz, k_c=3 * u.N/u.m,
                        Q=1000 * u.dimensionless)
    assert c1.f_c == 50 * u.kHz


class TestUnitCantilever(unittest.TestCase):

    def setUp(self):
        self.c = UnitCantilever(f_c=50*u.kHz, k_c=3*u.N/u.m,
                                Q=20000*u.dimensionless)

    def test_F_min(self):
        ex_F_min = 2.8125685411157023e-3 * u.pN
        assert_almost_equals(ex_F_min.magnitude,
                             self.c.F_min(300*u.K).magnitude)

    def test_Gamma_i(self):
        c = self.c
        ex_Gamma_i = 477.46482927568604 * u.pN * u.s / u.m
        assert_almost_equals(ex_Gamma_i.magnitude,
                             c.Gamma_i.magnitude)

# ---- Unit Experiment Tests --------------------------------------------

def test_UnitExperiment_init():
    e1 = UnitExperiment(jitter_f_f=4*u.Hz)
    e2 = UnitExperiment(d=0.3 * u.um, V_ts=10 * u.V)
    assert e1.jitter_f_f == 4*u.Hz
    assert e2.d == 0.3 * u.um
    assert_raises(DimensionalityError, UnitExperiment, **{'d': 4 * u.K})
    assert_raises(ValueError, UnitExperiment, **{'d':-1*u.nm})


# ----- Unit Transistor Tests --------------------------------------------

def test_UnitTransistor_init():
    samp1 = UnitTransistor(semiconductor='TPD',
                 h=70 * u.nm, h_trans=1 * u.nm, h_i=300 * u.nm,
                 E_s1=3.5, E_s2=-0.0005,
                 E_i1=4.65, E_i2=0,
                 mobility=3e-6 * u.cm ** 2 / u.V / u.s, T=298 * u.K,
                 V_g=10 * u.V, rho=None)
    
    assert samp1.h == 70 * u.nm
    assert samp1.mobility == 3e-6 * u.cm ** 2 / u.V / u.s
    # Try some things that shoud raise errors.
    assert_raises(ValueError, UnitTransistor, **{'T': -23 * u.K})
    assert_raises(DimensionalityError, UnitTransistor, **{'h': 70 * u.s})

