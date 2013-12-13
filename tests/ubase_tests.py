"""
Test UnitCantilever Class
2013-12-12

Ryan Dwyer

"""

from jittermodel.ubase import u, UnitCantilever, UnitExperiment
from nose.tools import assert_raises, assert_almost_equals
from pint import DimensionalityError
import unittest

# TO DO
# Pint Helper function!


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
    UnitCantilever(f_c=50 * u.kHz, k_c=3 * u.N/u.m,
                   Q=1000 * u.dimensionless)


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


def test_UnitExperiment():
    UnitExperiment(jitter_f_f=4*u.Hz)
    UnitExperiment(d=0.3 * u.um, V_ts=10 * u.V)