# -*- coding: utf-8 -*-
"""
Test UnitBase Classes
=========================
2013-12-12

Ryan Dwyer

"""

from jittermodel.base import (u, SimpleCantilever, Cantilever,
                              Experiment, Transistor)
from jittermodel.tests import pint_assert_almost_equal
from nose.tools import assert_raises
from pint import DimensionalityError
import unittest


# ----- SUCantilever Tests ----------------------------------------------

class TestSimpleCantilver(unittest.TestCase):

    def test_SimpleCantilever_bad_init(self):
        """Test that SUnit Cantilever returns an error when
        too few arguements are given."""
        to_test = [{'f_c': -10 * u.Hz, 'Q': 39 * u.dimensionless},
                   {'Q': 39 * u.dimensionless, 'k_c': 1 * u.N / u.m},
                   {'f_c': 40 * u.kHz, 'k_c': 1 * u.N / u.m}]
        for kwargs in to_test:
            assert_raises(TypeError, SimpleCantilever, **kwargs)

# ---- Unit Cantilever Tests --------------------------------------------

class TestCantilever(unittest.TestCase):

    def setUp(self):
        self.c = Cantilever(f_c=50*u.kHz, k_c=3*u.N/u.m,
                            Q=20000*u.dimensionless)

    def test_Cantilever_input(self):
        """Make sure that defining a Cantilever with an incorrect geometry, or
        negative number raises a ValueError."""
        to_test = [{'f_c': -10 * u.Hz},
               {'Q': -39 * u.dimensionless},
               {'f_c': 40 * u.kHz, 'R_tip': -0.023 * u.nm},
               {'theta_tip': -10 * u.degrees},
               {'theta_tip': 100 * u.degrees},
               {'geometry_c': 'not perpendicular or parallel'}]

        for kwargs in to_test:
            assert_raises(ValueError, Cantilever, **kwargs)

    def test_Cantilever_input_units(self):
        to_test = [{'f_c': -10 * u.s},
               {'Q': -39},
               {'f_c': 40 * u.kHz, 'R_tip': -0.023 * u.C},
               {'theta_tip': 12 * u.degK},
               {'theta_tip': 100}]

        for kwargs in to_test:
            assert_raises(DimensionalityError, Cantilever, **kwargs)

    def test_Cantilever_init(self):
        """Make sure the unit cantilever initializes properly."""
        c1 = Cantilever(f_c=50 * u.kHz, k_c=3 * u.N/u.m,
                        Q=1000 * u.dimensionless)
        assert c1.f_c == 50 * u.kHz

    def test_F_min(self):
        ex_F_min = 2.8125685411157023e-3 * u.pN
        pint_assert_almost_equal(ex_F_min, self.c.F_min(300*u.K))

    def test_Gamma_i(self):
        c = self.c
        ex_Gamma_i = 477.46482927568604 * u.pN * u.s / u.m
        pint_assert_almost_equal(ex_Gamma_i, c.Gamma_i)

# ---- Unit Experiment Tests --------------------------------------------


class TestExperiment(unittest.TestCase):

    def test_Experiment_init(self):
        e1 = Experiment(jitter_f_f=4*u.Hz)
        e2 = Experiment(d=0.3 * u.um, V_ts=10 * u.V)
        assert e1.jitter_f_f == 4*u.Hz
        assert e2.d == 0.3 * u.um
        assert_raises(DimensionalityError, Experiment, **{'d': 4 * u.K})
        assert_raises(ValueError, Experiment, **{'d': -1*u.nm})


# ----- Unit Transistor Tests --------------------------------------------

class TestTransistor(unittest.TestCase):

    def test_Transistor_init_success(self):
        samp1 = Transistor(semiconductor='TPD',
                           h=70 * u.nm, h_trans=1 * u.nm, h_i=300 * u.nm,
                           E_s1=3.5, E_s2=-0.0005, E_i1=4.65, E_i2=0,
                           mobility=3e-6 * u.cm ** 2 / u.V / u.s, T=298 * u.K,
                           V_g=10 * u.V, rho=None)

        assert samp1.h == 70 * u.nm
        assert samp1.mobility == 3e-6 * u.cm ** 2 / u.V / u.s

    def test_Transistor_bad_init(self):
        """Try some things that shoud raise errors."""
        assert_raises(ValueError, Transistor, **{'T': -23 * u.K})
        assert_raises(DimensionalityError, Transistor, **{'h': 70 * u.s})
