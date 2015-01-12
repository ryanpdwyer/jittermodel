# -*- coding: utf-8 -*-
from __future__ import division, print_function
import numpy as np

from jittermodel import u, q2unitless
from jittermodel.simulation import (Simulation, SphereCapacitance, _alpha, 
                                     sum_sinh, _eta, _lambda, _thetaI,
                                     _thetaII)
from jittermodel._sim import _thetaI_c
from jittermodel.base import Cantilever, Experiment, Transistor
from numpy.testing import assert_allclose
from nose.tools import eq_, assert_almost_equal, assert_raises
from bunch import Bunch
from jittermodel.tests import expected_failure
import unittest

u.d = u.dimensionless  # For brevity

import mpmath as mp


def mp_sum_sinh(alpha):
    """Implements the infinite sum using mpmath, at very high precision.
    Method 'r+s+e' was found to work accurately for all values of alpha,
    unlike most other alogithms in Mathematica, python, etc."""
    summand = lambda n: mp.sinh(alpha) / mp.sinh(alpha * n)
    return mp.nsum(summand, [1, mp.inf], method='r+s+e')


class Test_sum_sinh(unittest.TestCase):

    @staticmethod
    def test_sum_sinh():
        """Test that the sum is working properly for a range of alpha values.
        The mpmath module is used to verify that the sum meets error
        specifications.
        """
        alphas = [2 ** i for i in xrange(-12, 7)]
        results = [sum_sinh(alpha) for alpha in alphas]
        mp_results = [mp_sum_sinh(alpha) for alpha in alphas]
        for mp_result, test_result in zip(mp_results, results):
            assert_almost_equal(mp_result, test_result, 7)


class MockSimulationCapacitance(object):
    """A mock simulation object only containing the parameters necessary to
    test SphereCapacitance"""
    units =  {"[mass]": u.pg, "[length]": u.um, "[time]": u.ms,
             "[current]": u.aC / u.ms, "[temperature]": u.K, "[angle]": u.rad}
    E_0 = q2unitless(u.epsilon_0, units)
    q = q2unitless(u.elementary_charge, units)
    k_B = q2unitless(u.boltzmann_constant, units)

    Samp = Bunch(h=0.1, E_s1=3)
    Cant = Bunch(R_tip=0.05)
    Expt = Bunch(d=0.15)

    def __init__(self):
        self.sphere = SphereCapacitance(self)


# TODO: Where do these test cases come from?
class TestSphereCapacitance(unittest.TestCase):
    def setUp(self):
        self.sim = MockSimulationCapacitance()

    def test_C(self):
        assert_almost_equal(0.00623177, self.sim.sphere.C())

    def test_Cd(self):
        assert_almost_equal(-0.00322151, self.sim.sphere.Cd())

    def test_Cd2(self):
        assert_almost_equal(0.0311542, self.sim.sphere.Cd2())


class TestSimulation(unittest.TestCase):

    @staticmethod
    def test_init_Simulation():
        cant = Cantilever(f_c=50*u.kHz, k_c=3.5*u.N/u.m, Q=20000*u.d,
                              R_tip=40*u.nm, L_tip=15*u.um, theta_tip=16*u.degrees,
                              geometry_c='perpendicular')
        trans = Transistor(semiconductor='TPD', h=70 * u.nm, h_trans=1 * u.nm,
                               h_i=300 * u.nm, E_s1=3.5, E_s2=-0.0005, E_i1=4.65,
                               E_i2=0, mobility=3e-6 * u.cm ** 2 / u.V / u.s,
                               T=298 * u.K, V_g=10 * u.V, rho=None)
        expt = Experiment(d=100 * u.nm, V_ts=5 * u.V, jitter_f_i=0.2 * u.Hz,
                              jitter_f_f=3 * u.Hz)

        sim = Simulation(cant, trans, expt)
        # Test some properties are correct
        eq_(sim.Cant.f_c, 50)
        eq_(sim.Expt.d, 0.1)
        eq_(sim.Samp.h_i, 0.3)
        assert_almost_equal(sim.Samp.diff, 0.0077038955272097955)

# These tests are all generated by implementing sympy code for the functions in
# validate-im-dielectric.ipynb. That should be a good comparison; sympy
# uses mpmath as a backend for its infinite precision arithmatic, so this
# should be robust against ordinary floating point errors.

class TestImDielectricHelperFunctions(unittest.TestCase):

    @staticmethod
    def test__eta():
        k = np.array([1, 10, 100, 1000, 10000, 100000])
        kappa = 3500
        E_s = 3 - 0.001j
        D = 0.005
        omega = 300
        # Expected values were calculated using sympy,
        # to 15 digits of precision.
        # See test_verification/validate-im-dielectric
        expected_eta = np.array([2020.78311260126 + 15.182507854811j,
                                 2020.80760652432 + 15.182323829782j,
                                 2023.25550170076 + 15.163955048756j,
                                 2254.66583909462 + 13.607584302718j,
                                 10202.1243828582 + 3.007271263178j,
                                 100020.414581093 + 0.30674293451j])
        eta = _eta(k, kappa, E_s, D, omega)
        assert_allclose(eta, expected_eta)

    @staticmethod
    def test__lambda():
        k = np.array([1, 10, 100, 1000, 10000, 100000])
        eta = np.array([2020.78311260126 + 15.182507854811j,
                        2020.80760652432 + 15.182323829782j,
                        2023.25550170076 + 15.163955048756j,
                        2254.66583909462 + 13.607584302718j,
                        10202.1243828582 + 3.007271263178j,
                        100020.414581093 + 0.30674293451j])

        E_eff = 3 - 100j
        E_s = 3 - 0.001j
        expected_lambda = np.array([0.0001184255261724 + 0.0164941987549306j,
                                    0.00118421087011718 + 0.164939988752172j,
                                    0.0117978533636026 + 1.64740475175451j,
                                    0.0842948437214929 + 14.7834985873234j,
                                   -0.00125999301746353 + 32.672603689536j,
                                   -0.0110065260871034 + 33.3261929274151j])
        Lambda = _lambda(k, eta, E_eff, E_s)
        assert_allclose(Lambda, expected_lambda)

    @staticmethod
    def test_thetaI():
        k = np.array([1, 10, 100, 1000, 10000, 100000])
        eta = np.array([2020.78311260126 + 15.182507854811j,
                        2020.80760652432 + 15.182323829782j,
                        2023.25550170076 + 15.163955048756j,
                        2254.66583909462 + 13.607584302718j,
                        10202.1243828582 + 3.007271263178j,
                        100020.414581093 + 0.30674293451j])
        Lambda = np.array([0.0001184255261724 + 0.0164941987549306j,
                           0.00118421087011718 + 0.164939988752172j,
                           0.0117978533636026 + 1.64740475175451j,
                           0.0842948437214929 + 14.7834985873234j,
                          -0.00125999301746353 + 32.672603689536j,
                          -0.0110065260871034 + 33.3261929274151j])

        expected_thetaI = np.array([0.00157126996626562 + 0.0210682675809495j,
                                     0.00672782406000677 + 0.0281575198774334j,
                                     0.050275664263775 + 0.0281213204722464j,
                                     0.443934273416263 + 0.0140052914999941j,
                                     0.980197277948465 + 0.000305155415174606j,
                                     0.999795989512753 + 3.05416795636227e-6j])

        h_s = 0.1
        alpha = 0.65 - 0.0002j
        E_s = 3 - 0.001j
        E_eff = 3 - 100j
        thetaI = [_thetaI(_k, h_s, alpha, _Lambda, _eta, E_s, E_eff) for
                  _k, _Lambda, _eta in zip(k, Lambda, eta)]
        thetaI = np.array(thetaI)
        assert_allclose(expected_thetaI, thetaI)

    @staticmethod
    def test_thetaI_c():
        k = np.array([1, 10, 100, 1000, 10000, 100000])
        eta = np.array([2020.78311260126 + 15.182507854811j,
                        2020.80760652432 + 15.182323829782j,
                        2023.25550170076 + 15.163955048756j,
                        2254.66583909462 + 13.607584302718j,
                        10202.1243828582 + 3.007271263178j,
                        100020.414581093 + 0.30674293451j])
        Lambda = np.array([0.0001184255261724 + 0.0164941987549306j,
                           0.00118421087011718 + 0.164939988752172j,
                           0.0117978533636026 + 1.64740475175451j,
                           0.0842948437214929 + 14.7834985873234j,
                          -0.00125999301746353 + 32.672603689536j,
                          -0.0110065260871034 + 33.3261929274151j])

        expected_thetaI = np.array([0.00157126996626562 + 0.0210682675809495j,
                                     0.00672782406000677 + 0.0281575198774334j,
                                     0.050275664263775 + 0.0281213204722464j,
                                     0.443934273416263 + 0.0140052914999941j,
                                     0.980197277948465 + 0.000305155415174606j,
                                     0.999795989512753 + 3.05416795636227e-6j])

        h_s = 0.1
        alpha = 0.65 - 0.0002j
        E_s = 3 - 0.001j
        E_eff = 3 - 100j
        thetaI = [_thetaI_c(_k, h_s, alpha, _Lambda, _eta, E_s, E_eff) for
                  _k, _Lambda, _eta in zip(k, Lambda, eta)]
        thetaI = np.array(thetaI)
        assert_allclose(expected_thetaI, thetaI)

    @staticmethod
    def test_thetaII():
        k = np.array([1, 10, 100, 1000, 10000, 100000])
        Lambda = np.array([0.0001184255261724 + 0.0164941987549306j,
                           0.00118421087011718 + 0.164939988752172j,
                           0.0117978533636026 + 1.64740475175451j,
                           0.0842948437214929 + 14.7834985873234j,
                          -0.00125999301746353 + 32.672603689536j,
                          -0.0110065260871034 + 33.3261929274151j])

        expected_thetaII = np.array([0.101145810077246 + 0.0296480635666554j,
                                     0.764320753451023 + 0.0123928030520502j,
                                     0.999999996277978 + 2.1003332939236e-10j,
                                     1.0 + 8.470329472543e-22j,
                                     1.00000000000000,
                                     1.00000000000000])

        h = 0.1
        E_s = 3 - 0.001j
        E_d = 3 - 0.001j
        E_eff = 3 - 100j

        thetaII = _thetaII(k, h, E_s, E_d, E_eff, Lambda)
        assert_allclose(expected_thetaII, thetaII)


class TestSimulationObject(unittest.TestCase):
    def setUp(self):
        self.cant = Cantilever(
            f_c=50*u.kHz, k_c=3.5*u.N/u.m, Q=20000*u.d,
            R_tip=40*u.nm, L_tip=15*u.um, theta_tip=16*u.degrees,
            geometry_c='perpendicular')
        self.trans = Transistor(
            semiconductor='TPD', h=70 * u.nm, h_trans=1 * u.nm, h_i=300 * u.nm,
            E_s1=3.5, E_s2=-0.0005, E_i1=4.65, E_i2=0,
            mobility=3e-6*u.cm**2/u.V/u.s, T=298 * u.K, V_g=10 * u.V, rho=None)
        self.expt = Experiment(
            d=100 * u.nm, V_ts=5 * u.V, jitter_f_i=0.2 * u.Hz,
            jitter_f_f=3 * u.Hz)

        self.sim = Simulation(self.cant, self.trans, self.expt)

    def test_check_assignment(self):
        eq_(self.sim.Cant.f_c, 50)
        self.sim.assign('f_c', 75 * u.kHz)
        eq_(self.sim.Cant.f_c, 75)


# TODO: Where does this calculation come from?
class TestCapacitanceCalculations(unittest.TestCase):
    def setUp(self):
        cant = Cantilever(
            f_c=46*u.kHz, Q=2500*u.dimensionless, k_c=0.85*u.N/u.m,
            R_tip=40*u.nm, L_tip=15*u.um, theta_tip=16*u.degrees,
            geometry_c='perpendicular')
        trans = Transistor(
            semiconductor='TPD', h=63*u.nm, h_trans=1*u.nm, h_i=300*u.nm,
            E_s1=3.5, E_s2=-0.0005, E_i1=4.65, E_i2=0,
            mobility=3e-6*u.cm**2/u.V/u.s, T=298*u.K, V_g=10*u.V, rho=None)
        expt = Experiment(
            d=100*u.nm, V_ts=5*u.V, jitter_f_i=0.2*u.Hz, jitter_f_f=3*u.Hz)

        self.sim = Simulation(cant, trans, expt)

    def test_capacitance(self):
        assert_almost_equal(self.sim.Sphere.C(), 5.1e-3, places=2)

    def test_capacitance_2nd_derivative(self):
        assert_almost_equal(self.sim.Sphere.Cd2(), 6.9e-2, places=2)


class TestImDielectric(unittest.TestCase):
    def setUp(self):
        exp = Experiment(V_ts=3*u.V, d=300*u.nm)
        cant = Cantilever(
            f_c=65*u.kHz, Q=2500*u.dimensionless, k_c=3.5*u.N/u.m,
            R_tip=40*u.nm, L_tip=15*u.um, theta_tip=16*u.degrees,
            geometry_c='perpendicular')
        trans = Transistor(
            h=72*u.nm, V_g=40*u.V, E_s1=3.4, mobility=2.7e-6*u.cm**2/u.V/u.s,
            E_s2=-0.05, E_i1=4.65)
        self.sim = Simulation(cant, trans, exp)

    @expected_failure
    def test_im_dielectric(self):
        assert_almost_equal(self.sim._im_dielectric(1), -0.003635)
