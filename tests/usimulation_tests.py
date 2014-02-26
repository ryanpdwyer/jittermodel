from jittermodel import u
from jittermodel.usimulation import UnitSimulation
from jittermodel.ubase import UnitCantilever, UnitExperiment, UnitTransistor
from nose.tools import eq_, assert_almost_equal, assert_raises
import unittest

u.d = u.dimensionless  # For brevity


def test_init_UnitSimulation():
    cant = UnitCantilever(f_c=50*u.kHz, k_c=3.5*u.N/u.m, Q=20000*u.d,
                          R_tip=40*u.nm, L_tip=15*u.um, theta_tip=16*u.degrees,
                          geometry_c='perpendicular')
    trans = UnitTransistor(semiconductor='TPD', h=70 * u.nm, h_trans=1 * u.nm,
                           h_i=300 * u.nm, E_s1=3.5, E_s2=-0.0005, E_i1=4.65,
                           E_i2=0, mobility=3e-6 * u.cm ** 2 / u.V / u.s,
                           T=298 * u.K, V_g=10 * u.V, rho=None)
    expt = UnitExperiment(d=100 * u.nm, V_ts=5 * u.V, jitter_f_i=0.2 * u.Hz,
                          jitter_f_f=3 * u.Hz)

    sim = UnitSimulation(cant, trans, expt)
    # Test some properties are correct
    eq_(sim.Cant.f_c, 50)
    eq_(sim.Expt.d, 0.1)
    eq_(sim.Samp.h_i, 0.3)
    assert_almost_equal(sim.Samp.diff, 0.00770475093633)


class TestUnitSimulation(unittest.TestCase):
    def setUp(self):
        self.cant = UnitCantilever(
            f_c=50*u.kHz, k_c=3.5*u.N/u.m, Q=20000*u.d,
            R_tip=40*u.nm, L_tip=15*u.um, theta_tip=16*u.degrees,
            geometry_c='perpendicular')
        self.trans = UnitTransistor(
            semiconductor='TPD', h=70 * u.nm, h_trans=1 * u.nm, h_i=300 * u.nm,
            E_s1=3.5, E_s2=-0.0005, E_i1=4.65, E_i2=0,
            mobility=3e-6*u.cm**2/u.V/u.s, T=298 * u.K, V_g=10 * u.V, rho=None)
        self.expt = UnitExperiment(
            d=100 * u.nm, V_ts=5 * u.V, jitter_f_i=0.2 * u.Hz,
            jitter_f_f=3 * u.Hz)

        self.sim = UnitSimulation(self.cant, self.trans, self.expt)

    def test_check_assignment(self):
        eq_(self.sim.Cant.f_c, 50)
        self.sim.assign('f_c', 75 * u.kHz)
        eq_(self.sim.Cant.f_c, 75)


class TestCapacitanceCalculations(unittest.TestCase):
    def setUp(self):
        cant = UnitCantilever(
            f_c=46*u.kHz, Q=2500*u.dimensionless, k_c=0.85*u.N/u.m,
            R_tip=40*u.nm, L_tip=15*u.um, theta_tip=16*u.degrees,
            geometry_c='perpendicular')
        trans = UnitTransistor(
            semiconductor='TPD', h=63*u.nm, h_trans=1*u.nm, h_i=300*u.nm,
            E_s1=3.5, E_s2=-0.0005, E_i1=4.65, E_i2=0,
            mobility=3e-6*u.cm**2/u.V/u.s, T=298*u.K, V_g=10*u.V, rho=None)
        expt = UnitExperiment(
            d=100*u.nm, V_ts=5*u.V, jitter_f_i=0.2*u.Hz, jitter_f_f=3*u.Hz)

        self.sim = UnitSimulation(cant, trans, expt)

    def test_capacitance(self):
        assert_almost_equal(self.sim.C_sphere(), 5.1e-3, places=1)
        assert_almost_equal(self.sim.Cd2_sphere(), 6.9e-2, places=1) 


class TestImDielectric(unittest.TestCase):
    def setUp(self):
        pass
