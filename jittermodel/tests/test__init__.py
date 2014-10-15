# -*- coding: utf-8 -*-
"""
__init__ tests
2013-12-16
Ryan Dwyer

Test general helper functions in the package and the base classes used to
derive specific classes we implement our simulation with.
"""
from __future__ import division
import pint
from numpy import pi
from jittermodel import (get_defaults, get_default_units, u, Assigner,
                         UnitAssigner, NoUnitAssigner, q2unitless, make_units,
                         silentremove)
from jittermodel.tests import pint_assert_almost_equal, expected_failure
import unittest
from nose.tools import (eq_, assert_not_equal,
                        assert_raises, assert_almost_equal)
import cPickle as pickle


pa_eq = pint_assert_almost_equal


class Test_pint_assert_almost_equal(unittest.TestCase):

    def test_pint_assert_almost_equal(self):
        first = 4.77464829276e-07 * u.N / u.kHz / u.m
        second = 477.464829276 * u.pN * u.s / u.m
        pa_eq(first, second)

    def test_pint_assert_almost_equal_not(self):
        first = 4.78e-07 * u.N / u.kHz / u.m
        second = 477.464829276 * u.pN * u.s / u.m
        assert_raises(AssertionError, pa_eq, **{'first': first,
                      'second': second, 'places': 7, 'unit': u.pN * u.s / u.m})

class Test_get_default_funcs(unittest.TestCase):

    @staticmethod
    def test_get_defaults():
        def no_defaults(x):
            pass

        def both_defaults(x=3, y=5):
            pass

        def no_default_then_default(x, y='adfe'):
            pass

        funcs = [no_defaults, both_defaults, no_default_then_default]
        exp_default_dicts = [{'x': None}, {'x': 3, 'y': 5},
                             {'x': None, 'y': 'adfe'}]

        for func, exp_dict in zip(funcs, exp_default_dicts):
            eq_(get_defaults(func), exp_dict)

    @staticmethod
    def test_get_default_units():
        def f1(x=2*u.nm):
            pass

        def f2(x='TPD', y=10*u.V):
            pass

        def f3(x=10*u.m, y=3*u.s, z=298*u.K):
            pass

        funcs = [f1, f2, f3]
        exp_unit_dicts = [{'x': u.nm.units}, {'y': u.V.units},
                          {'x': u.m.units, 'y': u.s.units, 'z': u.K.units}]

        for func, exp_dict in zip(funcs, exp_unit_dicts):
            eq_(get_default_units(func), exp_dict)


class Test_q2unitless(unittest.TestCase):
    def setUp(self):
        self.units = {"[mass]": u.pg, "[length]": u.um, "[time]": u.ms,
                      "[current]": u.aC / u.ms, "[temperature]": u.K,
                      "[angle]": u.rad}

    def test_make_units(self):
        dim_dict = {'[length]': 1, '[time]': 1, '[temperature]': -1}
        unit = u.um * u.ms / u.K
        pa_eq(unit, make_units(dim_dict, self.units))

    def test_normal(self):
        speed = 10 * u.m / u.s
        friction = 100 * u.pN * u.s / u.m
        eq_(10000, q2unitless(speed, self.units))
        assert_almost_equal(q2unitless(friction, self.units), 100)

    def test_dimensionless(self):
        Q = 1000 * u.dimensionless
        eq_(q2unitless(Q, self.units), 1000)

    def test_radians(self):
        """If this does not work, change your unit definition file to
        to include

        radian = [angle] = rad

        so that angles put something in the unitcontainer."""
        angle = 15 * u.degrees
        assert_almost_equal(15 * pi / 180, q2unitless(angle, self.units))

    def test_bits(self):
        """This asserts that calculations using bits fail (in an expected way).
        Note: Bits / bytes now works correctly in pint 0.5dev, so the test has
        been updated accordingly.
        """
        bytes = 4 * u.bytes
        eq_(q2unitless(bytes, self.units), 32)


class TestAssigner(unittest.TestCase):
    """A test to verify that Assigner functions properly."""
    def setUp(self):
        self.a = Assigner()

    def test_assign(self):
        self.a.assign('att', 10)
        eq_(self.a.att, 10)

    def test_update_assign(self):
        """Test that we can update an existing value using assign,
        and return the updated value with 'dot' syntax or the
        lookup method."""
        self.a.a = 2
        self.a.assign('a', 10)
        eq_(self.a.a, 10)
        eq_(self.a.lookup('a'), 10)

    def test_lookup(self):
        self.a.att = 10
        eq_(self.a.lookup('att'), 10)

    def test_all_attributes(self):
        a = self.a
        a.a = 2
        a.b = 10
        a.c = "Hello"
        exp_all_attributes = {'a', 'b', 'c'}
        eq_(exp_all_attributes, a._all_attributes)

    def test_hidden_all_attributes(self):
        """Hide attributes prefaced with '_', so that property getter / setters
        can use these values without interference."""
        a = self.a
        a.a = "hello"
        a._a = "I am hidden"
        eq_({'a'}, a._all_attributes)

    def test_all_attributes_with_method(self):
        class AssignerPlusMethod(Assigner):
            def a_method(self):
                pass

            @property
            def a_property(self):
                pass

            @property
            def property_with_setter(self):
                return self._property_with_setter

            @property_with_setter.setter
            def property_with_setter(self, value):
                self._property_with_setter = value

        a = AssignerPlusMethod()
        a.a = 2
        a.b = 10
        a.property_with_setter = 1200
        eq_({'a', 'b', 'property_with_setter'}, a._all_attributes)

    @staticmethod
    def test_pickle_Assigner():
        filename = 'test-pickle-Assigner.pkl'
        a = Assigner()
        a.b = 2
        with open(filename, 'w') as f:
            pickle.dump(a, f)

        with open(filename, 'r') as f:
            b = pickle.load(f)

        eq_(a.b, b.b)
        silentremove(filename)


class TestEqAssigner(unittest.TestCase):
    def setUp(self):
        b = Assigner()
        a = Assigner()

        a.a = 100
        a.b = "Hello"
        a.assign("att", 233)

        b.a = 100
        b.assign('b', "Hello")
        b.att = 233

        self.a = a
        self.b = b

    def test__eq__(self):
        eq_(self.a, self.b)

    def test__not__eq__(self):
        self.a.d = 4
        self.b.d = 5
        assert_not_equal(self.a, self.b)


class TUAssigner(UnitAssigner):
    def __init__(self, x=2*u.m, t=1*u.s):
        self.UnitlessClass = UnitlessTUAssigner
        self.x = x
        self.t = t
        self._default_units = {'x': u.m, 't': u.s}


class UnitlessTUAssigner(NoUnitAssigner, TUAssigner):
    pass


class UAssignerNoDefault(UnitAssigner):
    def __init__(self, x, t):
        self.x = x
        self.t = t
        self._default_units = {'x': u.m, 't': u.s}


class TestUnitAssigner(unittest.TestCase):
    def setUp(self):
        self.ua = TUAssigner(3*u.nm, 1*u.ms)
        self.ua_bad = TUAssigner(3*u.s, 1*u.ms)  # Wrong units
        self.und = UAssignerNoDefault(3 * u.nm, 1 * u.ms)

    def test_check_dimensionality_units(self):
        self.ua._check_dimensionality_units()
        assert_raises(pint.DimensionalityError,
                      self.ua_bad._check_dimensionality_units)

    def test_dot_assign_wrong_units_does_nothing(self):
        self.ua.x = 5 * u.s

    def test_assign_wrong_units_raises_error(self):
        assert_raises(pint.DimensionalityError, self.ua.assign, 'x', 5*u.s)

    def test_get_units(self):
        ua = self.ua
        eq_(ua._default_units, {'x': u.m, 't': u.s})
        ua._get_default_units()
        eq_(ua._default_units, {'x': u.m.units, 't': u.s.units})

    def test_get_default_units_error_no_defaults(self):
        """Make sure _get_default_units raises a helpful error
        when called on an object with no default values."""
        assert_raises(AttributeError, self.und._get_default_units)


class TestUnit_to_unitless(unittest.TestCase):
    filename = 'test-pickle-unitless-Assigner.pkl'

    def setUp(self):
        self.unit_a = TUAssigner(3*u.nm, 1*u.ms)
        self.unit_a._unitless_units = {'[length]': u.um, '[time]': u.ms}
        self.no_unit_a = self.unit_a.to_unitless()

    def tearDown(self):
        silentremove(self.filename)

    def test_unit_to_unitless_magnitude_correct(self):
        eq_(self.no_unit_a.x, 0.003)
        eq_(self.no_unit_a.t, 1)

    def test_pickle(self):
        with open(self.filename, 'w') as f:
            pickle.dump(self.no_unit_a, f)

        with open(self.filename, 'r') as f:
            no_unit_unpickled = pickle.load(f)

        eq_(self.no_unit_a.x, no_unit_unpickled.x)


class TUAssigner2(UnitAssigner):
    def __init__(self, x=2*u.m, t=1*u.s):
        self.UnitlessClass = UnitlessTUAssigner2

        self.x = x
        self.t = t
        self._default_units = {'x': u.m, 't': u.s}

    def speed(self):
        return self.x / self.t


class UnitlessTUAssigner2(NoUnitAssigner, TUAssigner2):
    pass


class TestNoUnitAssigner(unittest.TestCase):
    def setUp(self):
        self.unitted = TUAssigner2(3*u.nm, 1*u.ms)
        self.unitted.what = "Car"
        self.unitted._unitless_units = {'[length]': u.um, '[time]': u.ms}
        self.nu = self.unitted.to_unitless()

    def test_verify_values_passed(self):
        eq_(self.nu.what, "Car")
        eq_(self.nu.speed(), 0.003)
