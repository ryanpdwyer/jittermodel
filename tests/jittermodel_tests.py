"""
__init__ tests
2013-12-16
Ryan Dwyer

Test general helper functions in the package and
"""
from __future__ import division
import pint
from numpy import pi
from jittermodel import (get_defaults, get_default_units, u, Assigner,
                         UnitAssigner, NoUnitAssigner, q2unitless)
from jittermodel.tests import pint_assert_almost_equal
import unittest
from nose.tools import eq_, assert_not_equal, assert_raises, assert_almost_equal

pa_eq = pint_assert_almost_equal


def test_get_defaults():
    def f1(x):
        pass

    def f2(x=3, y=5):
        pass

    def f3(x, y='adfe'):
        pass

    funcs = [f1, f2, f3]
    exp_default_dicts = [{'x': None}, {'x': 3, 'y': 5},
                         {'x': None, 'y': 'adfe'}]

    for func, exp_dict in zip(funcs, exp_default_dicts):
        eq_(get_defaults(func), exp_dict)


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
        self.x = x
        self.t = t
        self._default_units = {'x': u.m, 't': u.s}


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

    def test_get_units(self):
        ua = self.ua
        eq_(ua._default_units, {'x': u.m, 't': u.s})
        ua._get_default_units()
        eq_(ua._default_units, {'x': u.m.units, 't': u.s.units})

    def test_get_default_units_error_no_defaults(self):
        """Make sure _get_default_units raises a helpful error
        when called on an object with no default values."""
        assert_raises(AttributeError, self.und._get_default_units)


def test_unit_to_unitless():
    unit_a = TUAssigner(3*u.nm, 1*u.ms)
    unit_a._unitless_units = {'[length]': u.um, '[time]': u.ms}
    no_unit_a = unit_a.to_unitless()
    eq_(no_unit_a.x, 0.003)
    eq_(no_unit_a.t, 1)


class TUAssigner2(UnitAssigner):
    def __init__(self, x=2*u.m, t=1*u.s):
        self.x = x
        self.t = t
        self._default_units = {'x': u.m, 't': u.s}

    def speed(self):
        return self.x / self.t


class TestNoUnitAssigner(unittest.TestCase):
    def setUp(self):
        self.unitted = TUAssigner2(3*u.nm, 1*u.ms)
        self.unitted.what = "Car"
        self.unitted._unitless_units = {'[length]': u.um, '[time]': u.ms}
        self.nu = self.unitted.to_unitless()

    def test_verify_values_passed(self):
        eq_(self.nu.what, "Car")
        eq_(self.nu.speed(), 0.003)

    # def test_to_unitted(self):
    #     """Make sure that we can make an object with units again
    #     from a NoUnit object."""
    #     unitted = self.nu.to_unitted()
    #     pa_eq(unitted.x, self.unitted.x)
    #     pa_eq(unitted.t, self.unitted.t)
    #     pa_eq(unitted.speed(), self.unitted.speed())
    #     eq_(unitted.what, self.unitted.what)

    # def test_get_default_units(self):
    #     assert_raises(AttributeError, self.nu._get_default_units)

    # def test_check_dimensionality_units(self):
    #     assert_raises(AttributeError, self.nu._check_dimensionality_units)

    def test_to_unitless(self):
        assert_raises(AttributeError, self.nu.to_unitless)
