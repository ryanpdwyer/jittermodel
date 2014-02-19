"""
__init__ tests
2013-12-16
Ryan Dwyer

Test general helper functions in the package and
"""
from __future__ import division
import pint
from jittermodel import (get_defaults, get_units, u, Assigner, UnitAssigner,
                         NoUnitAssigner)
import unittest
from nose.tools import eq_, assert_not_equal, assert_raises


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


def test_get_units():
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
        eq_(get_units(func), exp_dict)


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

    def test_all_attributes_with_method(self):
        class AssignerPlusMethod(Assigner):
            def a_method():
                pass

        a = AssignerPlusMethod()
        a.a = 2
        a.b = 10

        eq_({'a', 'b'}, a._all_attributes)


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
        self.units = {'x': u.m, 't': u.s}


class TestUnitAssigner(unittest.TestCase):
    def setUp(self):
        self.ua = TUAssigner(3*u.nm, 1*u.ms)
        self.ua_bad = TUAssigner(3*u.s, 1*u.ms)

    def test_check_dimensionality_units(self):
        self.ua._check_dimensionality_units()
        assert_raises(pint.DimensionalityError,
                      self.ua_bad._check_dimensionality_units)


# def test_unit_to_unitless():
#     unit_a = TUAssigner(3*u.nm, 1*u.ms)
#     unit_a._unitless_units = {'x': u.um, 't': u.ms}
#     no_unit_a = unit_a.to_unitless()
#     eq_(no_unit_a.x, 3000)
#     eq_(no_unit_a.t, 1)


to_eval = """class No{class_name}(NoUnitAssigner, {class_name}):
    pass"""

instantiate = "No{class_name}()"


class TUAssigner2(UnitAssigner):
    def __init__(self, x=2*u.m, t=1*u.s):
        self.x = x
        self.t = t
        self.units = {'x': u.m, 't': u.s}

    def speed(self):
        return self.x / self.t


# def test_mult_inherit_from_eval():
#     eval(to_eval.format(class_name='TUAssigner2')
#     tua2 = eval(instantiate.format(class_name='TUAssigner2'))