"""
__init__ tests
2013-12-16
Ryan Dwyer

Test general helper functions in the package and
"""

from jittermodel import (get_defaults, get_units, u, Assigner)
import unittest
from nose.tools import eq_, assert_not_equal


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
