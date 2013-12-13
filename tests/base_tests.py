#!/usr/bin/env python
# encoding: utf-8
"""
sample_tests.py

Created by Ryan Dwyer on 2013-10-11.
Copyright (c) 2013 Cornell University. All rights reserved.

Look at http://stackoverflow.com/q/9613932/2823213 for expected failures.
"""

import unittest
from nose.tools import assert_raises, assert_almost_equals
from jittermodel.base import Sample, Cantilever
from numpy import pi


def test_Cantilever_input():
    """Make sure that defining a cantilever with an incorrect geometry, or
    negative number raises a ValueError."""
    to_test = [{'f_c': -10},
               {'Q': -39},
               {'f_c': 40, 'R_tip': -0.023},
               {'ThetaDegrees_tip': -10},
               {'ThetaDegrees_tip': 100},
               {'geometry_c': 'not perpendicular or parallel'}]

    for kwargs in to_test:
        assert_raises(ValueError, Cantilever, **kwargs)


def test_Cantilever_eq():
    """Test the cantilever equal method."""
    c1 = Cantilever()
    c2 = Cantilever(f_c=c1.f_c)
    c3 = Cantilever(Q=c1.Q * 2)
    assert c1 == c2
    assert c1 != c3


def test_Cantilever_repr():
    """Make sure that the __repr__method regenerates
    an equivalent cantilever object."""
    c1 = Cantilever()
    assert c1 == eval(repr(c1))


class TestCantilever(unittest.TestCase):
    """Test the Cantilever properties to make sure they update
    automatically and raise errors when necessary."""
    def setUp(self):
        """Set up a default sample."""
        self.cantilever = Cantilever()

    def test_property_update(self):
        """Test that properties update appropriately after
        changing parameters."""
        c = self.cantilever
        old_f_c = c.f_c
        c.f_c = old_f_c * 2
        self.assertEqual(c.f_c, old_f_c * 2)
        self.assertAlmostEqual(c.omega_c, 2 * pi * c.f_c)

    def test_property_setting(self):
        """Check that we cannot update properties such as omega_c."""
        c = self.cantilever
        omega_c = c.omega_c
        with self.assertRaises(AttributeError):
            c.omega_c = 2
        self.assertEqual(c.omega_c, omega_c)

    def test_repr(self):
        self.assertEqual


def test_sample_eq():
    s1 = Sample(rho=8600000)
    s2 = Sample(h=s1.h, rho=8600000)
    s3 = Sample(h=s1.h * 2, rho=8600000)
    assert s1 == s2
    assert s1 != s3


# def test_Sample_repr():
#     s1 = Sample(rho=1000000)
#     assert s1 == eval(repr(s1))


def test_Sample_V_g_rho_inconsistent():
    """Make sure that the Sample class gives a Value Error
    when invalid combinations of V_g and rho are passed to
    the function."""
    with assert_raises(ValueError):
        Sample(V_g=10, rho=122)


def test_Sample_V_g_rho():
    s1 = Sample()
    s2 = Sample(rho=s1.rho)
    assert_almost_equals(s1.V_g, s2.V_g)
    assert_almost_equals(s1.rho, s2.rho)


class TestSample(unittest.TestCase):
    """Test the sample class, specifically to make sure that
    its properties update automatically, and raise errors when
    necessary."""
    def setUp(self):
        """Set up a default sample."""
        self.sample = Sample()

    def test_property_update(self):
        """Test that properties update appropriately
        after change parameters."""
        s = self.sample
        old_h = s.h
        s.h = old_h * 2
        self.assertEqual(s.h, old_h * 2)
        self.assertAlmostEqual(s.h_diel, s.h - s.h_trans)

    def test_property_setting(self):
        """Check that we cannot update properties such as h_diel."""
        s = self.sample
        h_diel = s.h_diel
        with self.assertRaises(AttributeError):
            s.h_diel = 2
        self.assertEqual(s.h_diel, h_diel)


def test_Experiment():
    """Experiment currently has no tests; it is basically just
a container at this point."""
    pass


if __name__ == '__main__':
    unittest.main()
