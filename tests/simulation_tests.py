"""
simulation_tests.py

Tests on the Simulation object.

Created by Ryan Dwyer on 2013-10-11.
Copyright (c) 2013 Cornell University. All rights reserved.
"""

import unittest
from nose.tools import assert_almost_equals
from jittermodel.simulation import sum_sinh
import mpmath as mp


def mp_sum_sinh(alpha):
    """Implements the infinite sum using mpmath, at very high precision.
    Method 'r+s+e' was found to work accurately for all values of alpha,
    unlike most other alogithms in Mathematica, python, etc."""
    summand = lambda n: mp.sinh(alpha) / mp.sinh(alpha * n)
    return mp.nsum(summand, [1, mp.inf], method='r+s+e')


def test_sum_sinh():
    """Test that the sum is working properly for a range of alpha values.
    The mpmath module is used to verify that the sum meets error
    specifications.
    """
    alphas = [2 ** i for i in xrange(-12, 10)]
    results = [sum_sinh(alpha) for alpha in alphas]
    mp_results = [mp_sum_sinh(alpha) for alpha in alphas]
    for mp_result, test_result in zip(mp_results, results):
        assert_almost_equals(mp_result, test_result, 7)


class simulation_tests(unittest.TestCase):

    def setUp(self):
        pass


# def test_func(func):
#     alphas = (2 ** i for i in xrange(-12, 10))
#     [func(alpha) for alpha in alphas]
#
# def test_func2(func):
#     alphas = (2 ** i for i in xrange(-2, 3))
#     [func(alpha) for alpha in alphas]
#
# if __name__ == '__main__':
#     from timeit import timeit
#     # print(timeit("test_func(mp_sum_sinh)",
#     #              setup="from __main__ import test_func, mp_sum_sinh",
#     #              number=100))
#     # print(timeit("test_func(sum_sinh)",
#     #              setup="from __main__ import test_func, sum_sinh",
#     #              number=100))
#     print(timeit("test_func2(mp_sum_sinh)",
#                  setup="from __main__ import test_func2, mp_sum_sinh",
#                  number=100))
#     print(timeit("test_func2(sum_sinh)",
#                  setup="from __main__ import test_func2, sum_sinh",
#                  number=100))

"""This test shows that the mpmath implementation of the infinite sum is
approximately 30 times slower over the whole range of values (see test_func)
(63.5 s vs. 2.17 s evaluation time). Over the range relevant to experiment,
it is about 100 times slower (46.4 s vs. 0.20 s)."""
