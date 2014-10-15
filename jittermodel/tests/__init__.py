from jittermodel import u
import nose
import functools
from nose.tools import assert_almost_equal, assert_raises
import unittest


def pint_assert_almost_equal(first, second, unit=None, places=None,
                             msg=None, delta=None):
    """assert_almost_equal for pint quantities. Use unit to
    specify the unit to make the comparison in.

    The default behaviour converts the second quantity to the first unit."""
    if unit is None:
        second.ito(first.units)
    else:
        first.ito(unit)
        second.ito(unit)

    return assert_almost_equal(first.magnitude, second.magnitude,
                               places=places, msg=msg, delta=delta)


def expected_failure(test):
    """A wrapper so that tests which are expected to fail don't cause nose
    to crash.

    See http://stackoverflow.com/q/9613932/2823213."""
    @functools.wraps(test)
    def inner(*args, **kwargs):
        try:
            test(*args, **kwargs)
        except Exception:
            raise nose.SkipTest
        else:
            raise AssertionError('Failure expected')
    return inner
