from jittermodel import u
from nose.tools import assert_almost_equal, assert_raises


def pint_assert_almost_equal(first, second, unit=None, places=None, msg=None, delta=None):
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


def test_pint_assert_almost_equal():
    first = 4.77464829276e-07 * u.N / u.kHz / u.m
    second = 477.464829276 * u.pN * u.s / u.m
    pint_assert_almost_equal(first, second)


def test_pint_assert_almost_equal_not():
    first = 4.78e-07 * u.N / u.kHz / u.m
    second = 477.464829276 * u.pN * u.s / u.m
    assert_raises(AssertionError, pint_assert_almost_equal, **{'first': first,
                  'second': second, 'places': 7, 'unit': u.pN * u.s / u.m})
