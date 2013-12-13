"""
Test UnitCantilever Class
2013-12-12

Ryan Dwyer
"""

from jittermodel.ucant import UnitCantilever, u
from nose.tools import assert_raises, assert_almost_equals
from pint import DimensionalityError

def test_UnitCantilever_input():
    """Make sure that defining a UnitCantilever with an incorrect geometry, or
    negative number raises a ValueError."""
    to_test = [{'f_c': -10 * u.Hz},
               {'Q': -39 * u.dimensionless},
               {'f_c': 40 * u.kHz, 'R_tip': -0.023 * u.nm},
               {'theta_tip': -10 * u.degrees},
               {'theta_tip': 100 * u.degrees},
               {'geometry_c': 'not perpendicular or parallel'}]

    for kwargs in to_test:
        assert_raises(ValueError, UnitCantilever, **kwargs)

def test_UnitCantilever_input_units():
    to_test = [{'f_c': -10 * u.s},
               {'Q': -39},
               {'f_c': 40 * u.kHz, 'R_tip': -0.023 * u.C},
               {'theta_tip': 12 * u.degK},
               {'theta_tip': 100}]
    
    for kwargs in to_test:
        assert_raises(DimensionalityError, UnitCantilever, **kwargs)


def test_UnitCantilever_init():
    """Make sure the unit cantilever initializes properly."""
    c = UnitCantilever(f_c=50 * u.kHz, k_c=3 * u.N/u.m,
                       Q=1000 * u.dimensionless)
