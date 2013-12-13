"""
Test UnitCantilever Class
2013-12-12

Ryan Dwyer
"""

from jittermodel.ucant import UnitCantilever, u

def test_UnitCantilever_init():
    """Make sure the unit cantilever initializes properly."""
    c = UnitCantilever(f_c=50*u.kHz, k_c=3*u.N/u.m, Q=1000*u.dimensionless)
