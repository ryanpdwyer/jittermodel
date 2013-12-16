"""
__init__ tests
2013-12-16
Ryan Dwyer

Test general helper functions in the package and
"""

from jittermodel import get_defaults, get_units, u


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
        assert get_defaults(func) == exp_dict


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
        assert get_units(func) == exp_dict
