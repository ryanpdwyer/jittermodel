# -*- coding: utf-8 -*-
"""
jittermodel
Ryan Dwyer
2013-12-13

TO DO: Add tests for the base classes.

Should _units be _default_units? I think so. That way we can transfer those
properties between the UnitAssigner and NoUnitAssigner completely intact.

"""

import pint
import inspect
import os
import errno

u = pint.UnitRegistry()

# Universal Constants
E_0 = 1*u.epsilon_0
k_B = 1*u.boltzmann_constant
q = 1*u.elementary_charge


def get_defaults(func):
    """Return a dictionary containing argument names and defaults for a
    function. Dictionary contains None if an argument has no default value."""
    argspec = inspect.getargspec(func)
    names = argspec.args
    if argspec.defaults is not None:
        vals = argspec.defaults
    else:
        vals = ()  # Empty tuple, so it has length zero.

    # Handles missing default arguments by inserting None
    n_missing_default = len(names) - len(vals)
    vals_with_nones = [None for i in xrange(n_missing_default)]
    vals_with_nones.extend(vals)

    return dict(zip(names, vals_with_nones))


def get_default_units(func):
    """Return a dictionary of the units of a function
    from the default values."""
    default_dict = get_defaults(func)
    return {name: val.units for name, val in default_dict.viewitems()
            if isinstance(val, u.Quantity)}


def make_units(dim_dict, base_dict):
    """Takes a dimension dictionary in the form returned by the pint
    method `dimensionality`, along with a set of base_units, and creates
    appropriate unit."""
    units = u.dimensionless
    for dimension, power in dim_dict.items():
        units *= base_dict[dimension] ** power
    return units


def q2unitless(quant, base_dict):
    """Convert a pint quanitity to the unitless variables implied by
    base_dict. Usage:

    >>> speed = 10 * u.m/u.s
    >>> units = {'[length]': u.mm, '[time]': u.us}
    >>> q2unitless(speed, units)
    0.01

    """
    units = make_units(quant.dimensionality, base_dict)
    return quant.to(units).magnitude


class Assigner(object):
    """This class provides an update method, that allows Updater class
    properties to be updated using the syntax,

    >>> u1 = Assigner()
    >>> u1.assign('f_c', 50)
    >>> print(u1.f_c)
    50

    """

    def assign(self, attr, val):
        """This assign method allows an attribute 'attr' (given as a string)
        to be assigned to the value 'val'."""
        setattr(self, attr, val)

    def lookup(self, attr):
        """Lookup an attribute in the class with a string."""
        return getattr(self, attr)

    @property
    def _all_attributes(self):
        """Return a tuple of all the non-magic, non-hidden attributes
        of the class.

        See http://goo.gl/4juRRI for more information."""
        all_attrs = set([attr for attr in dir(self) if not attr.startswith('_')])

        # This prevents an infinite recursion when we run inspect.isroutine
        all_attrs.discard('_all_attributes')

        filtered = {attr for attr in all_attrs if not
                    inspect.isroutine(getattr(self, attr))}

        to_discard = set()
        for attr in filtered:
            try:
                value = getattr(self, attr)
                setattr(self, attr, value)
            except AttributeError:
                to_discard.add(attr)

        return filtered.difference(to_discard)

    def __eq__(self, other):
        """Define two Assigners to be equal if their internal
        dictionaries are the same. A relatively sane check for equality."""
        return self.__dict__ == other.__dict__


class UnitAssigner(Assigner):
    """An Assigner that uses units for numerical inputs.

    The unit behavior is specified by a dictionary self._default_units
    containing the name of the variable as the key,
    and default unit as the value."""

    def assign(self, attr, val):
        setattr(self, attr, val)
        self._check_dimensionality_units()

    def _get_default_units(self):
        """Get the units of the arguments to initialize the object, inferring
        them from the class's __init__ method."""
        _units = get_default_units(self.__init__)
        if _units == {}:
            raise AttributeError("No default values with units.")
        else:
            self._default_units = _units

    def _check_number_inputs_positive(self):
        """Return a ValueError if the number inputs are not positive."""
        greater_than_zero = self._default_units.viewkeys()

        for attr in greater_than_zero:
            if self.lookup(attr).magnitude <= 0:
                raise ValueError("The attribute '{attr}'\
must be positive.".format(attr=attr))

    def _check_dimensionality_units(self):
        """Return a DimensionalityError if unitted attributes
        (set in self._default_units  have the wrong dimensionality."""
        items = self._default_units.viewitems()
        for attr, unit in items:
            quant = self.lookup(attr)
            if type(quant) != u.Quantity:
                raise pint.DimensionalityError("No unit", unit)
            quant.to(unit)

    def to_unitless(self):
        """I want this to programatically generate a new, unitless,
        version of the current class. I'll need to reassign all values,
        and also use unitless_units to convert all the numbers.

        See http://www.webcitation.org/6NUjbigHH."""
        try:
            self._unitless
        except AttributeError:
            self._unitless = self.UnitlessClass()

        # Take all unitted quanities, and make them unitless before reassigning
        # them.
        for attr in self._all_attributes:
            val = self.lookup(attr)
            if isinstance(val, u.Quantity):
                unitless_val = q2unitless(val, self._unitless_units)
                self._unitless.assign(attr, unitless_val)
            else:
                # No units. Just assign the attr to its value with no changes.
                self._unitless.assign(attr, val)


        # We want the unitless version of the object to retain the default unit,
        # unitless unit variables.
        self._unitless._default_units = self._default_units
        self._unitless._unitless_units = self._unitless_units

        return self._unitless


class NoUnitAssigner(object):
    """A class with blank, dummy, _get_default_units and
    _check_number_inputs_positive. Meant to be used as a mix-in class, with
    some descendent of UnitAssigner."""

    def _get_default_units(self):
        pass

    def _check_dimensionality_units(self):
        pass

    def to_unitless(self):
        raise AttributeError


def silentremove(filename):
    """If a file exists, delete it. Otherwise, return nothing.
       See http://stackoverflow.com/q/10840533/2823213"""
    try:
        os.remove(filename)
    except OSError as e:  # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occured


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
