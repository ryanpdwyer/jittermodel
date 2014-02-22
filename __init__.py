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

u = pint.UnitRegistry()


def get_defaults(func):
    """Return a dictionary containing argument names and defaults.
    Dictionary contains None if an argument has no default value."""
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
    default_dict = get_defaults(func)
    return {name: val.units for name, val in default_dict.viewitems()
            if type(val) == u.Quantity}


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
        all = set([attr for attr in dir(self) if not attr.startswith('_')])

        # This prevents an infinite recursion when we run inspect.isroutine
        all.discard('_all_attributes')

        filtered = {attr for attr in all if not
                    inspect.isroutine(getattr(self, attr))}

        to_discard = set()
        for attr in filtered:
            try:
                value = getattr(self, attr)
                setattr(self, attr, value)
            except AttributeError:
                to_discard.add(attr)

        return filtered.difference(to_discard)

        return filtered

    def __eq__(self, other):
        """Define two Assigners to be equal if their internal
        dictionaries are the same. A relatively sane check for equality."""
        return self.__dict__ == other.__dict__


class UnitAssigner(Assigner):
    """An Assigner that uses units for numerical inputs.

    The unit behavior is specified by a dictionary self._default_units containing
    the name of the variable as the key, and default unit as the value."""

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
        class_name = self.__class__.__name__
        new_class_name = "No"+class_name
        unitless_class = type(new_class_name,
                              (NoUnitAssigner, self.__class__), {})

        unitless = unitless_class()

        for attr in self._all_attributes:
            val = self.lookup(attr)
            if attr in self._unitless_units:
                unit = self._unitless_units[attr]
                unitless.assign(attr, val.to(unit).magnitude)
            else:
                unitless.assign(attr, val)

        unitless._default_units = self._default_units
        unitless._unitless_units = self._unitless_units

        return unitless


class NoUnitAssigner(Assigner):
    """A class with blank, dummy, _get_default_units and
    _check_number_inputs_positive. Meant to be used as a mix-in class."""
    def to_unitted(self):
        """I need to look up how to make a member of a class without calling
        the class' __init__ method."""
        unit_class = self.__class__.__bases__[-1]
        init_keys = set(get_defaults(unit_class.__init__).keys())
        init_keys.remove('self')

        unit_attributes = {attr: self.lookup(attr) for attr
                           in self._all_attributes}

        for attr in self._all_attributes:
            unit_attributes[attr] = self.lookup(attr)
            if attr in self._unitless_units:
                unit_attributes[attr] *= self._unitless_units[attr]

        init_dict = {attr: val for attr, val in unit_attributes.viewitems()
                     if attr in init_keys}

        unitted = unit_class(**init_dict)

        for attr, val in unit_attributes.viewitems():
            unitted.assign(attr, val)

        return unitted

    def _get_default_units(self):
        raise AttributeError

    def _check_dimensionality_units(self):
        raise AttributeError

    def to_unitless(self):
        raise AttributeError

    def __init__(self):
        pass
