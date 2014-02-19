"""
jittermodel
Ryan Dwyer
2013-12-13

TO DO: Add tests for the base classes.

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


def get_units(func):
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
        """Return a tuple of all the non-magic attributes of the class.

        See http://goo.gl/4juRRI for more information."""
        all = set([attr for attr in dir(self) if not attr.startswith('__')])

        # This prevents an infinite recursion when we run inspect.isroutine
        all.discard('_all_attributes')

        return {attr for attr in all if not
                inspect.isroutine(getattr(self, attr))}

    def __eq__(self, other):
        """Define two Assigners to be equal if their internal
        dictionaries are the same. A relatively sane check for equality."""
        return self.__dict__ == other.__dict__


class UnitAssigner(Assigner):
    """An Assigner that uses units for numerical inputs.

    The unit behavior is specified by a dictionary self.units containing
    the name of the variable as the key, and default unit as the value."""

    def _get_units(self):
        """Get the units of the arguments to initialize the object, inferring
        them from the class's __init__ method."""
        self.units = get_units(self.__init__)

    def _check_number_inputs_positive(self):
        """Return a ValueError if the number inputs are not positive."""
        greater_than_zero = self.units.viewkeys()

        for attr in greater_than_zero:
            if self.lookup(attr).magnitude <= 0:
                raise ValueError("The attribute '{attr}'\
must be positive.".format(attr=attr))

    def _check_dimensionality_units(self):
        """Return a DimensionalityError if unitted attributes
        (set in self.units) have the wrong dimensionality."""
        items = self.units.viewitems()
        for attr, unit in items:
            quant = self.lookup(attr)
            if type(quant) != u.Quantity:
                raise pint.DimensionalityError("No unit", unit)
            quant.to(unit)

    def to_unitless(self):
        """I want this to programatically generate a new, unitless,
        version of the current class. I'll need to reassign all values,
        and also use unitless_units to convert all the numbers."""
        current_class = self.__class__.__name__
        to_eval = """class No{class_name}(NoUnitAssigner, {class_name}):
    pass
"""
        eval(to_eval.format(class_name=current_class))
        new_class = eval("No{class_name}()".format(class_name=current_class))

        for attr, unit in self.unitless_units:
            new_class.assign(attr, unit)

class NoUnitAssigner(Assigner):
    """A class with blank, dummy, _get_units and
    _check_number_inputs_positive"""
    def _get_units(self):
        pass

    def _check_number_inputs_positive(self):
        pass

    def _check_dimensionality_units(self):
        pass

    def to_unitless(self):
        pass

    def __init__(self):
        pass
