from __future__ import print_function, division, absolute_import

__all__ = ("Field", "RegionField", "IntField", "StrField", "DateTimeField",
           "ForeignKey", "ReverseForeignKey")


class Field(object):
    """An abstract base class for `Unit` descriptors whose values are stored
    in a database.

    `Field`s must only be used as attributes on subclasses of `Unit`, and they
    must be included in the `class` statement body directly when the class is
    first defined (they may not be assignd as class attributes later).

    The set of allowed `Fields` is fully defined in this file and may not be
    extended elsewhere, as any class in the polymorphic `Backend` hierarchy
    must have specialized code for dealing with each type of `Field`.

    Design Notes
    ------------
    Values for `Fields` are stored in a "_storage" dict in the unit instance,
    with keys corresponding to the attribute name of the descriptor.

    When constructed in a class body, a `Field` has its `name` attribute set
    to `None`, and is hence unusuable.  The `UnitMeta` class calls the `attach`
    method on each `Field` when a new `Unit` class is defined, passing in the
    name of the class attribute the `Field` occupies so the `Field` can
    initialize its `name` attribute and add itself to the `Unit` class's
    `fields` dictionary.  Subclasses may override `attach` to do additional
    initialization.

    Overriding the `finalize` method gives each `field` an opportunity to help
    initialize an *instance* of the `Unit` to which it is attached.  The base
    class implementation of this method does nothing, as `Fields` can generally
    assume that they will be initalized with the correct types in the
    `_storage` dict.
    """

    def __init__(self, optional=False):
        self.name = None
        self.optional = optional

    def __get__(self, instance, owner=None):
        if instance is not None:
            return instance._storage.get(self.name, None)
        return self

    def attach(self, cls, name):
        """Attach this `Field` to a new `Unit` type.

        Should only be called by `UnitMeta`.  Subclasses that override this
        method must call the base class implementation.

        Parameters
        ----------
        cls : `type` derived from `Unit`
            `Unit` class to which this field has been added.
        name : `str`
            Name of the class attribute used to store this field; will be used
            to set the `Field`'s internal name attribute'
        """
        self.name = name
        cls.fields[name] = self

    def finalize(self, instance, others):
        """Complete initialization of a `Unit` instance holding this field.

        Should only be called by `Backend` implementations when constructing a
        `Unit` instance from SQL query results.

        See `ForeignKey` for the only nontrivial implementation of this method.

        Parameters
        ----------
        instance : `Unit`
            An instance of the `Unit` subclass to which this `Field` is
            attaached.  The instance will have a `_storage` dict will all
            field keys present, and `finalize` may modify the dictonaries
            values to coerce them to the types needed by the `Field` or
            otherwise complete their initialization.
        others : nested `dict` of `{UnitClass: {`int`: `Unit`}}`
            A nested dictionary containing other partially-initialized `Units`
            instances, indexed first by class type and then by integer ID.
        """
        pass


class RegionField(Field):
    """A Field class that holds a spatial region on the sky.

    Regions in units are expected to be approximate and inclusive: they should
    cover *at least* the true of the entity the represent, but may cover more.
    Using a too-large regions may impact performance, but not correctness.

    Parameters
    ----------
    optional : `bool`
        If True, the field's value may be `None` in Python or NULL in a
        database.

    Design Notes
    ------------
    The actual Python type used by RegionField is TBD, and hence the class is
    essentially a placeholder right now.  It is envisioned to be something like
    a `sphgeom.RangeSet` associated with a `sphgeom.Pixelization` (with all
    regions in a database using the same pixelization).
    """
    pass


class IntField(Field):
    """A Field class that holds an integer.

    Parameters
    ----------
    optional : `bool`
        If True, the field's value may be `None` in Python or NULL in a
        database.
    """
    pass


class StrField(Field):
    """A Field class that holds an string.

    Parameters
    ----------
    optional : `bool`
        If True, the field's value may be `None` in Python or NULL in a
        database.
    """
    pass


class DateTimeField(Field):
    """A Field class that holds a date/time value.

    Parameters
    ----------
    optional : `bool`
        If True, the field's value may be `None` in Python or NULL in a
        database.
    """
    pass


class ForeignKey(Field):
    """A Field class that represents a link to another Unit type.

    Parameters
    ----------
    UnitClass : `type` derived from `Unit`
        The type of the instances this field holds.
    reverse : `str` or `None`
        If not `None`, the name of a `ReverseForeignKey` attribute in
        `UnitClass` that holds a set containing the `Unit` instances with
        `ForeignKey` fields that point to it.
    optional : `bool`
        If True, the field's value may be `None` in Python or NULL in a
        database.

    Design Notes
    ------------
    When being constructed by a `Backend` from a SQL query, a `ForeignKey`'s
    storage value may be set to the integer ID of the foreign `Unit` instance,
    rather than the instance itself, and `ReverseForeignKey` `sets` may be left
    unfilled.  When this is the case, `ForeignKey.finalize` may be called to
    transform the integer ID values to `Unit` instances and populate the
    `ReverseForeignKey` sets.
    """

    def __init__(self, UnitClass, reverse=None, optional=False):
        Field.__init__(self, optional=optional)
        self.UnitClass = UnitClass
        self._reverse = reverse

    def attach(self, cls, name):
        """Attach this `Field` to a new `Unit` type.

        Should only be called by `UnitMeta`.  Subclasses that override this
        method must call the base class implementation.

        The `attach` implementation for `ForeignKey` takes care of attaching
        its corresponding `ReverseForeignKey` to its class
        (`ReverseForeignKey` is not a `Field`, and is not attached directly by
        `UnitMeta`). Putting all of this logic in `ForeignKey` avoids circular
        dependency problems.

        Parameters
        ----------
        cls : `type` derived from `Unit`
            `Unit` class to which this field has been added.
        name : `str`
            Name of the class attribute used to store this field; will be used
            to set the `Field`'s internal name attribute'
        """
        Field.attach(self, cls, name)
        if self._reverse is not None:
            reverse = getattr(self.UnitClass, self._reverse)
            reverse.attach(self.UnitClass, self._reverse)

    def finalize(self, instance, others):
        """
        Complete initialization of a `Unit` instance holding this field.

        Should only be called by `Backend` implementations when constructing a
        `Unit` instance from SQL query results.

        When being constructed by a `Backend` from a SQL query, a
        `ForeignKey`'s storage value may be set to the integer ID of the
        foreign `Unit` instance, rather than the instance itself, and
        `ReverseForeignKey` `sets` may be left unfilled.  When this is the
        case, `ForeignKey.finalize` may be called to transform the integer ID
        values to `Unit` instances and populate the `ReverseForeignKey` sets.
        """
        id = instance._storage.get(self.name, None)
        if id is None:
            assert self.optional
            return
        target = others[self.UnitClass][id]
        instance._storage[self.name] = target
        rev = target._reversed.setdefault(self._reverse, set())
        rev.add(instance)


class ReverseForeignKey(object):
    """A descriptor for `Unit` classes that captures a one-to-many
    relationship.

    A `ReverseForeignKey` must be paired with a `ForeignKey` class in the
    `Unit` subclass to which it refers, and is not a true `Field` becaues it
    relies on the `ForeignKey` for its initialization and database storage.

    `ReverseForeignKey`s values are stored in a `_reversed` dictionary
    attribute on the `Unit` instance that holds it, and are typically populated
    by `ForeignKey.finalize`.

    `ReverseForeignKey`s may be empty and unusable when the database query that
    generates the `Unit`s does not extend to the `Unit` type referred to by
    the `ReverseForeignKey` (queries are always expanded to include the `Unit`s
    referred to by a `ForeignKey`, but this behavior is deliberately asymettric
    because a `Unit`'s `ReverseForeignKey` values cannot be used to uniquely
    define it).'
    """

    def __init__(self):
        self.name = None

    def __get__(self, instance, owner=None):
        if instance is not None:
            return instance._reversed.get(self.name, None)
        return self

    def attach(self, cls, name):
        """Attach this descriptor to a `Unit` type.

        Should only be called by `ForeignKey.attach`.

        Parameters
        ----------
        name : `str`
            Name of the class attribute used to store this field; will be used
            to set the descriptor's internal name attribute'
        """
        self.name = name

