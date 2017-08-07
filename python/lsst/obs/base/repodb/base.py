from __future__ import print_function, division, absolute_import

from future.utils import with_metaclass

__all__ = ("Field", "RegionField", "IntField", "StrField", "DateTimeField",
           "ForeignKey", "ReverseForeignKey", "Alias",
           "UnitMeta", "Unit", "categorizeUnit")


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


class Alias(object):
    """A descriptor for `Unit` classes that makes an attribute of a another
    `Unit` class (related via e.g. `ForeignKey`) available for convenience.

    `Alias` does not result in any extra database fields, and is hence not a
    true `Field`.

    Parameters
    ----------
    local : descriptor
        A `Field` or other `Unit` descriptor instance in the same class as
        the `Alias` (typically a`ForeignKey`) that can be used to retrieve an
        instance of a different `Unit` class.
    remote : descriptor
        A `Field` or other `Unit` descriptor instance that retrieves the
        a value from an instance of the `Unit` class retreived by `local`.

    TODO: point to an example in common.py when we have one.
    """

    def __init__(self, local, remote):
        self.local = local
        self.remote = remote

    def __get__(self, instance, owner=None):
        if instance is not None:
            return self.remote.__get__(self.local.__get__(instance))
        return self


class UnitMeta(type):
    """A custom metaclass for `Unit` and all types derived from it.

    `UnitMeta` is responsible for calling `attach` on all `Field`s present
    in a `Unit` class definition and making sure unique constraints are valid.

    Design Notes
    ------------
    `UnitMeta` probably could enforce the requirement that no `Fields` be added
    to a `Unit` type after it has been defined by overriding `setattr`.  It's
    not clear if that's necessary, given how rarely `Unit`s will be defined and
    the complexity involved in implementing `setattr`.
    """

    def __init__(self, name, bases, dct):
        self.fields = {}
        self.aliases = {}
        for k, v in dct.items():
            if isinstance(v, Field) or isinstance(v, ReverseForeignKey):
                v.attach(self, name=k)
            if isinstance(v, Alias):
                self.aliases[k] = v
        unique = dct.get("unique", None)
        if unique is not None:
            for f in unique:
                if not isinstance(f, Field):
                    raise ValueError("Unique constraints must be Fields")
        if len(bases) > 1:
            raise TypeError("Multiple inheritance is not supported for Units")
        if bases[0] is not object and bases[0] is not Unit and bases[0].__bases__[0] is not Unit:
            raise TypeError("Unit classes must inherit directly or at one remove from Unit.")


NO_COMPARE_MESSAGE = ("Cannot compare Units with no ID (to add an ID, "
                      "insert it into a RepoDatabase)")


class Unit(with_metaclass(UnitMeta, object)):
    """An abstract base class for types that represent units of data.

    `Unit`s are used to represent concepts that can be used to label
    `Dataset`s within a repository, such as a visit (an observation with a
    particular camera) or a patch (an area of sky defined by a skymap). There
    are also `Unit` classes to represent cameras and skymaps themselves,
    allowing multiple cameras and skymaps to coexist within the same
    repository.  The set of available `Unit` types is static: it may not be
    modified at runtime.

    `Units` are explicitly camera-generic: cameras may provide custom
    identifiers to label a `Unit` (implementation TBD), but these are always
    considered aliases to one or more camera-generic core identifiers that can
    fully define the `Unit`.  This allows `Unit`s to be used to define
    `Dataset`s in a camera-independent way, and it allows algorithmic code to
    perform `Unit`-dependent operations without having to work around camera-
    dependent concepts (e.g. aggregating all sensors in a visit).

    The creation of `Unit` *instances* is also restricted: while `SuperTask`s
    will almost always create new `Dataset`s in a repository, they cannot in
    general create new `Unit`s with which to label them. Instead, `Unit`s are
    added by *registration* commands, not pipeline processing commands.  This
    includes both the familiar process of ingesting raw data and the process
    of registering a new skymap for deep-sky data products.

    The `RepoDatabase` provides an interface to the `Unit`s in a repository
    (as well as the `Datasets`), and is responsible for loading `Unit`s from
    an on-disk SQL database to Python.  Usually `Units` loaded from the
    `Database` are part of a `RepoGraph`, which may be defined to represent
    only a subset of the `Unit`s in a repository.

    In Python, a `Unit` is typically defined by a set of *descriptors*:
    property-like objects that each define a particular attribute.  Descriptors
    that inherit from `Field` represent values that are stored in the SQL
    database, and these are used (by implementatons of the `Backend` interface)
    to define a database schema to store `Unit`s.

    Construction and Unique IDs
    ---------------------------
    Most `Unit` classes may be constructed simply by passing keyword arguments
    whose key-value pairs correspond to their `Field` attributes.

    A `Unit` that is created directly in Python (i.e. not loaded from the
    database) is considered "incomplete": it will have its `id` attribute set
    to `None`, and it cannot be compared with other `Unit`s or included in
    `set`s or `dict`s that require complete `Units`.  Adding the `Unit` to the
    database will fill in its `id` attribute and complete it.

    Because a complete `Unit` uses its `id` field for equality comparison, it
    is possible that `Unit` instances with some different content will compare
    as equal if they are associated with different `RepoGraph`s that placed
    different restrictions on associated `Unit`s or `Datasets`.  `Units` that
    compare as equal always refer to the same underlying label in the database,
    however.

    Associated Datasets
    -------------------
    In addition to the attributes provided by a `Unit`'s descriptors, all
    `Unit`s have a `datasets` dictionary that holds `Dataset`s
    that are labeled by that `Unit` (possibly in some restricted context; see
    `RepoGraph`).  The keys of the `datasets` dictionary are `Dataset` class
    types, and values are `set`s of `Dataset` instances.
    """

    unique = None

    def __init__(self, **kwds):
        self._storage = kwds
        self._storage.setdefault("id", None)
        self._reversed = {}
        self.datasets = {}

    @property
    def id(self):
        """Return an integer ID that uniquely identifies this unit in its data
        repository.
        """
        return self._storage["id"]

    def __eq__(self, other):
        if self.id is None:
            raise ValueError(NO_COMPARE_MESSAGE)
        return type(self) == type(other) and self.id == other.id

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        if self.id is None:
            raise ValueError(NO_COMPARE_MESSAGE)
        return self.id

    def __repr__(self):
        items = ["id={}".format(self.id)]
        if self.unique is not None:
            for f in self.unique:
                items.append("{}={}".format(f.name, repr(f.__get__(self))))
        return "{}({})".format(type(self).__name__, ", ".join(items))

    def getDataIdValue(self):
        """Return a value for a dictionary-style data ID.

        Design Notes
        ------------
        This interface is a temporary approach to making `Unit`-based datasets
        work with the `Butler` before the `RepoDatabase` system is fully
        integrated with the `Butler`/`Mapper`.  It should eventually be
        unnecessary.
        """
        raise NotImplementedError()


def categorizeUnit(UnitClass):
    """Categorize a Unit classes.

    Returns
    -------
    BaseClass : `type`
        The immediate base class of the given `Unit` type, or None if the
        immediate base class is `Unit` itself.
    hasTable : `bool`
        Whether the class requires its own table in a SQL representation of
        the Unit.

    Classes that inherit directly from Unit are "primary" units, and always
    have `hasTable==True`, even if they have no fields (besides the automatic
    primary key).  "Derived" units, which inherit from a primary unit, may or
    may not have a table; if they do, they must have a one-to-one relationship
    with the primary unit's table.  Units that inherit from a derived unit and
    multiple inheritance are not supported.
    """
    BaseClass = UnitClass.__bases__[0]
    if BaseClass is Unit:
        return (None, True)
    return (BaseClass, len(UnitClass.fields) > 0)
