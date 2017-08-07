from __future__ import print_function, division, absolute_import

from future.utils import with_metaclass

from .fields import Field, ReverseForeignKey

__all__ = ("UnitMeta", "Unit")


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
        if len(bases) > 1:
            raise TypeError("Multiple inheritance is not supported for Units")
        if bases[0] is not object and bases[0] is Unit:
            self.fields = {}
        for k, v in dct.items():
            if isinstance(v, Field) or isinstance(v, ReverseForeignKey):
                if bases[0] is Unit:
                    v.attach(self, name=k)
                else:
                    raise ValueError("Only direct subclasses of Unit may have fields")
        unique = dct.get("unique", None)
        if unique is not None:
            for f in unique:
                if not isinstance(f, Field):
                    raise ValueError("Unique constraints must be Fields")

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
