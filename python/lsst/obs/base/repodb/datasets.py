from __future__ import print_function, division, absolute_import

from future.utils import with_metaclass

__all__ = ("Dataset",)


class UnitProperty(object):
    """A descriptor class used to define `Unit` attributes for `Dataset`s

    All `Unit` attributes on a `Dataset` class must be defined by
    `UnitProperty` instances, and each `Dataset` class holds a `properties`
    dictionary containing all `UnitProperty` instances for easy iteration.

    Because `Dataset` types are not defined by users in the usual way, there
    should be no need for users to construct them directly, and no reason to
    interact with them other than for normal attribute access on instances.

    Design Notes
    ------------
    Because all `UnitProperty` instances are created and added to their
    `Dataset` subclasses by `DatasetMeta`, there's no need for the two-step
    construct-attach pattern used by the `Unit` descriptors (which is used to
    save the user from having to repeat the name of the attribute);
    `UnitProperty` can just take its name in `__init__`.
    """

    def __init__(self, name, UnitClass):
        self.name = name
        self.UnitClass = UnitClass

    def __get__(self, instance, owner=None):
        if instance is not None:
            return instance._storage.get(self.name, None)
        return self

    def finalize(self, instance, units):
        """Finish constructing a `Dataset` instance, translating integer ID
        values in a `Dataset` instances `_storage` dict to the `Unit` instances
        they corresond to.

        Should only be called by `Backend` implementations.

        Parameters
        ----------
        instance : `Dataset`
            Instance of the `Dataset` type to which this property belongs.
            Must have a `_storage` dict containing a integer IDs, not `Unit`
            instances.
        units : nested `dict` of `{UnitClass: {`int`: `Unit`}}
            Nested dictionary of `Unit` instances, indexed first by `Unit` type
            objects and then by integer ID.
        """
        id = instance._storage.get(self.name, None)
        assert isinstance(id, int)
        assert id is not None
        target = units[self.UnitClass][id]
        instance._storage[self.name] = target
        rev = target.datasets.setdefault(type(instance), set())
        rev.add(instance)


class DatasetMeta(type):
    """A custom metaclass for `Dataset` and types derived from it.

    `DatasetMeta` prevents construction of new `Dataset` types through the
    `class` statement, ensuring consistent construction through
    `Dataset.subclass`.  This includes registering all `Dataset` types in a
    global registry to ensure `Dataset` type names are unique; trying to
    define the same `Dataset` type in multiple places will return the same
    type object, and trying to define different `Dataset`s with the same name
    will raise `ValueError`.
    """

    registry = {}

    def __new__(cls, name, bases=None, dct=None, **UnitClasses):
        if len(UnitClasses) == 0:
            if name != "Dataset" or bases[0] != object:
                raise ValueError(
                    "Only Dataset.subclass can create new Dataset types."
                )
            return type.__new__(cls, name, bases, dct)
        r = cls.registry.get(name, None)
        if r is not None:
            assert {k: v.UnitClass for k, v in r.properties.items()} \
                == UnitClasses
            return r
        d = {k: UnitProperty(k, v) for k, v in UnitClasses.items()}
        d["properties"] = d.copy()
        d["name"] = name
        r = type.__new__(cls, name, (Dataset,), d)
        DatasetMeta.registry[name] = r
        return r

    def __init__(self, name, bases=None, dct=None, **UnitClasses):
        pass


class Dataset(with_metaclass(DatasetMeta, object)):
    """An abstract base class for types that define data products.

    A subclass of `Dataset` is conceptually just the combination of a
    unique string name with a dictionary of `Unit` classes, and a
    `Dataset` instance just provides values for those `Unit`s.

    `Dataset` types may *only* be defined using the `subclass` static method
    (not the usual Python `class` statement).  This reflects the fact that a
    `Dataset` type is usually created dynamically by the `SuperTask`s that
    produced it when the `SuperTask` is configured.

    The `Dataset` instances in a `RepoDatabase` represent the actual data
    products in that repository (or repositories chained to it).  It is
    expected that `Butler.put` will be responsible for updating the
    `RepoDatabase` when data products are added.  `Dataset` instances may be
    added to `RepoGraph` before being added to `RepoDatabase` to represent
    data products that *should* be produced by a `Pipeline` run.
    """

    @staticmethod
    def subclass(name, **UnitClasses):
        """Create a new type that inherits from `Dataset`.

        Parameters
        ----------
        name : `str`
            A unique name for this type.  If a `Dataset` type with this name
            has already been created, that type will be returned instread of
            creating a new one.

        Additional keyword arguments are interpreted as the names and `Unit`
        classes of the units of data that label the new `Dataset`.

        Exceptions
        ----------
        Raises `ValueError` if another `Dataset` with the same name but
        different units already exists.
        """
        return DatasetMeta(name, **UnitClasses)

    def __init__(self, **kwds):
        self._storage = kwds

    def __eq__(self, other):
        return self._storage == other._storage

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(tuple(self._storage.items()))

    def __repr__(self):
        items = []
        for p in self.properties.values():
            items.append("{}={}".format(p.name, repr(p.__get__(self))))
        return "({}: {})".format(self.name, ", ".join(items))

    def getDataId(self):
        """Return a dictionary-style data ID corresponding to this `Dataset`.

        This is a temporary approach that allows `Dataset` types to be used
        with the existing `Mapper`/`Butler` implementations before they have
        been integrated with the `RepoDatabase` system.  It should not be
        permanent.
        """
        return {p.name: p.__get__(self).getDataIdValue()
                for p in self.properties.values()}

    def get(self, butler, **kwds):
        """Read a `Dataset` using a `Butler`.

        This is a temporary interface that should be replaced by a version
        of `Butler.get` that accepts a `Dataset` instance in the future.
        """
        return butler.get(self.name, self.getDataId(), **kwds)

    def put(self, butler, value, **kwds):
        """Write a `Dataset` using a `Butler`.

        This is a temporary interface that should be replaced by a version
        of `Butler.put` that accepts a `Dataset` instance in the future.
        """
        return butler.put(value, self.name, self.getDataId(), **kwds)

    def exists(self, butler):
        """Test whether a `Dataset` exists on disk using a `Butler`.

        This is a temporary interface that delegates to `Butler.datasetExists`.
        In the future, `Butler` should guarantee that any `Dataset` in a
        `RepoDatabase` exists on disk, and this method may not need to exist
        at all.
        """
        return butler.datasetExists(self.name, self.getDataId())
