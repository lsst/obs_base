from __future__ import print_function, division, absolute_import

from future.utils import with_metaclass

__all__ = ("Dataset",)


class UnitProperty(object):

    def __init__(self, name, UnitClass):
        self.name = name
        self.UnitClass = UnitClass

    def __get__(self, instance, owner=None):
        if instance is not None:
            return instance._storage.get(self.name, None)
        return self

    def finalize(self, instance, units):
        id = instance._storage.get(self.name, None)
        assert isinstance(id, int)
        assert id is not None
        target = units[self.UnitClass][id]
        instance._storage[self.name] = target
        rev = target.datasets.setdefault(type(instance), set())
        rev.add(instance)


class DatasetMeta(type):

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
            assert {k: v.UnitClass for k, v in cls.properties} == UnitClasses
            return r
        d = {k: UnitProperty(k, v) for k, v in UnitClasses.items()}
        d["properties"] = d.copy()
        d["name"] = name
        return type.__new__(cls, name, (Dataset,), d)

    def __init__(self, name, bases=None, dct=None, **UnitClasses):
        pass


class Dataset(with_metaclass(DatasetMeta, object)):

    @staticmethod
    def subclass(name, **UnitClasses):
        return DatasetMeta(name, **UnitClasses)

    def __init__(self, storage=None):
        if storage is None:
            storage = {}
        self._storage = storage

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
