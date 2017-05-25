from __future__ import print_function, division, absolute_import

from future.utils import with_metaclass

__all__ = ("Field", "RegionField", "IntField", "StrField", "DateTimeField",
           "ForeignKey", "ReverseForeignKey", "Alias",
           "UnitMeta", "Unit", "SpatialUnit",)


class Field(object):

    def __init__(self, optional=False):
        self.name = None
        self.optional = optional

    def __get__(self, instance, owner=None):
        if instance is not None:
            return instance._storage.get(self.name, None)
        return self

    def attach(self, cls, name):
        self.name = name
        cls.fields[name] = self


class RegionField(Field):
    sqlType = "BLOB"


class IntField(Field):
    sqlType = "INTEGER"


class StrField(Field):
    sqlType = "TEXT"


class DateTimeField(Field):
    sqlType = "INTEGER"


class ForeignKey(Field):

    def __init__(self, UnitClass, reverse=None, optional=False):
        Field.__init__(self, optional=optional)
        self.UnitClass = UnitClass
        self._reverse = reverse

    def attach(self, cls, name):
        Field.attach(self, cls, name)
        if self._reverse is not None:
            reverse = getattr(self.UnitClass, self._reverse)
            reverse.attach(self.UnitClass, self._reverse)


class ReverseForeignKey(object):

    def __init__(self):
        self.name = None

    def __get__(self, instance, owner=None):
        if instance is not None:
            return instance._reversed.get(self.name, None)
        return self

    def attach(self, cls, name):
        self.name = name


class Alias(object):

    def __init__(self, local, remote):
        self.local = local
        self.remote = remote

    def __get__(self, instance, owner=None):
        if instance is not None:
            return self.remote.__get__(self.local.__get__(instance))
        return self


class UnitMeta(type):

    def __init__(self, name, bases, dct):
        self.fields = {}
        for k, v in dct.iteritems():
            if isinstance(v, Field):
                v.attach(self, name=k)
        unique = dct.get("unique", None)
        if unique is not None:
            for f in unique:
                if not isinstance(f, Field):
                    raise ValueError("Unique constraints must be Fields")


class Unit(with_metaclass(UnitMeta, object)):

    __metaclass__ = UnitMeta

    def __init__(self, storage=None):
        if storage is None:
            storage = {}
        self._storage = storage
        self._reversed = {}
        self.datasets = {}

    @property
    def id(self):
        return self._storage["id"]

    def __eq__(self, other):
        return type(self) == type(other) and self.id == other.id

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(type(self), self.id)


class SpatialUnit(Unit):

    region = RegionField()

    def __init__(self):
        Unit.__init__(self)
        self.overlapping = {}
