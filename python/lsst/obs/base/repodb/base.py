from __future__ import print_function, division, absolute_import

from future.utils import with_metaclass

__all__ = ("Field", "RegionField", "IntField", "StrField", "DateTimeField",
           "ForeignKey", "ReverseForeignKey", "Alias",
           "UnitMeta", "Unit",)


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

    def finalize(self, instance, others):
        pass


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

    def finalize(self, instance, others):
        id = instance._storage.get(self.name, None)
        if id is None:
            assert self.optional
            return
        target = others[self.UnitClass][id]
        instance._storage[self.name] = target
        rev = target._reversed.setdefault(self._reverse, set())
        rev.add(instance)


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


NO_COMPARE_MESSAGE = ("Cannot compare Units with no ID (to add an ID, "
                      "insert it into a RepoDatabase)")


class Unit(with_metaclass(UnitMeta, object)):

    __metaclass__ = UnitMeta

    def __init__(self, **kwds):
        self._storage = kwds
        self._storage.setdefault("id", None)
        self._reversed = {}
        self.datasets = {}

    @property
    def id(self):
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
