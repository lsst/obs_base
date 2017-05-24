from __future__ import print_function, division, absolute_import

from future.utils import with_metaclass

__all__ = ("Field", "RegionField", "IntField", "StrField", "DateTimeField",
           "PythonTypeField", "ForeignKey", "ReverseForeignKey", "Alias",
           "UnitMeta", "Unit", "SpatialUnit",
           "sqlCreateTable")


class Field(object):

    def __init__(self, optional=False):
        self.name = None
        self._parent = None
        self.optional = optional

    def __get__(self, instance, owner=None):
        if instance is not None:
            return instance._storage.get(self.name, None)
        return self

    def attach(self, cls, name):
        self.name = name
        self.parent = cls
        cls.fields[name] = self


class RegionField(Field):
    sqlType = "BLOB"


class IntField(Field):
    sqlType = "INTEGER"


class StrField(Field):
    sqlType = "TEXT"


class DateTimeField(Field):
    sqlType = "DATETIME"


class PythonTypeField(Field):
    sqlType = "TEXT"


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

    @property
    def sqlType(self):
        return "INTEGER REFERENCES {} (id)".format(self.UnitClass.__name__)


class ReverseForeignKey(Field):
    sqlType = None


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
                elif isinstance(f, ReverseForeignKey):
                    raise ValueError(
                        "Unique constraints must not be ReverseForeignKey"
                    )


class Unit(with_metaclass(UnitMeta, object)):

    __metaclass__ = UnitMeta

    def __init__(self, storage=None):
        if storage is None:
            storage = {}
        self._storage = storage
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


def sqlCreateTable(UnitClass):
    items = ["id INTEGER PRIMARY KEY"]
    for k, v in UnitClass.fields.items():
        if v.sqlType is None:
            continue
        t = [k, v.sqlType]
        if not v.optional:
            t.append("NOT NULL")
        items.append(" ".join(t))
    items.append(
        "UNIQUE ({})".format(
            ", ".join(f.name for f in UnitClass.unique)
        )
    )
    return "CREATE TABLE {} (\n    {}\n)".format(
        UnitClass.__name__, ",\n    ".join(items)
    )
