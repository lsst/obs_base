from __future__ import print_function, division, absolute_import

from .unit import Unit
from .fields import StrField, IntField, ForeignKey, ReverseForeignKey, DateTimeField

__all__ = ("CameraUnit", "SkyMapUnit", "TractUnit", "PatchUnit",
           "AbstractFilterUnit", "PhysicalFilterUnit",
           "VisitUnit", "SensorUnit")


class CameraUnit(Unit):
    name = StrField()
    filters = ReverseForeignKey()
    visits = ReverseForeignKey()
    sensors = ReverseForeignKey()
    unique = (name,)

    def register(self, repodb):
        raise NotImplementedError()


class SkyMapUnit(Unit):
    name = StrField()
    tracts = ReverseForeignKey()
    unique = (name,)


class TractUnit(Unit):
    number = IntField()
    skymap = ForeignKey(SkyMapUnit, reverse="tracts")
    unique = (skymap, number)

    def getDataIdValue(self):
        return self.number


class PatchUnit(Unit):
    skymap = ForeignKey(SkyMapUnit)
    x = IntField()
    y = IntField()
    unique = (skymap, x, y)

    def getDataIdValue(self):
        return "{},{}".format(self.x, self.y)


class AbstractFilterUnit(Unit):
    name = StrField()
    physical = ReverseForeignKey()
    unique = (name,)

    def getDataIdValue(self):
        return self.name


class PhysicalFilterUnit(Unit):
    name = StrField()
    abstract = ForeignKey(AbstractFilterUnit, reverse="physical")
    camera = ForeignKey(CameraUnit, reverse="filters")
    visits = ReverseForeignKey()
    unique = (camera, name)

    def getDataIdValue(self):
        return self.name


class VisitUnit(Unit):
    camera = ForeignKey(CameraUnit, reverse="visits")
    number = IntField()
    filter = ForeignKey(PhysicalFilterUnit, reverse="visits")
    sensors = ReverseForeignKey()
    dateobs = DateTimeField()


class SensorUnit(Unit):
    camera = ForeignKey(CameraUnit, reverse="sensors")
    number = IntField()
