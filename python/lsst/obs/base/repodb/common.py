from __future__ import print_function, division, absolute_import

from . import base

__all__ = ("CameraUnit", "SkyMapUnit", "TractUnit", "PatchUnit",
           "AbstractFilterUnit", "PhysicalFilterUnit",
           "VisitUnit", "SensorUnit")


class CameraUnit(base.Unit):
    name = base.StrField()
    filters = base.ReverseForeignKey()
    visits = base.ReverseForeignKey()
    sensors = base.ReverseForeignKey()
    unique = (name,)

    def register(self, repodb):
        raise NotImplementedError()


class SkyMapUnit(base.Unit):
    name = base.StrField()
    tracts = base.ReverseForeignKey()
    unique = (name,)


class TractUnit(base.Unit):
    number = base.IntField()
    skymap = base.ForeignKey(SkyMapUnit, reverse="tracts")
    unique = (skymap, number)

    def getDataIdValue(self):
        return self.number


class PatchUnit(base.Unit):
    skymap = base.ForeignKey(SkyMapUnit)
    x = base.IntField()
    y = base.IntField()
    unique = (skymap, x, y)

    def getDataIdValue(self):
        return "{},{}".format(self.x, self.y)


class AbstractFilterUnit(base.Unit):
    name = base.StrField()
    physical = base.ReverseForeignKey()
    unique = (name,)

    def getDataIdValue(self):
        return self.name


class PhysicalFilterUnit(base.Unit):
    name = base.StrField()
    abstract = base.ForeignKey(AbstractFilterUnit, reverse="physical")
    camera = base.ForeignKey(CameraUnit, reverse="filters")
    visits = base.ReverseForeignKey()
    unique = (camera, name)

    def getDataIdValue(self):
        return self.name


class VisitUnit(base.Unit):
    camera = base.ForeignKey(CameraUnit, reverse="visits")
    number = base.IntField()
    filter = base.ForeignKey(PhysicalFilterUnit, reverse="visits")
    sensors = base.ReverseForeignKey()
    dateobs = base.DateTimeField()


class SensorUnit(base.Unit):
    camera = base.ForeignKey(CameraUnit, reverse="sensors")
    number = base.IntField()
