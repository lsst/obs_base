from __future__ import print_function, division, absolute_import

from . import base

__all__ = ("CameraUnit", "SkyMapUnit", "TractUnit", "PatchUnit", "FilterUnit",
           "VisitUnit", "SensorUnit")


class CameraUnit(base.Unit):
    name = base.StrField()
    filters = base.ReverseForeignKey()
    visits = base.ReverseForeignKey()
    unique = (name,)


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


class FilterUnit(base.Unit):
    name = base.StrField()
    camera = base.ForeignKey(CameraUnit, reverse="filters")
    visits = base.ReverseForeignKey()
    unique = (camera, name)

    def getDataIdValue(self):
        return self.name


class VisitUnit(base.Unit):
    number = base.IntField()
    camera = base.ForeignKey(CameraUnit, reverse="visits")
    filter = base.ForeignKey(FilterUnit, reverse="visits")
    sensors = base.ReverseForeignKey()
    dateobs = base.DateTimeField()
    unique = (camera, number)


class SensorUnit(base.Unit):
    number = base.IntField()
    camera = base.ForeignKey(CameraUnit)
    unique = (camera, number)


if False:  # TODO: put these back in after relating them better to visit/sensor

    class RawUnit(base.Unit):
        sensor = base.ForeignKey(SensorUnit, reverse="raw")
        visit = base.Alias(sensor, SensorUnit.visit)
        camera = base.Alias(sensor, SensorUnit.camera)
        filter = base.Alias(sensor, SensorUnit.filter)
        dateobs = base.Alias(sensor, SensorUnit.dateobs)
        name = base.StrField()
        unique = (sensor, name)

    class MasterCalibrationUnit(base.Unit):
        begin = base.DateTimeField()
        end = base.DateTimeField()
        filter = base.ForeignKey(FilterUnit, reverse=None, optional=True)
        camera = base.StrField()
        unique = (begin, end, filter, camera)
