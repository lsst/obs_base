from __future__ import print_function, division, absolute_import

from .unit import Unit
from .fields import StrField, IntField, ForeignKey, ReverseForeignKey, DateTimeField, LabeledObjectField

__all__ = ("TractUnit", "PatchUnit",
           "AbstractFilterUnit", "PhysicalFilterUnit",
           "VisitUnit", "SensorUnit")


class TractUnit(Unit):
    number = IntField()
    skymap = LabeledObjectField("skymap")
    unique = (skymap, number)

    def getDataIdValue(self):
        return self.number


class PatchUnit(Unit):
    skymap = LabeledObjectField("skymap")
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
    camera = LabeledObjectField("camera")
    visits = ReverseForeignKey()
    unique = (camera, name)

    def getDataIdValue(self):
        return self.name


class VisitUnit(Unit):
    camera = LabeledObjectField("camera")
    number = IntField()
    filter = ForeignKey(PhysicalFilterUnit, reverse="visits")
    sensors = ReverseForeignKey()
    dateobs = DateTimeField()


class SensorUnit(Unit):
    camera = LabeledObjectField("camera")
    number = IntField()
