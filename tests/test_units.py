#!/usr/bin/env python

from __future__ import print_function, division, absolute_import

import unittest
import itertools
import datetime

from lsst.obs.base import repodb


def makeDateTime(hours):
    dt = datetime.datetime(2016, 5, 10, hours)
    return (dt - datetime.datetime(1970, 1, 1)).total_seconds()

EXAMPLE_DATA = [(
        repodb.TractUnit,
        ("id", "number", "skymap"),
        [(1, 8766, "RINGS_120"),
         (2, 8767, "RINGS_120")]
    ), (
        repodb.PatchUnit,
        ("id", "tract", "x", "y"),
        [(n, tract, x, y) for n, (tract, x, y) in enumerate(
            itertools.product((1, 2), range(8), range(8))
        )]
    ), (
        repodb.FilterUnit,
        ("id", "name", "camera"),
        [(n, f, "HSC") for n, f in enumerate("gri")]
    ), (
        repodb.VisitUnit,
        ("id", "number", "filter", "camera", "dateobs"),
        [(1, 1001, 0, "HSC", makeDateTime(10)),
         (2, 1002, 0, "HSC", makeDateTime(11)),
         (3, 1003, 1, "HSC", makeDateTime(12)),
         (4, 1004, 2, "HSC", makeDateTime(13))],
    ), (
        repodb.SensorUnit,
        ("id", "number", "visit"),
        [(n, number, visit) for n, (number, visit) in enumerate(
            itertools.product(range(45, 55), (1, 2, 3, 4))
         )]
    )]


class UnitsTestCase(unittest.TestCase):

    def setUp(self):
        self.backend = repodb.SqliteBackend(":memory:")
        for UnitClass in repodb.COMMON_UNITS:
            self.backend.createTable(UnitClass)
        for UnitClass, names, values in EXAMPLE_DATA:
            tableName = self.backend.getTableName(UnitClass)
            sql = "INSERT INTO {} ({}) VALUES ({})".format(
                tableName,
                ", ".join(names),
                ", ".join(["?"]*len(names))
            )
            self.backend.db.executemany(sql, values)
        self.backend.db.commit()

    def testCreateTables(self):
        pass

    def testUniverseGraph(self):
        graph = self.backend.makeGraph(repodb.COMMON_UNITS)
        self.assertItemsEqual(graph.units.keys(), repodb.COMMON_UNITS)
        for UnitClass, names, values in EXAMPLE_DATA:
            for unit, data in zip(graph.units[UnitClass], values):
                for i, name in enumerate(names):
                    attr = getattr(unit, name)
                    if isinstance(attr, repodb.Unit):
                        self.assertEqual(attr.id, data[i])
                    elif isinstance(attr, datetime.datetime):
                        self.assertEqual(
                            attr,
                            datetime.datetime.fromtimestamp(data[i])
                        )
                    else:
                        self.assertEqual(
                            attr,
                            data[i],
                            msg="{}.{}".format(UnitClass.__name__, name)
                        )


if __name__ == "__main__":
    unittest.main()
