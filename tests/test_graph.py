#!/usr/bin/env python

from __future__ import print_function, division, absolute_import


import unittest
import datetime

import lsst.obs.base.repodb.tests
from lsst.obs.base import repodb


class GraphTestCase(unittest.TestCase):

    def setUp(self):
        self.backend = repodb.tests.makeBackend()

    def testCreateTables(self):
        pass

    def testUnitUniverseGraph(self):
        graph = self.backend.makeGraph(repodb.COMMON_UNITS, ())
        self.assertItemsEqual(graph.units.keys(), repodb.COMMON_UNITS)
        for UnitClass, names, values in repodb.tests.EXAMPLE_UNITS:
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

    def testNeededDatasets(self):
        universe = self.backend.makeGraph(repodb.COMMON_UNITS)
        graph = self.backend.makeGraph(repodb.COMMON_UNITS,
                                       NeededDatasets=(repodb.tests.CalExp,))
        self.assertEqual(
            set(calexp.sensor.id for calexp
                in graph.datasets[repodb.tests.CalExp]),
            set(sensor.id for sensor in graph.units[repodb.SensorUnit])
        )
        self.assertLess(
            graph.units[repodb.SensorUnit],
            universe.units[repodb.SensorUnit]
        )
        # TODO: test various linkages

    def testFutureDatasets(self):
        universe = self.backend.makeGraph(repodb.COMMON_UNITS)
        graph = self.backend.makeGraph(repodb.COMMON_UNITS,
                                       FutureDatasets=(repodb.tests.CalExp,))
        self.assertLess(
            set(calexp.sensor.id for calexp
                in graph.datasets[repodb.tests.CalExp]),
            set(sensor.id for sensor in graph.units[repodb.SensorUnit])
        )
        self.assertEqual(graph.units[repodb.SensorUnit],
                         universe.units[repodb.SensorUnit])
        # TODO: test various linkages


if __name__ == "__main__":
    unittest.main()
