#!/usr/bin/env python

from __future__ import print_function, division, absolute_import


import unittest

import lsst.obs.base.repodb.tests
from lsst.obs.base import repodb


class GraphTestCase(unittest.TestCase):

    def setUp(self):
        self.db = repodb.tests.makeRepoDatabase()

    def testUnitUniverseGraph(self):
        graph = self.db.makeGraph([repodb.CameraUnit, repodb.SkyMapUnit, repodb.AbstractFilterUnit])
        self.assertItemsEqual(graph.units.keys(), [repodb.CameraUnit, repodb.SkyMapUnit,
                                                   repodb.AbstractFilterUnit])
        self.assertItemsEqual(
            [unit.name for unit in graph.units[repodb.CameraUnit]],
            ["HSC"]
        )
        self.assertItemsEqual(
            [unit.name for unit in graph.units[repodb.SkyMapUnit]],
            ["DISCRETE_2"]
        )
        self.assertItemsEqual(
            [unit.name for unit in graph.units[repodb.AbstractFilterUnit]],
            "ugrizy"
        )
        # Much more to test here

    def testDatasets(self):
        universe = self.db.makeGraph(FutureDatasets=(repodb.tests.Coadd,))
        restricted = self.db.makeGraph(NeededDatasets=(repodb.tests.Coadd,))
        self.assertItemsEqual(
            [unit.name for unit in restricted.units[repodb.AbstractFilterUnit]],
            ["r"]
        )
        self.assertLess(
            restricted.units[repodb.AbstractFilterUnit],
            universe.units[repodb.AbstractFilterUnit]
        )
        self.assertEqual(
            restricted.datasets[repodb.tests.Coadd],
            universe.datasets[repodb.tests.Coadd]
        )
        self.assertItemsEqual(
            [(d.tract, d.patch, d.filter.name)
             for d in restricted.datasets[repodb.tests.Coadd]],
            [(tract, patch, 'r')
             for tract in restricted.units[repodb.TractUnit]
             for patch in restricted.units[repodb.PatchUnit]]
        )
        self.assertItemsEqual(
            [(d.tract, d.patch, d.filter.name)
             for d in universe.datasets[repodb.tests.Coadd]],
            [(tract, patch, 'r')
             for tract in universe.units[repodb.TractUnit]
             for patch in universe.units[repodb.PatchUnit]]
        )

    def testUserWhere(self):
        universe = self.db.makeGraph([repodb.AbstractFilterUnit])
        restricted = self.db.makeGraph([repodb.AbstractFilterUnit],
                                       where="AbstractFilterUnit.name in ('g', 'r')")
        self.assertLess(
            restricted.units[repodb.AbstractFilterUnit],
            universe.units[repodb.AbstractFilterUnit],
        )
        self.assertEqual(
            set(filter.name for filter in restricted.units[repodb.AbstractFilterUnit]),
            set("gr")
        )

    def testPickle(self):
        with lsst.utils.tests.getTempFilePath(".sqlite") as tmpFile:
            import pickle
            db1 = repodb.tests.makeRepoDatabase(tmpFile)
            s = pickle.dumps(db1, protocol=pickle.HIGHEST_PROTOCOL)
            db2 = pickle.loads(s)
            self.assertEqual(db1._cameras.keys(), db2._cameras.keys())
            self.assertEqual(db1._skyMaps.keys(), db2._skyMaps.keys())
            graph1 = db1.makeGraph(db1.DEFAULT_UNIT_CLASSES)
            graph2 = db2.makeGraph(db2.DEFAULT_UNIT_CLASSES)
            self.assertEqual(graph1.units, graph2.units)


if __name__ == "__main__":
    unittest.main()
