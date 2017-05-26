#!/usr/bin/env python

from __future__ import print_function, division, absolute_import


import unittest

import lsst.obs.base.repodb.tests
from lsst.obs.base import repodb


class GraphTestCase(unittest.TestCase):

    def setUp(self):
        self.db = repodb.tests.makeRepoDatabase()

    def testUnitUniverseGraph(self):
        graph = self.db.makeGraph(self.db.UNIT_CLASSES)
        self.assertItemsEqual(graph.units.keys(), self.db.UNIT_CLASSES)
        self.assertItemsEqual(
            [unit.name for unit in graph.units[repodb.CameraUnit]],
            ["HSC"]
        )
        self.assertItemsEqual(
            [unit.name for unit in graph.units[repodb.SkyMapUnit]],
            ["DISCRETE_2"]
        )
        self.assertItemsEqual(
            [unit.name for unit in graph.units[repodb.FilterUnit]],
            "grizy"
        )
        # Much more to test here

    def tesDatasets(self):
        universe = self.db.makeGraph(FutureDatasets=(repodb.tests.Coadd,))
        restricted = self.db.makeGraph(NeededDatasets=(repodb.tests.Coadd,))
        self.assertItemsEqual(
            [unit.name for unit in restricted.units[repodb.FilterUnit]],
            ["r"]
        )
        self.assertLess(
            restricted.units[repodb.FilterUnit],
            universe.units[repodb.FilterUnit]
        )
        self.assertLess(
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
        universe = self.db.makeGraph(self.db.UNIT_CLASSES)
        restricted = self.db.makeGraph(self.db.UNIT_CLASSES,
                                       where="FilterUnit.name in ('g', 'r')")
        self.assertLess(
            restricted.units[repodb.FilterUnit],
            universe.units[repodb.FilterUnit],
        )
        self.assertEqual(
            set(filter.name for filter in restricted.units[repodb.FilterUnit]),
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
            graph1 = db1.makeGraph(db1.UNIT_CLASSES)
            graph2 = db2.makeGraph(db2.UNIT_CLASSES)
            self.assertEqual(graph1.units, graph2.units)


if __name__ == "__main__":
    unittest.main()
