#!/usr/bin/env python

from __future__ import print_function, division, absolute_import

import unittest
import sqlite3

from lsst.obs.base import repodb


class UnitsTestCase(unittest.TestCase):

    def setUp(self):
        self.db = sqlite3.connect(":memory:")

    def testCreateTables(self):
        for UnitClass in repodb.COMMON_UNITS:
            print(repodb.sqlCreateTable(UnitClass))
            self.db.execute(repodb.sqlCreateTable(UnitClass))
        self.db.commit()


if __name__ == "__main__":
    unittest.main()
