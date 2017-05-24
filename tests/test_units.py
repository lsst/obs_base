#!/usr/bin/env python

from __future__ import print_function, division, absolute_import

import unittest

from lsst.obs.base import repodb


class UnitsTestCase(unittest.TestCase):

    def setUp(self):
        self.backend = repodb.SqliteBackend(":memory:")

    def testCreateTables(self):
        for UnitClass in repodb.COMMON_UNITS:
            self.backend.createTable(UnitClass)
        self.backend.commit()


if __name__ == "__main__":
    unittest.main()
