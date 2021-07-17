# This file is part of obs_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import pickle
import shutil
import tempfile
import unittest

import lsst.daf.butler as dafButler
import lsst.daf.butler.tests as butlerTests

from lsst.obs.base import DefineVisitsTask


TESTDIR = os.path.dirname(__file__)


class DefineVisitsTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Create a new butler once only."""
        cls.root = tempfile.mkdtemp(dir=TESTDIR)

        dataIds = {
            "instrument": ["DummyCam"],
            "physical_filter": ["d-r"],
            "exposure": [42, 43, 44],
            "visit": [42, 43, 44],
        }

        cls.creatorButler = butlerTests.makeTestRepo(cls.root, dataIds)

        # Create dataset types used by the tests
        cls.storageClassFactory = dafButler.StorageClassFactory()
        for datasetTypeName, storageClassName in (("raw", "ExposureF"),
                                                  ):
            storageClass = cls.storageClassFactory.getStorageClass(storageClassName)
            butlerTests.addDatasetType(cls.creatorButler,
                                       datasetTypeName,
                                       {"instrument", "exposure"},
                                       storageClass)

    @classmethod
    def tearDownClass(cls):
        if cls.root is not None:
            shutil.rmtree(cls.root, ignore_errors=True)

    def setUp(self):
        self.butler = butlerTests.makeTestCollection(self.creatorButler)

        self.config = DefineVisitsTask.ConfigClass()
        self.config.computeVisitRegions.active.padding = 42  # non-default value
        self.task = DefineVisitsTask(config=self.config, butler=self.butler)

    def testPickleTask(self):
        stream = pickle.dumps(self.task)
        copy = pickle.loads(stream)
        self.assertEqual(self.task.getFullName(), copy.getFullName())
        self.assertEqual(self.task.log.name, copy.log.name)
        self.assertEqual(self.task.config, copy.config)
        self.assertEqual(self.task.butler._config, copy.butler._config)
        self.assertEqual(self.task.butler.collections, copy.butler.collections)
        self.assertEqual(self.task.butler.run, copy.butler.run)
        self.assertEqual(self.task.universe, copy.universe)


if __name__ == "__main__":
    unittest.main()
