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

import lsst.daf.butler.tests as butlerTests
from lsst.daf.butler import DatasetType, Butler, DataCoordinate
from lsst.daf.butler.core.utils import getFullTypeName

from lsst.obs.base.ingest_tests import IngestTestBase
from lsst.obs.base import RawIngestTask


TESTDIR = os.path.dirname(__file__)


class DummyCamRawIngestTask(RawIngestTask):
    """For DummyCam we ingest a different dataset type that can return
    a non-Exposure."""

    def getDatasetType(self):
        """Return the DatasetType of the datasets ingested by this Task.
        """
        return DatasetType("raw_dict", ("instrument", "detector", "exposure"), "StructuredDataDict",
                           universe=self.butler.registry.dimensions)


class RawIngestTestCase(IngestTestBase, unittest.TestCase):
    """Test ingest using JSON sidecar files."""

    ingestDatasetTypeName = "raw_dict"
    rawIngestTask = getFullTypeName(DummyCamRawIngestTask)
    curatedCalibrationDatasetTypes = ()
    ingestDir = TESTDIR
    instrumentClassName = "lsst.obs.base.instrument_tests.DummyCam"
    file = os.path.join(TESTDIR, "fakedata", "dataset_1.yaml")
    dataIds = [dict(instrument="DummyCam", exposure=100, detector=0)]

    @property
    def visits(self):
        butler = Butler(self.root, collections=[self.outputRun])
        return {
            DataCoordinate.standardize(
                instrument="DummyCam",
                visit=100,
                universe=butler.registry.dimensions
            ): [
                DataCoordinate.standardize(
                    instrument="DummyCam",
                    exposure=100,
                    universe=butler.registry.dimensions
                )
            ]
        }

    def testWriteCuratedCalibrations(self):
        """There are no curated calibrations in this test instrument"""
        pass


class RawIngestIndexTestCase(RawIngestTestCase):
    """Test ingest using JSON index files."""
    file = os.path.join(TESTDIR, "fakedata2", "dataset_1.yaml")


class TestRawIngestTaskPickle(unittest.TestCase):
    """Test that pickling of the RawIngestTask works properly."""

    @classmethod
    def setUpClass(cls):
        cls.root = tempfile.mkdtemp(dir=TESTDIR)
        cls.creatorButler = butlerTests.makeTestRepo(cls.root, {})

    @classmethod
    def tearDownClass(cls):
        if cls.root is not None:
            shutil.rmtree(cls.root, ignore_errors=True)

    def setUp(self):
        self.butler = butlerTests.makeTestCollection(self.creatorButler)

        self.config = RawIngestTask.ConfigClass()
        self.config.transfer = "copy"  # safe non-default value
        self.task = RawIngestTask(config=self.config, butler=self.butler)

    def testPickleTask(self):
        stream = pickle.dumps(self.task)
        copy = pickle.loads(stream)
        self.assertEqual(self.task.getFullName(), copy.getFullName())
        self.assertEqual(self.task.log.getName(), copy.log.getName())
        self.assertEqual(self.task.config, copy.config)
        self.assertEqual(self.task.butler._config, copy.butler._config)
        self.assertEqual(self.task.butler.collections, copy.butler.collections)
        self.assertEqual(self.task.butler.run, copy.butler.run)
        self.assertEqual(self.task.universe, copy.universe)
        self.assertEqual(self.task.datasetType, copy.datasetType)


if __name__ == "__main__":
    unittest.main()
