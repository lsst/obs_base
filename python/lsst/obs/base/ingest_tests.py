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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Base class for writing Gen3 raw data ingest tests.
"""

__all__ = ("IngestTestBase",)

import abc
import tempfile
import unittest
import os
import shutil

from lsst.daf.butler import Butler
import lsst.obs.base


class IngestTestBase(metaclass=abc.ABCMeta):
    """Base class for tests of gen3 ingest. Subclass from this, then
    `unittest.TestCase` to get a working test suite.
    """

    ingestDir = ""
    """Root path to ingest files into. Typically `obs_package/tests/`; the
    actual directory will be a tempdir under this one.
    """

    instrument = None
    """The instrument to be registered and tested."""

    dataIds = []
    """list of butler data IDs of files that should have been ingested."""

    file = ""
    """Full path to a file to ingest in tests."""

    RawIngestTask = lsst.obs.base.RawIngestTask
    """The task to use in the Ingest test."""

    curatedCalibrationDatasetTypes = None
    """List or tuple of Datasets types that should be present after calling
    writeCuratedCalibrations. If `None` writeCuratedCalibrations will
    not be called and the test will be skipped."""

    def setUp(self):
        # Use a temporary working directory
        self.root = tempfile.mkdtemp(dir=self.ingestDir)
        Butler.makeRepo(self.root)
        self.butler = Butler(self.root, run="raw")

        # Register the instrument and its static metadata
        self.instrument.register(self.butler.registry)

        # Make a default config for test methods to play with
        self.config = self.RawIngestTask.ConfigClass()
        self.config.instrument = \
            f"{self.instrument.__class__.__module__}.{self.instrument.__class__.__name__}"

    def tearDown(self):
        if os.path.exists(self.root):
            shutil.rmtree(self.root, ignore_errors=True)

    def runIngest(self, files=None):
        """
        Initialize and run RawIngestTask on a list of files.

        Parameters
        ----------
        files : `list` [`str`], or None
            List of files to be ingested, or None to use ``self.file``
        """
        if files is None:
            files = [self.file]
        task = self.RawIngestTask(config=self.config, butler=self.butler)
        task.log.setLevel(task.log.FATAL)  # silence logs, since we expect a lot of warnings
        task.run(files)

    def runIngestTest(self, files=None):
        """
        Test that RawIngestTask ingested the expected files.

        Parameters
        ----------
        files : `list` [`str`], or None
            List of files to be ingested, or None to use ``self.file``
        """
        self.runIngest(files)
        datasets = self.butler.registry.queryDatasets('raw', collections=...)
        self.assertEqual(len(list(datasets)), len(self.dataIds))
        for dataId in self.dataIds:
            exposure = self.butler.get("raw", dataId)
            metadata = self.butler.get("raw.metadata", dataId)
            self.assertEqual(metadata.toDict(), exposure.getMetadata().toDict())

            # Since components follow a different code path we check that
            # WCS match and also we check that at least the shape
            # of the image is the same (rather than doing per-pixel equality)
            # Check the observation type before trying to check WCS
            obsType = self.butler.registry.expandDataId(dataId).records["exposure"].observation_type
            if obsType == "science":
                wcs = self.butler.get("raw.wcs", dataId)
                self.assertEqual(wcs, exposure.getWcs())

            rawImage = self.butler.get("raw.image", dataId)
            self.assertEqual(rawImage.getBBox(), exposure.getBBox())

            self.checkRepo(files=files)

    def checkRepo(self, files=None):
        """Check the state of the repository after ingest.

        This is an optional hook provided for subclasses; by default it does
        nothing.

        Parameters
        ----------
        files : `list` [`str`], or None
            List of files to be ingested, or None to use ``self.file``
        """
        pass

    def testSymLink(self):
        self.config.transfer = "symlink"
        self.runIngestTest()

    def testCopy(self):
        self.config.transfer = "copy"
        self.runIngestTest()

    def testHardLink(self):
        self.config.transfer = "hardlink"
        try:
            self.runIngestTest()
        except PermissionError as err:
            raise unittest.SkipTest("Skipping hard-link test because input data"
                                    " is on a different filesystem.") from err

    def testInPlace(self):
        """Test that files already in the directory can be added to the
        registry in-place.
        """
        # symlink into repo root manually
        newPath = os.path.join(self.butler.datastore.root, os.path.basename(self.file))
        os.symlink(os.path.abspath(self.file), newPath)
        self.config.transfer = None
        self.runIngestTest([newPath])

    def testFailOnConflict(self):
        """Re-ingesting the same data into the repository should fail.
        """
        self.config.transfer = "symlink"
        self.runIngest()
        with self.assertRaises(Exception):
            self.runIngest()

    def testWriteCuratedCalibrations(self):
        """Test that we can ingest the curated calibrations"""
        if self.curatedCalibrationDatasetTypes is None:
            raise unittest.SkipTest("Class requests disabling of writeCuratedCalibrations test")

        self.instrument.writeCuratedCalibrations(self.butler)

        dataId = {"instrument": self.instrument.getName()}
        for datasetTypeName in self.curatedCalibrationDatasetTypes:
            with self.subTest(dtype=datasetTypeName, dataId=dataId):
                datasets = list(self.butler.registry.queryDatasets(datasetTypeName, collections=...,
                                                                   dataId=dataId))
                self.assertGreater(len(datasets), 0, f"Checking {datasetTypeName}")
