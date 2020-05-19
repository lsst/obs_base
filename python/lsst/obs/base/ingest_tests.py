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
import click.testing
import tempfile
import unittest
import os
import shutil

from lsst.daf.butler import Butler
from lsst.daf.butler.cli.butler import cli as butlerCli
import lsst.obs.base
from lsst.utils import doImport
from .utils import getInstrument


class IngestTestBase(metaclass=abc.ABCMeta):
    """Base class for tests of gen3 ingest. Subclass from this, then
    `unittest.TestCase` to get a working test suite.
    """

    ingestDir = ""
    """Root path to ingest files into. Typically `obs_package/tests/`; the
    actual directory will be a tempdir under this one.
    """

    dataIds = []
    """list of butler data IDs of files that should have been ingested."""

    file = ""
    """Full path to a file to ingest in tests."""

    rawIngestTask = "lsst.obs.base.RawIngestTask"
    """The task to use in the Ingest test."""

    curatedCalibrationDatasetTypes = None
    """List or tuple of Datasets types that should be present after calling
    writeCuratedCalibrations. If `None` writeCuratedCalibrations will
    not be called and the test will be skipped."""

    defineVisitsTask = lsst.obs.base.DefineVisitsTask
    """The task to use to define visits from groups of exposures.
    This is ignored if ``visits`` is `None`.
    """

    visits = {}
    """A dictionary mapping visit data IDs the lists of exposure data IDs that
    are associated with them.
    If this is empty (but not `None`), visit definition will be run but no
    visits will be expected (e.g. because no exposures are on-sky
    observations).
    """

    outputRun = "raw"
    """The name of the output run to use in tests.
    """

    @property
    @abc.abstractmethod
    def instrumentClassName(self):
        """The fully qualified instrument class name.

        Returns
        -------
        `str`
            The fully qualified instrument class name.
        """
        pass

    @property
    def instrumentClass(self):
        """The instrument class."""
        return doImport(self.instrumentClassName)

    @property
    def instrumentName(self):
        """The name of the instrument.

        Returns
        -------
        `str`
            The name of the instrument.
        """
        return self.instrumentClass.getName()

    def setUp(self):
        # Use a temporary working directory
        self.root = tempfile.mkdtemp(dir=self.ingestDir)
        self._createRepo()

        # Register the instrument and its static metadata
        self._registerInstrument()

    def tearDown(self):
        if os.path.exists(self.root):
            shutil.rmtree(self.root, ignore_errors=True)

    def verifyIngest(self, files=None, cli=False):
        """
        Test that RawIngestTask ingested the expected files.

        Parameters
        ----------
        files : `list` [`str`], or None
            List of files to be ingested, or None to use ``self.file``
        """
        butler = Butler(self.root, run=self.outputRun)
        datasets = butler.registry.queryDatasets(self.outputRun, collections=...)
        self.assertEqual(len(list(datasets)), len(self.dataIds))
        for dataId in self.dataIds:
            exposure = butler.get(self.outputRun, dataId)
            metadata = butler.get("raw.metadata", dataId)
            self.assertEqual(metadata.toDict(), exposure.getMetadata().toDict())

            # Since components follow a different code path we check that
            # WCS match and also we check that at least the shape
            # of the image is the same (rather than doing per-pixel equality)
            wcs = butler.get("raw.wcs", dataId)
            self.assertEqual(wcs, exposure.getWcs())

            rawImage = butler.get("raw.image", dataId)
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

    def _createRepo(self):
        """Use the Click `testing` module to call the butler command line api
        to create a repository."""
        runner = click.testing.CliRunner()
        result = runner.invoke(butlerCli, ["create", self.root])
        self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")

    def _ingestRaws(self, transfer):
        """Use the Click `testing` module to call the butler command line api
        to ingest raws.

        Parameters
        ----------
        transfer : `str`
            The external data transfer type.
        """
        runner = click.testing.CliRunner()
        result = runner.invoke(butlerCli, ["ingest-raws", self.root,
                                           "--output-run", self.outputRun,
                                           "--file", self.file,
                                           "--transfer", transfer,
                                           "--ingest-task", self.rawIngestTask])
        self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")

    def _registerInstrument(self):
        """Use the Click `testing` module to call the butler command line api
        to register the instrument."""
        runner = click.testing.CliRunner()
        result = runner.invoke(butlerCli, ["register-instrument", self.root,
                                           "--instrument", self.instrumentClassName])
        self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")

    def _writeCuratedCalibrations(self):
        """Use the Click `testing` module to call the butler command line api
        to write curated calibrations."""
        runner = click.testing.CliRunner()
        result = runner.invoke(butlerCli, ["write-curated-calibrations", self.root,
                                           "--instrument", self.instrumentName,
                                           "--output-run", self.outputRun])
        self.assertEqual(result.exit_code, 0, f"output: {result.output} exception: {result.exception}")

    def testLink(self):
        self._ingestRaws(transfer="link")
        self.verifyIngest()

    def testSymLink(self):
        self._ingestRaws(transfer="symlink")
        self.verifyIngest()

    def testCopy(self):
        self._ingestRaws(transfer="copy")
        self.verifyIngest()

    def testHardLink(self):
        try:
            self._ingestRaws(transfer="hardlink")
            self.verifyIngest()
        except PermissionError as err:
            raise unittest.SkipTest("Skipping hard-link test because input data"
                                    " is on a different filesystem.") from err

    def testInPlace(self):
        """Test that files already in the directory can be added to the
        registry in-place.
        """
        # symlink into repo root manually
        butler = Butler(self.root, run=self.outputRun)
        newPath = os.path.join(butler.datastore.root, os.path.basename(self.file))
        os.symlink(os.path.abspath(self.file), newPath)
        self._ingestRaws(transfer=None)
        self.verifyIngest()

    def testFailOnConflict(self):
        """Re-ingesting the same data into the repository should fail.
        """
        self._ingestRaws(transfer="symlink")
        with self.assertRaises(Exception):
            self._ingestRaws(transfer="symlink")

    def testWriteCuratedCalibrations(self):
        """Test that we can ingest the curated calibrations"""
        if self.curatedCalibrationDatasetTypes is None:
            raise unittest.SkipTest("Class requests disabling of writeCuratedCalibrations test")

        self._writeCuratedCalibrations()

        dataId = {"instrument": self.instrumentName}
        butler = Butler(self.root, run=self.outputRun)
        for datasetTypeName in self.curatedCalibrationDatasetTypes:
            with self.subTest(dtype=datasetTypeName, dataId=dataId):
                datasets = list(butler.registry.queryDatasets(datasetTypeName, collections=...,
                                                              dataId=dataId))
                self.assertGreater(len(datasets), 0, f"Checking {datasetTypeName}")

    def testDefineVisits(self):
        if self.visits is None:
            self.skipTest("Expected visits were not defined.")
        self._ingestRaws(transfer="link")

        config = self.defineVisitsTask.ConfigClass()
        butler = Butler(self.root, run=self.outputRun)
        instr = getInstrument(self.instrumentName, butler.registry)
        instr.applyConfigOverrides(self.defineVisitsTask._DefaultName, config)
        task = self.defineVisitsTask(config=config, butler=butler)
        task.run(self.dataIds)

        # Test that we got the visits we expected.
        visits = set(butler.registry.queryDimensions(["visit"], expand=True))
        self.assertCountEqual(visits, self.visits.keys())
        camera = instr.getCamera()
        for foundVisit, (expectedVisit, expectedExposures) in zip(visits, self.visits.items()):
            # Test that this visit is associated with the expected exposures.
            foundExposures = set(butler.registry.queryDimensions(["exposure"], dataId=expectedVisit,
                                                                 expand=True))
            self.assertCountEqual(foundExposures, expectedExposures)
            # Test that we have a visit region, and that it contains all of the
            # detector+visit regions.
            self.assertIsNotNone(foundVisit.region)
            detectorVisitDataIds = set(butler.registry.queryDimensions(["visit", "detector"],
                                                                       dataId=expectedVisit,
                                                                       expand=True))
            self.assertEqual(len(detectorVisitDataIds), len(camera))
            for dataId in detectorVisitDataIds:
                self.assertTrue(foundVisit.region.contains(dataId.region))
