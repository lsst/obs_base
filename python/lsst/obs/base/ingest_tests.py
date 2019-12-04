# This file is part of obs_subaru.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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

"""Base class for writing Gen3 raw data ingest tests.
"""

__all__ = ("IngestTestBase", "InstrumentSignatureDataIds")

import abc
import dataclasses
import tempfile
import unittest
import os
import shutil

import numpy as np

import lsst.afw.image
from lsst.daf.butler import Butler
from lsst.obs.base import RawIngestTask, InstrumentSignatureDataIngestTask


@dataclasses.dataclass
class InstrumentSignatureDataIds:
    """Additions to the dataIds required to check that instrument signature
    data were properly ingested.

    Set values to `None` to not test that specific type of ingestion.
    """
    camera: dict()
    """We always need to test that a Camera was ingested."""
    bfKernel: dict() = None
    defects: dict() = None
    transmission_filter: dict() = None
    transmission_detector: dict() = None
    transmission_optics: dict() = None
    transmission_atmosphere: dict() = None

    instrument: dataclasses.InitVar[str] = None
    """This is added to each of the above as the "instrument" dimension."""

    def __post_init__(self, instrument):
        for field in self.__dict__:
            if getattr(self, field) is not None:
                getattr(self, field)['instrument'] = instrument


class IngestTestBase(metaclass=abc.ABCMeta):
    """Base class for tests of gen3 `~lsst.obs.base.Instrument` registration
    and data ingestion.

    Subclass from this, then `unittest.TestCase` to get a working test suite.
    """

    instrument = None
    """The instrument to be registered and tested."""

    dataId = {}
    """Butler data ID of a file to ingest when testing."""

    file = ""
    """Full path to a file to ingest in tests."""

    instrumentSignatureDataIds = None
    """DataIds to test for correct ingestion."""

    def setUp(self):
        self.root = tempfile.mkdtemp()
        Butler.makeRepo(self.root)
        self.butler = Butler(self.root, run="raw")

        # Register the instrument and its static metadata
        self.instrument.register(self.butler.registry)

        # Make a default config for test methods to play with
        self.config = RawIngestTask.ConfigClass()
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
        task = RawIngestTask(config=self.config, butler=self.butler)
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
        exposure = self.butler.get("raw", self.dataId)
        metadata = self.butler.get("raw.metadata", self.dataId)
        image = self.butler.get("raw.image", self.dataId)
        self.assertImagesEqual(exposure.image, image)
        self.assertEqual(metadata.toDict(), exposure.getMetadata().toDict())
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
        os.symlink(self.file, newPath)
        self.config.transfer = None
        self.runIngestTest([newPath])

    def testFailOnConflict(self):
        """Re-ingesting the same data into the repository should fail.
        """
        self.config.transfer = "symlink"
        self.runIngest()
        with self.assertRaises(Exception):
            self.runIngest()

    def testInstrumentSignatureDataIngestTask(self):
        """Test ingesting instrument signature data."""
        with self.assertRaisesRegex(KeyError, "bfKernel"):
            self.butler.get('bfKernel', dataId=self.instrumentSignatureDataIds.bfKernel)
        with self.assertRaisesRegex(KeyError, "camera"):
            self.butler.get('camera', dataId=self.instrumentSignatureDataIds.camera)
        with self.assertRaisesRegex(KeyError, "defects"):
            self.butler.get('defects', dataId=self.instrumentSignatureDataIds.defects)

        task = InstrumentSignatureDataIngestTask(config=self.config, butler=self.butler)
        task.run()

        camera = self.butler.get('camera', dataId=self.instrumentSignatureDataIds.camera)
        self.assertEqual(camera.getName(), self.instrument.getName())

        if self.instrumentSignatureDataIds.defects is not None:
            defects = self.butler.get('defects', dataId=self.instrumentSignatureDataIds.defects)
            self.assertIn("lsst.meas.algorithms.defects.Defects", str(type(defects)))
            # Can't use assertIsInstance because obs_base cannot have a dependency on meas_algorithms
            # self.assertIsInstance(defects, lsst.meas.algorithms.defects.Defects)

        if self.instrumentSignatureDataIds.bfKernel is not None:
            bfKernel = self.butler.get('bfKernel', dataId=self.instrumentSignatureDataIds.bfKernel)
            self.assertIsInstance(bfKernel, np.ndarray)

        def checkTransmissionCurve(name, dataId):
            """Check that butler.get(dataId) returns a TransmissionCurve."""
            if dataId is not None:
                transmission = self.butler.get(name, dataId)
                self.assertIsInstance(transmission, lsst.afw.image.TransmissionCurve)

        checkTransmissionCurve("transmission_optics", self.instrumentSignatureDataIds.transmission_optics)
        checkTransmissionCurve("transmission_sensor", self.instrumentSignatureDataIds.transmission_detector)
        checkTransmissionCurve("transmission_filter", self.instrumentSignatureDataIds.transmission_filter)
        checkTransmissionCurve("transmission_atmosphere",
                               self.instrumentSignatureDataIds.transmission_atmosphere)
