# This file is part of obs_base.
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

"""Unit test base class for the gen2 to gen3 converter.
"""

import abc
import itertools
import shutil
import tempfile
import unittest


import lsst.afw.image
import lsst.afw.table
import lsst.daf.persistence
import lsst.daf.butler
import lsst.meas.algorithms
from lsst.obs.base.script import convert
import lsst.utils.tests
from lsst.utils import doImport


class ConvertGen2To3TestCase(metaclass=abc.ABCMeta):
    """Test the `butler convert` command.

    Subclass this, and then `lsst.utils.tests.TestCase` and set the below
    attributes.  Uses the `butler convert` command line command to do the
    conversion.
    """

    gen2root = ""
    """Root path to the gen2 repo to be converted."""

    gen2calib = None
    """Path to the gen2 calib repo to be converted."""

    @property
    @abc.abstractmethod
    def instrumentClassName(self):
        """Full path to the `Instrument` class of the data to be converted,
        e.g. ``lsst.obs.decam.DarkEnergyCamera``.

        Returns
        -------
        className : `str`
            The fully qualified instrument class name.
        """
        pass

    @property
    def instrumentClass(self):
        """The instrument class."""
        return doImport(self.instrumentClassName)

    @property
    def instrumentName(self):
        """Name of the instrument for the gen3 registry, e.g. "DECam".

        Returns
        -------
        name : `str`
            The name of the instrument.
        """
        return self.instrumentClass.getName()

    config = None
    """Full path to a config override for ConvertRepoTask, to be applied after
    the Instrument overrides when running the conversion function."""

    biases = []
    """List dataIds to use to load gen3 biases to test that they exist."""

    biasName = "bias"
    """Name of the dataset that the biases are loaded into."""

    flats = []
    """List dataIds to use to load gen3 flats to test that they exist."""

    flatName = "flat"
    """Name of the dataset that the flats are loaded into."""

    darks = []
    """List dataIds to use to load gen3 darks to test that they exist."""

    darkName = "dark"
    """Name of the dataset that the darks are loaded into."""

    kwargs = {}
    """Other keyword arguments to pass directly to the converter function,
    as a dict."""

    refcats = []
    """Names of the reference catalogs to query for the existence of in the
    converted gen3 repo."""

    collections = set()
    """Additional collections that should appear in the gen3 repo.

    This will automatically be populated by the base `setUp` to include
    ``"{instrumentName}/raw"``, ``"refcats"`` (if the ``refcats``
    class attribute is non-empty), and ``"skymaps"`` (if ``skymapName`` is
    not `None`).
    """

    detectorKey = "ccd"
    """Key to use in a gen2 dataId to refer to a detector."""

    exposureKey = "visit"
    """Key to use in a gen2 dataId to refer to a visit or exposure."""

    calibFilterType = "physical_filter"
    """Gen3 dimension that corresponds to Gen2 ``filter``. Should be
    physical_filter or band."""

    skymapName = None
    """Name of the Gen3 skymap."""

    skymapConfig = None
    """Path to skymap config file defining the new gen3 skymap."""

    def setUp(self):
        self.gen3root = tempfile.mkdtemp()
        self.gen2Butler = lsst.daf.persistence.Butler(root=self.gen2root, calibRoot=self.gen2calib)
        self.collections = set(type(self).collections)
        self.collections.add(self.instrumentClass.makeDefaultRawIngestRunName())
        if len(self.refcats) > 0:
            self.collections.add("refcats")
        if self.skymapName is not None:
            self.collections.add("skymaps")

        # We always write a default calibration collection
        # containing at least the camera
        self.collections.add(self.instrumentClass.makeCollectionName("calib"))

    def tearDown(self):
        shutil.rmtree(self.gen3root, ignore_errors=True)

    def _run_convert(self):
        """Convert a gen2 repo to gen3 for testing.
        """

        # Turn on logging
        log = lsst.log.Log.getLogger("convertRepo")
        log.setLevel(log.INFO)
        log.info("Converting %s to %s", self.gen2root, self.gen3root)

        convert(repo=self.gen3root,
                gen2root=self.gen2root,
                skymap_name=self.skymapName,
                skymap_config=self.skymapConfig,
                config_file=self.config,
                calibs=self.gen2calib,
                reruns=None,
                transfer="auto")

    def check_raw(self, gen3Butler, exposure, detector):
        """Check that a raw was converted correctly.

        Parameters
        ----------
        gen3Butler : `lsst.daf.butler.Butler`
            The Butler to be tested.
        exposure : `int`
            The exposure/vist identifier ``get`` from both butlers.
        detector : `int`
            The detector identifier to ``get`` from both butlers.
        """
        dataIdGen2 = {self.detectorKey: detector, self.exposureKey: exposure}
        try:
            gen2Exposure = self.gen2Butler.get("raw", dataId=dataIdGen2)
        except lsst.daf.persistence.butlerExceptions.NoResults:
            # ignore datasets that don't actually exist in the gen2 butler.
            return
        dataIdGen3 = dict(detector=detector, exposure=exposure, instrument=self.instrumentName)
        gen3Exposure = gen3Butler.get("raw", dataId=dataIdGen3)
        # Check that we got an Exposure, but not what type; there is
        # inconsistency between different obs packages.
        self.assertIsInstance(gen3Exposure, lsst.afw.image.Exposure)
        self.assertEqual(gen3Exposure.getInfo().getDetector().getId(), detector)
        self.assertMaskedImagesEqual(gen2Exposure.maskedImage, gen3Exposure.maskedImage)

    def check_calibs(self, calibName, calibIds, gen3Butler):
        """Test that we can get converted bias/dark/flat from the gen3 repo.

        Note: because there is no clear way to get calibrations from a gen2
        repo, we just test that the thing we got is an ExposureF here, and
        assume that formatter testing is handled properly elsewhere.

        Parameters
        ----------
        calibName : `str`
            The name of the calibration to attempt to get ("bias", "flat").
        calibIds : `list` of `dict`
            The list of calibration dataIds to get.
        gen3Butler : `lsst.daf.butler.Butler`
            The Butler to use to get the data.
        """
        for dataId in calibIds:
            with self.subTest(dtype=calibName, dataId=dataId):
                datasets = list(gen3Butler.registry.queryDatasets(calibName, collections=..., dataId=dataId))
                gen3Exposure = gen3Butler.getDirect(datasets[0])
                self.assertIsInstance(gen3Exposure, lsst.afw.image.ExposureF)

    def check_defects(self, gen3Butler, detectors):
        """Test that we can get converted defects from the gen3 repo.

        Parameters
        ----------
        gen3Butler : `lsst.daf.butler.Butler`
            The Butler to be tested.
        detectors : `list` of `int`
            The detector identifiers to ``get`` from the gen3 butler.
        """
        for detector in detectors:
            dataId = dict(detector=detector, instrument=self.instrumentName)
            # Fill out the missing parts of the dataId, as we don't a-priori
            # know e.g. the "calibration_label". Use the first element of the
            # result because we only need to check one.
            with self.subTest(dtype="defects", dataId=dataId):
                datasets = list(gen3Butler.registry.queryDatasets("defects", collections=..., dataId=dataId))
                if datasets:
                    gen3Defects = gen3Butler.getDirect(datasets[0])
                    self.assertIsInstance(gen3Defects, lsst.meas.algorithms.Defects)

    def check_refcat(self, gen3Butler):
        """Test that each expected refcat is in the gen3 repo.

        Parameters
        ----------
        gen3Butler : `lsst.daf.butler.Butler`
            The Butler to be tested.
        """
        if len(self.refcats) > 0:
            for refcat in self.refcats:
                query = gen3Butler.registry.queryDatasets(refcat, collections=["refcats"])
                self.assertGreater(len(list(query)), 0,
                                   msg=f"refcat={refcat} has no entries in collection 'refcats'.")

    def check_collections(self, gen3Butler):
        """Test that the correct set of collections is in the gen3 repo.

        Parameters
        ----------
        gen3Butler : `lsst.daf.butler.Butler`
            The Butler to be tested.
        """
        self.assertEqual(set(gen3Butler.registry.queryCollections()), self.collections,
                         f"Compare with expected collections ({self.collections})")

    def test_convert(self):
        """Test that all data are converted correctly.
        """
        self._run_convert()
        gen3Butler = lsst.daf.butler.Butler(self.gen3root,
                                            collections=self.instrumentClass.makeDefaultRawIngestRunName())
        self.check_collections(gen3Butler)

        # check every raw detector that the gen2 butler knows about
        detectors = self.gen2Butler.queryMetadata("raw", self.detectorKey)
        exposures = self.gen2Butler.queryMetadata("raw", self.exposureKey)
        for exposure, detector in itertools.product(exposures, detectors):
            with self.subTest(mode="raw", exposure=exposure, detector=detector):
                self.check_raw(gen3Butler, exposure, detector)

        self.check_refcat(gen3Butler)
        self.check_defects(gen3Butler, detectors)
        self.check_calibs(self.biasName, self.biases, gen3Butler)
        self.check_calibs(self.flatName, self.flats, gen3Butler)
        self.check_calibs(self.darkName, self.darks, gen3Butler)


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
