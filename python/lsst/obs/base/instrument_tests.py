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

"""Helpers for writing tests against subclassses of Instrument.

These are not tests themselves, but can be subclassed (plus unittest.TestCase)
to get a functional test of an Instrument.
"""

import abc
import dataclasses
import pkg_resources
from typing import Set

from lsst.obs.base import Instrument, FilterDefinitionCollection, FilterDefinition
from lsst.obs.base.gen2to3 import TranslatorFactory
from lsst.obs.base.yamlCamera import makeCamera
from lsst.daf.butler import Registry
from lsst.daf.butler import RegistryConfig
from lsst.daf.butler.core.utils import getFullTypeName
from lsst.daf.butler.formatters.yaml import YamlFormatter

from .utils import createInitialSkyWcsFromBoresight

DUMMY_FILTER_DEFINITIONS = FilterDefinitionCollection(
    FilterDefinition(physical_filter="dummy_u", band="u", lambdaEff=0),
    FilterDefinition(physical_filter="dummy_g", band="g", lambdaEff=0),
)


class DummyCamYamlWcsFormatter(YamlFormatter):
    """Specialist formatter for tests that can make a WCS."""

    @classmethod
    def makeRawSkyWcsFromBoresight(cls, boresight, orientation, detector):
        """Class method to make a raw sky WCS from boresight and detector.

        This uses the API expected by define-visits. A working example
        can be found in `FitsRawFormatterBase`.

        Notes
        -----
        This makes no attempt to create a proper WCS from geometry.
        """
        return createInitialSkyWcsFromBoresight(boresight, orientation, detector, flipX=False)


class DummyCam(Instrument):

    filterDefinitions = DUMMY_FILTER_DEFINITIONS

    @classmethod
    def getName(cls):
        return "DummyCam"

    def getCamera(self):
        # Return something that can be indexed by detector number
        # but also has to support getIdIter.
        filename = pkg_resources.resource_filename("lsst.obs.base", "test/dummycam.yaml")
        return makeCamera(filename)

    def register(self, registry):
        """Insert Instrument, physical_filter, and detector entries into a
        `Registry`.
        """
        detector_max = 2
        dataId = {"instrument": self.getName(), "class_name": getFullTypeName(DummyCam),
                  "detector_max": detector_max}
        with registry.transaction():
            registry.syncDimensionData("instrument", dataId)
            self._registerFilters(registry)
            for d in range(detector_max):
                registry.syncDimensionData("detector",
                                           dict(dataId, id=d, full_name=f"RXX_S0{d}"))

    def getRawFormatter(self, dataId):
        # Docstring inherited fromt Instrument.getRawFormatter.
        return DummyCamYamlWcsFormatter

    def writeCuratedCalibrations(self, butler):
        pass

    def applyConfigOverrides(self, name, config):
        pass

    def makeDataIdTranslatorFactory(self) -> TranslatorFactory:
        return TranslatorFactory()


@dataclasses.dataclass
class InstrumentTestData:
    """Values to test against in subclasses of `InstrumentTests`.
    """

    name: str
    """The name of the Camera this instrument describes."""

    nDetectors: int
    """The number of detectors in the Camera."""

    firstDetectorName: str
    """The name of the first detector in the Camera."""

    physical_filters: Set[str]
    """A subset of the physical filters should be registered."""


class InstrumentTests(metaclass=abc.ABCMeta):
    """Tests of sublcasses of Instrument.

    TestCase subclasses must derive from this, then `TestCase`, and override
    ``data`` and ``instrument``.
    """

    data = None
    """`InstrumentTestData` containing the values to test against."""

    instrument = None
    """The `~lsst.obs.base.Instrument` to be tested."""

    def test_name(self):
        self.assertEqual(self.instrument.getName(), self.data.name)

    def test_getCamera(self):
        """Test that getCamera() returns a reasonable Camera definition.
        """
        camera = self.instrument.getCamera()
        self.assertEqual(camera.getName(), self.instrument.getName())
        self.assertEqual(len(camera), self.data.nDetectors)
        self.assertEqual(next(iter(camera)).getName(), self.data.firstDetectorName)

    def test_register(self):
        """Test that register() sets appropriate Dimensions.
        """
        registryConfig = RegistryConfig()
        registryConfig["db"] = "sqlite://"
        registry = Registry.createFromConfig(registryConfig)
        # Check that the registry starts out empty.
        self.assertFalse(registry.queryDataIds(["instrument"]).toSequence())
        self.assertFalse(registry.queryDataIds(["detector"]).toSequence())
        self.assertFalse(registry.queryDataIds(["physical_filter"]).toSequence())

        # Register the instrument and check that certain dimensions appear.
        self.instrument.register(registry)
        instrumentDataIds = registry.queryDataIds(["instrument"]).toSequence()
        self.assertEqual(len(instrumentDataIds), 1)
        instrumentNames = {dataId["instrument"] for dataId in instrumentDataIds}
        self.assertEqual(instrumentNames, {self.data.name})
        detectorDataIds = registry.queryDataIds(["detector"]).expanded().toSequence()
        self.assertEqual(len(detectorDataIds), self.data.nDetectors)
        detectorNames = {dataId.records["detector"].full_name for dataId in detectorDataIds}
        self.assertIn(self.data.firstDetectorName, detectorNames)
        physicalFilterDataIds = registry.queryDataIds(["physical_filter"]).toSequence()
        filterNames = {dataId['physical_filter'] for dataId in physicalFilterDataIds}
        self.assertGreaterEqual(filterNames, self.data.physical_filters)

        # Check that the instrument class can be retrieved.
        registeredInstrument = Instrument.fromName(self.instrument.getName(), registry)
        self.assertEqual(type(registeredInstrument), type(self.instrument))

        # Check that re-registration is not an error.
        self.instrument.register(registry)

    def testMakeTranslatorFactory(self):
        factory = self.instrument.makeDataIdTranslatorFactory()
        self.assertIsInstance(factory, TranslatorFactory)
        str(factory)  # Just make sure this doesn't raise.
