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

from __future__ import annotations

__all__ = [
    "DummyCam",
    "InstrumentTestData",
    "InstrumentTests",
    "DummyCamYamlWcsFormatter",
    "CuratedCalibration",
]

import abc
import dataclasses
from collections.abc import Sequence
from functools import lru_cache
from typing import TYPE_CHECKING, Any, ClassVar

from lsst.daf.butler import CollectionType, DatasetType, RegistryConfig
from lsst.daf.butler.formatters.yaml import YamlFormatter
from lsst.daf.butler.registry.sql_registry import SqlRegistry
from lsst.obs.base import FilterDefinition, FilterDefinitionCollection, Instrument
from lsst.obs.base.yamlCamera import makeCamera
from lsst.resources import ResourcePath
from lsst.utils.introspection import get_full_type_name
from pydantic import BaseModel

from .utils import createInitialSkyWcsFromBoresight

if TYPE_CHECKING:
    from lsst.daf.butler import Butler

DUMMY_FILTER_DEFINITIONS = FilterDefinitionCollection(
    FilterDefinition(physical_filter="dummy_u", band="u"),
    FilterDefinition(physical_filter="dummy_g", band="g"),
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


class CuratedCalibration(BaseModel):
    """Class that implements minimal read/write interface needed to support
    curated calibration ingest.
    """

    metadata: dict[str, Any]
    values: list[int]

    @classmethod
    def readText(cls, path: str) -> CuratedCalibration:
        with open(path) as f:
            data = f.read()
        return cls.model_validate_json(data)

    def writeText(self, path: str) -> None:
        with open(path, "w") as f:
            print(self.model_dump_json(), file=f)

    def getMetadata(self) -> dict[str, Any]:
        return self.metadata


class DummyCam(Instrument):
    """Instrument class used for testing."""

    filterDefinitions = DUMMY_FILTER_DEFINITIONS
    additionalCuratedDatasetTypes = frozenset(["testCalib"])
    policyName = "dummycam"
    dataPackageDir: str | None = ""
    raw_definition = (
        "raw_dict",
        ("instrument", "detector", "exposure"),
        "StructuredDataDict",
    )

    @classmethod
    def getName(cls):
        return "DummyCam"

    @classmethod
    @lru_cache  # For mypy
    def getObsDataPackageDir(cls) -> str | None:
        return cls.dataPackageDir

    def getCamera(self):
        # Return something that can be indexed by detector number
        # but also has to support getIdIter.
        with ResourcePath("resource://lsst.obs.base/test/dummycam.yaml").as_local() as local_file:
            return makeCamera(local_file.ospath)

    def register(self, registry, update=False):
        """Insert Instrument, physical_filter, and detector entries into a
        `Registry`.
        """
        detector_max = 2
        dataId = {
            "instrument": self.getName(),
            "class_name": get_full_type_name(DummyCam),
            "detector_max": detector_max,
            "visit_max": 1_000_000,
            "exposure_max": 1_000_000,
        }

        with registry.transaction():
            registry.syncDimensionData("instrument", dataId, update=update)
            self._registerFilters(registry, update=update)
            for d in range(detector_max):
                registry.syncDimensionData(
                    "detector",
                    dict(
                        instrument=self.getName(),
                        id=d,
                        full_name=f"RXX_S0{d}",
                    ),
                    update=update,
                )

    def getRawFormatter(self, dataId):
        # Docstring inherited fromt Instrument.getRawFormatter.
        return DummyCamYamlWcsFormatter

    def applyConfigOverrides(self, name, config):
        pass

    def writeAdditionalCuratedCalibrations(
        self, butler: Butler, collection: str | None = None, labels: Sequence[str] = ()
    ) -> None:
        # We want to test the standard curated calibration ingest
        # but we do not have a standard class to use. There is no way
        # at the moment to inject a new class into the standard list
        # that is a package constant, so instead use this "Additional"
        # method but call the standard curated calibration code.
        if collection is None:
            collection = self.makeCalibrationCollectionName(*labels)
        butler.registry.registerCollection(collection, type=CollectionType.CALIBRATION)

        datasetType = DatasetType(
            "testCalib",
            universe=butler.dimensions,
            isCalibration=True,
            dimensions=("instrument", "detector"),
            storageClass="CuratedCalibration",
        )
        runs: set[str] = set()
        self._writeSpecificCuratedCalibrationDatasets(
            butler, datasetType, collection, runs=runs, labels=labels
        )


@dataclasses.dataclass
class InstrumentTestData:
    """Values to test against in subclasses of `InstrumentTests`."""

    name: str
    """The name of the Camera this instrument describes."""

    nDetectors: int
    """The number of detectors in the Camera."""

    firstDetectorName: str
    """The name of the first detector in the Camera."""

    physical_filters: set[str]
    """A subset of the physical filters should be registered."""


class InstrumentTests(metaclass=abc.ABCMeta):
    """Tests of sublcasses of Instrument.

    TestCase subclasses must derive from this, then `TestCase`, and override
    ``data`` and ``instrument``.
    """

    data: ClassVar[InstrumentTestData | None] = None
    """`InstrumentTestData` containing the values to test against."""

    instrument: ClassVar[Instrument | None] = None
    """The `~lsst.obs.base.Instrument` to be tested."""

    def test_name(self):
        self.assertEqual(self.instrument.getName(), self.data.name)

    def test_getCamera(self):
        """Test that getCamera() returns a reasonable Camera definition."""
        camera = self.instrument.getCamera()
        self.assertEqual(camera.getName(), self.instrument.getName())
        self.assertEqual(len(camera), self.data.nDetectors)
        self.assertEqual(next(iter(camera)).getName(), self.data.firstDetectorName)

    def test_register(self):
        """Test that register() sets appropriate Dimensions."""
        registryConfig = RegistryConfig()
        registryConfig["db"] = "sqlite://"
        registry = SqlRegistry.createFromConfig(registryConfig)
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
        filterNames = {dataId["physical_filter"] for dataId in physicalFilterDataIds}
        self.assertGreaterEqual(filterNames, self.data.physical_filters)

        # Check that the instrument class can be retrieved.
        registeredInstrument = Instrument.fromName(self.instrument.getName(), registry)
        self.assertEqual(type(registeredInstrument), type(self.instrument))

        # Check that re-registration is not an error.
        self.instrument.register(registry)
