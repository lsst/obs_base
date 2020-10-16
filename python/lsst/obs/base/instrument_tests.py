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

from lsst.obs.base import Instrument
from lsst.obs.base.gen2to3 import TranslatorFactory
from lsst.daf.butler import Registry
from lsst.daf.butler import ButlerConfig


@dataclasses.dataclass
class InstrumentTestData:
    """Values to test against in sublcasses of `InstrumentTests`.
    """

    name: str
    """The name of the Camera this instrument describes."""

    nDetectors: int
    """The number of detectors in the Camera."""

    firstDetectorName: str
    """The name of the first detector in the Camera."""

    physical_filters: {str}
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
        registry = Registry.fromConfig(ButlerConfig())
        # check that the registry starts out empty
        self.assertFalse(registry.queryDataIds(["instrument"]).toSequence())
        self.assertFalse(registry.queryDataIds(["detector"]).toSequence())
        self.assertFalse(registry.queryDataIds(["physical_filter"]).toSequence())

        # register the instrument and check that certain dimensions appear
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

        # Check that the instrument class can be retrieved
        registeredInstrument = Instrument.fromName(self.instrument.getName(), registry)
        self.assertEqual(type(registeredInstrument), type(self.instrument))

        # Check that re-registration is not an error.
        self.instrument.register(registry)

    def testMakeTranslatorFactory(self):
        factory = self.instrument.makeDataIdTranslatorFactory()
        self.assertIsInstance(factory, TranslatorFactory)
        str(factory)  # just make sure this doesn't raise.
