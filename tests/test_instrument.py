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

import unittest

from lsst.obs.base import Instrument, FilterDefinitionCollection
from lsst.obs.base.gen2to3 import TranslatorFactory
from lsst.obs.base.instrument_tests import InstrumentTests, InstrumentTestData
from lsst.daf.butler.core.utils import getFullTypeName

"""Tests of the Instrument class.
"""


class DummyCam(Instrument):

    filterDefinitions = FilterDefinitionCollection()

    @classmethod
    def getName(cls):
        return "DummyCam"

    def getCamera(self):
        return None

    def register(self, registry):
        """Insert Instrument, physical_filter, and detector entries into a
        `Registry`.
        """
        dataId = {"instrument": self.getName(), "class_name": getFullTypeName(DummyCam)}
        registry.insertDimensionData("instrument", dataId)
        for f in ("dummy_g", "dummy_u"):
            registry.insertDimensionData("physical_filter",
                                         dict(dataId, physical_filter=f, band=f[-1]))
        for d in (1, 2):
            registry.insertDimensionData("detector",
                                         dict(dataId, id=d, full_name=str(d)))

    def getRawFormatter(self, dataId):
        # Docstring inherited fromt Instrument.getRawFormatter.
        return None

    def writeCuratedCalibrations(self, butler):
        pass

    def applyConfigOverrides(self, name, config):
        pass

    def makeDataIdTranslatorFactory(self) -> TranslatorFactory:
        return TranslatorFactory()


class InstrumentTestCase(InstrumentTests, unittest.TestCase):
    """Test for Instrument.
    """

    instrument = DummyCam()

    data = InstrumentTestData(name="DummyCam",
                              nDetectors=2,
                              firstDetectorName="1",
                              physical_filters={"dummy_g", "dummy_u"})

    def test_getCamera(self):
        """No camera defined in DummyCam"""
        return


if __name__ == "__main__":
    unittest.main()
