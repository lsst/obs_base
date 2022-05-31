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

import unittest

import lsst.afw.image
import lsst.pex.exceptions
import lsst.utils.tests
from lsst.obs.base import FilterDefinition, FilterDefinitionCollection


class TestFilterDefinitionCollection(lsst.utils.tests.TestCase):
    def setUp(self):
        self.filters1 = FilterDefinitionCollection(
            FilterDefinition(physical_filter="abc"),
            FilterDefinition(physical_filter="def", band="d", doc="This is a test filter."),
        )
        self.filters2 = FilterDefinitionCollection(
            FilterDefinition(physical_filter="abc"),
            FilterDefinition(physical_filter="def", band="dd"),
        )

    def test_findAll(self):
        self.assertEqual(set(self.filters1.findAll("r")), set())
        matches = self.filters1.findAll("abc")
        self.assertEqual(len(matches), 1)
        match = list(matches)[0]
        self.assertEqual(match.physical_filter, "abc")

    def test_physical_to_band(self):
        """Test that the physical_to_band dict returns expected values."""
        self.assertIsNone(self.filters1.physical_to_band["abc"])
        self.assertEqual(self.filters1.physical_to_band["def"], "d")
        self.assertIsNone(self.filters2.physical_to_band["abc"])
        self.assertEqual(self.filters2.physical_to_band["def"], "dd")


class TestFilterDefinition(lsst.utils.tests.TestCase):
    def setUp(self):
        self.filter_g = FilterDefinition(band="g", physical_filter="HSC-G", alias={"ABCDEFG"})
        self.filter_g2 = FilterDefinition(band="g", physical_filter="HSC-G2", afw_name="g2", alias={"HIJK"})

        self.physical_only = FilterDefinition(physical_filter="physical")
        self.afw_name = FilterDefinition(physical_filter="afw_name", afw_name="afw only")
        self.abstract = FilterDefinition(physical_filter="abstract", band="abstract only")

    def test_physical_only(self):
        """physical_filter is the only name this filter has."""
        self.assertEqual(
            self.physical_only.makeFilterLabel(), lsst.afw.image.FilterLabel(physical="physical")
        )

    def test_afw_name(self):
        """afw_name is the Filter name, physical_filter is an alias."""
        self.assertEqual(self.afw_name.makeFilterLabel(), lsst.afw.image.FilterLabel(physical="afw_name"))

    def test_abstract_only(self):
        """band is the Filter name, physical_filter is an alias."""
        self.assertEqual(
            self.abstract.makeFilterLabel(),
            lsst.afw.image.FilterLabel(band="abstract only", physical="abstract"),
        )


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
