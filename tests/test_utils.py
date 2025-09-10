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
from datetime import datetime

import astropy.table
import numpy as np

import lsst.geom as geom
import lsst.obs.base as obsBase
import lsst.utils.tests
from lsst.daf.base import PropertyList
from lsst.obs.base.utils import TableVStack, _store_str_header
from lsst.pipe.base._dataset_handle import InMemoryDatasetHandle


class BboxFromIrafTestCase(lsst.utils.tests.TestCase):
    """Demonstrate that we can correctly parse IRAF-style BBOXes."""

    def testValid(self):
        test_data = {
            "[1:1084,1:1024]": geom.BoxI(geom.PointI(0, 0), geom.PointI(1083, 1023)),
            "[0:0,0:0]": geom.BoxI(geom.PointI(-1, -1), geom.PointI(-1, -1)),
        }
        for val, result in test_data.items():
            self.assertEqual(obsBase.bboxFromIraf(val), result)

    def testInvalid(self):
        test_data = {
            "1:1084,1:1024": RuntimeError,
            "(1:1084,1:1024)": RuntimeError,
            ("1:1084", "1:1024"): TypeError,
        }
        for val, err in test_data.items():
            self.assertRaises(err, obsBase.bboxFromIraf, val)


class TestProvenanceAdd(unittest.TestCase):
    """Tests relating to provenance infrastructure."""

    def test_truncation(self):
        """Test that long headers can be truncated."""
        pl = PropertyList()

        _store_str_header(pl, "LSST BUTLER RUN", "short")
        self.assertEqual(pl["LSST BUTLER RUN"], "short")

        _store_str_header(pl, "LSST BUTLER RUN", "short", allow_long_headers=False)
        self.assertEqual(pl["LSST BUTLER RUN"], "short")

        long = "a123456789b123456789c123456789d123456789e123456789f123456789"
        _store_str_header(pl, "LSST BUTLER INPUT 0 RUN", long)
        self.assertEqual(pl["LSST BUTLER INPUT 0 RUN"], long)

        _store_str_header(pl, "LSST BUTLER INPUT 1 RUN", long, allow_long_headers=False)
        self.assertEqual(pl["LSST BUTLER INPUT 1 RUN"], "a123456789b123456789...e123456789f123456789")

        key = "LSSTX BUTLER VERY LONG KEYWORD THAT IS AT THE LIMIT OF ALL"
        _store_str_header(pl, key, "abc", allow_long_headers=False)
        self.assertEqual(pl[key], "abc")
        _store_str_header(pl, key, "abcdefghi", allow_long_headers=False)
        self.assertEqual(pl[key], "ab...ghi")

        with self.assertRaises(ValueError):
            _store_str_header(pl, key + "0", "abcdef", allow_long_headers=False)


class TestCalibDates(unittest.TestCase):
    """Tests relating to curated calibration dates."""

    def test_calib_dates(self):
        """Test file root creation and parsing."""
        for valid_start, file_root, standard_iso in (
            ("2025-04-30T12:25", "20250430T122500", "2025-04-30T12:25:00"),
            ("2025-04-30 12:25:50", "20250430T122550", "2025-04-30T12:25:50"),
            ("2025-04-30 12:25:50.123", "20250430T122550", "2025-04-30T12:25:50"),
            ("1970-01-01", "19700101T000000", "1970-01-01T00:00:00"),
            ("20100501T1223", "20100501T122300", "2010-05-01T12:23:00"),
        ):
            derived = obsBase.iso_date_to_curated_calib_file_root(valid_start)
            self.assertEqual(derived, file_root)
            self.assertEqual(datetime.fromisoformat(derived).isoformat(), standard_iso)
            # We know that fromisoformat is used internally by read_one_calib
            # so have explicit test here to make sure the date formats we
            # expect will work.
            self.assertEqual(datetime.fromisoformat(valid_start).isoformat(timespec="seconds"), standard_iso)


class TableVStackTestCase(unittest.TestCase):
    """Unit tests for lsst.pipe.base.utils.TableVStack."""

    def setUp(self):
        self.len_table1 = 5
        self.len_table2 = 6
        self.table1 = astropy.table.Table(
            {
                "x": np.arange(self.len_table1),
                "y": np.ma.masked_array(
                    np.arange(1.0, 1.0 + self.len_table1), mask=np.arange(self.len_table1) < 3
                ),
            }
        )
        self.table2 = astropy.table.Table(
            {
                "x": np.arange(self.len_table2),
                "y": np.ma.masked_array(
                    np.arange(-2.0, -2.0 + self.len_table2),
                    mask=(np.arange(self.len_table2) % 2) == 0,
                ),
            }
        )
        self.extra_values1 = {"z": -self.table1["x"]}
        self.extra_values2 = {"z": -self.table2["x"]}

    def test_vstack(self):
        handles = (
            InMemoryDatasetHandle(self.table1, storageClass="ArrowAstropy"),
            InMemoryDatasetHandle(self.table2, storageClass="ArrowAstropy"),
        )
        u_stack = TableVStack.vstack_handles(
            handles=handles,
            extra_values={
                idx: extra_values for idx, extra_values in enumerate((self.extra_values1, self.extra_values2))
            },
        )
        ap_stack = astropy.table.vstack((self.table1, self.table2))
        self.assertEqual(len(u_stack), self.len_table1 + self.len_table2)
        self.assertEqual(len(u_stack), len(ap_stack))

        ap_stack["z"] = np.concatenate((self.extra_values1["z"], self.extra_values2["z"]))
        self.assertEqual(u_stack.colnames, ap_stack.colnames)
        for colname in u_stack.colnames:
            col_u = u_stack[colname]
            col_ap = ap_stack[colname]
            if (mask := getattr(col_u, "mask", None)) is not None:
                np.testing.assert_array_equal(mask, col_ap.mask)
                col_u = col_u.data
                col_ap = col_ap.data
            np.testing.assert_array_equal(col_u, col_ap)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    """Test for file leaks."""


def setup_module(module):
    """Initialize pytest."""
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
