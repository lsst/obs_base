# This file is part of daf_butler.
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

"""Unit tests for the daf_butler CliLog utility. Code is implemented in
daf_butler but some only runs if lsst.log.Log can be imported so these parts of
it can't be tested there because daf_butler does not directly depend on
lsst.log, and only uses it if it has been setup by another package."""

import unittest

import lsst.log
from lsst.daf.butler.cli.cliLog import CliLog
from lsst.daf.butler.tests import CliLogTestBase


class CliLogTestCase(CliLogTestBase, unittest.TestCase):
    """Test log initialization, reset, and setting log levels on python
    `logging` and `lsst.log`.

    This test also runs in daf_butler but will not test `lsst.log` in CI
    because daf_butler does not directly depend on that package."""

    pass


class ConvertLsstLogLevelTestCase(unittest.TestCase):
    def test_convertToLsstLogLevel(self):
        """Test that the log levels accepted by the log_level_option are
        translated to lsst.log levels correctly."""
        self.assertEqual(lsst.log.Log.FATAL, CliLog._getLsstLogLevel("CRITICAL"))
        self.assertEqual(lsst.log.ERROR, CliLog._getLsstLogLevel("ERROR"))
        self.assertEqual(lsst.log.WARN, CliLog._getLsstLogLevel("WARNING"))
        self.assertEqual(lsst.log.INFO, CliLog._getLsstLogLevel("INFO"))
        self.assertEqual(lsst.log.DEBUG, CliLog._getLsstLogLevel("DEBUG"))


if __name__ == "__main__":
    unittest.main()
