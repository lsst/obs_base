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

"""Unit tests for the daf_butler cli.Log utility that test code depenent on
lsst.log. Code is implemented in daf_butler only runs if lsst.log.Log can be
imported. It can't be tested there because daf_butler does not directly depend
on lsst.log, and only uses it if it has been setup by another package."""

import logging
import unittest

from lsst.daf.butler.cli.cmd import create
from lsst.daf.butler.cli.log import Log as butlerCliLog
from lsst.daf.butler.tests import CliCmdTestBase
from lsst.log import Log as lsstLog


class ConvertLogLevelTestCase(CliCmdTestBase,
                              unittest.TestCase):
    """Test that the log levels accepted by the log_level_option are translated
    to lsst.log levels correctly.

    """

    command = create
    defaultExpected = dict(repo=None,
                           seed_config=None,
                           standalone=False,
                           override=False,
                           outfile=None)

    def test_logInit(self):
        """Initialize logging by running a butler command. Verify a handler
        has been added"""
        rootLogger = logging.getLogger()
        self.assertEqual(len(rootLogger.handlers), 0,
                         msg="The root logger should not have any handlers when this test is started.")
        self.run_test(["--log-level", "WARNING", "create", "foo"],
                      self.makeExpected(repo="foo"))
        self.assertEqual(len(rootLogger.handlers), 1,
                         msg="After running a butler command the root loger should have been initialized "
                             "with the lsst logger handler.")
        butlerCliLog.uninitLog()
        self.assertEqual(len(rootLogger.handlers), 0,
                         msg="After uninit the lsst logger handler should have been removed from the root "
                             "logger.")


    def test_convertToLsstLogLevel(self):
        self.assertEqual(lsstLog.FATAL, butlerCliLog.getLsstLogLevel("CRITICAL"))
        self.assertEqual(lsstLog.ERROR, butlerCliLog.getLsstLogLevel("ERROR"))
        self.assertEqual(lsstLog.WARN, butlerCliLog.getLsstLogLevel("WARNING"))
        self.assertEqual(lsstLog.INFO, butlerCliLog.getLsstLogLevel("INFO"))
        self.assertEqual(lsstLog.DEBUG, butlerCliLog.getLsstLogLevel("DEBUG"))


if __name__ == "__main__":
    unittest.main()
