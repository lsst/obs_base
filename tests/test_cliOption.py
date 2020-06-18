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

"""Unit tests for the butler instrument CLI option.
"""

import unittest

from lsst.daf.butler.tests import (OptHelpTest,
                                   OptMultipleTest,
                                   OptRequiredTest)
from lsst.obs.base.cli.opt import instrument_parameter


class InstrumentArgumentTestCase(OptHelpTest,
                                 OptMultipleTest,
                                 OptRequiredTest,
                                 unittest.TestCase):
    optionClass = instrument_parameter
    optionName = "instrument"
    shortOptionName = "i"
    isArgument = True
    isParameter = True


class InstrumentOptionTestCase(OptHelpTest,
                               OptMultipleTest,
                               OptRequiredTest,
                               unittest.TestCase):
    optionClass = instrument_parameter
    optionName = "instrument"
    shortOptionName = "i"
    isParameter = True


if __name__ == "__main__":
    unittest.main()
