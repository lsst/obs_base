#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2015 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import unittest
import lsst.utils.tests as utilsTests
from lsst.daf.butlerUtils.fsScanner import FsScanner


class FsScannerTestCase(unittest.TestCase):

    def test1(self):
        template = '%(visit)d%(state)1s.fits.fz[%(extension)d]'
        scanner = FsScanner(template)
        res = scanner.processPath('tests/testFsScanner')
        self.assertEqual(res, {'1038843o.fits.fz': {'state': 'o', 'visit': 1038843}})

    def test2(self):
        template = 'raw_v%(visit)d_f%(filter)1s.fits.gz'
        scanner = FsScanner(template)
        res = scanner.processPath('tests/testFsScanner')
        self.assertEqual(res, {'raw_v1_fg.fits.gz': {'visit': 1, 'filter': 'g'}})

def suite():
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(FsScannerTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)
