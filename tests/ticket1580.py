#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
from lsst.afw.image import PARENT

import os
import sys

import eups

if not eups.Eups().isSetup("obs_lsstSim") or \
        not eups.Eups().isSetup("afwdata"):
    print "obs_lsstSim or afwdata not setup; not running test."
    sys.exit(0)

import lsst.daf.persistence as dafPersist
from lsst.obs.lsstSim import LsstSimMapper

class Ticket1580TestCase(unittest.TestCase):

    def setUp(self):
        self.dataRoot = os.path.join(eups.productDir("afwdata"), "ImSim")

    def testTrimmedCalexp(self):
        #This was originally written with ticket/1580 in mind, but there is no
        #isTrimmed method any more.  We can still check that the image is the
        #trimmed size (getBBox())
        butler = dafPersist.ButlerFactory(
                mapper=LsstSimMapper(root=self.dataRoot)).create()
        kwargs = {'visit': 85408556, 'sensor': '1,1', 'raft': '2,3'}
        calexp = butler.get("calexp", **kwargs)
        self.assertEqual(calexp.getDetector().getBBox(), calexp.getBBox(PARENT))

def suite():
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(Ticket1580TestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)
