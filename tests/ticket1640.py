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

import os
import sys

import eups

if not eups.Eups().isSetup("obs_lsstSim") or \
        not eups.Eups().isSetup("afwdata"):
    print "obs_lsstSim or afwdata not setup; not running test."
    sys.exit(0)

import lsst.daf.persistence as dafPersist
from lsst.obs.lsstSim import LsstSimMapper

class Ticket1640TestCase(unittest.TestCase):

    def setUp(self):
        self.dataRoot = os.path.join(eups.productDir("afwdata"), "ImSim")

    def testTicket1640(self):
        butler = dafPersist.ButlerFactory(
                mapper=LsstSimMapper(root=self.dataRoot)).create()
        kwargs = {'visit': 85408556, 'sensor': '1,1', 'raft': '2,3'}
        self.assertEqual(
                butler.queryMetadata("raw", "visit", **kwargs)[0:1],
                [85408556])
        self.assertEqual(
                butler.queryMetadata("raw", "visit", visit=kwargs["visit"],
                    raft=kwargs["raft"], sensor=kwargs["sensor"])[0:1],
                [85408556])

def suite():
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(Ticket1640TestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)
