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

import collections
import unittest
import lsst.utils.tests as utilsTests

import lsst.daf.butlerUtils as butlerUtils


class PosixRegistryTestCase(unittest.TestCase):

    def test1(self):
        """Basic test case to verify parameter extraction from template + filename.

        :return:
        """
        testData = collections.namedtuple('testData', 'root template returnFields dataId expectedLookup')
        td = (\
              testData('tests/posixRegistry/repo01',
                       'foo-%(ccd)02d.fits',
                       ('ccd',),
                       {},
                       [(1,)]),
              testData('tests/posixRegistry/repo02',
                       'foo-%(ccd)02d-%(filter)s.fits',
                       ('ccd', 'filter'),
                       {},
                       [(1, 'g'), (1, 'h'), (2, 'g'), (2, 'h'), (3, 'i'),]),
              testData('tests/posixRegistry/repo02',
                       'foo-%(ccd)02d-%(filter)s.fits',
                       ('filter'), # intentionally no comma; it becomes a string not a tuple. This is handled.
                       {'ccd':1},
                       [('g',), ('h',),]),
              testData('tests/posixRegistry/repo02',
                       'foo-%(ccd)02d-%(filter)s.fits',
                       ('ccd',),
                       {'filter':'i'},
                       [(3,),]),
              )

        policyTables = None
        storage = 'fits'
        values = []
        for root, template, returnFields, dataId, expectedLookup in td:
            registry = butlerUtils.PosixRegistry(root)
            lookups = registry.lookup(returnFields, policyTables, dataId, template=template, storage=storage)
            lookups.sort()
            expectedLookup.sort()
            self.assertEqual(lookups, expectedLookup)

def suite():
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(PosixRegistryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)
