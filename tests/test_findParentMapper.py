#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2016 LSST Corporation.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import os
import shutil
import unittest
import tempfile

import lsst.utils.tests
import lsst.daf.persistence as dafPersist

ROOT = os.path.abspath(os.path.dirname(__file__))


class TestFindParentMapperV1Butler(unittest.TestCase):
    """Test that in an 'old' Butler repository the mapper can be found in a parent repo via the _parent
    symlink.
    """

    def setUp(self):
        self.testDir = tempfile.mkdtemp(dir=ROOT, prefix='TestFindParentMapperV1Butler-')
        self.parentRepoDir = os.path.join(self.testDir, 'parentRepo')
        os.makedirs(self.parentRepoDir)

        # note: do not put a _mapper file in the parent here, it will be added later when needed.

        self.childRepoDir = os.path.join(self.testDir, 'childRepo')
        os.makedirs(self.childRepoDir)
        os.symlink(self.parentRepoDir, os.path.join(self.childRepoDir, '_parent'))

    def tearDown(self):
        if os.path.exists(self.testDir):
            shutil.rmtree(self.testDir)

    def test(self):
        # put the mapper file in the parent repo, and see if the child repo can find it via the _parent
        # symlink.
        with open(os.path.join(self.parentRepoDir, '_mapper'), 'w') as f:
            f.write('lsst.obs.base.test.CompositeMapper')

        # this should not raise an error, no error indicates that the mapper can not be found.
        dafPersist.Butler(root=self.childRepoDir)

    def testNoMapper(self):
        # Since in this case the _mapper file does not exist in the parent repo, this should raise an error
        # indicating that the mapper can not be found.
        with self.assertRaises(RuntimeError):
            dafPersist.Butler(self.childRepoDir, policyDir=self.testDir)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
