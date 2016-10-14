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

import filecmp
import gc
import os
import shutil
import unittest
import lsst.utils.tests
from lsst.utils import getPackageDir

import lsst.daf.persistence as dafPersist
import lsst.daf.persistence.test as dpTest
import lsst.utils.tests
from lsst.utils import getPackageDir


class TestCompositeTestCase(unittest.TestCase):
    """A test case for composite object i/o."""

    def setUp(self):
        packageDir = getPackageDir('obs_base')
        self.testData = os.path.join(packageDir, 'tests', 'composite')
        self.firstRepoPath = os.path.join(self.testData, 'repo1')
        self.objA = dpTest.TestObject("abc")
        self.objB = dpTest.TestObject("def")
        policy = dafPersist.Policy({'camera': 'lsst.afw.cameraGeom.Camera',
                                    'datasets': {
                                        'basicObject1': {
                                            'python': 'lsst.daf.persistence.test.TestObject',
                                            'template': 'basic/id%(id)s.pickle',
                                            'storage': 'PickleStorage'},
                                        'basicObject2': {
                                            'python': 'lsst.daf.persistence.test.TestObject',
                                            'template': 'basic/name%(name)s.pickle',
                                            'storage': 'PickleStorage'},
                                        'basicPair': {
                                            'python': 'lsst.daf.persistence.test.TestObjectPair',
                                            'composite': {
                                                'a': {'datasetType': 'basicObject1'},
                                                'b': {'datasetType': 'basicObject2'}
                                            },
                                            'assembler': 'lsst.daf.persistence.test.TestObjectPair.assembler',
                                            'disassembler': 'lsst.daf.persistence.test.TestObjectPair.disassembler'

                                        }
                                    }})

        # We need a way to put policy into a repo. Butler does not support it yet. This is a cheat.
        # The ticket to fix it is DM-7777
        if not os.path.exists(self.firstRepoPath):
            os.makedirs(self.firstRepoPath)
        policy.dumpToFile(os.path.join(self.testData, 'policy.yaml'))
        del policy

        repoArgs = dafPersist.RepositoryArgs(root=self.firstRepoPath,
                                             mapper='lsst.obs.base.test.CompositeMapper',
                                             mapperArgs={'policyDir': self.testData})
        butler = dafPersist.Butler(outputs=repoArgs)
        butler.put(self.objA, 'basicObject1', dataId={'id': 'foo'})
        butler.put(self.objB, 'basicObject2', dataId={'name': 'bar'})
        del butler
        del repoArgs



    def tearDown(self):
        if os.path.exists(self.testData):
            shutil.rmtree(self.testData)

    def testType3GetAndPut(self):
        """1. Verify that a composite can be loaded and that its components are the same as when the type1
        components are loaded individually (verifies correct lookup in this case).
        2. Verify that when the individual components are put and when the composite is put (which
        disassembles into individual components) that the objects that are written are the same.
        3. Verify object caching & reuse - that 2 components of the same datasetType and the same used dataId
        are loaded from the same dataset the object is shared instead of loaded 2 times.
        4. Verify release & garbage collection of the cached objects when they are no longer used.
        """
        secondRepoPath = os.path.join(self.testData, 'repo2')
        # child repositories do not look up in-repo policies. We need to fix that.
        # The ticket to fix this is DM-7778
        repoArgs = dafPersist.RepositoryArgs(root=secondRepoPath,
                                             mapperArgs={'policyDir': self.testData})
        butler = dafPersist.Butler(inputs=self.firstRepoPath, outputs=repoArgs)
        objABPair = butler.get('basicPair', dataId={'id': 'foo', 'name': 'bar'})

        self.assertEqual(self.objA, objABPair.objA)
        self.assertEqual(self.objB, objABPair.objB)

        # For now also test that the type 1 and type 3 components are not the same object.
        # These objects are not yet in the butler cache because they have not been gotten yet (they have only
        # only been put)
        self.assertIsNot(self.objA, objABPair.objA)
        self.assertIsNot(self.objB, objABPair.objB)

        # Now, get a type 1 copy of objA and objB, and they should be the same instance as in the composite.
        objA = butler.get('basicObject1', dataId={'id': 'foo'}, immediate=True)
        objB = butler.get('basicObject2', dataId={'name': 'bar'}, immediate=True)
        self.assertIs(objA, objABPair.objA)
        self.assertIs(objB, objABPair.objB)

        butler.put(objABPair, 'basicPair', dataId={'id': 'foo', 'name': 'bar'})
        # comparing the output files directly works so long as the storage is posix:
        self.assertTrue(filecmp.cmp(os.path.join(self.firstRepoPath, 'basic', 'idfoo.pickle'),
                                    os.path.join(secondRepoPath, 'basic', 'idfoo.pickle')))
        self.assertTrue(filecmp.cmp(os.path.join(self.firstRepoPath, 'basic', 'namebar.pickle'),
                                    os.path.join(secondRepoPath, 'basic', 'namebar.pickle')))

        # check that objA and objB are cached
        self.assertIn(objA, butler.objectCache.itervalues())
        self.assertIn(objB, butler.objectCache.itervalues())
        del objA
        del objABPair
        gc.collect()
        # check that A is not cached but B is cached by verifying B is the only object in the cache
        self.assertIs(len(butler.objectCache), 1)
        self.assertIn(objB, butler.objectCache.itervalues())
        del objB
        gc.collect()
        # check that B is not cached
        self.assertIs(len(butler.objectCache), 0)

        del butler



    def testDottedDatasetType(self):
        """Verify that components of a composite can be loaded by dotted name in the form
        DatasetType.componentName
        """
        thirdRepoPath = os.path.join(self.testData, 'repo3')
        # child repositories do not look up in-repo policies. We need to fix that.
        repoArgs = dafPersist.RepositoryArgs(root=thirdRepoPath,
                                             mapperArgs={'policyDir': self.testData})
        butler = dafPersist.Butler(inputs=self.firstRepoPath, outputs=repoArgs)
        componentObjA = butler.get('basicPair.a', dataId={'id': 'foo', 'name': 'bar'})
        componentObjB = butler.get('basicPair.b', dataId={'id': 'foo', 'name': 'bar'})
        self.assertEqual(self.objA, componentObjA)
        self.assertEqual(self.objB, componentObjB)
        butler.put(componentObjA, 'basicPair.a', dataId={'id': 'foo', 'name': 'bar'})
        butler.put(componentObjB, 'basicPair.b', dataId={'id': 'foo', 'name': 'bar'})
        # comparing the output files directly works so long as the storage is posix:
        self.assertTrue(filecmp.cmp(os.path.join(self.firstRepoPath, 'basic', 'idfoo.pickle'),
                                    os.path.join(thirdRepoPath, 'basic', 'idfoo.pickle')))
        self.assertTrue(filecmp.cmp(os.path.join(self.firstRepoPath, 'basic', 'namebar.pickle'),
                                    os.path.join(thirdRepoPath, 'basic', 'namebar.pickle')))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
