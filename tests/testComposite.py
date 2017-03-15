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
        self.tearDown()
        self.firstRepoPath = os.path.join(self.testData, 'repo1')
        self.objA = dpTest.TestObject("abc")
        self.objB = dpTest.TestObject("def")
        self.policy = dafPersist.Policy(
                                   {'camera': 'lsst.afw.cameraGeom.Camera',
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

        repoArgs = dafPersist.RepositoryArgs(root=self.firstRepoPath, policy=self.policy,
                                             mapper='lsst.obs.base.test.CompositeMapper')
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
        repoArgs = dafPersist.RepositoryArgs(root=secondRepoPath, policy=self.policy)
        butler = dafPersist.Butler(inputs=self.firstRepoPath, outputs=repoArgs)
        verificationButler = dafPersist.Butler(inputs=secondRepoPath)
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
        verObjA = verificationButler.get('basicObject1', {'id': 'foo'})
        self.assertEqual(verObjA, objABPair.objA)
        verObjB = verificationButler.get('basicObject2', {'name': 'bar'})
        self.assertEqual(verObjB, objABPair.objB)

        # check that objA and objB are cached
        self.assertIn(objA, butler.objectCache.values())
        self.assertIn(objB, butler.objectCache.values())
        del objA
        del objABPair
        gc.collect()
        # check that A is not cached but B is cached by verifying B is the only object in the cache
        self.assertIs(len(butler.objectCache), 1)
        self.assertIn(objB, butler.objectCache.values())
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
        repoArgs = dafPersist.RepositoryArgs(root=thirdRepoPath, policy=self.policy)
        butler = dafPersist.Butler(inputs=self.firstRepoPath, outputs=repoArgs)
        verificationButler = dafPersist.Butler(inputs=thirdRepoPath)
        componentObjA = butler.get('basicPair.a', dataId={'id': 'foo', 'name': 'bar'})
        componentObjB = butler.get('basicPair.b', dataId={'id': 'foo', 'name': 'bar'})
        self.assertEqual(self.objA, componentObjA)
        self.assertEqual(self.objB, componentObjB)
        butler.put(componentObjA, 'basicPair.a', dataId={'id': 'foo', 'name': 'bar'})
        butler.put(componentObjB, 'basicPair.b', dataId={'id': 'foo', 'name': 'bar'})
        verObjA = verificationButler.get('basicObject1', {'id': 'foo'})
        self.assertEqual(verObjA, componentObjA)
        verObjB = verificationButler.get('basicObject2', {'name': 'bar'})
        self.assertEqual(verObjB, componentObjB)


class TestGenericAssembler(unittest.TestCase):
    """A test case for the generic assembler feature of composite datasets."""

    def setUp(self):
        packageDir = getPackageDir('obs_base')
        self.testData = os.path.join(packageDir, 'tests', 'genericAssembler')
        self.tearDown()
        self.firstRepoPath = os.path.join(self.testData, 'repo1')
        self.secondRepoPath = os.path.join(self.testData, 'repo2')
        self.objA = dpTest.TestObject("abc")
        self.objB = dpTest.TestObject("def")
        self.policy = dafPersist.Policy(
                                   {'camera': 'lsst.afw.cameraGeom.Camera',
                                    'datasets': {
                                        'basicObject1': {
                                            'python': 'lsst.daf.persistence.test.TestObject',
                                            'template': 'basic/id%(id)s.pickle',
                                            'storage': 'PickleStorage'
                                        },
                                        'basicObject2': {
                                            'python': 'lsst.daf.persistence.test.TestObject',
                                            'template': 'basic/name%(name)s.pickle',
                                            'storage': 'PickleStorage'
                                        },
                                        'basicPair': {
                                            'python': 'lsst.daf.persistence.test.TestObjectPair',
                                            'composite': {
                                                'a': {'datasetType': 'basicObject1'},
                                                'b': {'datasetType': 'basicObject2'}
                                            },
                                            # note, no assembler or disassembler specified here, will use
                                            # setter names inferred by component name.
                                        },
                                        # "generic assembler default constructor pair"
                                        'gaDefCtorPair': { # dataset defition that uses the default ctor
                                            'python': 'lsst.daf.persistence.test.TestObjectPair',
                                            'composite': {
                                                # note that the component names are the same as the argument
                                                # names in the TestObjectPair.__init__ func.
                                                'objA': {'datasetType': 'basicObject1',
                                                         'getter': 'get_a'},
                                                'objB': {'datasetType': 'basicObject2',
                                                         'getter': 'get_b'}
                                            },
                                            # note, no assembler or disassembler specified here.
                                        },
                                        # "generic assembler default "
                                        'gaPairWithSetter': {
                                            'python': 'lsst.daf.persistence.test.TestObjectPair',
                                            'composite': {
                                                # note that the component names do not match argument names
                                                # in the TestObjectPair.__init__ func or the set functions
                                                # in the python object.
                                                'z': {'datasetType': 'basicObject1',
                                                      'setter': 'set_a',
                                                      'getter': 'get_a'
                                                },
                                                'x': {'datasetType': 'basicObject2',
                                                      'setter': 'set_b',
                                                      'getter': 'get_b'
                                                }
                                            }
                                        },
                                        # simple object where setter and getter is named with underscore
                                        # separator
                                        'underscoreSetter': {
                                            'python': 'lsst.daf.persistence.test.TestObjectUnderscoreSetter',
                                            'composite': {
                                                'foo': {'datasetType': 'basicObject1'}
                                            }
                                        },
                                        # simple object where setter and getter is named with camelcase
                                        'camelCaseSetter': {
                                            'python': 'lsst.daf.persistence.test.TestObjectCamelCaseSetter',
                                            'composite': {
                                                'foo': {'datasetType': 'basicObject1'}
                                            }
                                        }
                                    }})

        repoArgs = dafPersist.RepositoryArgs(root=self.firstRepoPath, policy=self.policy,
                                             mapper='lsst.obs.base.test.CompositeMapper')
        butler = dafPersist.Butler(outputs=repoArgs)
        butler.put(self.objA, 'basicObject1', dataId={'id': 'foo'})
        butler.put(self.objB, 'basicObject2', dataId={'name': 'bar'})
        del butler
        del repoArgs

    def tearDown(self):
        if os.path.exists(self.testData):
            shutil.rmtree(self.testData)

    def testConstructor(self):
        """Test the case where the arguments to the default constructor match the component names and so the
        default constructor can be used by the generic assembler to assemble the object
        Uses getters named by the policy to disassemble the object.
        """
        repoArgs = dafPersist.RepositoryArgs(root=self.secondRepoPath, policy=self.policy)
        butler = dafPersist.Butler(inputs=self.firstRepoPath, outputs=repoArgs)
        verificationButler = dafPersist.Butler(inputs=self.secondRepoPath)

        objABPair = butler.get('basicPair', dataId={'id': 'foo', 'name': 'bar'})
        self.assertEqual(self.objA, objABPair.objA)
        self.assertEqual(self.objB, objABPair.objB)

        butler.put(objABPair, 'basicPair', dataId={'id': 'foo', 'name': 'bar'})
        # comparing the output files directly works so long as the storage is posix:

        # put the composite object and verify it disassembled into the right locations by loading the
        # components directly
        objA = verificationButler.get('basicObject1', dataId={'id': 'foo'})
        self.assertEqual(objA, objABPair.objA)
        objB = verificationButler.get('basicObject2', dataId={'name': 'bar'})
        self.assertEqual(objB, objABPair.objB)

        del objABPair

        objABPair = butler.get('gaDefCtorPair', dataId={'id': 'foo', 'name': 'bar'})
        self.assertEqual(self.objA, objABPair.objA)
        self.assertEqual(self.objB, objABPair.objB)
        self.assertTrue(objABPair.usedInitSetter)
        self.assertFalse(objABPair.usedASetter)
        self.assertFalse(objABPair.usedBSetter)

        # put the composite object and verify it disassembled into the right locations by loading the
        # components directly
        butler.put(objABPair, 'gaDefCtorPair', dataId={'id': 'baz', 'name': 'qux'})
        verObjA = verificationButler.get('basicObject1', dataId={'id': 'baz'})
        self.assertEqual(objABPair.objA, verObjA)
        verObjB = verificationButler.get('basicObject2', dataId={'name': 'qux'})
        self.assertEqual(objABPair.objB, verObjB)

    def testGenericAssemblerPolicySpecifiedSetterGetter(self):
        """Test the case where the component names do not have anything to do with the setter/getter names
        or the init function parameter names, and instead the component policy entry specifies the setter and
        getter names.
        """
        repoArgs = dafPersist.RepositoryArgs(root=self.secondRepoPath, policy=self.policy)
        butler = dafPersist.Butler(inputs=self.firstRepoPath, outputs=repoArgs)
        objABPair = butler.get('gaPairWithSetter', dataId={'id': 'foo', 'name': 'bar'})
        self.assertEqual(self.objA, objABPair.objA)
        self.assertEqual(self.objB, objABPair.objB)
        self.assertFalse(objABPair.usedInitSetter)
        self.assertTrue(objABPair.usedASetter)
        self.assertTrue(objABPair.usedBSetter)

        butler.put(objABPair, 'gaPairWithSetter', dataId={'id': 'foo', 'name': 'bar'})
        # comparing the output files directly works so long as the storage is posix:

        verificationButler = dafPersist.Butler(inputs=self.secondRepoPath)
        verObjA = verificationButler.get('basicObject1', dataId={'id': 'foo'})
        verObjB = verificationButler.get('basicObject2', dataId={'name': 'bar'})
        self.assertEqual(objABPair.objA, verObjA)
        self.assertEqual(objABPair.objB, verObjB)

    def testInferredNameUnderscoreSeparator(self):
        """Test the case where the name of the setter & getter is inferred by the policy name by prepending
        'set_' and get_
        """
        repoArgs = dafPersist.RepositoryArgs(root=self.secondRepoPath, policy=self.policy)
        butler = dafPersist.Butler(inputs=self.firstRepoPath, outputs=repoArgs)
        obj = butler.get('underscoreSetter', dataId={'id': 'foo'})
        self.assertEqual(self.objA, obj.get_foo())
        butler.put(obj, 'underscoreSetter', dataId={'id': 'foo'})

        verificationButler = dafPersist.Butler(inputs=self.secondRepoPath)
        componentObj = verificationButler.get('basicObject1', dataId={'id': 'foo'})
        self.assertEqual(componentObj, obj.get_foo())


    def testInferredNameCamelcase(self):
        """Test the case where the name of the setter & getter is inferred by the policy name by prepending
        'set' or 'get', to the capitalized component name. E.g. for component name 'foo' the setter and getter
        will be named setFoo and getFoo.
        """
        repoArgs = dafPersist.RepositoryArgs(root=self.secondRepoPath, policy=self.policy)
        butler = dafPersist.Butler(inputs=self.firstRepoPath, outputs=repoArgs)
        obj = butler.get('camelCaseSetter', dataId={'id': 'foo'})
        self.assertEqual(self.objA, obj.getFoo())
        butler.put(obj, 'camelCaseSetter', dataId={'id': 'foo'})

        verificationButler = dafPersist.Butler(inputs=self.secondRepoPath)
        componentObj = verificationButler.get('basicObject1', dataId={'id': 'foo'})
        self.assertEqual(componentObj, obj.getFoo())


def subsetAssembler(dataId, componentInfo, cls):
    obj = cls()
    obj.set_a(componentInfo['a'].obj)
    return obj

class TestSubset(unittest.TestCase):
    """A test case for composite object subset keyword."""

    def setUp(self):
        packageDir = getPackageDir('obs_base')
        self.testData = os.path.join(packageDir, 'tests', 'compositeSubset')
        self.firstRepoPath = os.path.join(self.testData, 'repo1')
        self.objA1 = dpTest.TestObject("abc")
        self.objA2 = dpTest.TestObject("ABC")
        self.objB = dpTest.TestObject("def")
        self.policy = dafPersist.Policy(
                                   {'camera': 'lsst.afw.cameraGeom.Camera',
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
                                                'a': {
                                                    'datasetType': 'basicObject1',
                                                    'subset': True
                                                }
                                            },
                                            'assembler': subsetAssembler
                                        },
                                    }})

        repoArgs = dafPersist.RepositoryArgs(root=self.firstRepoPath, policy=self.policy,
                                             mapper='lsst.obs.base.test.CompositeMapper')
        butler = dafPersist.Butler(outputs=repoArgs)
        butler.put(self.objA1, 'basicObject1', dataId={'id': 'foo1'})
        butler.put(self.objA2, 'basicObject1', dataId={'id': 'foo2'})
        butler.put(self.objB, 'basicObject2', dataId={'name': 'bar'})
        del butler
        del repoArgs


    def tearDown(self):
        if os.path.exists(self.testData):
            shutil.rmtree(self.testData)


    def test(self):
        """Verify that the generic assembler and disassembler work for objects that conform to the generic
        set/get API.
        """
        secondRepoPath = os.path.join(self.testData, 'repo2')
        repoArgs = dafPersist.RepositoryArgs(root=secondRepoPath, policy=self.policy)
        butler = dafPersist.Butler(inputs=self.firstRepoPath, outputs=repoArgs)
        # the name 'bar' will find the obj that was put as obj b. It expects to find n objects of dataset
        # type basicObject1. Since we don't specify any dataId that relates to basicObject1 (its only dataId
        # key is 'id'), it will return everything it finds according to its policy. In this case that should
        # be self.objA1 and self.objA2 that we put above. They will be in a list at objABPair.objA.
        objABPair = butler.get('basicPair', dataId={'name': 'bar'})
        objABPair.objA.sort()
        self.assertEqual(self.objA2, objABPair.objA[0])
        self.assertEqual(self.objA1, objABPair.objA[1])
        # subset is a get-only operation. To put, the dataId must be specified, so there's no put to test
        # here.


class TestInputOnly(unittest.TestCase):
    """A test case for composite input keyword."""

    def setUp(self):
        packageDir = getPackageDir('obs_base')
        self.testData = os.path.join(packageDir, 'tests', 'composite')
        self.firstRepoPath = os.path.join(self.testData, 'repo1')
        self.objA = dpTest.TestObject("abc")
        self.objB = dpTest.TestObject("def")
        self.policy = dafPersist.Policy(
                                   {'camera': 'lsst.afw.cameraGeom.Camera',
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
                                                'a': {
                                                    'datasetType': 'basicObject1'
                                                },
                                                'b': {
                                                    'datasetType': 'basicObject2',
                                                    'inputOnly' : True
                                                }
                                            },
                                            'assembler': 'lsst.daf.persistence.test.TestObjectPair.assembler',
                                            'disassembler': 'lsst.daf.persistence.test.TestObjectPair.disassembler'

                                        }
                                    }})

        repoArgs = dafPersist.RepositoryArgs(root=self.firstRepoPath,
                                             mapper='lsst.obs.base.test.CompositeMapper',
                                             policy=self.policy)
        butler = dafPersist.Butler(outputs=repoArgs)
        butler.put(self.objA, 'basicObject1', dataId={'id': 'foo'})
        butler.put(self.objB, 'basicObject2', dataId={'name': 'bar'})
        del butler
        del repoArgs

    def tearDown(self):
        if os.path.exists(self.testData):
            shutil.rmtree(self.testData)

    def test(self):
        """ Verify that when a type 3 dataset is put and one of its components is marked 'inputOnly' by the
        policy that the inputOnly comonent is not written.
        """
        secondRepoPath = os.path.join(self.testData, 'repo2')
        repoArgs = dafPersist.RepositoryArgs(root=secondRepoPath, policy=self.policy)
        butler = dafPersist.Butler(inputs=self.firstRepoPath, outputs=repoArgs)
        objABPair = butler.get('basicPair', dataId={'id': 'foo', 'name': 'bar'})

        verificationButler = dafPersist.Butler(outputs={'root': os.path.join(self.testData, 'repo3'),
                                                        'policy': self.policy,
                                                        'mapper': 'lsst.obs.base.test.CompositeMapper',
                                                        'mode': 'rw'})
        verificationButler.put(objABPair, 'basicPair', dataId={'id': 'foo', 'name': 'bar'})

        objA = verificationButler.get('basicObject1', {'id': 'foo'})
        self.assertEqual(objA, objABPair.objA)
        with self.assertRaises(RuntimeError):
            verificationButler.get('basicObject2', {'name': 'bar'}, immediate=True)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
