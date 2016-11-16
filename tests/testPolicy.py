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

import unittest
import lsst.utils.tests

import lsst.daf.persistence as dafPersist
import lsst.daf.persistence.test as dpTest
from lsst.utils import getPackageDir
import os
import shutil
import yaml


class TestPolicyInRepo(unittest.TestCase):
    def setUp(self):
        packageDir = getPackageDir('obs_base')
        self.testData = os.path.join(packageDir, 'tests', 'TestPolicyInRepo')

    def tearDown(self):
        if os.path.exists(self.testData):
            shutil.rmtree(self.testData)

    def test(self):
        """Verify that when specifying a repo policy that the policy gets written & loaded correctly.
        """

        objA = dpTest.TestObject("abc")
        objB = dpTest.TestObject("def")

        firstRepoPath = os.path.join(self.testData, 'repo1')
        secondRepoPath = os.path.join(self.testData, 'repo2')

        policy = dafPersist.Policy({'camera': 'lsst.afw.cameraGeom.Camera',
                                    'datasets': {
                                        'basicObject1': {
                                            'python': 'lsst.daf.persistence.test.TestObject',
                                            'template': 'basic/id%(id)s.pickle',
                                            'storage': 'PickleStorage'},
                                        }
                                    })

        firstRepoPath = os.path.join(self.testData, 'repo1')
        repoArgs = dafPersist.RepositoryArgs(root=firstRepoPath,
                                             mapper='lsst.obs.base.test.CompositeMapper',
                                             policy=policy)
        butler = dafPersist.Butler(outputs=repoArgs)

        with open(os.path.join(firstRepoPath, 'repositoryCfg.yaml')) as f:
            cfg = yaml.load(f)
        self.assertEqual(cfg.policy, policy)
        butler.put(objA, 'basicObject1', {'id':1})

        del butler
        del repoArgs

        # Test that a newly-initialized butler can find the policy in the repositoryCfg.
        repoArgs = dafPersist.RepositoryArgs(root=firstRepoPath)
        butler = dafPersist.Butler(inputs=repoArgs)
        reloadedObjA = butler.get('basicObject1', {'id':1})
        self.assertEqual(reloadedObjA, objA)



class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
