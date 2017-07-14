#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2017 LSST Corporation.
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

import cerberus
import unittest
import yaml

import lsst.utils.tests
import lsst.daf.persistence as dafPersist


class Test(unittest.TestCase):
    """A test case for composite object i/o."""

    def setUp(self):
        schemaFile = dafPersist.Policy.defaultPolicyFile(
            'obs_base', 'DatasetMappingDictionary.yaml', 'policy')
        with open(schemaFile) as f:
            self.schema = yaml.load(f)

    @staticmethod
    def _validate(schema, policy):
        v = cerberus.Validator(schema)
        return v.validate(policy)

    def testValidDatasetsDatasetParameters(self):
        """Verify that the datasets defined in dataset may contain specified parameters."""
        policy = {
            'policy': {
                'datasets': {
                    'processCcd_config': {
                        'persistable': 'Config',
                        'python': 'lsst.pipe.tasks.processCcd.ProcessCcdConfig',
                        'storage': 'ConfigStorage',
                        'template': 'config/processCcd.py'
                    }
                },
            }
        }
        v = cerberus.Validator(self.schema)
        self.assertTrue(v.validate(policy), v.errors)

    def testExtraDatasetsDatasetParameters(self):
        """Verify that the datasets defined in dataset may not contain unspecified parameters."""
        policy = {
            'policy': {
                'datasets': {
                    'processCcd_config': {
                        'template': 'config/processCcd.py',
                        'foo': 'bar',  # NB foo is not declared in the datasets' dataset schema.
                    }
                },
            }
        }
        v = cerberus.Validator(self.schema)
        self.assertFalse(v.validate(policy), v.errors)

    def allowUnknownDatasetCategory(self):
        """Test that the policy may contain top-level categories that are not declared in the schema.
        TODO do we want to allow this?"""
        policy = {
            'policy': {
                'foo': 123
            }
        }
        v = cerberus.Validator(self.schema)
        self.assertTrue(v.validate(policy), v.errors)

    def testEmptyDatasetTemplate(self):
        """Verify that the datasets template may not be empty."""
        policy = {
            'policy': {
                'datasets': {
                    'processCcd_config': {
                        'template': '',
                    }
                },
            }
        }
        v = cerberus.Validator(self.schema)
        self.assertFalse(v.validate(policy), v.errors)

    def testMissingDatasetTemplate(self):
        """Verify that the datasets template may not be missing."""
        policy = {
            'policy': {
                'datasets': {
                    'processCcd_config': {}
                },
            }
        }
        v = cerberus.Validator(self.schema)
        self.assertFalse(v.validate(policy), v.errors)

    def testDatasetIsNotDictFails(self):
        """Test that policy.dataset must be a dict"""
        policy = {
            'policy': {
                'datasets': 123
            }
        }
        v = cerberus.Validator(self.schema)
        self.assertFalse(v.validate(policy), v.errors)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
