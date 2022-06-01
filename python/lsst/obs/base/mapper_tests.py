# This file is part of obs_base.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import abc
import collections
import inspect
import os

import lsst.afw.geom
import lsst.afw.image
import lsst.daf.persistence
import lsst.utils.tests

__all__ = ["MapperTests"]


class MapperTests(metaclass=abc.ABCMeta):
    """
    Generic tests of obs_* package mapper functionality.

    In the subclasses's setUp():
        * Call setUp_mapper() to fill in required parameters.
    """

    def setUp_mapper(
        self,
        output=None,
        path_to_raw=None,
        keys=None,
        query_format=None,
        queryMetadata=None,
        metadata_output_path=None,
        map_python_type=None,
        map_python_std_type=None,
        map_cpp_type=None,
        map_storage_name=None,
        raw_filename=None,
        default_level=None,
        raw_levels=None,
        test_config_metadata=True,
    ):
        """
        Set up the necessary variables for mapper tests.

        Parameters
        ----------

        output : `str`
            full path to output repository (can be the same as data_dir input
            repository)
        path_to_raw : `str`
            full path to the raw file referenced by dataIds['raw']
        keys : `set`
            dictionary keys that this mapper should contain
        query_format : `list`
            format list for the results portion of queryMetadata
        queryMetadata : `tuple` of (`dict`, `tuple`)
            dataIds and the results of calling them in queryMetadata
        metadata_output_path : `str`
            path to metadata output associated with dataIds['raw']
        map_python_type : `type`
            full python type specification returned by the mapper for
            dataIds['raw']
        map_python_std_type : `type`
            full python type specification returned by the mapper for
            dataIds['raw'] after standardization
        map_cpp_type : `str`
            C++ type specification returned by the mapper for dataIds['raw']
        map_storage_name : `str`
            butler name for the storage type dataIds['raw']
        raw_filename : `str`
            Name of the raw files returned by the mapper for dataIds['raw']
        default_level : `str`
            value returned from mapper.getDefaultLevel
        raw_levels : `tuple` of (`str`, `set` of `str`)
            (level, expect) level and expected mapper return for
            mapper.getKeys('raw', level)
        test_config_metadata : `bool`
            Test persisted config and metadata?  These tests may not be
            appropriate for test stand data. Defaults to True.
        """
        fields = [
            "output",
            "path_to_raw",
            "keys",
            "query_format",
            "queryMetadata",
            "metadata_output_path",
            "map_python_type",
            "map_python_std_type",
            "map_cpp_type",
            "map_storage_name",
            "raw_filename",
            "default_level",
            "raw_levels",
            "test_config_metadata",
        ]
        MapperData = collections.namedtuple("MapperData", fields)
        self.mapper_data = MapperData(
            output=output,
            path_to_raw=path_to_raw,
            keys=keys,
            query_format=query_format,
            queryMetadata=queryMetadata,
            metadata_output_path=metadata_output_path,
            map_python_type=map_python_type,
            map_python_std_type=map_python_std_type,
            map_cpp_type=map_cpp_type,
            map_storage_name=map_storage_name,
            raw_filename=raw_filename,
            default_level=default_level,
            raw_levels=raw_levels,
            test_config_metadata=test_config_metadata,
        )

    def _check_mapper(self):
        if self.mapper is None:
            self.skipTest("No mapper defined. Likely a gen3 test environment.")

    def test_map_config_data(self):
        self._check_mapper()
        if not self.mapper_data.test_config_metadata:
            self.skipTest("Skipping %s as requested" % (inspect.currentframe().f_code.co_name))
        dataId = self.dataIds["raw"]
        butlerLocation = self.mapper.map("processCcd_config_filename", dataId)
        self.assertEqual(butlerLocation.getPythonType(), "lsst.pipe.tasks.processCcd.ProcessCcdConfig")
        self.assertEqual(butlerLocation.getCppType(), "Config")
        self.assertEqual(butlerLocation.getStorageName(), "ConfigStorage")
        processCcd_path = os.path.join("config", "processCcd.py")
        self.assertEqual(self.mapper.root, butlerLocation.getStorage().root)
        self.assertEqual(butlerLocation.getLocations(), [processCcd_path])
        for k, v in dataId.items():
            self.assertEqual(
                butlerLocation.getAdditionalData().getScalar(k), v, msg="Failed for key={}".format(k)
            )

    def test_map_metadata_data(self):
        self._check_mapper()
        if not self.mapper_data.test_config_metadata:
            self.skipTest("Skipping %s as requested" % (inspect.currentframe().f_code.co_name))
        dataId = self.dataIds["raw"]
        butlerLocation = self.mapper.map_processCcd_metadata(dataId)
        self.assertEqual(butlerLocation.getPythonType(), "lsst.daf.base.PropertySet")
        self.assertEqual(butlerLocation.getCppType(), "PropertySet")
        self.assertEqual(butlerLocation.getStorageName(), "YamlStorage")
        self.assertEqual(butlerLocation.getLocations(), [self.mapper_data.metadata_output_path])
        for k, v in dataId.items():
            self.assertEqual(
                butlerLocation.getAdditionalData().getScalar(k), v, msg="Failed for key={}".format(k)
            )

    def test_keys(self):
        self._check_mapper()
        self.assertEqual(set(self.mapper.keys()), self.mapper_data.keys)

    def test_get_dataset_types(self):
        self._check_mapper()
        if not self.mapper_data.test_config_metadata:
            self.skipTest("Skipping %s as requested" % (inspect.currentframe().f_code.co_name))
        someKeys = set(["raw", "processCcd_config", "processCcd_metadata"])
        self.assertTrue(set(self.mapper.getDatasetTypes()).issuperset(someKeys))

    def test_get_keys_raw(self):
        self._check_mapper()
        for level, expect in self.mapper_data.raw_levels:
            result = self.mapper.getKeys("raw", level)
            self.assertEqual(set(result), expect, msg="Failed for level={}".format(level))

    def test_get_default_level(self):
        self._check_mapper()
        self.assertEqual(self.mapper.getDefaultLevel(), self.mapper_data.default_level)

    def _test_map(self, butlerLocation, dataId):
        self._check_mapper()
        self.assertEqual(butlerLocation.getPythonType(), self.mapper_data.map_python_type)
        self.assertEqual(butlerLocation.getCppType(), self.mapper_data.map_cpp_type)
        self.assertEqual(butlerLocation.getStorageName(), self.mapper_data.map_storage_name)
        locationList = butlerLocation.getLocations()
        self.assertEqual(len(locationList), 1)
        fileName = os.path.basename(locationList[0])
        self.assertEqual(fileName, self.mapper_data.raw_filename)
        for k, v in dataId.items():
            self.assertEqual(
                butlerLocation.getAdditionalData().getScalar(k), v, msg="Failed for key={}".format(k)
            )

    def test_map(self):
        self._check_mapper()
        dataId = self.dataIds["raw"]
        location = self.mapper.map_raw(dataId)
        if not isinstance(location, lsst.daf.persistence.butlerLocation.ButlerComposite):
            self._test_map(location, dataId)
        else:
            self.log.warning(
                """ButlerComposite datasets are not tested for mapper functions.  Though
ButlerComposites duck type as ButlerLocations in some ways, they do not
share enough methods to be usefully tested by the same function.  Note
there are tests of the objects in the package in which they are implemented."""
            )
        #  This should be the same as above.  Testing that both the generic
        # and specific interface work for mapping the raw.
        location = self.mapper.map("raw", dataId)
        if not isinstance(location, lsst.daf.persistence.butlerLocation.ButlerComposite):
            self._test_map(location, dataId)

    def test_query_metadata(self):
        """
        Test expansion of incomplete information of the available data in this
        obs package's testdata repo.
        """
        self._check_mapper()
        for query, expect in self.mapper_data.queryMetadata:
            # queryMetadata returns tuples of available items of the 2nd
            # parameter.
            result = self.mapper.queryMetadata("raw", self.mapper_data.query_format, query)
            self.assertEqual(sorted(result), sorted(expect), msg="Failed for query={}".format(query))

    def test_can_standardize(self):
        self._check_mapper()
        self.assertTrue(self.mapper.canStandardize("raw"))
        self.assertFalse(self.mapper.canStandardize("camera"))
        if not self.mapper_data.test_config_metadata:
            self.assertFalse(self.mapper.canStandardize("processCcd_config"))
            self.assertFalse(self.mapper.canStandardize("processCcd_metadata"))

    def _test_validate(self, dataId):
        self._check_mapper()
        self.assertEqual(self.mapper.validate(dataId), dataId)

    def test_validate(self):
        self._check_mapper()
        self._test_validate({"visit": 1, "filter": "g"})
        self._test_validate({"visit": 2, "filter": "r"})
        self._test_validate({"visit": 3, "filter": "g", "tract": 4})
        # NOTE: when DM-7909 is completed, add assertRaises test here.
        # visit must be an integers
