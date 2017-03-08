from __future__ import absolute_import, division, print_function
from future.utils import with_metaclass
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
import abc
import inspect
import unittest
import collections


class ButlerGetTests(with_metaclass(abc.ABCMeta)):
    """
    Tests of obs_* Butler get() functionality.

    In the subclasses's setUp():
        * Call setUp_butler_get() to fill in required parameters.
    """

    def setUp_butler_get(self,
                         ccdExposureId_bits=None,
                         exposureIds=None,
                         filters=None,
                         exptimes=None,
                         detectorIds=None,
                         detector_names=None,
                         detector_serials=None,
                         dimensions=None,
                         sky_origin=None,
                         raw_subsets=None,
                         good_detectorIds=None,
                         bad_detectorIds=None,
                         linearizer_type=None
                         ):
        """
        Set up the necessary variables for butlerGet tests.

        All "exposure name" entries below should correspond to an entry in
        self.dataIds.

        Parameters
        ----------

        ccdExposureId_bits : `int`
            expected value of ccdExposureId_bits
        exposureIds : `dict`
            dict of exposure name : ccdExposureId (the number as returned by the butler)
        filters : `dict`
            dict of exposure name : filter name
        exptimes : `dict`
            dict of exposure name : exposure time
        detector_names : `dict`
            dict of exposure name : detector name
        detectorIds : `dict`
            dict of exposure name : detectorId
        detector_serials : `dict`
            dict of exposure name : detector serial
        dimensions : `dict`
            dict of exposure name : dimensions (as a geom.Extent2I)
        sky_origin : `tuple` of `float`
            Longitude, Latitude of 'raw' exposure
        raw_subsets : `tuple` of (kwargs, `int`)
            keyword args and expected number of subsets for butler.subset('raw', **kwargs)
        good_detectorIds : `list` of `int`
            list of valid ccd numbers
        bad_detectorIds : `list` of `int`
            list of invalid ccd numbers
        linearizer_type : `dict`
            dict of detectorId (usually `int`): LinearizerType
            (e.g. lsst.ip.isr.LinearizeLookupTable.LinearityType),
            or unittest.SkipTest to skip all linearizer tests.
        """

        fields = ['ccdExposureId_bits',
                  'exposureIds',
                  'filters',
                  'exptimes',
                  'detector_names',
                  'detectorIds',
                  'detector_serials',
                  'dimensions',
                  'sky_origin',
                  'raw_subsets',
                  'good_detectorIds',
                  'bad_detectorIds',
                  'linearizer_type'
                  ]
        ButlerGet = collections.namedtuple("ButlerGetData", fields)

        self.butler_get_data = ButlerGet(ccdExposureId_bits=ccdExposureId_bits,
                                         exposureIds=exposureIds,
                                         filters=filters,
                                         exptimes=exptimes,
                                         detectorIds=detectorIds,
                                         detector_names=detector_names,
                                         detector_serials=detector_serials,
                                         dimensions=dimensions,
                                         sky_origin=sky_origin,
                                         raw_subsets=raw_subsets,
                                         good_detectorIds=good_detectorIds,
                                         bad_detectorIds=bad_detectorIds,
                                         linearizer_type=linearizer_type
                                         )

    def test_exposureId_bits(self):
        bits = self.butler.get('ccdExposureId_bits')
        self.assertEqual(bits, self.butler_get_data.ccdExposureId_bits)

    def _test_exposure(self, name):
        if self.dataIds[name] is unittest.SkipTest:
            self.skipTest('Skipping %s as requested' % (inspect.currentframe().f_code.co_name))

        exp = self.butler.get(name, self.dataIds[name])
        self.assertEqual(exp.getDimensions(), self.butler_get_data.dimensions[name])
        self.assertEqual(exp.getDetector().getId(), self.butler_get_data.detectorIds[name])
        self.assertEqual(exp.getDetector().getName(), self.butler_get_data.detector_names[name])
        self.assertEqual(exp.getDetector().getSerial(), self.butler_get_data.detector_serials[name])
        self.assertEqual(exp.getFilter().getName(), self.butler_get_data.filters[name])
        exposureId = self.butler.get('ccdExposureId', dataId=self.dataIds[name])
        self.assertEqual(exposureId, self.butler_get_data.exposureIds[name])
        return exp

    def test_raw(self):
        exp = self._test_exposure('raw')
        # We only test the existence of WCS in the raw files, since it's only well-defined
        # for raw, and other exposure types could have or not have a WCS depending
        # on various implementation details.
        self.assertEqual(exp.hasWcs(), True)
        origin = exp.getWcs().getSkyOrigin()
        self.assertEqual(exp.getInfo().getVisitInfo().getExposureTime(), self.butler_get_data.exptimes['raw'])
        self.assertCoordsNearlyEqual(origin, self.butler_get_data.sky_origin)

    def test_bias(self):
        self._test_exposure('bias')

    def test_dark(self):
        self._test_exposure('dark')

    def test_flat(self):
        self._test_exposure('flat')

    @unittest.skip('Cannot test this, as there is a bug in the butler! DM-8097')
    def test_raw_sub_bbox(self):
        exp = self.butler.get('raw', self.dataIds['raw'], immediate=True)
        bbox = exp.getBBox()
        bbox.grow(-1)
        sub = self.butler.get("raw_sub", self.dataIds['raw'], bbox=bbox, immediate=True)
        self.assertEqual(sub.getImage().getBBox(), bbox)
        self.assertImagesEqual(sub, exp.Factory(exp, bbox))

    def test_subset_raw(self):
        for kwargs, expect in self.butler_get_data.raw_subsets:
            subset = self.butler.subset("raw", **kwargs)
            self.assertEqual(len(subset), expect, msg="Failed for kwargs: {}".format(kwargs))

    def test_get_linearizer(self):
        """Test that we can get a linearizer for good detectorIds."""
        if self.butler_get_data.linearizer_type is unittest.SkipTest:
            self.skipTest('Skipping %s as requested' % (inspect.currentframe().f_code.co_name))

        camera = self.butler.get("camera")
        for detectorId in self.butler_get_data.good_detectorIds:
            detector = camera[detectorId]
            linearizer = self.butler.get("linearizer", dataId=dict(ccd=detectorId), immediate=True)
            self.assertEqual(linearizer.LinearityType, self.butler_get_data.linearizer_type)
            linearizer.checkDetector(detector)

    def test_get_linearizer_bad_detectorIds(self):
        """Do bad detectorIds raise?"""
        if self.butler_get_data.linearizer_type is unittest.SkipTest:
            self.skipTest('Skipping %s as requested' % (inspect.currentframe().f_code.co_name))

        for badccd in self.butler_get_data.bad_detectorIds:
            with self.assertRaises(RuntimeError):
                self.butler.get("linearizer", dataId=dict(ccd=badccd), immediate=True)
