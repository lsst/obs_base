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
import collections

import lsst.afw.geom


class CameraTests(with_metaclass(abc.ABCMeta)):
    """
    Tests that the butler returns a useable Camera.

    In the subclasses's setUp():
        * Call setUp_camera() to fill in required parameters.
    """

    def setUp_camera(self,
                     camera_name=None,
                     n_detectors=None,
                     first_detector_name=None
                     ):
        """
        Set up the necessary variables for camera tests.

        Parameters
        ----------

         camera_name : `str`
             name of this camera
         n_detectors : `int`
            number of detectors in this camera
         first_detector_name : `str`
            name of the first detector in this camera
        """
        fields = ['camera_name',
                  'n_detectors',
                  'first_detector_name'
                  ]
        CameraData = collections.namedtuple("CameraData", fields)
        self.camera_data = CameraData(camera_name=camera_name,
                                      n_detectors=n_detectors,
                                      first_detector_name=first_detector_name
                                      )

    def test_iterable(self):
        """Simplest camera test: can we get a Camera instance, and does iterating return Detectors?"""
        camera = self.butler.get('camera', immediate=True)
        self.assertIsInstance(camera, lsst.afw.cameraGeom.Camera)
        for detector in camera:
            msg = "Failed for detector={}".format(detector)
            self.assertIsInstance(detector, lsst.afw.cameraGeom.Detector, msg=msg)

    def test_camera_butler(self):
        """Check that the butler returns the right type of camera."""
        camera = self.butler.get('camera', immediate=True)
        self.assertEqual(camera.getName(), self.camera_data.camera_name)
        self.assertEqual(len(camera), self.camera_data.n_detectors)
        self.assertEqual(next(iter(camera)).getName(), self.camera_data.first_detector_name)
