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
import math

import lsst.geom
from lsst.afw.cameraGeom import FIELD_ANGLE, FOCAL_PLANE
from lsst.daf.butler import Butler

__all__ = ["CameraTests"]


class CameraTests(metaclass=abc.ABCMeta):
    """Tests that the butler returns a useable Camera.

    In the subclasses's setUp():
        * Call setUp_camera() to fill in required parameters.
    """

    def setUp_camera(
        self,
        camera_name=None,
        n_detectors=None,
        first_detector_name=None,
        plate_scale=None,
    ):
        """Set up the necessary variables for camera tests.

        Parameters
        ----------

        camera_name : `str`
            name of this camera
        n_detectors : `int`
            number of detectors in this camera
        first_detector_name : `str`
            name of the first detector in this camera
        plate_scale : `lsst.geom.Angle`
            plate scale at center of focal plane, as angle-on-sky/mm
        """
        fields = [
            "camera_name",
            "n_detectors",
            "first_detector_name",
            "plate_scale",
        ]
        CameraData = collections.namedtuple("CameraData", fields)
        self.camera_data = CameraData(
            camera_name=camera_name,
            n_detectors=n_detectors,
            first_detector_name=first_detector_name,
            plate_scale=plate_scale,
        )

    def _butler_args(self):
        """Arguments to pass to butler get in addition to camera dataset type.

        Returns
        -------
        kwargs : `dict`
            The arguments to add.
        """
        kwargs = {}
        if isinstance(self.butler, Butler):
            kwargs = dict(instrument=self.camera_data.camera_name)
        return kwargs

    def test_iterable(self):
        """Simplest camera test: can we get a Camera instance, and does
        iterating return Detectors?"""
        camera = self.butler.get("camera", **self._butler_args())
        self.assertIsInstance(camera, lsst.afw.cameraGeom.Camera)
        for detector in camera:
            msg = "Failed for detector={}".format(detector)
            self.assertIsInstance(detector, lsst.afw.cameraGeom.Detector, msg=msg)

    def test_camera_butler(self):
        """Check that the butler returns the right type of camera."""
        camera = self.butler.get("camera", **self._butler_args())
        self.assertEqual(camera.getName(), self.camera_data.camera_name)
        self.assertEqual(len(camera), self.camera_data.n_detectors)
        self.assertEqual(next(iter(camera)).getName(), self.camera_data.first_detector_name)

    def test_plate_scale(self):
        """Check the plate scale at center of focal plane

        Check plate_scale using the FOCAL_PLANE to FIELD_ANGLE transform
        from the camera.
        """
        plate_scale = self.camera_data.plate_scale
        self.assertIsNotNone(plate_scale)
        camera = self.butler.get("camera", **self._butler_args())
        focalPlaneToFieldAngle = camera.getTransformMap().getTransform(FOCAL_PLANE, FIELD_ANGLE)
        focalPlaneRadiusMm = 0.001  # an offset small enough to be in the linear regime
        for offsetAngleRad in (0.0, 0.65, 1.3):  # direction of offset; a few arbitrary angles
            cosAng = math.cos(offsetAngleRad)
            sinAng = math.sin(offsetAngleRad)
            fieldAngleRadians = focalPlaneToFieldAngle.applyForward(
                lsst.geom.Point2D(cosAng * focalPlaneRadiusMm, sinAng * focalPlaneRadiusMm)
            )
            fieldAngleRadius = math.hypot(*fieldAngleRadians) * lsst.geom.radians
            measuredScale1 = fieldAngleRadius / focalPlaneRadiusMm
            self.assertAnglesAlmostEqual(measuredScale1, plate_scale)

            focalPlanePos = focalPlaneToFieldAngle.applyInverse(
                lsst.geom.Point2D(
                    fieldAngleRadius.asRadians() * cosAng, fieldAngleRadius.asRadians() * sinAng
                )
            )
            focalPlaneRadiusMm2 = math.hypot(*focalPlanePos)
            measureScale2 = fieldAngleRadius / focalPlaneRadiusMm2
            self.assertAnglesAlmostEqual(measureScale2, plate_scale)
