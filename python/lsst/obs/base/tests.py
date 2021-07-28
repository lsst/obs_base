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

"""
Test utilities for obs_base and concrete obs* packages.
"""

__all__ = (
    "ObsTests",
    "make_ramp_array",
    "make_ramp_exposure_trimmed",
    "make_ramp_exposure_untrimmed",
)

import logging
import numpy as np

from lsst.afw.image import Exposure
from lsst.afw.cameraGeom.utils import calcRawCcdBBox

from . import butler_tests
from . import mapper_tests
from . import camera_tests


class ObsTests(butler_tests.ButlerGetTests, mapper_tests.MapperTests,
               camera_tests.CameraTests):
    """Aggregator class for all of the obs_* test classes.

    Inherit from this class, then lsst.utils.tests.TestCase, in that order.

    Example subclass:

    .. code-block:: python

        class TestObs(lsst.obs.base.tests.ObsTests, lsst.utils.tests.TestCase):
            def setUp(self):
                self.setUp_tests(...)
                self.setUp_butler_get(...)
                self.setUp_mapper(...)
                self.setUp_camera(...)

    Notes
    -----
    The intention is for each obs package to have a single test class that
    inherits from this collector class, thus "automatically" getting all new
    tests. If those tests require setup that isn't defined in a given obs
    package, that obs package will be broken until updated. This is
    intentional, as a way to prevent obs packages from falling behind out of
    neglect.
    """

    def setUp_tests(self, butler, mapper, dataIds):
        """Set up the necessary shared variables used by multiple tests.

        Parameters
        ----------
        butler: lsst.daf.persistence.Butler
            A butler object, instantiated on the testdata repository for the
            obs package being tested.
        mapper: lsst.obs.CameraMapper
            A CameraMapper object for your camera, instantiated on the testdata
            repository the obs package being tested.
        dataIds: dict
            dictionary of (exposure name): (dataId of that exposure in the
            testdata repository), with unittest.SkipTest as the value for any
            exposures you do not have/do not want to test. It must contain a
            valid 'raw' dataId, in addition to 'bias','flat','dark', which may
            be set to SkipTest. For example::

                  self.dataIds = {'raw': {'visit': 1, 'filter': 'g'},
                                  'bias': {'visit': 1},
                                  'flat': {'visit': 1},
                                  'dark': unittest.SkipTest
                                 }
        """
        self.butler = butler
        self.mapper = mapper
        self.dataIds = dataIds
        self.log = logging.getLogger('ObsTests')

    def tearDown(self):
        del self.butler
        del self.mapper
        super(ObsTests, self).tearDown()


def make_ramp_array(bbox, pedestal):
    """Make a 2-d ramp array.

    Parameters
    ----------
    bbox : `lsst.geom.Box2I`
        Bounding box for the array.
    pedestal : `int`
        Minimum value for the ramp.

    Returns
    -------
    ramp : `np.ndarray`
        A 2-d array with shape ``(bbox.getHeight(), bbox.getWidth())``.
    end : `int`
        One past the maximum value in the ramp (for use as the
        pedestal for another box).
    """
    end = pedestal + bbox.getArea()
    return np.arange(pedestal, end).reshape(bbox.getHeight(), bbox.getWidth()), end


def make_ramp_exposure_untrimmed(detector, dtype=None):
    """Create an untrimmed, assembled exposure with different ramps for
    each sub-amplifier region.

    Parameters
    ----------
    detector : `lsst.afw.cameraGeom.Detector`
        Detector object that the new exposure should match.  Must have all amp
        flips and offsets set to False/zero (i.e. represent an already-
        assembled image).
    dtype : `np.dtype`, optional
        Type of the new exposure.  Defaults to ``int32``.

    Returns
    -------
    exposure : `lsst.afw.image.Exposure`
        New exposure with the given detector attached.
    """
    if dtype is None:
        dtype = np.dtype(np.int32)
    ramp_exposure = Exposure(calcRawCcdBBox(detector), dtype=np.dtype(dtype))
    ramp_exposure.setDetector(detector)
    pedestal = 0
    for amp in detector:
        for name in ("HorizontalOverscan", "VerticalOverscan", "Prescan", "Data"):
            bbox = getattr(amp, f"getRaw{name}BBox")()
            ramp, pedestal = make_ramp_array(bbox, pedestal)
            ramp_exposure.image[bbox].array[:, :] = ramp
    return ramp_exposure


def make_ramp_exposure_trimmed(detector, dtype=None):
    """Create a trimmed, assembled exposure with different ramps for
    each amplifier region.

    Parameters
    ----------
    detector : `lsst.afw.cameraGeom.Detector`
        Detector object that the new exposure should match.
    dtype : `np.dtype`, optional
        Type of the new exposure.  Defaults to ``int32``.

    Returns
    -------
    exposure : `lsst.afw.image.Exposure`
        New exposure with the given detector attached.
    """
    if dtype is None:
        dtype = np.dtype(np.int32)
    ramp_exposure = Exposure(detector.getBBox(), dtype=np.dtype(dtype))
    ramp_exposure.setDetector(detector)
    pedestal = 0
    for amp in detector:
        ramp, pedestal = make_ramp_array(amp.getBBox(), pedestal)
        ramp_exposure.image[amp.getBBox()].array[:, :] = ramp
    return ramp_exposure
