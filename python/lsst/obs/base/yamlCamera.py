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

from functools import lru_cache

import numpy as np
import yaml

import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.geom as afwGeom
import lsst.geom as geom
from lsst.afw.cameraGeom import Amplifier, Camera, ReadoutCorner

__all__ = ["makeCamera"]


@lru_cache
def makeCamera(cameraFile):
    """Construct an imaging camera (e.g. the LSST 3Gpix camera).

    Parameters
    ----------
    cameraFile : `str`
        Camera description YAML file.

    Returns
    -------
    camera : `lsst.afw.cameraGeom.Camera`
        The desired Camera
    """
    with open(cameraFile) as fd:
        cameraParams = yaml.load(fd, Loader=yaml.CLoader)

    cameraName = cameraParams["name"]

    #
    # Handle distortion models.
    #
    plateScale = geom.Angle(cameraParams["plateScale"], geom.arcseconds)
    nativeSys = cameraGeom.CameraSys(cameraParams["transforms"].pop("nativeSys"))
    transforms = makeTransformDict(nativeSys, cameraParams["transforms"], plateScale)

    ccdParams = cameraParams["CCDs"]
    detectorConfigList = makeDetectorConfigList(ccdParams)

    amplifierDict = {}
    for ccdName, ccdValues in ccdParams.items():
        amplifierDict[ccdName] = makeAmplifierList(ccdValues)

    return makeCameraFromCatalogs(cameraName, detectorConfigList, nativeSys, transforms, amplifierDict)


def makeDetectorConfigList(ccdParams):
    """Make a list of detector configs.

    Returns
    -------
    detectorConfig : `list` of `lsst.afw.cameraGeom.DetectorConfig`
        A list of detector configs.
    """
    detectorConfigs = []
    for name, ccd in ccdParams.items():
        detectorConfig = cameraGeom.DetectorConfig()
        detectorConfigs.append(detectorConfig)

        detectorConfig.name = name
        detectorConfig.id = ccd["id"]
        detectorConfig.serial = ccd["serial"]
        detectorConfig.detectorType = ccd["detectorType"]
        if "physicalType" in ccd:
            detectorConfig.physicalType = ccd["physicalType"]
        # This is the orientation we need to put the serial direction along
        # the x-axis
        detectorConfig.bbox_x0, detectorConfig.bbox_y0 = ccd["bbox"][0]
        detectorConfig.bbox_x1, detectorConfig.bbox_y1 = ccd["bbox"][1]
        detectorConfig.pixelSize_x, detectorConfig.pixelSize_y = ccd["pixelSize"]
        detectorConfig.transformDict.nativeSys = ccd["transformDict"]["nativeSys"]
        transforms = ccd["transformDict"]["transforms"]
        detectorConfig.transformDict.transforms = None if transforms == "None" else transforms
        detectorConfig.refpos_x, detectorConfig.refpos_y = ccd["refpos"]
        if len(ccd["offset"]) == 2:
            detectorConfig.offset_x, detectorConfig.offset_y = ccd["offset"]
            detectorConfig.offset_z = 0.0
        else:
            detectorConfig.offset_x, detectorConfig.offset_y, detectorConfig.offset_z = ccd["offset"]
        detectorConfig.transposeDetector = ccd["transposeDetector"]
        detectorConfig.pitchDeg = ccd["pitch"]
        detectorConfig.yawDeg = ccd["yaw"]
        detectorConfig.rollDeg = ccd["roll"]
        if "crosstalk" in ccd:
            detectorConfig.crosstalk = ccd["crosstalk"]

    return detectorConfigs


def makeAmplifierList(ccd):
    """Construct a list of AmplifierBuilder objects."""
    # Much of this will need to be filled in when we know it.
    assert len(ccd) > 0
    amp = list(ccd["amplifiers"].values())[0]

    rawBBox = makeBBoxFromList(amp["rawBBox"])  # total in file
    xRawExtent, yRawExtent = rawBBox.getDimensions()

    readCorners = {
        "LL": ReadoutCorner.LL,
        "LR": ReadoutCorner.LR,
        "UL": ReadoutCorner.UL,
        "UR": ReadoutCorner.UR,
    }

    amplifierList = []
    for name, amp in sorted(ccd["amplifiers"].items(), key=lambda x: x[1]["hdu"]):
        amplifier = Amplifier.Builder()
        amplifier.setName(name)

        ix, iy = amp["ixy"]
        perAmpData = amp["perAmpData"]
        if perAmpData:
            x0, y0 = 0, 0  # origin of data within each amp image
        else:
            x0, y0 = ix * xRawExtent, iy * yRawExtent

        rawDataBBox = makeBBoxFromList(amp["rawDataBBox"])  # Photosensitive area
        xDataExtent, yDataExtent = rawDataBBox.getDimensions()
        amplifier.setBBox(
            geom.BoxI(geom.PointI(ix * xDataExtent, iy * yDataExtent), rawDataBBox.getDimensions())
        )

        rawBBox = makeBBoxFromList(amp["rawBBox"])
        rawBBox.shift(geom.ExtentI(x0, y0))
        amplifier.setRawBBox(rawBBox)

        rawDataBBox = makeBBoxFromList(amp["rawDataBBox"])
        rawDataBBox.shift(geom.ExtentI(x0, y0))
        amplifier.setRawDataBBox(rawDataBBox)

        rawSerialOverscanBBox = makeBBoxFromList(amp["rawSerialOverscanBBox"])
        rawSerialOverscanBBox.shift(geom.ExtentI(x0, y0))
        amplifier.setRawHorizontalOverscanBBox(rawSerialOverscanBBox)

        rawParallelOverscanBBox = makeBBoxFromList(amp["rawParallelOverscanBBox"])
        rawParallelOverscanBBox.shift(geom.ExtentI(x0, y0))
        amplifier.setRawVerticalOverscanBBox(rawParallelOverscanBBox)

        rawSerialPrescanBBox = makeBBoxFromList(amp["rawSerialPrescanBBox"])
        rawSerialPrescanBBox.shift(geom.ExtentI(x0, y0))
        amplifier.setRawPrescanBBox(rawSerialPrescanBBox)

        if perAmpData:
            amplifier.setRawXYOffset(geom.Extent2I(ix * xRawExtent, iy * yRawExtent))
        else:
            amplifier.setRawXYOffset(geom.Extent2I(0, 0))

        amplifier.setReadoutCorner(readCorners[amp["readCorner"]])
        amplifier.setGain(amp["gain"])
        amplifier.setReadNoise(amp["readNoise"])
        amplifier.setSaturation(amp["saturation"])
        amplifier.setSuspectLevel(amp.get("suspect", np.nan))

        # flip data when assembling if needs be (e.g. data from the serial at
        # the top of a CCD)
        flipX, flipY = amp.get("flipXY")

        amplifier.setRawFlipX(flipX)
        amplifier.setRawFlipY(flipY)
        # linearity placeholder stuff
        amplifier.setLinearityCoeffs([float(val) for val in amp["linearityCoeffs"]])
        amplifier.setLinearityType(amp["linearityType"])
        amplifier.setLinearityThreshold(float(amp["linearityThreshold"]))
        amplifier.setLinearityMaximum(float(amp["linearityMax"]))
        amplifier.setLinearityUnits("DN")
        amplifierList.append(amplifier)
    return amplifierList


def makeAmpInfoCatalog(ccd):
    """Backward compatible name."""
    return makeAmplifierList(ccd)


def makeBBoxFromList(ylist):
    """Given a list [(x0, y0), (xsize, ysize)], probably from a yaml file,
    return a BoxI.
    """
    (x0, y0), (xsize, ysize) = ylist
    return geom.BoxI(geom.PointI(x0, y0), geom.ExtentI(xsize, ysize))


def makeTransformDict(nativeSys, transformDict, plateScale):
    """Make a dictionary of TransformPoint2ToPoint2s from yaml, mapping from
    nativeSys.

    Parameters
    ----------
    nativeSys : `lsst.afw.cameraGeom.CameraSys`
    transformDict : `dict`
        A dict specifying parameters of transforms; keys are camera system
        names.
    plateScale : `lsst.geom.Angle`
        The size of a pixel in angular units/mm (e.g. 20 arcsec/mm for LSST)

    Returns
    -------
    transforms : `dict`
        A dict of `lsst.afw.cameraGeom.CameraSys` :
        `lsst.afw.geom.TransformPoint2ToPoint2`

    The resulting dict's keys are `~lsst.afw.cameraGeom.CameraSys`,
    and the values are Transforms *from* NativeSys *to* CameraSys
    """
    # As other comments note this is required, and this is one function where
    # it's assumed
    assert nativeSys == cameraGeom.FOCAL_PLANE, "Cameras with nativeSys != FOCAL_PLANE are not supported."

    resMap = {}

    for key, transform in transformDict.items():
        transformType = transform["transformType"]
        knownTransformTypes = ["affine", "radial"]
        if transformType not in knownTransformTypes:
            raise RuntimeError(
                "Saw unknown transform type for {}: {} (known types are: [{}])".format(
                    key, transform["transformType"], ", ".join(knownTransformTypes)
                )
            )

        if transformType == "affine":
            affine = geom.AffineTransform(np.array(transform["linear"]), np.array(transform["translation"]))

            transform = afwGeom.makeTransform(affine)
        elif transformType == "radial":
            # radial coefficients of the form
            #   [0, 1 (no units), C2 (mm^-1), usually 0, C4 (mm^-3), ...]
            # Radial distortion is modeled as a radial polynomial that converts
            # from focal plane radius (in mm) to field angle (in radians).
            # The provided coefficients are divided by the plate
            # scale (in radians/mm) meaning that C1 is always 1.
            radialCoeffs = np.array(transform["coeffs"])

            radialCoeffs *= plateScale.asRadians()
            transform = afwGeom.makeRadialTransform(radialCoeffs)
        else:
            raise RuntimeError(
                'Impossible condition "{}" is not in: [{}])'.format(
                    transform["transformType"], ", ".join(knownTransformTypes)
                )
            )

        resMap[cameraGeom.CameraSys(key)] = transform

    return resMap


def makeCameraFromCatalogs(
    cameraName,
    detectorConfigList,
    nativeSys,
    transformDict,
    amplifierDict,
    pupilFactoryClass=cameraGeom.pupil.PupilFactory,
):
    """Construct a Camera instance from a dictionary of
    detector name and `lsst.afw.cameraGeom.Amplifier`.

    Parameters
    ----------
    cameraName : `str`
        The name of the camera
    detectorConfigList : `list`
        A list of `lsst.afw.cameraGeom.cameraConfig.DetectorConfig`
    nativeSys : `lsst.afw.cameraGeom.CameraSys`
        The native transformation type; must be
        `lsst.afw.cameraGeom.FOCAL_PLANE`
    transformDict : `dict`
        A dict of lsst.afw.cameraGeom.CameraSys :
        `lsst.afw.geom.TransformPoint2ToPoint2`
    amplifierDict : `dict` [`str`, `lsst.afw.cameraGeom.Amplifier.Builder` ]
        A dictionary of detector name and amplifier builders.
    pupilFactoryClass : `type` [ `lsst.default afw.cameraGeom.PupilFactory`], \
            optional
        Class to attach to camera.

    Returns
    -------
    camera : `lsst.afw.cameraGeom.Camera`
        New Camera instance.

    Notes
    -----
    Copied from `lsst.afw.cameraGeom.cameraFactory` with permission and
    encouragement from Jim Bosch.
    """
    # nativeSys=FOCAL_PLANE seems to be assumed in various places in this file
    # (e.g. the definition of TAN_PIXELS), despite CameraConfig providing the
    # illusion that it's configurable.
    # Note that we can't actually get rid of the nativeSys config option
    # without breaking lots of on-disk camera configs.
    assert nativeSys == cameraGeom.FOCAL_PLANE, "Cameras with nativeSys != FOCAL_PLANE are not supported."

    focalPlaneToField = transformDict[cameraGeom.FIELD_ANGLE]

    cameraBuilder = Camera.Builder(cameraName)
    cameraBuilder.setPupilFactoryClass(pupilFactoryClass)

    # Ensure all transforms in the camera transform dict are included.
    for toSys, transform in transformDict.items():
        cameraBuilder.setTransformFromFocalPlaneTo(toSys, transform)

    for detectorConfig in detectorConfigList:
        # This should build all detector pixel -> focalPlane transforms.
        cameraGeom.addDetectorBuilderFromConfig(
            cameraBuilder, detectorConfig, amplifierDict[detectorConfig.name], focalPlaneToField
        )

        # For reasons I don't understand, some obs_ packages (e.g. HSC) set
        # nativeSys to None for their detectors (which doesn't seem to be
        # permitted by the config class!), but they really mean PIXELS. For
        # backwards compatibility we use that as the default...
        detectorNativeSys = detectorConfig.transformDict.nativeSys
        detectorNativeSys = (
            cameraGeom.PIXELS if detectorNativeSys is None else cameraGeom.CameraSysPrefix(detectorNativeSys)
        )

        # ...well, actually, it seems that we've always assumed down in C++
        # that the answer is always PIXELS without ever checking that it is.
        # So let's assert that it is, since there are hints all over this file
        # (e.g. the definition of TAN_PIXELS) that other parts of the codebase
        # have regularly made that assumption as well.  Note that we can't
        # actually get rid of the nativeSys config option without breaking
        # lots of on-disk camera configs.
        assert detectorNativeSys == cameraGeom.PIXELS, "Detectors with nativeSys != PIXELS are not supported."
        detectorNativeSys = cameraGeom.CameraSys(detectorNativeSys, detectorConfig.name)

    return cameraBuilder.finish()
