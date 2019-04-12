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
import yaml

import numpy as np
import lsst.afw.cameraGeom as cameraGeom
import lsst.geom as geom
import lsst.afw.geom as afwGeom
from lsst.afw.table import AmpInfoCatalog, AmpInfoTable
from lsst.afw.cameraGeom.cameraFactory import makeDetectorData

__all__ = ["makeCamera"]


def makeCamera(cameraFile):
    """An imaging camera (e.g. the LSST 3Gpix camera)

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

    ampInfoCatDict = {}
    for ccdName, ccdValues in ccdParams.items():
        ampInfoCatDict[ccdName] = makeAmpInfoCatalog(ccdValues)

    return makeCameraFromCatalogs(cameraName, detectorConfigList, nativeSys, transforms, ampInfoCatDict)


def makeDetectorConfigList(ccdParams):
    """Make a list of detector configs

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
        detectorConfig.id = ccd['id']
        detectorConfig.serial = ccd['serial']
        detectorConfig.detectorType = ccd['detectorType']
        if 'physicalType' in ccd:
            detectorConfig.physicalType = ccd['physicalType']
        # This is the orientation we need to put the serial direction along the x-axis
        detectorConfig.bbox_x0, detectorConfig.bbox_y0 = ccd['bbox'][0]
        detectorConfig.bbox_x1, detectorConfig.bbox_y1 = ccd['bbox'][1]
        detectorConfig.pixelSize_x, detectorConfig.pixelSize_y = ccd['pixelSize']
        detectorConfig.transformDict.nativeSys = ccd['transformDict']['nativeSys']
        transforms = ccd['transformDict']['transforms']
        detectorConfig.transformDict.transforms = None if transforms == 'None' else transforms
        detectorConfig.refpos_x, detectorConfig.refpos_y = ccd['refpos']
        detectorConfig.offset_x, detectorConfig.offset_y = ccd['offset']
        detectorConfig.transposeDetector = ccd['transposeDetector']
        detectorConfig.pitchDeg = ccd['pitch']
        detectorConfig.yawDeg = ccd['yaw']
        detectorConfig.rollDeg = ccd['roll']
        if 'crosstalk' in ccd:
            detectorConfig.crosstalk = ccd['crosstalk']

    return detectorConfigs


def makeAmpInfoCatalog(ccd):
    """Construct an amplifier info catalog
    """
    # Much of this will need to be filled in when we know it.
    assert len(ccd['amplifiers']) > 0
    amp = list(ccd['amplifiers'].values())[0]

    rawBBox = makeBBoxFromList(amp['rawBBox'])  # total in file
    xRawExtent, yRawExtent = rawBBox.getDimensions()

    from lsst.afw.table import LL, LR, UL, UR
    readCorners = dict(LL=LL, LR=LR, UL=UL, UR=UR)

    schema = AmpInfoTable.makeMinimalSchema()

    linThreshKey = schema.addField('linearityThreshold', type=float)
    linMaxKey = schema.addField('linearityMaximum', type=float)
    linUnitsKey = schema.addField('linearityUnits', type=str, size=9)
    hduKey = schema.addField('hdu', type=np.int32)
    # end placeholder

    ampCatalog = AmpInfoCatalog(schema)
    for name, amp in sorted(ccd['amplifiers'].items(), key=lambda x: x[1]['hdu']):
        record = ampCatalog.addNew()
        record.setName(name)
        record.set(hduKey, amp['hdu'])

        ix, iy = amp['ixy']
        perAmpData = amp['perAmpData']
        if perAmpData:
            x0, y0 = 0, 0           # origin of data within each amp image
        else:
            x0, y0 = ix*xRawExtent, iy*yRawExtent

        rawDataBBox = makeBBoxFromList(amp['rawDataBBox'])  # Photosensitive area
        xDataExtent, yDataExtent = rawDataBBox.getDimensions()
        record.setBBox(afwGeom.BoxI(
            afwGeom.PointI(ix*xDataExtent, iy*yDataExtent), rawDataBBox.getDimensions()))

        rawBBox = makeBBoxFromList(amp['rawBBox'])
        rawBBox.shift(afwGeom.ExtentI(x0, y0))
        record.setRawBBox(rawBBox)

        rawDataBBox = makeBBoxFromList(amp['rawDataBBox'])
        rawDataBBox.shift(afwGeom.ExtentI(x0, y0))
        record.setRawDataBBox(rawDataBBox)

        rawSerialOverscanBBox = makeBBoxFromList(amp['rawSerialOverscanBBox'])
        rawSerialOverscanBBox.shift(afwGeom.ExtentI(x0, y0))
        record.setRawHorizontalOverscanBBox(rawSerialOverscanBBox)

        rawParallelOverscanBBox = makeBBoxFromList(amp['rawParallelOverscanBBox'])
        rawParallelOverscanBBox.shift(afwGeom.ExtentI(x0, y0))
        record.setRawVerticalOverscanBBox(rawParallelOverscanBBox)

        rawSerialPrescanBBox = makeBBoxFromList(amp['rawSerialPrescanBBox'])
        rawSerialPrescanBBox.shift(afwGeom.ExtentI(x0, y0))
        record.setRawPrescanBBox(rawSerialPrescanBBox)

        if perAmpData:
            record.setRawXYOffset(afwGeom.Extent2I(ix*xRawExtent, iy*yRawExtent))
        else:
            record.setRawXYOffset(afwGeom.Extent2I(0, 0))

        record.setReadoutCorner(readCorners[amp['readCorner']])
        record.setGain(amp['gain'])
        record.setReadNoise(amp['readNoise'])
        record.setSaturation(amp['saturation'])
        record.setHasRawInfo(True)
        # flip data when assembling if needs be (e.g. data from the serial at the top of a CCD)
        flipX, flipY = amp.get("flipXY")

        record.setRawFlipX(flipX)
        record.setRawFlipY(flipY)
        # linearity placeholder stuff
        record.setLinearityCoeffs([float(val) for val in amp['linearityCoeffs']])
        record.setLinearityType(amp['linearityType'])
        record.set(linThreshKey, float(amp['linearityThreshold']))
        record.set(linMaxKey, float(amp['linearityMax']))
        record.set(linUnitsKey, "DN")
    return ampCatalog


def makeBBoxFromList(ylist):
    """Given a list [(x0, y0), (xsize, ysize)], probably from a yaml file,
    return a BoxI
    """
    (x0, y0), (xsize, ysize) = ylist
    return afwGeom.BoxI(afwGeom.PointI(x0, y0), afwGeom.ExtentI(xsize, ysize))


def makeTransformDict(nativeSys, transformDict, plateScale):
    """Make a dictionary of TransformPoint2ToPoint2s from yaml, mapping from nativeSys

    Parameters
    ----------
    nativeSys : `lsst.afw.cameraGeom.CameraSys`
    transformDict : `dict`
        A dict specifying parameters of transforms; keys are camera system names.
    plateScale : `lsst.geom.Angle`
        The size of a pixel in angular units/mm (e.g. 20 arcsec/mm for LSST)

    Returns
    -------
    transforms : `dict`
        A dict of `lsst.afw.cameraGeom.CameraSys` : `lsst.afw.geom.TransformPoint2ToPoint2`

    The resulting dict's keys are `~lsst.afw.cameraGeom.CameraSys`,
    and the values are Transforms *from* NativeSys *to* CameraSys
    """
    # As other comments note this is required, and this is one function where it's assumed
    assert nativeSys == cameraGeom.FOCAL_PLANE, "Cameras with nativeSys != FOCAL_PLANE are not supported."

    resMap = dict()

    for key, transform in transformDict.items():
        transformType = transform["transformType"]
        knownTransformTypes = ["affine", "radial"]
        if transformType not in knownTransformTypes:
            raise RuntimeError("Saw unknown transform type for %s: %s (known types are: [%s])" % (
                key, transform["transformType"], ", ".join(knownTransformTypes)))

        if transformType == "affine":
            affine = geom.AffineTransform(np.array(transform["linear"]),
                                          np.array(transform["translation"]))

            transform = afwGeom.makeTransform(affine)
        elif transformType == "radial":
            # radial coefficients of the form [0, 1 (no units), C2 (rad), usually 0, C3 (rad^2), ...]
            # Radial distortion is modeled as a radial polynomial that converts from focal plane radius
            # (in mm) to field angle (in radians). The provided coefficients are divided by the plate
            # scale (in radians/mm) meaning that C1 is always 1.
            radialCoeffs = np.array(transform["coeffs"])

            radialCoeffs *= plateScale.asRadians()
            transform = afwGeom.makeRadialTransform(radialCoeffs)
        else:
            raise RuntimeError("Impossible condition \"%s\" is not in: [%s])" % (
                transform["transformType"], ", ".join(knownTransformTypes)))

        resMap[cameraGeom.CameraSys(key)] = transform

    return resMap


def makeCameraFromCatalogs(cameraName, detectorConfigList, nativeSys, transformDict, ampInfoCatDict,
                           pupilFactoryClass=cameraGeom.pupil.PupilFactory):
    """Construct a Camera instance from a dictionary of
       detector name : `lsst.afw.table.ampInfo.AmpInfoCatalog`

    Parameters
    ----------
    cameraName : `str`
        The name of the camera
    detectorConfig : `list`
        A list of `lsst.afw.cameraGeom.cameraConfig.DetectorConfig`
    nativeSys : `lsst.afw.cameraGeom.CameraSys`
        The native transformation type; must be `lsst.afw.cameraGeom.FOCAL_PLANE`
    transformDict : `dict`
        A dict of lsst.afw.cameraGeom.CameraSys : `lsst.afw.geom.TransformPoint2ToPoint2`
    ampInfoCatDict : `dict`
        A dictionary of detector name :
                           `lsst.afw.table.ampInfo.AmpInfoCatalog`
    pupilFactoryClass : `type`, optional
        Class to attach to camera;
             `lsst.default afw.cameraGeom.PupilFactory`

    Returns
    -------
    camera : `lsst.afw.cameraGeom.Camera`
        New Camera instance.

    Notes
    ------
    Copied from `lsst.afw.cameraGeom.cameraFactory` with permission and encouragement
    from Jim Bosch
    """

    # nativeSys=FOCAL_PLANE seems to be assumed in various places in this file
    # (e.g. the definition of TAN_PIXELS), despite CameraConfig providing the
    # illusion that it's configurable.
    # Note that we can't actually get rid of the nativeSys config option
    # without breaking lots of on-disk camera configs.
    assert nativeSys == cameraGeom.FOCAL_PLANE, "Cameras with nativeSys != FOCAL_PLANE are not supported."

    focalPlaneToField = transformDict[cameraGeom.FIELD_ANGLE]
    transformMapBuilder = cameraGeom.TransformMap.Builder(nativeSys)
    transformMapBuilder.connect(transformDict)

    # First pass: build a list of all Detector ctor kwargs, minus the
    # transformMap (which needs information from all Detectors).
    detectorData = []
    for detectorConfig in detectorConfigList:

        # Get kwargs that could be used to construct each Detector
        # if we didn't care about giving each of them access to
        # all of the transforms.
        thisDetectorData = makeDetectorData(
            detectorConfig=detectorConfig,
            ampInfoCatalog=ampInfoCatDict[detectorConfig.name],
            focalPlaneToField=focalPlaneToField,
        )

        # Pull the transforms dictionary out of the data dict; we'll replace
        # it with a TransformMap argument later.
        thisDetectorTransforms = thisDetectorData.pop("transforms")

        # Save the rest of the Detector data dictionary for later
        detectorData.append(thisDetectorData)

        # For reasons I don't understand, some obs_ packages (e.g. HSC) set
        # nativeSys to None for their detectors (which doesn't seem to be
        # permitted by the config class!), but they really mean PIXELS. For
        # backwards compatibility we use that as the default...
        detectorNativeSys = detectorConfig.transformDict.nativeSys
        detectorNativeSys = (cameraGeom.PIXELS if detectorNativeSys is None else
                             cameraGeom.CameraSysPrefix(detectorNativeSys))

        # ...well, actually, it seems that we've always assumed down in C++
        # that the answer is always PIXELS without ever checking that it is.
        # So let's assert that it is, since there are hints all over this file
        # (e.g. the definition of TAN_PIXELS) that other parts of the codebase
        # have regularly made that assumption as well.  Note that we can't
        # actually get rid of the nativeSys config option without breaking
        # lots of on-disk camera configs.
        assert detectorNativeSys == cameraGeom.PIXELS, \
            "Detectors with nativeSys != PIXELS are not supported."
        detectorNativeSys = cameraGeom.CameraSys(detectorNativeSys, detectorConfig.name)

        # Add this detector's transform dict to the shared TransformMapBuilder
        transformMapBuilder.connect(detectorNativeSys, thisDetectorTransforms)

    # Now that we've collected all of the Transforms, we can finally build the
    # (immutable) TransformMap.
    transformMap = transformMapBuilder.build()

    # Second pass through the detectorConfigs: actually make Detector instances
    detectorList = [cameraGeom.Detector(transformMap=transformMap, **kw) for kw in detectorData]

    return cameraGeom.Camera(cameraName, detectorList, transformMap, pupilFactoryClass)
