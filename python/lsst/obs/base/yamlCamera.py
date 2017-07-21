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
from builtins import range
import yaml

import numpy as np
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.geom as afwGeom
from lsst.afw.table import AmpInfoCatalog, AmpInfoTable
from lsst.afw.cameraGeom.cameraFactory import makeDetector

class YamlCamera(cameraGeom.Camera):
    """The Commissioning Camera (comCam)
    """

    def __init__(self, cameraYamlFile):
        """Construct a Camera
        """
        with file(cameraYamlFile) as fd:
            cameraParams = yaml.load(fd, Loader=yaml.Loader)

        plateScale = afwGeom.Angle(cameraParams["plateScale"], afwGeom.arcseconds)
        radialCoeffs = np.array(cameraParams["radialCoeffs"])/plateScale.asRadians()
        focalPlaneToPupil = afwGeom.RadialXYTransform(radialCoeffs)
        pupilToFocalPlane = afwGeom.InvertedXYTransform(focalPlaneToPupil)
        cameraTransformMap = cameraGeom.CameraTransformMap(cameraGeom.FOCAL_PLANE,
                                                           {cameraGeom.PUPIL: pupilToFocalPlane})
        detectorList = self._makeDetectorList(cameraParams["CCDs"], pupilToFocalPlane, plateScale)
        cameraGeom.Camera.__init__(self, cameraParams["name"], detectorList, cameraTransformMap)

    def _makeDetectorList(self, ccdParams, focalPlaneToPupil, plateScale):
        """!Make a list of detectors
        @param[in] ccdParams  Dict of YAML descriptions of CCDs
        @param[in] focalPlaneToPupil  lsst.afw.geom.XYTransform from FOCAL_PLANE to PUPIL coordinates
        @param[in] plateScale  plate scale, in angle on sky/mm
        @return a list of detectors (lsst.afw.cameraGeom.Detector)
        """
        detectorList = []
        detectorConfigList = self._makeDetectorConfigList(ccdParams)
        for ccd, detectorConfig in zip(ccdParams.values(), detectorConfigList):
            ampInfoCatalog = self._makeAmpInfoCatalog(ccd)
            detector = makeDetector(detectorConfig, ampInfoCatalog, focalPlaneToPupil)
            detectorList.append(detector)
        return detectorList

    def _makeDetectorConfigList(self, ccdParams):
        """!Make a list of detector configs

        @return a list of detector configs (lsst.afw.cameraGeom.DetectorConfig)
        """
        detectorConfigs = []
        for name, ccd in ccdParams.items():
            detectorConfig = cameraGeom.DetectorConfig()
            detectorConfigs.append(detectorConfig)

            detectorConfig.name = name
            detectorConfig.id = ccd['id']
            detectorConfig.serial = ccd['serial']
            detectorConfig.detectorType = ccd['detectorType']
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
        
        return detectorConfigs

    @staticmethod
    def _makeBBoxFromList(ylist):
        """Given a list [(x0, y0), (xsize, ysize)], probably from a yaml file, return a BoxI
            """
        (x0, y0), (xsize, ysize) = ylist
        return  afwGeom.BoxI(afwGeom.PointI(x0, y0), afwGeom.ExtentI(xsize, ysize))

    def _makeAmpInfoCatalog(self, ccd):
        """Construct an amplifier info catalog
        """
        # Much of this will need to be filled in when we know it.
        assert len(ccd['amplifiers']) > 0
        amp = ccd['amplifiers'].values()[0]

        rawBBox = self._makeBBoxFromList(amp['rawBBox']) # total in file
        xRawExtent, yRawExtent = rawBBox.getDimensions()
        
        from lsst.afw.table import LL, LR, UL, UR
        readCorners = dict(LL = LL, LR = LR, UL = UL, UR = UR)

        schema = AmpInfoTable.makeMinimalSchema()

        linThreshKey = schema.addField('linearityThreshold', type=float)
        linMaxKey = schema.addField('linearityMaximum', type=float)
        linUnitsKey = schema.addField('linearityUnits', type=str, size=9)
        hduKey = schema.addField('hdu', type=np.int32)
        # end placeholder
        self.ampInfoDict = {}
        ampCatalog = AmpInfoCatalog(schema)
        for name, amp in sorted(ccd['amplifiers'].items(), key=lambda x : x[1]['hdu']):
            record = ampCatalog.addNew()
            record.setName(name)
            record.set(hduKey, amp['hdu'])

            ix, iy = amp['ixy']
            perAmpData = amp['perAmpData']
            if perAmpData:
                x0, y0 = 0, 0           # origin of data within each amp image
            else:
                x0, y0 = ix*xRawExtent, iy*yRawExtent

            rawDataBBox = self._makeBBoxFromList(amp['rawDataBBox']) # Photosensitive area
            xDataExtent, yDataExtent = rawDataBBox.getDimensions()
            record.setBBox(afwGeom.BoxI(
                afwGeom.PointI(ix*xDataExtent, iy*yDataExtent), rawDataBBox.getDimensions()))

            rawBBox = self._makeBBoxFromList(amp['rawBBox'])
            rawBBox.shift(afwGeom.ExtentI(x0, y0))
            record.setRawBBox(rawBBox)
            
            rawDataBBox = self._makeBBoxFromList(amp['rawDataBBox'])
            rawDataBBox.shift(afwGeom.ExtentI(x0, y0))
            record.setRawDataBBox(rawDataBBox)

            rawSerialOverscanBBox = self._makeBBoxFromList(amp['rawSerialOverscanBBox'])
            rawSerialOverscanBBox.shift(afwGeom.ExtentI(x0, y0))
            record.setRawHorizontalOverscanBBox(rawSerialOverscanBBox)

            rawParallelOverscanBBox = self._makeBBoxFromList(amp['rawParallelOverscanBBox'])
            rawParallelOverscanBBox.shift(afwGeom.ExtentI(x0, y0))
            record.setRawVerticalOverscanBBox(rawParallelOverscanBBox)

            rawSerialPrescanBBox = self._makeBBoxFromList(amp['rawSerialPrescanBBox'])
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
