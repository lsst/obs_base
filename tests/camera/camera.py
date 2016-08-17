# based on obs_cfht's camera, but just one sensor
import lsst.afw.cameraGeom.cameraConfig
assert type(config) == lsst.afw.cameraGeom.cameraConfig.CameraConfig, 'config is of type %s.%s instead of lsst.afw.cameraGeom.cameraConfig.CameraConfig' % (
    type(config).__module__, type(config).__name__)
config.plateScale = 13.7
config.transformDict.nativeSys = 'FocalPlane'
config.transformDict.transforms = {}
config.transformDict.transforms['Pupil'] = lsst.afw.geom.transformConfig.TransformConfig()
config.transformDict.transforms['Pupil'].transform['multi'].transformDict = None
config.transformDict.transforms['Pupil'].transform['affine'].translation = [0.0, 0.0]
config.transformDict.transforms['Pupil'].transform['affine'].linear = [1.0, 0.0, 0.0, 1.0]
config.transformDict.transforms['Pupil'].transform['radial'].coeffs = None
import lsst.afw.geom.xyTransformFactory
config.transformDict.transforms['Pupil'].transform['inverted'].transform.retarget(
    target=lsst.afw.geom.xyTransformFactory.makeRadialXYTransform, ConfigClass=lsst.afw.geom.xyTransformFactory.RadialXYTransformConfig)
config.transformDict.transforms['Pupil'].transform[
    'inverted'].transform.coeffs = [0.0, 14805.4, 13619.3, 426637.0]
config.transformDict.transforms['Pupil'].transform.name = 'inverted'
config.detectorList = {}
config.detectorList[0] = lsst.afw.cameraGeom.cameraConfig.DetectorConfig()
config.detectorList[0].bbox_y0 = 0
config.detectorList[0].bbox_y1 = 4611
config.detectorList[0].bbox_x1 = 2047
config.detectorList[0].bbox_x0 = 0
config.detectorList[0].name = 'ccd00'
config.detectorList[0].pixelSize_x = 0.0135
config.detectorList[0].transformDict.nativeSys = None
config.detectorList[0].transformDict.transforms = None
config.detectorList[0].refpos_x = 1023.5
config.detectorList[0].refpos_y = 2305.5
config.detectorList[0].pixelSize_y = 0.0135
config.detectorList[0].detectorType = 0
config.detectorList[0].offset_x = -114.399
config.detectorList[0].offset_y = 99.46125
config.detectorList[0].transposeDetector = None
config.detectorList[0].yawDeg = 180.0
config.detectorList[0].rollDeg = 0.0
config.detectorList[0].serial = '834175'
config.detectorList[0].pitchDeg = 0.0
config.detectorList[0].id = 0
config.radialCoeffs = None
config.name = 'Trivial Camera'
