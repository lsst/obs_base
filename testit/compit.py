import argparse
from shutil import copyfile
#   dataset which have tract or patch replacements
keys = ['goodSeeingCoadd_src', 'goodSeeing_srcMatch', 'deepCoadd_ref', 'deep_assembleCoadd_metadata', 'goodSeeingCoadd_icMatch', 'deepCoadd_diffsrc', 'calibrated_src', 'deepCoadd_multibandReprocessing', 'deepCoadd_meas', 'deepCoadd_multiModelfits', 'goodSeeingCoadd_depth', 'goodSeeingCoadd_apCorr', 'deepCoadd_icSrc', 'diffsources', 'deepCoadd_forced_src', 'deepCoadd_tempExp_diffsrc', 'deepCoadd_calexp_background', 'chiSquared_coadd_metadata', 'goodSeeing_makeCoaddTempExp_metadata', 'chiSquaredCoadd_apCorr', 'measureCoaddSources_metadata', 'mergeCoaddMeasurements_metadata', 'chiSquared_assembleCoadd_metadata', 'mergeCoaddDetections_metadata', 'goodSeeing_assembleCoadd_metadata', 'deepCoadd_extract', 'detectCoaddSources_metadata', 'deepCoadd_calexp_detBackground', 'deepCoadd_src', 'goodSeeingCoadd_calexpBackground', 'brightObjectMask', 'deep_safeClipAssembleCoadd_metadata', 'goodSeeing_coadd_metadata', 'deepCoadd_modelfits', 'goodSeeingCoadd_icSrc', 'goodSeeing_processCoadd_metadata', 'chiSquaredCoadd_icSrc', 'forcedPhotCcd_metadata', 'processStack_metadata', 'deepCoadd_icMatch', 'deep_measureMulti_metadata', 'deepCoadd_calexpBackground', 'deepCoadd_depth', 'forced_src', 'chiSquaredCoadd_icMatch', 'deepCoadd_srcMatchFull', 'deepCoadd_srcMatch', 'deepCoadd_forced_metadata', 'deepCoadd_mergeDet', 'deepCoadd_apCorr', 'deepCoad_calexpBackground', 'deep_coadd_metadata', 'forcedPhot_metadata', 'deep_forcedPhotCoadd_metadata', 'deepCoadd_initPsf', 'deep_makeCoaddTempExp_metadata', 'deepCoadd_psf', 'diffpsf', 'deep_measureCoadd_metadata', 'chiSquared_processCoadd_metadata', 'chiSquaredCoadd_src', 'deepCoadd_det', 'deep_processCoadd_metadata','forcedPhotCcd_config', 'src_schema', 'deepDiff_kernelSrc', 'badSource_schema', 'plotSeeingRough', 'diffsources_schema', 'deep_assembleCoadd_config', 'goodSeeingDiff_diaSrc_schema', 'sourceAssoc_config', 'deepCoadd_peak_schema', 'deepDiff_diaSrc', 'tsField', 'eups_versions', 'plotEllipseMap', 'goodSeeingDiff_metadata', 'source', 'deep_safeClipAssembleCoadd_config', 'flat_config', 'chiSquared_makeSkyMap_config', 'plotPsfModelGrid', 'chiSquared_makeSkyMap_metadata', 'focusSweepPlot', 'psf', 'deepCoadd_icSrc_schema', 'fringe_config', 'ccdExposureId', 'logDir', 'singleFrameDriver_metadata', 'deepCoadd_modelfits_schema', 'isr_config', 'ossThumb', 'chiSquared_processCoadd_config', 'deepCoadd_multiModelfits_schema', 'deepCoadd_forced_src_schema', 'goodSeeing_makeCoaddTempExp_config', 'plotEllipticityMap', 'goodSeeing_coadd_config', 'deepCoadd_ref_schema', 'deepCoadd_src_schema', 'deepCoadd_meas_schema', 'mergeCoaddMeasurements_config', 'invalidSource_schema', 'eimage', 'stack_config', 'transformSrcMeasurement_metadata', 'multiband_config', 'isr_metadata', 'stackExposureId', 'transformed_src', 'cal_ref_cat', 'goodSeeingCoadd_skyMap', 'goodSeeingCoadd_src_schema', 'deep_forcedPhotCoadd_config', 'deepCoadd_forced_config', 'deepDiff_diaSrc_schema', 'mergeCoaddDetections_config', 'deep_makeSkyMap_metadata', 'test_config', 'goodSeeingCoaddId', 'Mosaic_metadata', 'deepCoadd_mergeDet_schema', 'forced_metadata', 'transformed_src_schema', 'goodSeeingCoaddId_bits', 'forcedPhotCoadd_config', 'badSourceHist', 'chiSquaredCoaddId', 'forced_schema', 'processStack_config', 'deepMergedCoaddId', 'stackExposureId_bits', 'goodSeeing_makeSkyMap_config', 'fitsFwhmGrid', 'src', 'processFocus_metadata', 'apCorr', 'deep_makeCoaddTempExp_config', 'goodSeeingDiff_kernelSrc', 'calexpBackground', 'goodSeeingDiff_differenceExp', 'plotEllipticityGrid', 'srcMatch', 'source_schema', 'deepCoadd_forced_schema', 'plotMagHist', 'refcat', 'deepCoaddId_bits', 'forced_src_schema', 'icSrc', 'goodSeeing_makeSkyMap_metadata', 'processCcdDecam_config', 'calibrated_src_schema', 'deepCoadd_det_schema', 'transformSrcMeasurement_config', 'deepDiff_differenceExp', 'solvetansip_config', 'dark_config', 'processEimage_metadata', 'log', 'test_metadata', 'icMatchFull', 'characterizeImage_metadata', 'deepDiff_config', 'modelfits_schema', 'processCcd_metadata', 'deep_measureMulti_config', 'deep_processCoadd_config', 'deep_makeSkyMap_config', 'ampExposureId', 'psField', 'tableSeeingGrid', 'multiBandDriver_config', 'measureCoaddSources_config', 'sourceHist', 'deep_makeDiscreteSkyMap_metadata', 'processCcdDecam_metadata', 'calexpThumb', 'flattenedThumb', 'bias_config', 'focusPlot', 'chiSquaredCoadd_src_schema', 'forcedCcd_config', 'deepCoaddId', 'goodSeeingCoadd_icSrc_schema', 'invalidSource', 'processEimage_config', 'deepCoadd_skyMap', 'fitsEllipticityGrid', 'crDiffimSrc', 'chiSquaredCoaddId_bits', 'object_schema', 'processFocus_config', 'badSource', 'detectCoaddSources_config', 'coaddDriver_config', 'chiSquared_assembleCoadd_config', 'deepDiff_matchedExp', 'fitsEllPaGrid', 'fitsPsfSrcGrid', 'Mosaic_config', 'deepMergedCoaddId_bits', 'forcedCoadd_config', 'icExpBackground', 'deep_coadd_config', 'asTrans', 'goodSeeingDiff_diaSrc', 'crDiffimSrc_schema', 'deepDiff_metadata', 'calibrate_metadata', 'other_photo_astro_ref', 'tableSeeingMap', 'fitsPsfModelGrid', 'goodSeeing_processCoadd_config', 'goodSeeingDiff_config', 'icSrc_schema', 'plotEllipseGrid', 'ccdExposureId_bits', 'chiSquaredCoadd_icSrc_schema', 'warppsf', 'object', 'processFocusSweep_config', 'plotSeeingRobust', 'measureCcd_metadata', 'sourceAssoc_metadata', 'goodSeeingDiff_matchedExp', 'measureCcd_config', 'deep_measureCoadd_config', 'plotPsfSrcGrid', 'srcMatchFull', 'icMatch', 'singleFrameDriver_config', 'ampExposureId_bits', 'modelfits', 'chiSquared_coadd_config', 'goodSeeing_assembleCoadd_config', 'forced_config', 'plotFwhmGrid', 'chiSquaredCoadd_skyMap', 'plotSeeingMap']
pafs = ("HscMapper.paf","SuprimecamMapper.paf", "DecamMapper.paf", "LsstSimMapper.paf","MegacamMapper.paf","SdssMapper.paf","testMapper.paf")

if __name__ == "__main__":
  fout2 = open("all.summary", "w")
  for key in keys:
    fout = open(key + ".summary", "w")
    bfirst = None
    bdiffers = False
    which = []
    for file in pafs:
        fout.write("--------" + file + "--------\n")
        fin = open(file, "r")
        print "Doing ", file
        lines = fin.readlines()
        bset = None
        for line in lines:
            if bset:
                #  In comparing, don't include lines which have raw_skyType but don't use it
                if not (line.find("raw_skyTile") >= 0 and not bset.find("(skyTile)") >= 0):
                    bset = bset + line 
                if len(line) > 0 and  line.find(' }') > 0:
                    if name == key:
                        if bfirst is None:
                            bfirst = bset
                        if not bfirst.strip() == bset.strip():
                            bdiffers = True
                            fout.write("Differs\n");
                            which.append(file[0].lower())
                        else:
                            which.append(file[0])
                        fout.write(bset)
                    bset = None
            if line.find('{') > 0 and line.find(":") > 0:
                bset = line
                name = line[4: line.find(":")]
                bhasparam = False
            if bset:
                if line.find("tract") > 0 or line.find("patch") > 0:
                     bhasparam = True
        fin.close()           
    fout.close()
    if bdiffers:
        copyfile(key + ".summary", key + ".differs") 
    else:
        copyfile(key + ".summary", key + ".nodiffs")
    fout2.write("%s: %s\n"%(key, "".join(which)))
  fout2.close() 
