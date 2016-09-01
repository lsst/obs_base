import argparse

rmKeys = ('deepCoadd_calexp_background', 'deepCoadd_det', 'deepCoadd_mergeDet', 'deepCoadd_meas', 'deepCoadd_ref', 'deepCoadd_forced_src', 'detectCoaddSources_metadata', 'mergeCoaddDetections_metadata', 'measureCoaddSources_metadata', 'mergeCoaddMeasurements_metadata', 'deepCoadd_forced_metadata', 'isr_config', 'transformSrcMeasurement_config', 'singleFrameDriver_config', 'Mosaic_config', 'processStack_config', 'deep_makeCoaddTempExp_config', 'deep_assembleCoadd_config', 'deep_safeClipAssembleCoadd_config', 'deep_coadd_config', 'deep_processCoadd_config', 'bias_config', 'dark_config', 'flat_config', 'fringe_config', 'solvetansip_config', 'coaddDriver_config', 'forcedCoadd_config', 'forcedCcd_config', 'deepCoadd_forced_config', 'forcedPhotCcd_config', 'processFocus_config', 'processFocusSweep_config', 'detectCoaddSources_config', 'mergeCoaddDetections_config', 'measureCoaddSources_config', 'mergeCoaddMeasurements_config', 'multiBandDriver_config',)

pafs = ("DecamMapper.paf","HscMapper.paf","LsstSimMapper.paf","MegacamMapper.paf","SdssMapper.paf","SuprimecamMapper.paf","testMapper.paf")

if __name__ == "__main__":
    keys = {}
    for file in pafs:
        fin = open(file, "r")
        fout = open(file + ".rem", "w")
        lines = fin.readlines()
        bset = None
        for line in lines:
            if bset:
                bset = bset + line 
                if len(line) > 0 and  line.find(' }') > 0:
                    if not name in rmKeys:
                        fout.write(bset)
                    keys[name] = None
                    bset = None
                continue
             
            if line.find('{') > 0 and line.find(":") > 0:
                bset = line
                name = line[4: line.find(":")]
            else:
                fout.write(line)
        fin.close()           
        fout.close()           
