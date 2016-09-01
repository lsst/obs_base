import argparse

keys_added = {
    'datasets': ('forcedPhotCcd_config', 'src_schema', 'diffsources_schema', 'deep_assembleCoadd_config', 'deepCoadd_calexp_background', 'deepCoadd_peak_schema', 'eups_versions', 'deepCoadd_src', 'deep_safeClipAssembleCoadd_config', 'flat_config', 'characterizeImage_config', 'deepCoadd_icSrc_schema', 'fringe_config', 'bias_config', 'deepCoadd_srcMatchFull', 'deepCoadd_srcMatch', 'deepCoadd_forced_metadata', 'isr_config', 'deepCoadd_forced_src_schema', 'deep_makeCoaddTempExp_metadata', 'deepCoadd_psf', 'deepCoadd_ref_schema', 'deepCoadd_src_schema', 'deepCoaddId', 'mergeCoaddMeasurements_config', 'deep_assembleCoadd_metadata', 'deepCoadd_meas', 'deepCoadd_forced_config', 'deepCoadd_forced_src', 'mergeCoaddDetections_config', 'deep_makeSkyMap_metadata', 'measureCoaddSources_metadata', 'mergeCoaddMeasurements_metadata', 'transformSrcMeasurement_config', 'processCcd_config', 'deepCoadd_mergeDet_schema', 'brightObjectMask', 'deep_safeClipAssembleCoadd_metadata', 'transformed_src_schema', 'processStack_config', 'deepMergedCoaddId', 'stackExposureId_bits', 'packages', 'detectCoaddSources_metadata', 'deep_makeCoaddTempExp_config', 'deepCoaddId_bits', 'deepCoadd_ref', 'forced_src_schema', 'deepCoadd_diffsrc', 'calibrated_src_schema', 'deepCoadd_det_schema', 'stackExposureId', 'deepCoadd_apCorr', 'solvetansip_config', 'dark_config', 'IngestIndexedReferenceTask_config', 'calibrate_config', 'deep_processCoadd_config', 'multiBandDriver_config', 'measureCoaddSources_config', 'deepCoadd_icMatch', 'deep_makeDiscreteSkyMap_metadata', 'deepCoadd_mergeDet', 'ccdExposureId', 'deep_coadd_metadata', 'forcedCcd_config', 'deepCoadd_meas_schema', 'deepCoadd_skyMap', 'processFocus_config', 'deepCoadd_det', 'detectCoaddSources_config', 'coaddDriver_config', 'deepCoadd_multibandReprocessing', 'deepCoadd_icSrc', 'Mosaic_config', 'deepMergedCoaddId_bits', 'forcedCoadd_config', 'singleFrameDriver_config', 'deep_coadd_config', 'mergeCoaddDetections_metadata', 'icSrc_schema', 'ccdExposureId_bits', 'processFocusSweep_config', 'deepCoadd_initPsf', 'deep_processCoadd_metadata'),
    'exposures': ('deepCoadd_diff','deepCoadd','deepCoadd_depth','deepCoadd_calexp','deepCoadd_bg', 'deepCoadd_bgRef'),
    'calibrations': (),
    'images': (),
}

def readSection(filename):
    try:
        fin = open(filename, "r")
    except:
        return None
    result = {}
    lines = fin.readlines()
    bset = None
    for line in lines:
        if len(line) > 0 and not line[0] == ' ' and line.find(":") > 0:
            if bset:
                result[name] = bset
            bset = line
            name = line[:line.find(":")]
        else:
            bset = bset + line 
    fin.close()           
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--name", type=str, help="mapper name", default=None) 
    parser.add_argument("-r", "--reftype", type=str, help="path to ref", default=None) 
    args = parser.parse_args()
    name = args.name
    reftype = args.reftype
    if reftype == None:
        reftype = name
        fout = open(name + ".selfdiffs", "w")
    else:
        fout = open(name + ".%sdiffs"%(reftype,), "w")
    errors = 0
    types = ("datasets", "images", "exposures", "calibrations")
    for type in types:
        ref = readSection("refs/%s.%s"%(reftype, type,))
        comp = readSection("ones/%s.%s"%(name, type))
        if ref is None and comp is None:
            continue
        if ref is None and not comp is None:
            fout.write("Section %s doesn't appear in reference of %s"%(type, reftype))
            continue
        if not ref is None and comp is None:
            fout.write("Section %s doesn't appear in current of %s"%(type,name))
            continue
        fout.write("%s.%s\n"%(name,type))
        for key in ref.keys():
            if not key in comp.keys():
                fout.write("    %s.%s key of reference %s not present.\n"%(type, key, reftype))
                errors = errors + 1
            
            elif not ref[key] == comp[key] and (key in keys_added[type]):
                fout.write("    %s.%s different\n"%(type, key,))
                refs = ref[key].split("\n")
                comps = comp[key].split("\n")
                
                for i in range(max(len(refs), len(comps))):
                    if i > (len(comps)-1):
                        comps.append("")
                    if i > (len(refs)-1):
                        refs.append("")
                    if not refs[i] == comps[i]:
                        fout.write("    different:" + comps[i] + "\n")
                        fout.write("    different:" + refs[i] + " in %s reference\n"%(reftype))
                errors = errors + 1
        # see if key in current was in the reference
        for key in comp.keys():
            if not key in ref.keys():
                if name == reftype:
                    if not key in keys_added[type]:
                        fout.write("   %s.%s in %s was added, and not by butlerUtils\n"%(type, key, name))
                else:
                    fout.write("    %s.%s in %s is not in %s\n"%(type, key, name, reftype))
