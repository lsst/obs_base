import argparse

rmKeys = ('diffsources_schema', 'eups_versions', 'deepCoadd_src', 'deepCoadd_src_schema', 'deepCoadd_icSrc', 'deepCoadd_icSrc_schema', 'deepCoadd_srcMatch', 'deepCoadd_srcMatchFull', 'deepCoadd_icMatch', 'deepCoadd_icMatchFull', 'deepCoadd_psf', 'processStack_config', 'stackExposureId', 'stackExposureId_bits', 'deepCoadd_diffsrc', 'deepCoadd_apCorr', 'solvetansip_config', 'deep_processCoadd_config', 'deep_coadd_metadata', 'deep_coadd_config', 'deepCoadd_multibandReprocessing', 'deepCoadd_icSrc', 'forcedCoadd_config', 'deepCoadd_initPsf', 'deep_processCoadd_metadata', 'deepCoadd_diff', 'deepCoadd_depth', 'deepCoadd_bg', 'deepCoadd_bgRef', 'deepCoadd_calexp')

def getIndent(line):
    for i in range(len(line)):
        if not line[i].isspace():
            break
    return i

def readPafSection(filename, type=None):
    result = {}
    try:
        fin = open(filename, "r")
    except:
        return result
    lines = fin.readlines()
    if type is None:
        indent = -4
    else:
        indent = None    
    for line in lines:
        line = line.rstrip()
        if line.find("#") >= 0:
            continue
        if len(line.strip()) == 0:
            continue
        if indent is None:
            if line.find(type + ":") >= 0:
                indent = getIndent(line)               
            continue
        thisindent = getIndent(line)
        if thisindent <= indent:
                break 
        if line.find(':') > 0:
            if thisindent == indent + 4:
                name = line[:line.find(":")].strip()
                result[name] = None
    fin.close()           
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file1", type=str, help="first file", default=None) 
    parser.add_argument("file2", type=str, help="second file", default=None) 
    parser.add_argument("-c", "--change", type=str, help="move, obs, or both", default=None) 
    args = parser.parse_args()
    types = ("datasets", "images", "exposures", "calibrations")
    excludekeys = []
    if args.change == "both" or args.change == "move":
        excludekeys.extend(readPafSection("movekeys.yaml").keys())
    if args.change == "both" or args.change == "obs":
        excludekeys.extend(readPafSection("obskeys.yaml").keys())
    for type in types:
        rmKeys = readPafSection(type+".yaml.new").keys()
        ref = readPafSection(args.file1, type)
        comp = readPafSection(args.file2, type)
        print type, len(ref.keys()), len(comp.keys())
        for key in ref.keys():
            if not key in comp.keys() and not key in excludekeys:
                print ("    %s.%s key not present in file2.\n"%(type, key))
        # see if key in current was in the reference
        for key in comp.keys():
            if not key in ref.keys() and not key in excludekeys:
                print ("    %s.%s key not present in file1.\n"%(type, key))
