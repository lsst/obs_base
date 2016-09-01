import argparse
import sys
from lsst.obs.cfht import MegacamMapper
from lsst.obs.decam import DecamMapper
from lsst.obs.lsstSim import LsstSimMapper
from lsst.obs.test import TestMapper
from lsst.obs.sdss import SdssMapper
from lsst.obs.hsc import HscMapper
from lsst.obs.suprimecam import SuprimecamMapper

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--name", type=str, help="mapper name", default=None) 
    parser.add_argument("-t", "--type", type=str, help="policy section", default=None) 
    args = parser.parse_args()
    name = args.name
    if name == "megacam": mapper = MegacamMapper()
    elif name == "decam": mapper = DecamMapper()
    elif name == "lsstSim": mapper = LsstSimMapper()
    elif name == "test": mapper = TestMapper()
    elif name == "sdss": mapper = SdssMapper()
    elif name == "hsc": 
        mapper = HscMapper(root="/sandbox/lsstshared/pgee/mylsst12/Linux64/ci_hsc/DATA/")
    elif name == "suprimecam": mapper = SuprimecamMapper(root="/sandbox/lsstshared/pgee/mylsst12/Linux64/ci_hsc/DATA/")
    else:
        print "Name of mapper (-n name) must be specified as an existing mapper."
        sys.exit() 
    camName = mapper.getCameraName()
    if args.type == None:
        types = ("exposures", "datasets", "calibrations", "images")
    else:
        types = (args.type,)
    print types
    for type in types:
        if not type in mapper.__dict__.keys():
            continue
        fout = open(camName + "." + type, "w")
        sets = mapper.__dict__[type]
        for key in sets.keys():
            template = str(sets[key].template)
            
            #if not key in keys:
            #    continue
            fout.write(key + ":\n")
            fout.write("    persistable: " + str(sets[key].persistable) + "\n")
            fout.write("    storage: " + str(sets[key].storage) + "\n")
            fout.write("    python: " + str(sets[key].python) + "\n")
            tables = sets[key].tables
            if not tables is None:
                fout.write("    tables: " + " ".join(tables) + "\n")
            fout.write("    template: " + str(sets[key].template))
            fout.write("\n")
        fout.close()
