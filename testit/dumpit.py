import argparse
import sys
from lsst.obs.suprimecam import SuprimecamMapper
from utils import getMapperSet, getMapper

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--name", type=str, help="mapper name", default=None) 
    parser.add_argument("-t", "--type", type=str, help="policy section", default=None) 
    parser.add_argument("-d", "--dir", type=str, help="output directory", default=".") 
    args = parser.parse_args()
    name = args.name
    for type in ("datasets", "exposures", "calibrations"):
        sets = getMapperSet(name, type)
        if sets is None:
            print ("Name of mapper (-n name) must be specified:")
            print ("    hsc, suprimecam, decam, test, sdss, megacam, lsstSim")
            sys.exit(1)
        fout = open(args.dir + "/" + name + "." + type, "w")
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
