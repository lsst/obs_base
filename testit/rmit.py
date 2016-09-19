import argparse
import sys
from utils import readPafSection, compare

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=str, help="first file", default=None) 
    parser.add_argument("type", type=str, help="remove obsolete, or moved", default="obs") 
    args = parser.parse_args()
    keys = {}
    movedict = readPafSection("datasets.yaml.new")
    movedict.update(readPafSection("exposures.yaml.new"))
    obsdict = readPafSection("obskeys.yaml")
    pafdict = readPafSection(args.file + ".paf", "datasets")
    pafdict.update(readPafSection(args.file + ".paf", "exposures"))
    fin = open(args.file + ".paf", "r")
    if args.type.startswith("mov"):
        fout = open(args.file + ".rem_mov.paf", "w")
        fout2 = open(args.file + ".notmoved", "w")
    elif args.type.startswith("obs"):
        fout = open(args.file + ".rem_obs.paf", "w")
    else:
        print "type of removal must be specified"
        sys.exit()
    lines = fin.readlines()
    bset = None
    name = None
    for line in lines:
        if bset:
            if line.find("tables:") < 0 or not args.type.startswith("mov") or not name in movedict.keys(): 
                bset = bset + line
            # This is the ending condition
            if line.find(' }') > 0:
                if args.type.startswith("obs"):
                    if not name in obsdict.keys():
                        fout.write(bset)
                if args.type.startswith("mov"):
                    if name == "deepCoadd" and name == "deepCoadd_calexp":
                        pdb.set_trace()

                    if not name in movedict.keys() or compare(name, pafdict, movedict, ignoretables=True):
                        fout.write(bset)
                        if name in movedict.keys():
                            fout2.write("%s: not moved in %s because of mismatch:\n"%(name, args.file))
                            compare(name, pafdict, movedict, ignoretables=True, fout=fout2)
                keys[name] = None
                bset = None
            continue
        parts = line.split() 
        if len(line) > 4 and line.startswith("    ") and line[4].isalpha() and len(parts) == 2 and parts[0].endswith(":") and parts[1] == '{':
            bset = line
            template = ""
            name = line[4: line.find(":")]
        else:
            fout.write(line)
