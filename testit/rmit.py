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
    movedict.update(readPafSection("datasets.yaml"))
    movedict.update(readPafSection("exposures.yaml"))
    obsdict = readPafSection("obskeys.yaml")
    pafdict = readPafSection(args.file + ".paf", "datasets")
    pafdict.update(readPafSection(args.file + ".paf", "exposures"))
    fin = open(args.file + ".paf", "r")
    if args.type.startswith("mov"):
        fout = open(args.file + ".rem_mov.paf", "w")
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
            bset = bset + line
            # This is the ending condition
            if line.find(' }') > 0:
                if args.type.startswith("obs"):
                    if not name in obsdict.keys():
                        fout.write(bset)
                if args.type.startswith("mov"):
                    if not name in movedict.keys():
                        fout.write(bset)
                    else:
                        differ = compare(name, pafdict, movedict, ignoretables=False)
                        if "coaddtemplate" in differ:
                            differ.append("template")
                        if len(differ) > 0:
                            fout.write("    # dataset defined in obs_base modifield in this mapper as follows:\n")
                            for item in bset.split("\n"):
                                if item.find(":") >= 0:
                                    key = item[:item.find(":")].strip()
                                    if not key == name and not key in differ:
                                        continue
                                if len(item) > 0:
                                    fout.write(item + "\n")
                        else:
                            fout.write(bset)
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
