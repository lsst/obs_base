import argparse
import sys

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
    parser.add_argument("file", type=str, help="first file", default=None) 
    parser.add_argument("type", type=str, help="remove obsolete or moved", default="obs") 
    args = parser.parse_args()
    keys = {}
    movekeys = readPafSection("movekeys.yaml")
    obskeys = readPafSection("obskeys.yaml")
    fin = open(args.file + ".paf", "r")
    if args.type.startswith("mov"):
        rmKeys = movekeys
        fout = open(args.file + ".rem_mov.paf", "w")
    elif args.type.startswith("obs"):
        rmKeys = obskeys
        fout = open(args.file + ".rem_obs.paf", "w")
    else:
        print "type of removal must be specified"
        sys.exit()
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
