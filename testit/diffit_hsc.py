import argparse

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
    if bset:
        result[name] = bset
    fin.close()           
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--name", type=str, help="mapper name", default=None) 
    parser.add_argument("-r", "--refmapper", type=str, help="path to ref", default=None) 
    parser.add_argument("-c", "--cpath", type=str, help="path to comps", default=None) 
    parser.add_argument("-p", "--rpath", type=str, help="path to comps", default="refs") 
    obskeys = readPafSection("obskeys.yaml")
    movekeys = readPafSection("movekeys.yaml")
    args = parser.parse_args()
    name = args.name
    # refmapper is the mappers we are comparing against.  
    refmapper = args.refmapper
    if refmapper == None:
        refmapper = name
        outname = name + ".selfdiffs"
    else:
        outname = name + ".%sdiffs"%(refmapper,)
    if not args.cpath is None:
        outname = args.cpath + "/" + outname
    cpath = args.cpath
    rpath = args.rpath
    fout = open(outname, "w")
    types = ("datasets", "images", "exposures", "calibrations")
    #  "ref" is the original typeset, and "comp" is the new one
    for type in types:
        ref = readSection("%s/%s.%s"%(args.rpath, refmapper, type))
        comp = readSection("%s/%s.%s"%(args.cpath, name, type))

        if ref is None and comp is None:
            continue
        if ref is None and not comp is None:
            fout.write("Section %s doesn't appear in reference of %s"%(type, refmapper))
            continue
        if not ref is None and comp is None:
            fout.write("Section %s doesn't appear in current of %s"%(type,name))
            continue
      
        # Compare the keys in the reference with the keys in the compare
        # The focus is on how the compare differs from the reference, which
        # is the way the output is written.  For example, if the comp is missing
        # something in the reference, it is a "missing" key.  If the comp has
        # somthing which the reference does not, it is an "added" key
        count = 0
        for key in ref.keys():
            if not key in comp.keys():
                if count == 0:
                    count = count + 1
                    fout.write("Keys missing from mapper %s:%s which appear in %s:%s\n"%(cpath, name, rpath, refmapper))
                fout.write("    %s.%s: "%(type, key))
                fout.write("(M)" if key in movekeys else ("(O)" if key in obskeys else "(-)"))
                fout.write("\n")

        count = 0
        for key in comp.keys():
            # just look at the ones we have not already looked at
            if not key in ref.keys():
                if count == 0:
                    count = count + 1
                    fout.write("Keys not in %s:%s which exist in reference %s:%s\n"%(cpath, name, rpath, refmapper))
                fout.write("    %s.%s: "%(type, key))
                fout.write("(M)" if key in movekeys else ("(O)" if key in obskeys else "(-)"))
                fout.write("\n")

        count = 0
        for key in ref.keys():
            if key in comp.keys() and not ref[key] == comp[key]:
                if count == 0:
                    count = count + 1
                    fout.write("Keys which appear in both but differ:\n")
                fout.write("    %s.%s: "%(type, key))
                fout.write("(M)" if key in movekeys else ("(O)" if key in obskeys else "(-)"))
                fout.write("\n")
                refs = ref[key].split("\n")
                comps = comp[key].split("\n")
                
                for i in range(max(len(refs), len(comps))):
                    if i > (len(comps)-1):
                        comps.append("")
                    if i > (len(refs)-1):
                        refs.append("")
                    if not refs[i] == comps[i]:
                        fout.write("        " + comps[i] + "\n")
                        fout.write("        " + refs[i] + " in %s reference\n"%(refmapper))
        # see if all the keys in the comp are also in the reference
