import argparse
from utils import readPafSection, compare, parseEntry

def getItem(reg, key, subkey):
    if key in reg.keys():
        items = parseEntry(reg[key])
        if subkey in items.keys():
            return items[subkey]
    return None
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--name", type=str, help="mapper name", default=None) 
    parser.add_argument("-r", "--refmapper", type=str, help="path to ref", default=None) 
    parser.add_argument("-c", "--cpath", type=str, help="path to comps", default='.') 
    parser.add_argument("-p", "--rpath", type=str, help="path to refs", default=".") 
    parser.add_argument("-x", "--xtags", type=bool, help="print xtagged items", default=True) 
    parser.add_argument("-a", "--added", type=bool, help="print items added in comp", default=True) 
    parser.add_argument("-m", "--missing", type=bool, help="print items missing from comp", default=False) 
    obsdict = readPafSection("obskeys.yaml")
    movedict = readPafSection("datasets.yaml.new")
    movedict.update(readPafSection("exposures.yaml.new"))
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
    types = ("datasets", "exposures", "calibrations")
    #  "ref" is the original typeset, and "comp" is the new one
    for type in types:
        ref = readPafSection("%s/%s.%s"%(args.rpath, refmapper, type))
        comp = readPafSection("%s/%s.%s"%(args.cpath, name, type))

        if ref is None and comp is None:
            continue
        if ref is None and not comp is None:
            fout.write("Section %s doesn't appear in reference of %s"%(type, refmapper))
            continue
        if not ref is None and comp is None:
            fout.write("Section %s doesn't appear in current of %s"%(type,name))
            continue

    #  Do the keys which appear in both, but differ
    count = 0
    xtagged =  {}
    for type in types:
        ref = readPafSection("%s/%s.%s"%(args.rpath, refmapper, type))
        comp = readPafSection("%s/%s.%s"%(args.cpath, name, type))
        for key in ref.keys():
            ignoretables = key in movedict.keys()
            if key in comp.keys():
              diffs = compare(key, ref, comp)
              if diffs:
                tag = ''
                if key in movedict.keys():
                    tag = tag + 'M'
                    if count == 0:
                        fout.write("\n1. Datasets should have been moved (M), but retain a definition in %s.%s:\n"%(args.cpath, name))
                    count = count + 1
                    fout.write("    %s.%s\n"%(type, key))
                    compare(key, ref, comp, refName=refmapper, indent = "      ", fout=fout)
                else:
                    if 'template' in diffs:
                        tag = tag + 'x'
                    if 'coaddtemplate' in diffs:
                        tag = tag + 'X'
                    if 'python' in diffs:
                        tag = tag + 'P'
                    if 'persistable' in diffs:
                        tag = tag + 'p'
                    if 'tables' in diffs:
                        tag = tag + 'T'
                    if 'storage' in diffs:
                        tag = tag + 'S'
                    if 'level' in diffs:
                        tag = tag + 'L'
                    outstr = "    %s %s.%s\n"%(tag.ljust(4), type, key)
                    xtagged[key] = outstr
    #  Do the keys in the reference which do not appear in the comp    
    if args.added:
      count = 0
      for type in types:
        ref = readPafSection("%s/%s.%s"%(args.rpath, refmapper, type))
        comp = readPafSection("%s/%s.%s"%(args.cpath, name, type))
        for key in comp.keys():
            # just look at the ones we have not already looked at
            if not key in ref.keys():
                if count == 0:
                    count = count + 1
                    fout.write("\n2. Datasets in %s:%s which are not in reference %s:%s\n"%(cpath, name, rpath, refmapper))
                    fout.write("-- If these are keys unique to this mapper, keeping them is fine.\n")
                    fout.write("-- But if they correspond to one of the HSC keys, the name should be changed.\n")
                    fout.write("-- And if they are no longer in use, they should be deleted.\n")
                fout.write("    %s.%s: "%(type, key))
                fout.write("(M)" if key in movedict.keys() else ("(O)" if key in obsdict.keys() else ""))
                fout.write("\n")

    #  Do the keys missing from comp which appear in the reference 
    count = 0
    for type in types:
        ref = readPafSection("%s/%s.%s"%(args.rpath, refmapper, type))
        comp = readPafSection("%s/%s.%s"%(args.cpath, name, type))
        # Compare the keys in the reference with the keys in the compare
        # The focus is on how the compare differs from the reference, which
        # is the way the output is written.  For example, if the comp is missing
        # something in the reference, it is a "missing" key.  If the comp has
        # somthing which the reference does not, it is an "added" key
        for key in ref.keys():
            if not key in comp.keys():
                if count == 0:
                    if args.missing:    
                        fout.write("\n3. Datasets missing from mapper %s:%s which appear in %s:%s\n"%(cpath, name, rpath, refmapper))
                        fout.write("-- If these are keys are not appropriate for this mapper, that's fine.\n")
                        fout.write("-- But if they belong in this mapper and  have a different name, consider modifying.\n")
                count = count + 1
                if args.missing:    
                    fout.write("    %s.%s: "%(type, key))
                    fout.write("(M)" if key in movedict.keys() else ("(O)" if key in obsdict.keys() else ""))
                    fout.write("\n")
    if count > 0 and not args.missing:
        fout.write("\n3. %d datasets missing from mapper %s:%s which appear in %s:%s\n"%(count, cpath, name, rpath, refmapper))


    if len(xtagged.keys()) > 0:
        if not args.xtags:
            fout.write("\n4. There were %d datasets common names, but differed from %s to %s."%(len(xtagged), name, refmapper))
            fout.write("\nThese are not defined in daf_butlerUtils, however.\n")
        else:
            fout.write("\n4. Datasets with common names in %s and %s, but which differ.\n"%(name, refmapper))
            fout.write("These are not defined in daf_butlerUtils, however.\n")
            fout.write("-- Difference is marked by X=template, T=tables, S=storage, L=level, p=persistable, and P=Python.\n")
            fout.write("-- x=template differs but uses columns besides filter, patch, or tract:\n")
            for key in xtagged.keys():
                fout.write(xtagged[key])

