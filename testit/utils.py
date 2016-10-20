import argparse
import sys

mappers = {}

def getIndent(line):
    for i in range(len(line)):
        if not line[i].isspace():
            break
    return i

def getMapper(name):
    if name in mappers.keys():
        return mappers[name]
    mapper = None
    if name == "megacam":
        from lsst.obs.cfht import MegacamMapper
        mapper = MegacamMapper()
    elif name == "decam":
        from lsst.obs.decam import DecamMapper
        mapper = DecamMapper()
    elif name == "lsstSim":
        from lsst.obs.lsstSim import LsstSimMapper
        mapper = LsstSimMapper()
    elif name == "test":
        from lsst.obs.test import TestMapper
        mapper = TestMapper()
    elif name == "sdss":
        from lsst.obs.sdss import SdssMapper
        mapper = SdssMapper()
    elif name == "hsc": 
        from lsst.obs.hsc import HscMapper
        mapper = HscMapper(root="/sandbox/lsstshared/pgee/mylsst12/Linux64/ci_hsc/DATA/")
    elif name == "suprimecam":
        from lsst.obs.suprimecam import SuprimecamMapper
        mapper = SuprimecamMapper(root="/sandbox/lsstshared/pgee/mylsst12/Linux64/ci_hsc/DATA/")
    if not mapper is None:
        mappers[name] = mapper
    return mapper

def getMapperSet(name, type):
    mapper = getMapper(name)
    return mapper.__dict__[type]

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
        if line.lstrip().find("#") >= 0:
            line = line[:line.find("#")].strip()
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
                result[name] = line 
            if thisindent == indent + 8:
                result[name] = result[name] + line
    fin.close()           
    return result

def parseEntry(text):
    result = {}
    lines = text.split("\n")
    for line in lines:
        line = line.strip()
        if len(line) == 0:
            continue
        if line == '}':
            continue
        name = line[0: line.find(":")]
        entry = line[line.find(":")+1:]
        entry = entry.strip()
        if entry.find("{") >= 0:
            entry = entry[entry.find("{") + 1:]
        if entry == '':
            continue
        if entry.startswith('"') and entry.endswith('"'):
            entry = entry[1: len(entry) -1]
        if entry.startswith("'") and entry.endswith("'"):
            entry = entry[1: len(entry) -1]
        if name in result.keys():
            result[name] = result[name] + " " + entry
        else:
            result[name] = entry
    return result

def onlyCoadd(template):
    coaddKeys = ["filter", "patch", "tract"]
    parts = template.split("%(")
    parts.pop(0)
    for part in parts:
        part = part[0: part.find(')')]
        if not part in coaddKeys:
           return False
    return True

def compare(name, refdict, compdict, refName="ref", compName="comp", ignoretables=False, fout=None, indent="    "):
    result = list() 
    if not name in refdict.keys() or not name in compdict.keys():
        return result
    items1 = parseEntry(refdict[name])
    items2 = parseEntry(compdict[name])
    for key in items1.keys():
        if key == 'tables' and ignoretables:
            continue
        if not key in items2.keys():
            result.append(key)
            if fout: fout.write(indent + "%s key not in comp\n"%key)
    for key in items2.keys():
        if key == 'tables' and ignoretables:
            continue
        if not key in items1.keys():
            result.append(key)
            if fout: fout.write(indent + "%s key not in ref\n"%key)
    for key in items1.keys():
        if key == 'tables' and ignoretables:
            continue
        if key in items2.keys() and not items1[key] == items2[key] and not key==name:
            if key == 'template' and (onlyCoadd(items2[key]) or onlyCoadd(items1[key])):
                result.append("coaddtemplate")
            else:
                result.append(key)
            if fout:
                line = items2[key]
                line2 = None
                if len(line) > 50:
                    line2 = line[50:]
                    line = line[0:50]
                fout.write(indent + "    " + key + "\t" + line + "\n")
                if not line2 is None:
                    fout.write(indent + "            " + "\t" + line2 + "\n")
            if fout:
                line = items1[key]
                line2 = None
                if len(line) > 50:
                    line2 = line[50:]
                    line = line[0:50]
                fout.write(indent + "ref:" + key + "\t" + line + "\n")
                if not line2 is None:
                    fout.write(indent + "            " + "\t" + line2 + "\n")
    result.sort()
    return result
