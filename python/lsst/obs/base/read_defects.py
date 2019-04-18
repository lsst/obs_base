from lsst.meas.algorithms import Defect
from lsst.afw.geom import Box2I, Point2I, Extent2I
import numpy
import os
import glob
from dateutil import parser
def defect_reader(filename):
    defects = []
    defect_array = numpy.genfromtxt(filename, dtype=[('x0', '<i8'), ('y0', '<i8'), 
                                                     ('x_extent', '<i8'), ('y_extent', '<i8')])
    for row in defect_array:
        pt = Point2I(row['x0'], row['y0'])
        ext = ExtentI(row['x_extent'], row['y_extent'])
        box = Box2I(pt, ext)
        defects.append(box)
    return defects

def read_defects_one_chip(root, chip_name):
    files = glob.glob(os.path.join(root, chip_name, '*.dat'))
    defect_dict = {}
    for f in files:
        date_str = os.path.splitext(os.path.basename(f))[0]
        valid_start = parser.parse(date_str)
        defect_dict[valid_start] = defect_reader(f)
    return defect_dict

def read_all_defects(root):
    dirs = os.listdir(root)  # assumes all directories contain defects
    defects_by_chip = {}
    for d in dirs:
        chip_name = os.path.basename(d)
        defects_by_chip[chip_name] = read_defects_one_chip(root, chip_name)
    return defects_by_chip   

