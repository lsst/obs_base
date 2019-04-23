from lsst.meas.algorithms import Defects
import os
import glob
from dateutil import parser


def read_defects_one_chip(root, chip_name):
    files = glob.glob(os.path.join(root, chip_name, '*.dat'))
    defect_dict = {}
    for f in files:
        date_str = os.path.splitext(os.path.basename(f))[0]
        valid_start = parser.parse(date_str)
        defect_dict[valid_start] = Defects.readLsstDefectsFile(f)
    return defect_dict


def read_all_defects(root):
    dirs = os.listdir(root)  # assumes all directories contain defects
    defects_by_chip = {}
    for d in dirs:
        chip_name = os.path.basename(d)
        defects_by_chip[chip_name] = read_defects_one_chip(root, chip_name)
    return defects_by_chip
