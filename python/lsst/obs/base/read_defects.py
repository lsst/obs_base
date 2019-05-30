from lsst.meas.algorithms import Defects
import os
import glob
from dateutil import parser


def read_defects_one_chip(root, chip_name, camera):
    files = glob.glob(os.path.join(root, chip_name.lower(), '*.ecsv'))
    parts = os.path.split(root)
    instrument = os.path.split(parts[0])[1]  # convention is that these reside at <instrument>/defects
    defect_dict = {}
    for f in files:
        date_str = os.path.splitext(os.path.basename(f))[0]
        valid_start = parser.parse(date_str)
        defect_dict[valid_start] = Defects.readText(f)
        check_metadata(defect_dict[valid_start], valid_start, instrument, camera[chip_name].getId())
    return defect_dict


def check_metadata(defects, valid_start, instrument, chip_id):
    md = defects.getMetadata()
    finst = md.get('INSTRUME')
    fchip_id = md.get('DETECTOR')
    fcalib_date = md.get('CALIBDATE')
    if not (finst, int(fchip_id), fcalib_date) == (instrument, chip_id, valid_start.isoformat()):
        raise ValueError("Path and file metadata do not agree:\n" +
                         "Path metadata: %s, %s, %s\n"%(instrument, chip_id, valid_start.isoformat()) +
                         "File metadata: %s, %s, %s\n"%(finst, fchip_id, fcalib_date))


def read_all_defects(root, camera):
    dirs = os.listdir(root)  # assumes all directories contain defects
    dirs = [d for d in dirs if os.path.isdir(os.path.join(root, d))]
    defects_by_chip = {}
    name_map = {det.getName().lower(): det.getName() for
                det in camera}  # we assume the directories have been lowered
    for d in dirs:
        chip_name = name_map[os.path.basename(d)]
        defects_by_chip[chip_name] = read_defects_one_chip(root, chip_name, camera)
    return defects_by_chip
