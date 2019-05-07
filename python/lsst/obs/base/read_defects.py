from lsst.meas.algorithms import Defects
import os
import glob
from dateutil import parser


def read_defects_one_chip(root, chip_name, camera):
    files = glob.glob(os.path.join(root, chip_name.lower(), '*.dat'))
    parts = os.path.split(root)
    instrument = os.path.split(parts[0])[1]  # convention is that these reside at <instrument>/defects
    defect_dict = {}
    for f in files:
        date_str = os.path.splitext(os.path.basename(f))[0]
        valid_start = parser.parse(date_str)
        defect_dict[valid_start] = Defects.readLsstDefectsFile(f)
        fix_up_metadata(defect_dict[valid_start], valid_start, instrument, camera[chip_name].getId())
    return defect_dict


def fix_up_metadata(defects, valid_start, instrument, chip_id):
    md = defects.getMetadata()
    md['INSTRUME'] = instrument
    md['DETECTOR'] = chip_id
    md['CALIBDATE'] = valid_start.isoformat()
    md['FILTER'] = None
    md['CALIB_ID'] = f'detector={chip_id} calibDate={valid_start.isoformat()} ccd={chip_id} ccdnum={chip_id} filter=None'

def read_all_defects(root, camera):
    dirs = os.listdir(root)  # assumes all directories contain defects
    defects_by_chip = {}
    name_map = {det.getName().lower():det.getName() for det in camera}  # we assume the directories have been lowered
    for d in dirs:
        chip_name = name_map[os.path.basename(d)]
        defects_by_chip[chip_name] = read_defects_one_chip(root, chip_name, camera)
    return defects_by_chip
