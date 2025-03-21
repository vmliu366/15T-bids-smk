import sys
import pprint
import json
import pydicom
from pybruker import jcamp
import numpy as np
from glob import glob


# adapted from conversion scripts written by Naila Rahman (Corey Baron Lab)

#get series number from json
with open(snakemake.params.jsons[0]) as f:
   metadata = json.load(f)

series_num = metadata['SeriesNumber']


dcm = glob(snakemake.input.dcm_folder + f'/*/*/*/*/*/{series_num}/*.dcm')[-1]


# read dicom header
H = pydicom.read_file(dcm, stop_before_pixels=True)


# read bruker headers
method = jcamp.jcamp_parse(
    H[0x0177, 0x1100].value.decode('utf-8').splitlines()
)
print(metadata)

#visu_pars = jcamp.jcamp_parse(
#    H[0x0177, 0x1101].value.decode('utf-8').splitlines()
#)

with open(snakemake.output.method_json, 'w') as fp:
    json.dump(method,  fp, indent=4)