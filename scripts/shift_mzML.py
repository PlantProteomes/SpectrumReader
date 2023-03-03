#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"

# Uses PPM Shift as a parsed variable
# to run the program (use terminal): python shift_mzML.py --input_filename ..\data\HFX_9850_GVA_DLD1_2_180719_subset.mzML --output_filename ..\data\HFX_9850_GVA_DLD1_2_180719_subset_calibrated.mzML --ppm_shift 5
# to run the program (use terminal): python shift_mzML.py --input_filename ..\data\HFX_9850_GVA_DLD1_2_180719.mzML.gz --output_filename ..\data\HFX_9850_GVA_DLD1_2_180719_calibrated.mzML.gz --ppm_shift 5
# to run the program (use terminal): python shift_mzML.py --input_filename ..\data\QEX03_210305_PTMscan_wt_IP_63.mzML.gz --output_filename ..\data\QEX03_210305_PTMscan_wt_IP_63_calibrated.mzML.gz --ppm_shift 5
# to run the program (use terminal): python shift_mzML.py --input_filename ..\data\06CPTAC_BCprospective_W_BI_20161116_BA_f17.mzML.gz --output_filename ..\data\06CPTAC_BCprospective_W_BI_20161116_BA_f17_calibrated.mzML.gz --ppm_shift 5
# to run the program (use terminal): python shift_mzML.py --input_filename ..\data\Q20181210_06.mzML.gz --output_filename ..\data\Q20181210_06_calibrated.mzML.gz --ppm_shift 5
# to run the program (use terminal): python shift_mzML.py --input_filename ..\data\OR13_A2_20161014_CMP_Ecoli_glycerol_exp_2mg_IMAC_R2.mzML.gz --output_filename ..\data\OR13_A2_20161014_CMP_Ecoli_glycerol_exp_2mg_IMAC_R2_calibrated.mzML.gz --ppm_shift 5

# Uses json file to shift the file
# to run the program (use terminal): python shift_mzML.py --input_filename ..\data\HFX_9850_GVA_DLD1_2_180719_subset.mzML --output_filename ..\data\HFX_9850_GVA_DLD1_2_180719_subset_calibrated.mzML --json_filename ..\data\HFX_9850_GVA_DLD1_2_180719_subset.calibration.json
# to run the program (use terminal): python shift_mzML.py --input_filename ..\data\HFX_9850_GVA_DLD1_2_180719.mzML.gz --output_filename ..\data\HFX_9850_GVA_DLD1_2_180719_calibrated.mzML.gz --json_filename ..\data\HFX_9850_GVA_DLD1_2_180719.calibration.json
# to run the program (use terminal): python shift_mzML.py --input_filename ..\data\QEX03_210305_PTMscan_wt_IP_63.mzML.gz --output_filename ..\data\QEX03_210305_PTMscan_wt_IP_63_calibrated.mzML.gz --json_filename ..\data\QEX03_210305_PTMscan_wt_IP_63.calibration.json
# to run the program (use terminal): python shift_mzML.py --input_filename ..\data\06CPTAC_BCprospective_W_BI_20161116_BA_f17.mzML.gz --output_filename ..\data\06CPTAC_BCprospective_W_BI_20161116_BA_f17_calibrated.mzML.gz --json_filename ..\data\06CPTAC_BCprospective_W_BI_20161116_BA_f17.calibration.json
# to run the program (use terminal): python shift_mzML.py --input_filename ..\data\Q20181210_06.mzML.gz --output_filename ..\data\Q20181210_06_calibrated.mzML.gz --json_filename ..\data\Q20181210_06.calibration.json
# to run the program (use terminal): python shift_mzML.py --input_filename ..\data\OR13_A2_20161014_CMP_Ecoli_glycerol_exp_2mg_IMAC_R2.mzML.gz --output_filename ..\data\OR13_A2_20161014_CMP_Ecoli_glycerol_exp_2mg_IMAC_R2_calibrated.mzML.gz --json_filename ..\data\OR13_A2_20161014_CMP_Ecoli_glycerol_exp_2mg_IMAC_R2.calibration.json

import sys
import os
import argparse
import os.path
import timeit
import re
import json
import numpy
import pickle
import gzip
import numpy as np

from pyteomics import mzml, auxiliary
from scipy import interpolate
from psims.transform.mzml import MzMLTransformer

def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)


class MyMzMLTransformer(MzMLTransformer):
    def __init__(self, input_stream, output_stream, transform=None, transform_description=None,
                 sort_by_scan_time=False, ppm_shift=None, correction_values=None):
        super().__init__(input_stream, output_stream, transform=transform, transform_description=transform_description,
                 sort_by_scan_time=sort_by_scan_time)

        if ppm_shift is not None:
            self.ppm_shift = float(ppm_shift)
            self.has_spline_correction = False
        else:
            self.ppm_shift = float(correction_values['crude correction'])
            self.has_spline_correction = True
            self.t = correction_values['t']
            self.c = correction_values['c']
            self.k = int(correction_values['k'])

    def format_spectrum(self, spectrum):
        new_spectrum = super().format_spectrum(spectrum)

        #### Shift all the m/z values
        new_spectrum['mz_array'] += self.ppm_shift * new_spectrum['mz_array'] / 1e6
        if self.has_spline_correction:
            x_values = np.array([new_spectrum['mz_array']])
            y_values = interpolate.BSpline(self.t, self.c, self.k)(x_values)
            new_spectrum['mz_array'] += y_values[0] * new_spectrum['mz_array'] / 1e6

        #### Get the MS level
        ms_level = None
        for param in new_spectrum['params']:
            if param['name'] == 'MS:1000511':
                ms_level = param['value']

        #### Correct the precursor m/z values by the requested shift
        if ms_level is not None and ms_level > 1:
            precursor_mz = new_spectrum['precursor_information'][0]['mz']
            precursor_mz += self.ppm_shift * precursor_mz / 1e6
            if self.has_spline_correction:
                if precursor_mz <= 400:
                    x_values = np.array([precursor_mz])
                else:
                    x_values = np.array([400])
                y_values = interpolate.BSpline(self.t, self.c, self.k)(x_values)
                precursor_mz += y_values[0] * precursor_mz / 1e6
            new_spectrum['precursor_information'][0]['mz'] = precursor_mz

        return new_spectrum

def main():

    argparser = argparse.ArgumentParser(description='Read an input mzML file and write out a new one with the m/zs shifted')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--input_filename', type=str, action='store', required=True, help='Name of the input mzML file')
    argparser.add_argument('--output_filename', type=str, action='store', required=True, help='Name of the output mzML file shifted m/zs')
    argparser.add_argument('--ppm_shift', type=str, action='store', help='Offset to shift all spectra in units of PPM')
    argparser.add_argument('--json_filename', type=str, action='store', help='Name of the json file with calibration values')
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    #### Check that the filenames are defined and exist
    if not os.path.isfile(params.input_filename):
        print(f"ERROR: --input_filename '{params.input_filename}' not found or not a file")
        return

    #### Ensure the output is not the same as the input
    if params.output_filename == params.input_filename:
        print(f"ERROR: --output_filename '{params.input_filename}' may not be the same as the --input_filename")
        return

    #### Open in the input file
    print(f"INFO: Opening {params.input_filename} for reading")
    if params.input_filename.endswith('.gz'):
        infile = gzip.open(params.input_filename)
    else:
        infile = open(params.input_filename, 'rb')

    #### Open the output file
    print(f"INFO: Opening {params.output_filename} for writing")
    try:
        outfile = open(params.output_filename, 'wb')
    except:
        print(f"ERROR: --output_filename '{params.output_filename}' is not writable. Permissions problem?")
        return
    
    # Determines if json file is used or ppm shift value
    if params.json_filename is None or params.json_filename == "":
        if params.ppm_shift is None:
            print("ERROR: --json_filename and --ppm_shift were not found")
            return
        MyMzMLTransformer(infile, outfile, ppm_shift=params.ppm_shift,
            transform_description=f"Shifted all spectra by {params.ppm_shift}").write()
    else:
        if not os.path.isfile(params.json_filename):
            print(f"ERROR: --json_filename '{params.json_filename}' not found or not a file")
            return
            
        correction_file = open(f'{params.json_filename}')
        correction_values = json.load(correction_file)
        MyMzMLTransformer(infile, outfile, correction_values=correction_values, 
            transform_description=f"Shifted all spectra by values {correction_values['crude correction']} and spline fit").write()

#### For command line usage
if __name__ == "__main__": main()
