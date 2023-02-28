#!/usr/bin/env python3

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

from pyteomics import mzml, auxiliary
from psims.transform.mzml import MzMLTransformer

def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)


####################################################################################################
class MyMzMLTransformer(MzMLTransformer):
    def __init__(self, input_stream, output_stream, transform=None, transform_description=None,
                 sort_by_scan_time=False, ppm_shift=None):
        super().__init__(input_stream, output_stream, transform=transform, transform_description=transform_description,
                 sort_by_scan_time=sort_by_scan_time)

        self.ppm_shift = float(ppm_shift)

    def format_spectrum(self, spectrum):
        new_spectrum = super().format_spectrum(spectrum)

        #### Shift all the m/z values
        new_spectrum['mz_array'] += self.ppm_shift * new_spectrum['mz_array'] / 1e6

        #### Get the MS level
        ms_level = None
        for param in new_spectrum['params']:
            if param['name'] == 'MS:1000511':
                ms_level = param['value']

        #### Correct the precursor m/z values by the requested shift
        if ms_level is not None and ms_level > 1:
            precursor_mz = new_spectrum['precursor_information'][0]['mz']
            precursor_mz += self.ppm_shift * precursor_mz / 1e6
            new_spectrum['precursor_information'][0]['mz'] = precursor_mz

        return new_spectrum


####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Read an input mzML file and write out a new one with the m/zs shifted')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--input_filename', type=str, action='store', required=True, help='Name of the input mzML file')
    argparser.add_argument('--output_filename', type=str, action='store', required=True, help='Name of the output mzML file shifted m/zs')
    argparser.add_argument('--ppm_shift', type=str, action='store', required=True, help='Offset to shift all spectra in units of PPM')
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

    MyMzMLTransformer(infile, outfile, ppm_shift=params.ppm_shift,
        transform_description=f"Shifted all spectra by {params.ppm_shift}").write()


#### For command line usage
if __name__ == "__main__": main()
