#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# python generate_input_files.py --unknown_mz 105.0659  

import os
import argparse
import os.path
import timeit
import matplotlib.pyplot as plt
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
import csv

from tabulate import tabulate
from pyteomics import mzml, auxiliary

class inputGenerator:

    def __init__(self):
        # sets up a module with a description
        argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
        argparser.add_argument('--unknown_mz', action='store', type=float, help='Precursor m/z to select')

        # runs the parse_args method to store all the results into the params variable
        params = argparser.parse_args()

        if params.unknown_mz is None or params.unknown_mz == 0.0:
            print('ERROR: Parameter --unknown_mz must be provided. See --help for more information')
            return
        else:
            self.unknown_mz = params.unknown_mz

        self.known_ions = [
            ['IQ', 101.0709],
            ['IK', 101.1073],
            ['IM', 104.0529],
            ['IH', 110.0712],
            ['IF', 120.0809],
            ['b-Q', 129.0658],
            ['b-K', 129.1022],
            ['IC[Carbamidomethyl]', 133.0430],
            ['IY', 136.0757],
            ['yK', 147.1128],
            ['IW', 159.0917],
            ['y-R', 175.1190]
        ]

    def generate_file(self):
        with open(f'correlations_with_{self.unknown_mz}.tsv', 'w') as file:
            writer = csv.writer(file, delimiter='\t', lineterminator='\n')

            for ion in self.known_ions:
                line = [self.unknown_mz, ion[0], ion[1]]
                writer.writerow(line)

def main():

    input_file = inputGenerator()
    input_file.generate_file()


if __name__ == "__main__": main()
