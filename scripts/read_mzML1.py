#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python read_mzML1.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML

import os
import argparse
import os.path
import timeit

from pyteomics import mzml, auxiliary

def main():

    argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
    argparser.add_argument('--mzml_file', action='store', help='Name of the mzML file to read')
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    params = argparser.parse_args()

    #### Ensure that mzml_file was passed
    if params.mzml_file is None or params.mzml_file == "":
        print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
        return

    if not os.path.isfile(params.mzml_file):
        print(f"ERROR: File '{params.mzml_file}' not found or not a file")
        return


    #### Read spectra from the file
    t0 = timeit.default_timer()
    stats = { 'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0 }
    with mzml.read(params.mzml_file) as reader:
        for spectrum in reader:
            if stats['counter'] == 0:
                auxiliary.print_tree(spectrum)

            #### Update counters and print progress
            stats['counter'] += 1

            # check to see how i can access the ms level
            # if stats['counter'] == 1:
                # print(spectrum)

            if spectrum['ms level'] == 1:
                stats['ms1spectra'] += 1
            elif spectrum['ms level'] == 2:
                stats['ms2spectra'] += 1

            if stats['counter']/1000 == int(stats['counter']/1000):
                print(f"  {stats['counter']}")


    #### Print final timing information
    t1 = timeit.default_timer()
    print(f"INFO: Read {stats['counter']} spectra from {params.mzml_file}")
    print(f"The number of ms1spectra is {stats['ms1spectra']}")
    print(f"The number of ms2spectra is {stats['ms2spectra']}")
    print(f"INFO: Elapsed time: {t1-t0}")
    print(f"INFO: Processed {stats['counter']/(t1-t0)} spectra per second")


if __name__ == "__main__": main()
