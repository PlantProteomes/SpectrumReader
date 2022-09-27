#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python read_mzML1.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML
# can do python read_mzML1.py --help and look at the usage statements for what I can do

import os
import argparse
import os.path
import timeit
import matplotlib.pyplot as plt
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup

from pyteomics import mzml, auxiliary

def main():

    # sets up a module with a description
    argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
    argparser.add_argument('--mzml_file', action='store', help='Name of the mzML file to read')
    argparser.add_argument('--precursor_mz', action='store', type=float, help='Precursor m/z to select')
    argparser.add_argument('--tolerance_mz', action='store', type=float, help='Tolerance of a m/z and the measured peak m/z')
    argparser.add_argument('--required_peaks', action='store', help='Tolerance of a specified peak m/z and a measured peak m/z')
    argparser.add_argument('--tolerance_peaks', action='store', type=float, help='Tolerance of a m/z and the measured peak m/z')
    # runs the parse_args method to store all the results into the params variable
    params = argparser.parse_args()

    #### Ensure that mzml_file was passed, look in the params variable for the specified mzml_file. if you want to open precursor, it would be params.precursor_mz
    if params.mzml_file is None or params.mzml_file == "":
        print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
        return

    # tries to open the params file specified
    if not os.path.isfile(params.mzml_file):
        print(f"ERROR: File '{params.mzml_file}' not found or not a file")
        return

    hasMZCheck = False
    if not (params.precursor_mz is None or params.precursor_mz == 0.0):
        hasMZCheck = True

    # ensures there is always a tolerance
    if params.tolerance_mz is None or params.tolerance_mz == 0.0:
        params.tolerance_mz = 0.01
    print(f"tolerance_mz is {params.tolerance_mz}")

    if params.tolerance_peaks is None or params.tolerance_peaks == 0.0:
        params.tolerance_peaks = 0.01
    print(f"tolerance_peaks is {params.tolerance_peaks}")

    if (hasMZCheck):
        lowerBound = params.precursor_mz - params.tolerance_mz
        upperBound = params.precursor_mz + params.tolerance_mz

    needToCheckPeaks = False
    # if not (params.required_peaks is None or params.required_peaks == ""):
    if params.required_peaks is not None and params.required_peaks != "":
        requiredPeaksList = [float(i) for i in params.required_peaks.split(',')]
        needToCheckPeaks = True
        print(f"required peaks are {requiredPeaksList}")

    validSpectra = []

    #### Read spectra from the file
    t0 = timeit.default_timer()
    stats = { 'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0 }
    with mzml.read(params.mzml_file) as reader:
        for spectrum in reader:
            if stats['counter'] == 0:
                auxiliary.print_tree(spectrum)

            #### Update counters and print progress
            stats['counter'] += 1

            satisfyRequirements = True
            if spectrum['ms level'] == 1:
                stats['ms1spectra'] += 1
            elif spectrum['ms level'] == 2:
                stats['ms2spectra'] += 1
                precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                if hasMZCheck and (precursor_mz < lowerBound or precursor_mz > upperBound):
                    satisfyRequirements = False
                if needToCheckPeaks and satisfyRequirements:
                    for peak in requiredPeaksList:
                        if spectrum['m/z array'][len(spectrum['m/z array']) - 1] < lowerPeakBound:
                            satisfyRequirements = False
                        else:
                            lowerPeakBound = peak - params.tolerance_peak
                            upperPeakBound = peak + params.tolerance_peak
                            for mz in spectrum['m/z array']:
                                if mz < lowerPeakBound:
                                    continue
                                if mz > upperPeakBound:
                                    satisfyRequirements = False
                                    break
                                break
                        
                if satisfyRequirements:
                    validSpectra.append(spectrum)

            if stats['counter']/1000 == int(stats['counter']/1000):
                print(f"  {stats['counter']}")


    #### Print final timing information
    t1 = timeit.default_timer()
    print(f"INFO: Read {stats['counter']} spectra from {params.mzml_file}")
    print(f"The number of ms1spectra is {stats['ms1spectra']}")
    print(f"The number of ms2spectra is {stats['ms2spectra']}")
    print(f"Number of valid spectra: {len(validSpectra)}")
    # write out mz and intensity array
    # print(f"Valid Spectra: {validSpectra}")
    print(f"INFO: Elapsed time: {t1-t0}")
    print(f"INFO: Processed {stats['counter']/(t1-t0)} spectra per second")

if __name__ == "__main__": main()
