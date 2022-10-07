#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python read_mzML1.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML
# can do python read_mzML1.py --help and look at the usage statements for what I can do

# test case: 506.888, 5 valids
# python find_spectra.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --precursor_mz 506.888 --tolerance_mz 0.1
# test case: 506.888, required peaks 244.166 and 245.168, 2 valids
# python find_spectra.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --precursor_mz 506.888 --tolerance_mz 0.1 --required_peaks "244.166, 245.168" --tolerance_peaks 0.1
# test case: 395.1828, required peaks 217.08191, 217.54974, 217.5828, 7 valids
# python find_spectra.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --precursor_mz 395.1828 --tolerance_mz 0.1 --required_peaks "217.08191, 217.54974, 217.5828" --tolerance_peaks 0.1
# test case: 387.1953, required peaks 129.10266, 402.17227, 7 valids
# python find_spectra.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --precursor_mz 387.1953 --tolerance_mz 0.1 --required_peaks "129.10266, 402.17227" --tolerance_peaks 0.1
# test case: 387.1953, required peaks 129.10266, optional peak 402.17227, 29 valids
# python find_spectra.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --precursor_mz 387.1953 --tolerance_mz 0.1 --required_peaks "129.10266" --tolerance_peaks 0.1 --optional_peaks "402.17227"

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

def main():

    # sets up a module with a description
    argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
    argparser.add_argument('--mzml_file', action='store', help='Name of the mzML file to read')
    argparser.add_argument('--precursor_mz', action='store', type=float, help='Precursor m/z to select')
    argparser.add_argument('--tolerance_mz', action='store', type=float, help='Tolerance of a m/z and the measured peak m/z')
    argparser.add_argument('--required_peaks', action='store', help='Required peaks that must show up')
    argparser.add_argument('--tolerance_peaks', action='store', type=float, help='Tolerance of a m/z and the measured peak m/z')
    argparser.add_argument('--optional_peaks', action='store', help='Optional peaks that can show up')

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

    has_MZ_check = False
    if not (params.precursor_mz is None or params.precursor_mz == 0.0):
        has_MZ_check = True

    # ensures there is always a tolerance
    if params.tolerance_mz is None or params.tolerance_mz == 0.0:
        params.tolerance_mz = 0.01

    if params.tolerance_peaks is None or params.tolerance_peaks == 0.0:
        params.tolerance_peaks = 0.01

    if has_MZ_check:
        lower_bound = params.precursor_mz - params.tolerance_mz
        upper_bound = params.precursor_mz + params.tolerance_mz

    need_to_check_peaks = False
    # if not (params.required_peaks is None or params.required_peaks == ""):
    if params.required_peaks is not None and params.required_peaks != "":
        required_peaks_list = [float(i) for i in params.required_peaks.split(',')]
        need_to_check_peaks = True
        required_peaks = len(required_peaks_list)

    if params.optional_peaks is not None and params.optional_peaks != "":
        optional_peaks_list = [float(i) for i in params.optional_peaks.split(',')]
        optional_peaks = len(optional_peaks_list)
    else:
        optional_peaks = 0

    valid_spectra = []

    #### Read spectra from the file
    t0 = timeit.default_timer()
    stats = { 'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0 }
    with mzml.read(params.mzml_file) as reader:
        for spectrum in reader:
            if stats['counter'] == 0:
                auxiliary.print_tree(spectrum)

            #### Update counters and print progress
            stats['counter'] += 1

            satisfy_requirements = True
            if spectrum['ms level'] == 1:
                stats['ms1spectra'] += 1
            elif spectrum['ms level'] == 2:
                stats['ms2spectra'] += 1
                precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                if has_MZ_check and (precursor_mz < lower_bound or precursor_mz > upper_bound):
                    satisfy_requirements = False
                if need_to_check_peaks and satisfy_requirements:
                    matching_peaks = []
                    for peak in required_peaks_list:
                        lower_peak_bound = peak - params.tolerance_peaks
                        upper_peak_bound = peak + params.tolerance_peaks
                        if spectrum['m/z array'][len(spectrum['m/z array']) - 1] < lower_peak_bound:
                            satisfy_requirements = False
                        else:
                            current_peak = []
                            valid = False
                            for index in range(len(spectrum['m/z array'])):
                                mz = spectrum['m/z array'][index]
                                if mz < lower_peak_bound:
                                    continue
                                elif mz > upper_peak_bound:
                                    if not valid:
                                        satisfy_requirements = False
                                    break
                                else:
                                    intensity = spectrum['intensity array'][index]
                                    current_peak.append([mz, intensity])
                                    valid = True
                                    continue
                            matching_peaks.append(current_peak)

                    for peak in optional_peaks_list:
                        lower_peak_bound = peak - params.tolerance_peaks
                        upper_peak_bound = peak + params.tolerance_peaks
                        if spectrum['m/z array'][len(spectrum['m/z array']) - 1] < lower_peak_bound:
                            satisfy_requirements = False
                        else:
                            current_peak = []
                            for index in range(len(spectrum['m/z array'])):
                                mz = spectrum['m/z array'][index]
                                if mz < lower_peak_bound:
                                    continue
                                elif mz > upper_peak_bound:
                                    break
                                else:
                                    intensity = spectrum['intensity array'][index]
                                    current_peak.append([mz, intensity])
                                    continue
                            matching_peaks.append(current_peak)
                        
                # compute: how many peaks are in the bin? -> len(matchingPeaks)
                # what is the average m/z? -> take the average of matchingPeaks
                # what is the sum of the intensities?
                if satisfy_requirements:
                    spectra_data = ['HFX_9850_GVA_DLD1_2_180719', spectrum['index'], spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']]
                    if need_to_check_peaks:
                        spectra_data.append(matching_peaks)
                    valid_spectra.append(spectra_data)

            if stats['counter']/1000 == int(stats['counter']/1000):
                print(f"  {stats['counter']}")

    #### Print final timing information
    t1 = timeit.default_timer()
    print(f"INFO: Read {stats['counter']} spectra from {params.mzml_file}")
    print(f"The number of ms1spectra is {stats['ms1spectra']}")
    print(f"The number of ms2spectra is {stats['ms2spectra']}")
    print(f"Number of valid spectra: {len(valid_spectra)}")
    # write out mz and intensity array
    if need_to_check_peaks:
        print_valid_spectra_peaks(valid_spectra, required_peaks, optional_peaks)
    else: 
        print_valid_spectra(valid_spectra)
    print(f"INFO: Elapsed time: {t1-t0}")
    print(f"INFO: Processed {stats['counter']/(t1-t0)} spectra per second")

def print_valid_spectra_peaks(valid_spectra, required_peaks, optional_peaks):
    with open('valid_spectra.csv', 'w') as file:
        writer = csv.writer(file, delimiter='\t', lineterminator='\n') # separated by tabs
        # writer = csv.writer(file) # separated by commas

        spectra_table = []
        header = ['MS run name', 'scan number', 'precursor m/z']
        for index in range(required_peaks + optional_peaks):
            header.append(f'peak {index + 1} n matches')
            header.append(f'peak {index + 1} tallest m/z')
            header.append(f'peak {index + 1} intensity')
        writer.writerow(header)
        # at the end, add the number of peaks in each bin, the average m/z, and the sum of the intensities

        spectra_table.append(header)
        for spectra in valid_spectra:
            print(spectra)
            spectra_data = [spectra[0]]
            spectra_data.append(spectra[1])
            spectra_data.append(spectra[2])

            for current_peak in spectra[3]:
                spectra_data.append(len(current_peak))
                if len(current_peak) != 0:
                    tallest_peak = current_peak[0]
                    for matched_peak in current_peak:
                        if tallest_peak[1] < matched_peak[1]: # change to [0] if comparing m/z instead
                            tallest_peak = matched_peak
                    spectra_data.append(tallest_peak[0])
                    spectra_data.append(tallest_peak[1])
                    writer.writerow(spectra_data)
                else:
                    spectra_data.append('')
                    spectra_data.append('')
                    writer.writerow(spectra_data)
                spectra_table.append(spectra_data)

        print(tabulate(spectra_table, headers='firstrow', tablefmt='grid'))

def print_valid_spectra(valid_spectra):
    with open('valid_spectra.csv', 'w') as file:
        writer = csv.writer(file, delimiter='\t', lineterminator='\n') # separated by tabs
        # writer = csv.writer(file) # separated by commas
        spectra_table = []
        header = ['MS run name', 'scan number', 'precursor m/z']
        writer.writerow(header)

        spectra_table.append(header)
        for spectra in valid_spectra:
            spectra_data = [spectra[0]]
            spectra_data.append(spectra[1])
            spectra_data.append(spectra[2])
            writer.writerow(spectra_data)
            spectra_table.append(spectra_data)

    print(tabulate(spectra_table, headers='firstrow', tablefmt='grid'))

if __name__ == "__main__": main()