#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python scan_ms2s.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML


import os
import argparse
import os.path
import timeit
import matplotlib.pyplot as plt
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
import numpy as np
import sys
import csv

from pyteomics import mzml, auxiliary

# focus on the range between 0 and 400 (peak)
# as you read the 10000 spectra, keep track of all the peaks between 0 and 400 - create an array of bins of 0.0001, when you find
# a peak, put it in the closest bin. after cycling through all the spectra, bin counts all the peaks i see, then sort them by
# intensity and find the most common bin
# see what constitutes one specific peak
# top 50 ions, how many there are, and what fraction of spectra contain them
# peak at 129.1026, the average of the peaks, what the standard deviation is of the measurements

def main():

    # sets up a module with a description
    argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
    argparser.add_argument('--mzml_file', action='store', help='Name of the mzML file to read')

    # runs the parse_args method to store all the results into the params variable
    params = argparser.parse_args()
    
    if params.mzml_file is None or params.mzml_file == "":
        print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
        return

    # tries to open the params file specified
    if not os.path.isfile(params.mzml_file):
        print(f"ERROR: File '{params.mzml_file}' not found or not a file")
        return

    root_mzml = params.mzml_file[params.mzml_file.rfind("\\")+1: params.mzml_file.rfind(".")]

    # create three arrays of 4 million for bins
    by_count = np.zeros((4000000,), dtype=int)
    by_intensity = np.zeros((4000000,), dtype=int)
    by_strength = np.zeros((4000000,), dtype=int)

    all_peaks = []
    #### Read spectra from the file
    t0 = timeit.default_timer()
    stats = { 'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0 }
    with mzml.read(params.mzml_file) as reader:
        for spectrum in reader:
            # if stats['counter'] == 0:
                # auxiliary.print_tree(spectrum)

            #### Update counters and print progress
            stats['counter'] += 1

            # if stats['counter'] == 101:
                # break

            satisfy_requirements = True
            if spectrum['ms level'] == 1:
                stats['ms1spectra'] += 1
            elif spectrum['ms level'] == 2:
                stats['ms2spectra'] += 1
                smallest_peak_intensity = sys.maxsize
                for index in range(len(spectrum['m/z array'])):
                    peak = spectrum['m/z array'][index]
                    if peak > 400:
                        break
                    else:
                        intensity = spectrum['intensity array'][index]
                        all_peaks.append(intensity)
                        by_count[int(10000 * peak)] += 1
                        by_intensity[int(10000 * peak)] += intensity
                        smallest_peak_intensity = min(smallest_peak_intensity, intensity)

                for index in range(len(spectrum['m/z array'])):
                    peak = spectrum['m/z array'][index]
                    if peak > 400:
                        break
                    else:
                        intensity = spectrum['intensity array'][index]
                        by_strength[int(10000 * peak)] += get_strength(intensity, smallest_peak_intensity)

    # sort them by size
    # find the average
    # find the standard deviation

            if stats['counter']/1000 == int(stats['counter']/1000):
                print(f"  {stats['counter']}")

    #### Print final timing information
    t1 = timeit.default_timer()
    print(f"INFO: Read {stats['counter']} spectra from {params.mzml_file}")
    print(f"The number of ms1spectra is {stats['ms1spectra']}")
    print(f"The number of ms2spectra is {stats['ms2spectra']}")
    # print out top spectra
    # GET THE HISTOGRAM WORKING!!
    plt.hist(all_peaks, bins=200, range=[0, 20000])
    plt.show()
    # populate tallest_bins with only the tallest ones
    # ex: 101.0714 1577 99204530 2857 <- disregard
    # 101.0715 1789 225788258 3557 <- keep
    # 101.0716 1054 57483130 1875 <- disregard
    with open('popular_spectra.tsv', 'a') as file:
        writer = csv.writer(file, delimiter='\t', lineterminator='\n')
        tallest_bins = []
        for i in range(len(by_count)):
            if by_count[i] > 100:
                print(i/10000, by_count[i], by_intensity[i], by_strength[i])
                data = [i/10000, by_count[i]]
                writer.writerow(data)
    
    print(f"INFO: Elapsed time: {t1-t0}")
    print(f"INFO: Processed {stats['counter']/(t1-t0)} spectra per second")

def get_strength(intensity, smallest_peak_intensity):
    if intensity <= 3 * smallest_peak_intensity:
        return 1
    elif intensity <= 10 * smallest_peak_intensity:
        return 2
    elif intensity <= 30 * smallest_peak_intensity:
        return 3
    else:
        return 4

if __name__ == "__main__": main()
