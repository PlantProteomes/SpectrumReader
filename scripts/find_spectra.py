#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python read_mzML1.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML
# can do python read_mzML1.py --help and look at the usage statements for what I can do

# pxd is PXD019252
# http://proteomecentral.proteomexchange.org/usi/

# test case: 506.888, 5 valids
# python find_spectra.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --precursor_mz 506.888 --tolerance_mz 0.1
# test case: 506.888, required peaks 244.166 and 245.168, 2 valids
# python find_spectra.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --precursor_mz 506.888 --tolerance_mz 0.1 --required_windows "244.166, 245.168" --tolerance_peaks 0.1
# test case: 395.1828, required peaks 217.08191, 217.54974, 217.5828, 7 valids
# python find_spectra.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --precursor_mz 395.1828 --tolerance_mz 0.1 --required_windows "217.08191, 217.54974, 217.5828" --tolerance_peaks 0.1
# test case: 387.1953, required peaks 129.10266, 402.17227, 7 valids
# python find_spectra.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --precursor_mz 387.1953 --tolerance_mz 0.1 --required_windows "129.10266, 402.17227" --tolerance_peaks 0.1
# test case: 387.1953, required peaks 129.10266, minimum optional peak 0, optional peak 402.17227, 29 valids
# python find_spectra.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --precursor_mz 387.1953 --tolerance_mz 0.1 --required_windows "129.10266" --tolerance_peaks 0.1 --minimum_optional_windows 0 --optional_windows "402.17227"
# test case: 387.1953, required peaks 129.10266, minimum optional peak 1, optional peak 402.17227, 594.26514, 7 valids
# python find_spectra.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --precursor_mz 387.1953 --tolerance_mz 0.1 --required_windows "129.10266" --tolerance_peaks 0.1 --minimum_optional_windows 1 --optional_windows "402.17227, 594.26514"

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

class MSRunSpectraFinder:

    def __init__(self):
        # sets up a module with a description
        argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
        argparser.add_argument('--pxd', action='store', help='ProteomeXchange Dataset identifier')
        argparser.add_argument('--mzml_file', action='store', help='Name of the mzML file to read')
        argparser.add_argument('--precursor_mz', action='store', type=float, help='Precursor m/z to select')
        argparser.add_argument('--tolerance_mz', action='store', type=float, help='Tolerance of a m/z and the measured peak m/z')
        argparser.add_argument('--required_windows', action='store', help='Required windows that must show up')
        argparser.add_argument('--tolerance_peaks', action='store', type=float, help='Tolerance of a m/z and the measured peak m/z')
        argparser.add_argument('--minimum_optional_windows', action='store', type=int, help='Minimum number of optional peaks that show up')
        argparser.add_argument('--optional_windows', action='store', help='Optional peaks that can show up')

        # runs the parse_args method to store all the results into the params variable
        params = argparser.parse_args()

        if params.mzml_file is None or params.mzml_file == "":
            print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
            return

        # tries to open the params file specified
        if not os.path.isfile(params.mzml_file):
            print(f"ERROR: File '{params.mzml_file}' not found or not a file")
            return
        
        self.mzml_file = params.mzml_file
        self.root_mzml = params.mzml_file[params.mzml_file.rfind("\\")+1: params.mzml_file.rfind(".")]

        self.has_pxd = False
        if params.pxd is not None and params.pxd != "":
            self.has_pxd = True

        self.has_MZ_check = False
        if not (params.precursor_mz is None or params.precursor_mz == 0.0):
            self.has_MZ_check = True
            self.precursor_mz = params.precursor_mz

        # ensures there is always a tolerance
        if params.tolerance_mz is None or params.tolerance_mz == 0.0:
            self.tolerance_mz = 0.01
        else:
            self.tolerance_mz = params.tolerance_mz

        if params.tolerance_peaks is None or params.tolerance_peaks == 0.0:
            self.tolerance_peaks = 0.01
        else:
            self.tolerance_peaks = params.tolerance_peaks

        if self.has_MZ_check:
            self.lower_bound = self.precursor_mz - self.tolerance_mz
            self.upper_bound = self.precursor_mz + self.tolerance_mz

        self.need_to_check_windows = False
        if params.required_windows is not None and params.required_windows != "":
            self.required_windows_list = [float(i) for i in params.required_windows.split(',')]
            self.need_to_check_windows = True
            self.required_windows = len(self.required_windows_list)

        if params.minimum_optional_windows is None or params.minimum_optional_windows == "":
            self.minimum_optional_windows = 0
        else:
            self.minimum_optional_windows = params.minimum_optional_windows

        if params.optional_windows is not None and params.optional_windows != "":
            self.optional_windows_list = [float(i) for i in params.optional_windows.split(',')]
            self.optional_windows = len(self.optional_windows_list)
            if self.optional_windows < self.minimum_optional_windows:
                print(f"ERROR: Minimum optional windows is too small")
                return
        else:
            self.optional_windows = 0

        self.valid_spectra = []

        #### Read spectra from the file
        self.t0 = timeit.default_timer()
        self.stats = { 'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0 }

    def find_spectra(self):
        with mzml.read(self.mzml_file) as reader:
            for spectrum in reader:
                # if stats['counter'] == 0:
                    # auxiliary.print_tree(spectrum)

                #### Update counters and print progress
                self.stats['counter'] += 1

                satisfy_requirements = True
                if spectrum['ms level'] == 1:
                    self.stats['ms1spectra'] += 1
                elif spectrum['ms level'] == 2:
                    self.stats['ms2spectra'] += 1
                    precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                    if self.has_MZ_check and (precursor_mz < self.lower_bound or precursor_mz > self.upper_bound):
                        satisfy_requirements = False
                    if self.need_to_check_windows and satisfy_requirements:
                        matching_peaks = []
                        for peak in self.required_windows_list:
                            lower_peak_bound = peak - self.tolerance_peaks
                            upper_peak_bound = peak + self.tolerance_peaks
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

                        if self.optional_windows > 0:
                            optional_windows_found = []
                            optionals_found = 0
                            for peak in self.optional_windows_list:
                                optional_found = False
                                lower_peak_bound = peak - self.tolerance_peaks
                                upper_peak_bound = peak + self.tolerance_peaks
                                current_peak = []
                                if spectrum['m/z array'][len(spectrum['m/z array']) - 1] < lower_peak_bound:
                                    optional_windows_found.append(current_peak)
                                else:
                                    for index in range(len(spectrum['m/z array'])):
                                        mz = spectrum['m/z array'][index]
                                        if mz < lower_peak_bound:
                                            continue
                                        elif mz > upper_peak_bound:
                                            break
                                        else:
                                            intensity = spectrum['intensity array'][index]
                                            if not optional_found:
                                                optionals_found += 1
                                                optional_found = True
                                            current_peak.append([mz, intensity])
                                            continue
                                    optional_windows_found.append(current_peak)

                            if optionals_found < self.minimum_optional_windows:
                                satisfy_requirements = False
                            else:
                                for peak in optional_windows_found:
                                    matching_peaks.append(peak)
                            
                    # find a way to use metadata of the scan, pull out the true scan number and not by adding 10000
                    scan_number = 10000 + spectrum['index']
                    if satisfy_requirements and self.has_pxd:
                        spectra_data = [f'mzspec:{self.pxd}:{self.root_mzml}:scan:{scan_number}', spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']]
                        if self.need_to_check_windows:
                            spectra_data.append(matching_peaks)
                        self.valid_spectra.append(spectra_data)
                    elif satisfy_requirements:
                        spectra_data = [self.root_mzml, scan_number, spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']]
                        if self.need_to_check_windows:
                            spectra_data.append(matching_peaks)
                        self.valid_spectra.append(spectra_data)

                if self.stats['counter']/1000 == int(self.stats['counter']/1000):
                    print(f"  {self.stats['counter']}")

    def show_stats(self):
        #### Print final timing information
        t1 = timeit.default_timer()
        print(f"INFO: Read {self.stats['counter']} spectra from {self.mzml_file}")
        print(f"The number of ms1spectra is {self.stats['ms1spectra']}")
        print(f"The number of ms2spectra is {self.stats['ms2spectra']}")
        print(f"Number of valid spectra: {len(self.valid_spectra)}")
        print(f"INFO: Elapsed time: {t1-self.t0}")
        print(f"INFO: Processed {self.stats['counter']/(t1-self.t0)} spectra per second")

    def print_valid_spectra_peaks_pxd(self):
        with open('valid_spectra.tsv', 'a') as file:
            writer = csv.writer(file, delimiter='\t', lineterminator='\n') # separated by tabs
            # writer = csv.writer(file) # separated by commas
            # usi - of the letters mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555
            # if PXD is specified, then make column 1 in the format above ^

            spectra_table = []
            header = ['USI', 'precursor m/z']
            for index in range(self.required_windows + self.optional_windows):
                header.append(f'peak {index + 1} n matches')
                header.append(f'peak {index + 1} tallest m/z')
                header.append(f'peak {index + 1} tallest intensity')
            writer.writerow(header)
            # at the end, add the number of peaks in each bin, the average m/z, and the sum of the intensities

            spectra_table.append(header)
            for spectra in self.valid_spectra:
                # print(spectra)
                spectra_data = [spectra[0]]
                spectra_data.append(round(spectra[1], 4))

                for current_peak in spectra[2]:
                    spectra_data.append(len(current_peak))
                    if len(current_peak) != 0:
                        tallest_peak = current_peak[0]
                        for matched_peak in current_peak:
                            if tallest_peak[1] < matched_peak[1]: # change to [0] if comparing m/z instead
                                tallest_peak = matched_peak
                        spectra_data.append(round(tallest_peak[0], 4))
                        spectra_data.append(int(tallest_peak[1]))
                    else:
                        spectra_data.append('')
                        spectra_data.append('')
                
                # print(tabulate(spectra_table, headers='firstrow', tablefmt='grid'))
                writer.writerow(spectra_data)
                spectra_table.append(spectra_data)

    def print_valid_spectra_peaks(self):
        with open('valid_spectra.tsv', 'a') as file:
            writer = csv.writer(file, delimiter='\t', lineterminator='\n') # separated by tabs
            # writer = csv.writer(file) # separated by commas
            # usi - of the letters mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555
            # if PXD is specified, then make column 1 in the format above ^

            spectra_table = []
            header = ['MS run name', 'scan number', 'precursor m/z']
            for index in range(self.required_windows + self.optional_windows):
                header.append(f'peak {index + 1} n matches')
                header.append(f'peak {index + 1} tallest m/z')
                header.append(f'peak {index + 1} tallest intensity')
            writer.writerow(header)
            # at the end, add the number of peaks in each bin, the average m/z, and the sum of the intensities

            spectra_table.append(header)
            for spectra in self.valid_spectra:
                # print(spectra)
                spectra_data = [spectra[0]]
                spectra_data.append(spectra[1])
                spectra_data.append(round(spectra[2], 4))

                for current_peak in spectra[3]:
                    spectra_data.append(len(current_peak))
                    if len(current_peak) != 0:
                        tallest_peak = current_peak[0]
                        for matched_peak in current_peak:
                            if tallest_peak[1] < matched_peak[1]: # change to [0] if comparing m/z instead
                                tallest_peak = matched_peak
                        spectra_data.append(round(tallest_peak[0], 4))
                        spectra_data.append(int(tallest_peak[1]))
                    else:
                        spectra_data.append('')
                        spectra_data.append('')
                
                # print(tabulate(spectra_table, headers='firstrow', tablefmt='grid'))
                writer.writerow(spectra_data)
                spectra_table.append(spectra_data)

    def print_valid_spectra_pxd(self):
        with open('valid_spectra.tsv', 'a') as file:
            writer = csv.writer(file, delimiter='\t', lineterminator='\n') # separated by tabs
            # writer = csv.writer(file) # separated by commas
            spectra_table = []
            header = ['USI', 'precursor m/z']
            writer.writerow(header)

            spectra_table.append(header)
            for spectra in self.valid_spectra:
                spectra_data = [spectra[0]]
                spectra_data.append(round(spectra[1], 4))
                writer.writerow(spectra_data)
                spectra_table.append(spectra_data)
        # print(tabulate(spectra_table, headers='firstrow', tablefmt='grid'))

    def print_valid_spectra(self):
        with open('valid_spectra.tsv', 'w') as file:
            writer = csv.writer(file, delimiter='\t', lineterminator='\n') # separated by tabs
            # writer = csv.writer(file) # separated by commas
            spectra_table = []
            header = ['MS run name', 'scan number', 'precursor m/z']
            writer.writerow(header)

            spectra_table.append(header)
            for spectra in self.valid_spectra:
                spectra_data = [spectra[0]]
                spectra_data.append(spectra[1])
                spectra_data.append(round(spectra[2], 4))
                writer.writerow(spectra_data)
                spectra_table.append(spectra_data)
        # print(tabulate(spectra_table, headers='firstrow', tablefmt='grid'))

def main():

    spectra_finder = MSRunSpectraFinder()

    spectra_finder.find_spectra()

    # write out mz and intensity array
    if spectra_finder.need_to_check_windows and spectra_finder.has_pxd:
        spectra_finder.print_valid_spectra_peaks_pxd()
    elif spectra_finder.need_to_check_windows and not spectra_finder.has_pxd:
        spectra_finder.print_valid_spectra_peaks()
    elif spectra_finder.has_pxd: 
        spectra_finder.print_valid_spectra_pxd()
    else:
        spectra_finder.print_valid_spectra()

    spectra_finder.show_stats()

if __name__ == "__main__": main()
