#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# python correlate_peaks.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --peak_file peak_pairs.txt
# python correlate_peaks.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --peak_file correlations_with_105.0659.tsv
# first 3, 5 are correlated, 4th and last one is not

import os
import argparse
import os.path
import timeit
import matplotlib.pyplot as plt
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup

from pyteomics import mzml, auxiliary
from matplotlib.backends.backend_pdf import PdfPages

class MSRunSpectraFinder:

    def __init__(self):
        # sets up a module with a description
        argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
        argparser.add_argument('--pxd', action='store', help='ProteomeXchange Dataset identifier')
        argparser.add_argument('--mzml_file', action='store', help='Name of the mzML file to read')
        argparser.add_argument('--peak_file', action='store', help='peaks that must show up in a file')
        argparser.add_argument('--peak_tolerance', action='store', type=float, help='Tolerance of a m/z and the measured peak m/z')

        # runs the parse_args method to store all the results into the params variable
        params = argparser.parse_args()

        if params.mzml_file is None or params.mzml_file == "":
            print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
            return

        # tries to open the params file specified
        if not os.path.isfile(params.mzml_file):
            print(f"ERROR: File '{params.mzml_file}' not found or not a file")
            return

        if params.peak_file is None or params.peak_file == "":
            print('ERROR: Parameter --peak_file must be provided. See --help for more information')
            return

        # tries to open the params file specified
        if not os.path.isfile(params.peak_file):
            print(f"ERROR: File '{params.peak_file}' not found or not a file")
            return
        
        self.mzml_file = params.mzml_file
        self.peak_file = params.peak_file[0:len(params.peak_file) - 4]

        self.correlate_pairs = []
        with open(params.peak_file) as file:
            lines = [line.rstrip() for line in file]
            for pairs in lines:
                line_split = [i for i in pairs.split("\t")]
                pair = [float(line_split[0]), line_split[1], float(line_split[2])]
                self.correlate_pairs.append(pair)

        self.has_pxd = False
        if params.pxd is not None and params.pxd != "":
            self.has_pxd = True

        # ensures there is always a tolerance
        if params.peak_tolerance is None or params.peak_tolerance == 0.0:
            self.peak_tolerance = 0.01
        else:
            self.peak_tolerance = params.peak_tolerance

        self.peak_known_intensity = []
        self.peak_unknown_intensity = []

        self.lower_bound = 0
        self.upper_bound = 0

        #### Read spectra from the file
        self.t0 = timeit.default_timer()
        self.stats_counter = 0

    def find_spectra(self, pair):
        peak_known_intensity = []
        peak_unknown_intensity = []
        with mzml.read(self.mzml_file) as reader:
            for spectrum in reader:
                #### Update counters and print progress
                self.stats_counter += 1

                if spectrum['ms level'] == 2:
                    
                    self.lower_bound = pair[0] - self.peak_tolerance
                    self.upper_bound = pair[0] + self.peak_tolerance
                    if spectrum['m/z array'][len(spectrum['m/z array']) - 1] < self.lower_bound:
                        continue
                    else:
                        for index in range(len(spectrum['m/z array'])):
                            mz = spectrum['m/z array'][index]
                            if mz < self.lower_bound:
                                continue
                            elif mz > self.upper_bound:
                                break
                            else:
                                current_intensity = spectrum['intensity array'][index]

                                lower_peak_bound = pair[2] - self.peak_tolerance
                                upper_peak_bound = pair[2] + self.peak_tolerance
                                if spectrum['m/z array'][len(spectrum['m/z array']) - 1] < lower_peak_bound:
                                        continue
                                else:
                                    for index in range(len(spectrum['m/z array'])):
                                        mz = spectrum['m/z array'][index]
                                        if mz < lower_peak_bound:
                                            continue
                                        elif mz > upper_peak_bound:
                                            break
                                        else:
                                            known_intensity = spectrum['intensity array'][index]
                                            peak_unknown_intensity.append(current_intensity)
                                            peak_known_intensity.append(known_intensity)
                                            continue

                    if self.stats_counter/1000 == int(self.stats_counter/1000):
                        print(f"  {self.stats_counter}")
        self.peak_unknown_intensity.append(peak_unknown_intensity)
        self.peak_known_intensity.append(peak_known_intensity)
               
    def show_stats(self):
        #### Print final timing information
        t1 = timeit.default_timer()
        print(f"INFO: Read {self.stats['counter']} spectra from {self.mzml_file}")
        print(f"INFO: Elapsed time: {t1-self.t0}")
        print(f"INFO: Processed {self.stats_counter/(t1-self.t0)} spectra per second")

    def correlate_peaks(self):
        x = 0
        y = 0
        pdf = PdfPages(f'{self.peak_file}.pdf')
        fig, ax = plt.subplots(5, 3, figsize=(8.5,11))
        for index in range(len(self.correlate_pairs)):
            fig.subplots_adjust(top=0.98, bottom=0.05, left=0.12, right=0.98, hspace=0.2, wspace=0.38)
            ax[x, y].scatter(self.peak_unknown_intensity[index], self.peak_known_intensity[index], 0.5)
            ax[x, y].set_xscale("log")
            ax[x, y].set_yscale("log")
            ax[x, y].tick_params(axis='x', labelsize='xx-small')
            ax[x, y].tick_params(axis='y', labelsize='xx-small')
            ax[x, y].text(min(self.peak_unknown_intensity[index]), max(self.peak_known_intensity[index]) * 1.03, f'unknown: \n    {self.correlate_pairs[index][0]}\nknown: \n    {self.correlate_pairs[index][1]}, {self.correlate_pairs[index][2]}', fontsize='xx-small', ha='left', va='top',)

            y += 1
            y %= 3
            if y == 0:
                x += 1
                x %= 5
                if x == 0:
                    pdf.savefig(fig)
                    plt.close('all')
                    fig, ax = plt.subplots(self.rows, self.columns,figsize=(8.5,11))
        
        if x != 0:
            pdf.savefig(fig)
        pdf.close()
        plt.close()

def main():

    spectra_finder = MSRunSpectraFinder()

    for pair in spectra_finder.correlate_pairs:
        # find all scans with both peaks
        spectra_finder.find_spectra(pair)
        
    # print out intensities as a scatterplot
    spectra_finder.correlate_peaks()

    spectra_finder.show_stats()

if __name__ == "__main__": main()
