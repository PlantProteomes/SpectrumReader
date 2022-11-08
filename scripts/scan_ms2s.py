#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python scan_ms2s.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --rows 3 --columns 5

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

from pyteomics import mzml, auxiliary, mass
from collections import OrderedDict
from matplotlib.backends.backend_pdf import PdfPages

class MSRunPeakFinder:

    def __init__(self, file_name, rows, columns):

        # variables declaration
        self.file_name = file_name
        self.by_count = np.zeros((4000000,), dtype=int)
        # self.by_intensity = np.zeros((4000000,), dtype=int)
        # self.by_strength = np.zeros((4000000,), dtype=int)
        self.all_peaks_intensities = [] # keeps track of all intensities
        self.observed_peaks = [] # keeps track of all peaks over a certain intensities, identified and unidentified
        self.known_ions = []
        self.peak_correction_factor = 0.0006 # eventually make this not a constant
        self.t0 = timeit.default_timer()
        self.stats = { 'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0 }

        # for pdf graph output
        # rename to plot_output_rows/columns
        if rows is None or rows <= 0:
            self.rows = 5
        else:
            self.rows = rows
        if columns is None or columns <= 0:
            self.columns = 3
        else:
            self.columns = rows
        # if tolerance is None or tolerance <= 0:
            # self.tolerance = 0.002
        # else:
            # self.olerance = tolerance

        # verify valid file to read
        if file_name is None or file_name == "":
            print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
            return

        # tries to open the params file specified
        if not os.path.isfile(file_name):
            print(f"ERROR: File '{file_name}' not found or not a file")
            return

    def aggregate_spectra(self):   

        with mzml.read(self.file_name) as reader:
            for spectrum in reader:

                #### Update counters and print progress
                self.stats['counter'] += 1

                # start with a smaller sample size
                # if self.stats['counter'] == 101:
                    # break

                # add most popular mz values to be printed out
                if spectrum['ms level'] == 1:
                    self.stats['ms1spectra'] += 1
                elif spectrum['ms level'] == 2:
                    self.stats['ms2spectra'] += 1
                    # smallest_peak_intensity = sys.maxsize
                    for index in range(len(spectrum['m/z array'])):
                        peak = spectrum['m/z array'][index]
                        peak -= self.peak_correction_factor
                        if peak > 400:
                            break
                        else:
                            intensity = spectrum['intensity array'][index]
                            self.all_peaks_intensities.append(intensity)
                            self.by_count[int(10000 * peak + 0.5)] += 1
                            # self.by_intensity[int(10000 * peak + 0.5)] += intensity
                            # self.smallest_peak_intensity = min(smallest_peak_intensity, intensity)

                    # compare intensities to the smallest intensity
                    # for index in range(len(spectrum['m/z array'])):
                        # peak = spectrum['m/z array'][index]
                        # if peak > 400:
                            # break
                        # else:
                            # intensity = spectrum['intensity array'][index]
                            # self.by_strength[int(10000 * peak + 0.5)] += get_strength(intensity, smallest_peak_intensity)

                # updates terminal with the progress of reading peaks
                if self.stats['counter']/1000 == int(self.stats['counter']/1000):
                    print(f"  {self.stats['counter']}")

    def show_intensity_histogram(self):

        plt.hist(self.all_peaks_intensities, bins=200, range=[0, 20000])
        plt.show()

    def get_theoretical_ions(self):

        ion_types = ['a', 'b', 'y']
        amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        amino_acid_modifications = {
            'C[Carbamidomethyl]': {'mz': 57.021464, 'amino acid': 'C'},
            'M[Oxidation]': {'mz': 15.994915, 'amino acid': 'M'}
        } # append this to amino_acids using the key, treat it as another amino acid, everything I do
        # downstream (such as pairs of amino acids), do it with Carbamidomethyl again
        aa_immonium_losses = {
                'G': [],
                'A': [],
                'S': [],
                'P': [],
                'V': [ '-CH2-NH3', '-NH3', '+CO-NH3-CH2' ],
                'T': [ '+CO-NH3'],
                'C': [],
                'L': [ '-C3H6', '-CH2' ],
                'I': [ '-C3H6', '-CH2' ],
                'N': [ '-NH3' ],
                'D': [ '-H2O' ],
                'Q': [ '-CO-NH3', '-NH3', '+CO'],
                'K': [ '+CO-NH3', '-NH3', '+CO', '-C2H4-NH3', '+CO+H2ON2', '-NH3', '-C4H7N', '+CO+CO-C2H3N3', '+CO+H2O'],
                'E': [],
                'M': [ '-C2H2-NH3'],
                'H': [ '-CH2N', '+CO-NH2', '+CO-NH3', '+CO-NH', '+CO+H2O' ],
                'F': [ '-CH3N'],
                'R': [ '-C3H6N2', '-CH5N3', '-CH6N2', '-C2H4N2', '-CH2N2', '-CH3N', '-NH3', '-C4H7N', '+H2O+H2O-N3H7', '+CO+H2O' ],
                'Y': [ '-CO-NH3', '-CH3N' ],
                'W': [ '+CO', '-C4H6N2', '-C2H4N', '-CH3N', '-CHN', '+CO-NH3', '-NH3'],
            }
        modification_deltas = {}

        # calculate the additional deltas from possible explanations in aa_immonium_losses
        # turns formula into numbers
        for modification in aa_immonium_losses:
            delta_masses = {}
            for variant in aa_immonium_losses[modification]:
                delta_mass = 0
                variant_add = variant.split('-')[0].split('+')[1:]
                variant_subtract = variant.split('-')[1:]
                for variant_formula in variant_add:
                    delta_mass += mass.calculate_mass(formula=variant_formula)
                for variant_formula in variant_subtract:
                    delta_mass -= mass.calculate_mass(formula=variant_formula)
                delta_masses[variant] = delta_mass
            modification_deltas[modification] = delta_masses
        
        # populate amino_acid_mass with theoretical ions and their mass
        # when giving the sequence = acid, if possible on line 162, if you give it C[Carbamidomethyl], it will work but most likely it will give an error
        # if it has a bracket, do something special, then add the m/z (57.x) to the end
        for acid in amino_acids:
            for ion in ion_types:
                acid_mass = mass.calculate_mass(sequence=acid, ion_type=ion, charge=1)
                acid_mass = int(acid_mass * 10000 + 0.5) / 10000.0

                # add possible modifications, but only for a ions
                if ion == 'a':
                    # add base a ion
                    self.known_ions.append([acid_mass, 'I' + acid])
                    # add base a ion's isotope
                    self.known_ions.append([int(10000 * (acid_mass + 1.003355)) / 10000, 'I' + acid + '+i'])
                    # check if there are any possible modifications, if there are, add them too
                    for modification in modification_deltas[acid]:
                        acid_mass = mass.calculate_mass(sequence=acid, ion_type=ion, charge=1)
                        acid_mass += modification_deltas[acid][modification]
                        acid_mass = int(acid_mass * 10000 + 0.5) / 10000.0
                        self.known_ions.append([acid_mass, 'I' + acid + modification])
                # add b and y ions
                else:
                    self.known_ions.append([acid_mass, 'I' + acid + ion])

        # another loop to do pairs - do it for a, b, and y for all combinations of 2 amino acids
        # always sort them so the lower one is first (alphanumerically one is first)
        # AS, don't bother with SA, but do include AA or CC
        # have a list of amino_acid_modifications - goes into known_ions list and this pair list
        # treat C[carbamidomethyl] as another amino acid
        # http://www.unimod.org/modifications_list.php? - login as guest, use monoisotopic mass

        sorted(self.known_ions)

    def find_peaks(self):
        # scan through peaks to find all peaks with over 50
        previous_peak = [0, 0]
        counted = False
        for i in range(len(self.by_count)):
            if self.by_count[i] >= 50:
                if self.by_count[i] < previous_peak[1] and not counted:
                    self.observed_peaks.append(previous_peak)
                    counted = True
                previous_peak = [i/10000, self.by_count[i]]
            else:
                if not counted and previous_peak[0] != 0:
                    self.observed_peaks.append(previous_peak)
                counted = False
                previous_peak = [0, 0]

    def identify_peaks(self):
        # out of the intense peaks, see which ones are identifiable
        # includes all possible identifications
        for index in range(len(self.observed_peaks)):
            peak = self.observed_peaks[index]
            for amino_acid in self.known_ions:
                if peak[0] > amino_acid[0] - 0.002 and peak[0] < amino_acid[0] + 0.002:
                    self.observed_peaks[index].append(amino_acid[0])
                    self.observed_peaks[index].append(int(10000*(amino_acid[0] - peak[0]) + 0.5) / 10000)
                    self.observed_peaks[index].append(amino_acid[1])

    def write_output(self):
        with open('popular_spectra.tsv', 'w') as file:
            writer = csv.writer(file, delimiter='\t', lineterminator='\n')
            writer.writerows(self.observed_peaks)

    def show_intense_peaks(self):
        # plots each step graph of an identified peak, going 15 m/z to the left and right
        for peak in self.observed_peaks:
            mz_values = []
            intensity_values = []
            for index in range(31):
                add_index = int(peak[0] * 10000 + index - 16)
                mz_values.append(add_index / 10000)
                intensity_values.append(self.by_count[add_index])
            plt.step(mz_values, intensity_values, where='mid')
            plt.show()

    def plot_intense_peaks(self):
        x = 0
        y = 0
        pdf = PdfPages('common_peaks.pdf')
        fig, ax = plt.subplots(self.rows, self.columns, figsize=(8.5,11))
        for peak in self.observed_peaks:
            fig.subplots_adjust(top=0.98, bottom=0.05, left=0.12, right=0.98, hspace=0.2, wspace=0.38)
            mz_values = []
            intensity_values = []
            for index in range(31):
                # change it to be 0
                add_index = int(peak[0] * 10000 + index - 16)
                mz_values.append(add_index / 10000 - peak[0])
                intensity_values.append(self.by_count[add_index])
            # change the axis labels to be a smaller text
            ax[x,y].step(mz_values, intensity_values, where='mid')
            ax[x,y].tick_params(axis='x', labelsize='xx-small')
            ax[x,y].locator_params(axis='x', nbins=5)
            # ax[x,y].set_xticks([-0.0015, -0.0007, 0, 0.0007, 0.0015])
            ax[x,y].tick_params(axis='y', labelsize='xx-small')
            if len(peak) == 5:
                ax[x,y].axvline(x=peak[3], color='black', lw=1, linestyle='--')
                # show all possible identifications in the text, change the way it is written "IQ"
                ax[x,y].text(-0.0017, 45*max(intensity_values)/64, f'peak center: \n{peak[0]}\nion m/z: \n{peak[2]}\nion id: \n{peak[4][0]} type {peak[4][1]}', fontsize='xx-small')
                # ax[x,y].set_title(f"Peak {peak[0]}, Amino Acid {peak[4][0]} type {peak[4][1]}", fontsize="xx-small")
            else:
                ax[x,y].text(-0.0017, 15*max(intensity_values)/16, f'peak center: \n{peak[0]}', fontsize='xx-small')
                # ax[x,y].set_title(f"Peak {peak[0]}", fontsize="xx-small")

            # creates a new figure if full
            y += 1
            y %= self.columns
            if y == 0:
                x += 1
                x %= self.rows
                if x == 0:
                    pdf.savefig(fig)
                    fig, ax = plt.subplots(self.rows, self.columns,figsize=(8.5,11))
        
        pdf.savefig(fig)
        pdf.close()

    def show_stats(self):
        t1 = timeit.default_timer()
        print(f"INFO: Read {self.stats['counter']} spectra from {self.file_name}")
        print(f"The number of ms1spectra is {self.stats['ms1spectra']}")
        print(f"The number of ms2spectra is {self.stats['ms2spectra']}")

        print(f"INFO: Elapsed time: {t1-self.t0}")
        print(f"INFO: Processed {self.stats['counter']/(t1-self.t0)} spectra per second")

def get_strength(intensity, smallest_peak_intensity):
    if intensity <= 3 * smallest_peak_intensity:
        return 1
    elif intensity <= 10 * smallest_peak_intensity:
        return 2
    elif intensity <= 30 * smallest_peak_intensity:
        return 3
    else:
        return 4

# put main at the end of the program, define identify peaks method first
def main():

    argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
    argparser.add_argument('--mzml_file', action='store', help='Name of the mzML file to read')
    argparser.add_argument('--rows', action='store', type=int, help='Number of rows for graphs')
    argparser.add_argument('--columns', action='store', type=int, help='Number of columns for output')

    params = argparser.parse_args()

    peak_finder = MSRunPeakFinder(params.mzml_file, params.rows, params.columns)
    # find the intensity of each m/z value
    peak_finder.aggregate_spectra()
    # create a histogram of all the intensities between 0 and 20000
    # peak_finder.show_intensity_histogram()
    # create an array of all identifiable theoretical and their masses
    peak_finder.get_theoretical_ions()
    # identify all the spectra with a certain intensity
    peak_finder.find_peaks()
    peak_finder.identify_peaks()
    # save the identified peaks to a file
    peak_finder.write_output()
    # plot graphs of each peak as a pop-up
    # peak_finder.show_intense_peaks()
    # save graphs of each peak to a pdf
    peak_finder.plot_intense_peaks()
    # print out data, including run time, number of peaks found, etc
    peak_finder.show_stats()

if __name__ == "__main__": main()
