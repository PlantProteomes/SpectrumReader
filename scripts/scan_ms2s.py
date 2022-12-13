#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python scan_ms2s.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --rows 3 --columns 5

# bad outputs for some of these, ex: peak 228.1349

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
from scipy.optimize import curve_fit
from numpy import exp

class MSRunPeakFinder:

    def __init__(self, file_name, rows, columns, tolerance):

        # variables declaration
        self.file_name = file_name
        self.by_count = np.zeros((4000000,), dtype=int)
        self.by_intensity = np.zeros((4000000,), dtype=int)
        self.by_strength = np.zeros((4000000,), dtype=int)
        self.all_peaks_intensities = [] # keeps track of all intensities
        self.observed_peaks = [] # keeps track of all peaks over a certain intensities, identified and unidentified
        self.known_ions = []
        self.peak_correction_factor = 0.0004 # eventually make this not a constant
        self.t0 = timeit.default_timer()
        self.stats = { 'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0 }

        # for pdf graph output
        # rename to plot_output_rows/columns
        if rows is None or rows <= 0:
            self.rows = 5
        else:
            self.rows = rows
        # if columns is None or columns <= 0:
        self.columns = 3
        # else:
            # self.columns = columns
        if tolerance is None or tolerance <= 0:
            self.tolerance = 5
        else:
            self.tolerance = tolerance

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
                    self.smallest_peak_intensity = sys.maxsize
                    for index in range(len(spectrum['m/z array'])):
                        peak_mz = spectrum['m/z array'][index]
                        if peak_mz < 120:
                            peak_mz -= (120 - peak_mz) * 0.00001
                        peak_mz -= self.peak_correction_factor
                        if peak_mz > 400:
                            break
                        else:
                            intensity = spectrum['intensity array'][index]
                            self.all_peaks_intensities.append(intensity)
                            self.by_count[int(10000 * peak_mz + 0.5)] += 1
                            self.by_intensity[int(10000 * peak_mz + 0.5)] += intensity
                            self.smallest_peak_intensity = min(self.smallest_peak_intensity, intensity)

                    # compare intensities to the smallest intensity
                    for index in range(len(spectrum['m/z array'])):
                        peak_mz = spectrum['m/z array'][index]
                        peak_mz -= self.peak_correction_factor
                        if peak_mz < 120:
                            peak_mz -= (120 - peak_mz) * 0.00001
                        if peak_mz > 400:
                            break
                        else:
                            intensity = spectrum['intensity array'][index]
                            self.by_strength[int(10000 * peak_mz + 0.5)] += get_strength(intensity, self.smallest_peak_intensity)

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
        }

        for key in amino_acid_modifications:
            amino_acids.append(key)

        amino_acids.sort()

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
                'K': [ '+CO-NH3', '-NH3', '+CO', '-C2H4-NH3', '+CO+H2ON2', '-NH3', '-C4H7N', '+CO+CO-C2H3N3', '+CO+H2O', '+CO+H2O-NH3'],
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
        for acid in amino_acids:
            for ion in ion_types:
                if len(acid) == 1:
                    base_acid_mass = mass.calculate_mass(sequence=acid, ion_type=ion, charge=1)
                else:
                    base_acid_mass = mass.calculate_mass(sequence=acid[0], ion_type=ion, charge=1) + amino_acid_modifications[acid]['mz']

                base_acid_mass = int(base_acid_mass * 100000 + 0.5) / 100000.0

                # add possible modifications, but only for a ions
                if ion == 'a':
                    # add base a ion
                    self.known_ions.append([base_acid_mass, 'I' + acid, False])
                    # add base a ion's isotope
                    self.known_ions.append([int(100000 * (base_acid_mass + 1.0034)) / 100000, 'I' + acid + '+i', False])
                    # check if there are any possible modifications, if there are, add them too
                    for modification in modification_deltas[acid[0]]:
                        acid_mass = base_acid_mass
                        acid_mass += modification_deltas[acid[0]][modification]
                        acid_mass = int(acid_mass * 100000 + 0.5) / 100000.0
                        self.known_ions.append([acid_mass, 'I' + acid + modification, False])
                # add b and y ions
                else:
                    self.known_ions.append([base_acid_mass, ion + '-' + acid, False])
                    self.known_ions.append([int(100000 * (base_acid_mass + 1.0034)) / 100000, ion + '-' + acid + '+i', False])

        # double nested for loops to identify pairs of amino acids (i.e. A and A)
        # first amino acid
        for identifier_index_1 in range(len(amino_acids)):
            for ion in ion_types:
                pair_acid_1 = amino_acids[identifier_index_1]
                # checks if it is Carbamidomethyl or Oxidation
                if len(pair_acid_1) == 1:
                    pair_mass_1 = mass.calculate_mass(sequence=pair_acid_1, ion_type=ion, charge=1)
                else:
                    pair_mass_1 = amino_acid_modifications[pair_acid_1]['mz'] + mass.calculate_mass(sequence=pair_acid_1[0], ion_type=ion, charge=1)
                
                # second amino acid
                for identifier_index_2 in range(len(amino_acids) - identifier_index_1):
                    identifier_index_2 += identifier_index_1
                    pair_acid_2 = amino_acids[identifier_index_2]
                    if len(pair_acid_2) == 1:
                        pair_mass_2 = mass.calculate_mass(sequence=pair_acid_2, ion_type='b', charge=0)
                    else:
                        pair_mass_2 = amino_acid_modifications[pair_acid_2]['mz'] + mass.calculate_mass(sequence=pair_acid_2[0], ion_type='b', charge=0)
                    
                    # add two amino acids together
                    pair_mass = int(100000 * (pair_mass_1 + pair_mass_2) + 0.5) / 100000
                    self.known_ions.append([pair_mass, ion + '-' + pair_acid_1 + pair_acid_2, False])
                    # if pair_acid_1 + pair_acid_2 == 'RR':
                        # print(f'1: {pair_mass_1}, 2: {pair_mass_2}, ion: {ion}')

        # triple nested for loops to identify trios of amino acids (i.e. A and A)
        # first amino acid
        for identifier_index_1 in range(len(amino_acids)):
            for ion in ion_types:
                trio_acid_1 = amino_acids[identifier_index_1]
                # checks if it is Carbamidomethyl or Oxidation
                if len(trio_acid_1) == 1:
                    trio_mass_1 = mass.calculate_mass(sequence=trio_acid_1, ion_type=ion, charge=1)
                else:
                    trio_mass_1 = amino_acid_modifications[trio_acid_1]['mz'] + mass.calculate_mass(sequence=trio_acid_1[0], ion_type=ion, charge=1)
                
                # second amino acid
                for identifier_index_2 in range(len(amino_acids) - identifier_index_1):
                    identifier_index_2 += identifier_index_1
                    trio_acid_2 = amino_acids[identifier_index_2]
                    if len(trio_acid_2) == 1:
                        trio_mass_2 = mass.calculate_mass(sequence=trio_acid_2, ion_type='b', charge=0)
                    else:
                        trio_mass_2 = amino_acid_modifications[trio_acid_2]['mz'] + mass.calculate_mass(sequence=trio_acid_2[0], ion_type='b', charge=0)
                    
                    # third amino acid
                    for identifier_index_3 in range(len(amino_acids) - identifier_index_2):
                        identifier_index_3 += identifier_index_2
                        trio_acid_3 = amino_acids[identifier_index_3]
                        if len(trio_acid_3) == 1:
                            trio_mass_3 = mass.calculate_mass(sequence=trio_acid_3, ion_type='b', charge=0)
                        else:
                            trio_mass_3 = amino_acid_modifications[trio_acid_3]['mz'] + mass.calculate_mass(sequence=trio_acid_3[0], ion_type='b', charge=0)

                        # add all 3 together
                        trio_mass = int(100000 * (trio_mass_1 + trio_mass_2 + trio_mass_3) + 0.5) / 100000
                        if trio_mass <= 400:
                            self.known_ions.append([trio_mass, ion + '-' + trio_acid_1 + trio_acid_2 + trio_acid_3, False])
                        else:
                            continue

        # consider water and ammonia losses
        water = mass.calculate_mass(formula='H2O')
        ammonia = mass.calculate_mass(formula='NH3')
        water_ammonia = water + ammonia
        for index in range(len(self.known_ions)):
            self.known_ions.append([self.known_ions[index][0] - water, self.known_ions[index][1] + '-H2O', False])
            self.known_ions.append([self.known_ions[index][0] - ammonia, self.known_ions[index][1] + '-NH3', False])
            self.known_ions.append([self.known_ions[index][0] - water_ammonia, self.known_ions[index][1] + '-H2O-NH3', False])

        # print(len(self.known_ions))
        # print(self.known_ions[800 : 850])
        self.known_ions.sort()

    def find_peaks(self):
        # scan through peaks to find all peaks with over 50
        previous_peak = [0, 0]
        counted = False
        for i in range(len(self.by_count)):
            if self.by_count[i] >= 30:
                if self.by_count[i] < previous_peak[1] and not counted:
                    if self.by_count[i + 1] < self.by_count[i]:
                        self.observed_peaks.append(previous_peak)
                        counted = True
                previous_peak = [i/10000, self.by_count[i]]
            else:
                if not counted and previous_peak[0] != 0 and self.by_count[i + 1] < self.by_count[i]:
                    self.observed_peaks.append(previous_peak)
                counted = False
                previous_peak = [0, 0]

    def identify_peaks(self):
        # out of the intense peaks, see which ones are identifiable
        # includes all possible identifications
        removable_peaks_index = []
        for index in range(len(self.observed_peaks)):
            peak = self.observed_peaks[index]
            for ion_index in range(len(self.known_ions)):
                amino_acid = self.known_ions[ion_index]
                # check tolerance based on ppm
                peak_tolerance = self.tolerance * peak[0] / 1e6
                if peak[0] > amino_acid[0] - peak_tolerance and peak[0] < amino_acid[0] + peak_tolerance and not amino_acid[2]:
                    identified_peak = [amino_acid[0], int(100000*(amino_acid[0] - peak[0]) + 0.5) / 100000, amino_acid[1]]
                    self.observed_peaks[index].append(identified_peak)
                    self.known_ions[ion_index][2] = True
                elif peak[0] > amino_acid[0] - peak_tolerance and peak[0] < amino_acid[0] + peak_tolerance and amino_acid[2]:
                    removable_peaks_index.append(index)

        for index in range(len(removable_peaks_index)):
            if len(self.observed_peaks[removable_peaks_index[index] - index]) == 2:
                del self.observed_peaks[removable_peaks_index[index] - index]

    def write_output(self):
        with open('common_peaks.tsv', 'w') as file:
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
            
            # add gaussian fitting
            n = len(mz_values)
            center = int(n/2)
            binsize = mz_values[center]-mz_values[center-1]

            try:
            #if 1:
                popt,pcov = curve_fit(gaussian_function,mz_values,intensity_values,p0=[intensity_values[center],mz_values[center],binsize])
                ax[x,y].plot(mz_values,gaussian_function(mz_values,*popt),'r:')
            except:
                continue
            
            peak_fit_center = int(100000*(peak[0] + popt[1])) / 100000

            if len(peak) >= 3:
                ax[x,y].axvline(x=peak[2][0] - peak[0], color='black', lw=1, linestyle='--')
                index = 2
                identified_ion_name = ''
                while index <= len(peak) - 1:
                    identified_ion_name += peak[index][2] + '\n    '
                    index += 1
                ax[x,y].text(-0.0017, max(intensity_values) * 1.03, f'peak fit center: \n    {peak_fit_center}\nion m/z: \n    {peak[2][0]}\nion id: \n    {identified_ion_name}', fontsize='xx-small', ha='left', va='top')
                # ax[x,y].set_title(f"Peak {peak[0]}, Amino Acid {peak[4]}", fontsize="xx-small")
            else:
                ax[x,y].text(-0.0017, max(intensity_values) * 1.03, f'peak fit center: \n    {peak_fit_center}', fontsize='xx-small', ha='left', va='top',)
                # ax[x,y].set_title(f"Peak {peak[0]}", fontsize="xx-small")

            # creates a new figure if full
            y += 1
            y %= self.columns
            if y == 0:
                x += 1
                x %= self.rows
                if x == 0:
                    pdf.savefig(fig)
                    plt.close('all')
                    fig, ax = plt.subplots(self.rows, self.columns,figsize=(8.5,11))
        
        if x != 0:
            pdf.savefig(fig)
        pdf.close()

    def plot_three_histograms(self):
        x = 0
        pdf = PdfPages('common_peaks.pdf')
        fig, ax = plt.subplots(self.rows, self.columns, figsize=(8.5,11))
        for peak in self.observed_peaks:
            fig.subplots_adjust(top=0.98, bottom=0.05, left=0.12, right=0.98, hspace=0.2, wspace=0.38)
            mz_values = []
            # 0 is bycount, 1 is byintensity, 2 is bystrength
            total_intensity_values = [[], [], []]
            for index in range(31):
                # change it to be 0
                add_index = int(peak[0] * 10000 + index - 16)
                mz_values.append(add_index / 10000 - peak[0])
                total_intensity_values[0].append(self.by_count[add_index])
                total_intensity_values[1].append(self.by_intensity[add_index])
                total_intensity_values[2].append(self.by_strength[add_index])
            # change the axis labels to be a smaller text
            for y in range(3):
                ax[x,y].step(mz_values, total_intensity_values[y], where='mid')
                ax[x,y].tick_params(axis='x', labelsize='xx-small')
                ax[x,y].locator_params(axis='x', nbins=5)
                # ax[x,y].set_xticks([-0.0015, -0.0007, 0, 0.0007, 0.0015])
                ax[x,y].tick_params(axis='y', labelsize='xx-small')

                # add gaussian fitting
                n = len(mz_values)
                center = int(n/2)
                binsize = mz_values[center]-mz_values[center-1]

                try:
                    popt,pcov = curve_fit(gaussian_function,mz_values,total_intensity_values[y],p0=[total_intensity_values[y][center],mz_values[center],binsize])
                    ax[x,y].plot(mz_values,gaussian_function(mz_values,*popt),'r:')
                except:
                    continue

                peak_fit_center = int(100000*(peak[0] + popt[1])) / 100000

                if len(peak) >= 3:
                    ax[x,y].axvline(x=peak[2][0] - peak[0], color='black', lw=1, linestyle='--')
                    index = 2
                    identified_ion_name = ''
                    while index <= len(peak) - 1:
                        identified_ion_name += peak[index][2] + '\n    '
                        index += 1
                    ax[x,y].text(-0.0017, max(total_intensity_values[y]) * 1.03, f'peak center: \n    {peak_fit_center}\nion m/z: \n    {peak[2][0]}\nion id: \n    {identified_ion_name}', fontsize='xx-small', ha='left', va='top')
                    # ax[x,y].set_title(f"Peak {peak[0]}, Amino Acid {peak[4]}", fontsize="xx-small")
                else:
                    ax[x,y].text(-0.0017, max(total_intensity_values[y]) * 1.03, f'peak center: \n    {peak_fit_center}', fontsize='xx-small', ha='left', va='top',)
                    # ax[x,y].set_title(f"Peak {peak[0]}", fontsize="xx-small")

            # creates a new figure if full
            x += 1
            x %= self.rows
            if x == 0:
                pdf.savefig(fig)
                plt.close('all')
                fig, ax = plt.subplots(self.rows, self.columns,figsize=(8.5,11))
        
        if x != 0:
            pdf.savefig(fig)
        pdf.close()

    def plot_peaks_strength(self):
        x = 0
        y = 0
        pdf = PdfPages('common_peaks.pdf')
        fig, ax = plt.subplots(self.rows, self.columns, figsize=(8.5,11))
        for peak in self.observed_peaks:
            fig.subplots_adjust(top=0.98, bottom=0.05, left=0.12, right=0.98, hspace=0.2, wspace=0.38)
            mz_values = []
            intensity_values = []
            ppm_values = []
            ppm_delta = 15
            mz_delta = int((ppm_delta / 1e6) * peak[0] * 10000)
            for index in range(mz_delta * 2 + 1):
                # change it to be 0
                add_index = int(peak[0] * 10000 + index - mz_delta - 1)
                mz_values.append(add_index / 10000 - peak[0])
                ppm_values.append((add_index / 10000 - peak[0]) * 1e6 / peak[0])
                intensity_values.append(self.by_strength[add_index])
            # change the axis labels to be a smaller text

            ax[x,y].step(ppm_values, intensity_values, where='mid')
            ax[x,y].tick_params(axis='x', labelsize='xx-small')
            ax[x,y].locator_params(axis='x', nbins=5)
            # ax[x,y].set_xticks([-0.0015, -0.0007, 0, 0.0007, 0.0015])
            ax[x,y].tick_params(axis='y', labelsize='xx-small')

            # add gaussian fitting
            n = len(ppm_values)
            center = int(n/2)
            binsize = ppm_values[center]-ppm_values[center-1]

            try:
            #if 1:
                popt,pcov = curve_fit(gaussian_function,ppm_values,intensity_values,p0=[intensity_values[center],ppm_values[center],binsize])
                ax[x,y].plot(ppm_values,gaussian_function(ppm_values,*popt),'r:')
            except:
                continue

            peak_fit_center = int(100000*(peak[0] + (popt[1] * peak[0] / 1e6))) / 100000

            if len(peak) >= 3:
                line_x = (peak[2][0] - peak[0]) * 1e6 / peak[0]
                ax[x,y].axvline(x=line_x, color='black', lw=1, linestyle='--')
                count = 2
                identified_ion_name = ''
                ion_mz = int(100000 * peak[2][0] + 0.5) / 100000
                while count <= len(peak) - 1 and count <= 10:
                    identified_ion_name += peak[count][2] + '\n    '
                    count += 1
                    if count == 11 and count <= len(peak) - 1:
                        identified_ion_name += "..."
                ax[x,y].text(-16, max(intensity_values) * 1.03, f'peak fit center: \n    {peak_fit_center}\nion m/z: \n    {ion_mz}\nion id: \n    {identified_ion_name}', fontsize='xx-small', ha='left', va='top')
                # ax[x,y].set_title(f"Peak {peak[0]}, Amino Acid {peak[4]}", fontsize="xx-small")
            else:
                ax[x,y].text(-16, max(intensity_values) * 1.03, f'peak fit center: \n    {peak_fit_center}', fontsize='xx-small', ha='left', va='top',)
                # ax[x,y].set_title(f"Peak {peak[0]}", fontsize="xx-small")

            # creates a new figure if full
            y += 1
            y %= self.columns
            if y == 0:
                x += 1
                x %= self.rows
                if x == 0:
                    pdf.savefig(fig)
                    plt.close('all')
                    fig, ax = plt.subplots(self.rows, self.columns,figsize=(8.5,11))
        
        if x != 0:
            pdf.savefig(fig)
        pdf.close()
        plt.close()

    def delta_scatterplot(self):
        mz_values = []
        intensity_values = []
        for peak in self.observed_peaks:
            if len(peak) >= 3:
                mz_values.append(peak[0])
                intensity_values.append(peak[2][1])
            else:
                continue
        plt.scatter(mz_values, intensity_values)
        plt.savefig('delta_scatterplot.pdf')

    def show_stats(self):
        t1 = timeit.default_timer()
        print(f"INFO: Read {self.stats['counter']} spectra from {self.file_name}")
        print(f"The number of ms1spectra is {self.stats['ms1spectra']}")
        print(f"The number of ms2spectra is {self.stats['ms2spectra']}")

        print(f"INFO: Elapsed time: {t1-self.t0}")
        print(f"INFO: Processed {self.stats['counter']/(t1-self.t0)} spectra per second")

    # Gaussian function used for curve fitting
def gaussian_function(x,a,x0,sigma):

    return a*exp(-(x-x0)**2/(2*sigma**2))

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
    argparser.add_argument('--tolerance', action='store', type=float, help='Tolerance for identifying peaks in ppm')

    params = argparser.parse_args()

    peak_finder = MSRunPeakFinder(params.mzml_file, params.rows, params.columns, params.tolerance)
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
    peak_finder.plot_peaks_strength()
    peak_finder.delta_scatterplot()
    # print out data, including run time, number of peaks found, etc
    peak_finder.show_stats()

if __name__ == "__main__": main()
