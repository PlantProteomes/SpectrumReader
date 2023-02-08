#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python scan_ms2s.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719_subset.mzML --rows 3 --columns 5
# to run the program (use terminal): python scan_ms2s.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML.gz
# to run the program (use terminal): python scan_ms2s.py --mzml_file ..\data\QEX03_210305_PTMscan_wt_IP_63.mzML.gz
# to run the program (use terminal): python scan_ms2s.py --mzml_file ..\data\06CPTAC_BCprospective_W_BI_20161116_BA_f17.mzML.gz
# to run the program (use terminal): python scan_ms2s.py --mzml_file ..\data\Q20181210_06.mzML.gz
# to run the program (use terminal): python scan_ms2s.py --mzml_file ..\data\OR13_A2_20161014_CMP_Ecoli_glycerol_exp_2mg_IMAC_R2.mzML.gz

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
import gzip
import statistics
import warnings

from scipy.optimize import OptimizeWarning
from scipy import interpolate
from pyteomics import mzml, auxiliary, mass
from collections import OrderedDict
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from numpy import exp

warnings.simplefilter("ignore", OptimizeWarning)

class MSRunPeakFinder:

    def __init__(self):

        argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
        argparser.add_argument('--mzml_file', action='store', help='Name of the mzML file to read')
        argparser.add_argument('--rows', action='store', type=int, help='Number of rows for graphs')
        argparser.add_argument('--columns', action='store', type=int, help='Number of columns for output')
        argparser.add_argument('--tolerance', action='store', type=float, help='Tolerance for identifying peaks in ppm')
        argparser.add_argument('--correction_factor', action='store', help='Correction Factor for every peak in ppm')

        params = argparser.parse_args()

        # variables declaration
        self.by_count = np.zeros((4000000,), dtype=int)
        self.by_intensity = np.zeros((4000000,), dtype=float)
        self.by_strength = np.zeros((4000000,), dtype=int)
        self.all_peaks_intensities = [] # keeps track of all intensities
        self.triggered_peaks = []
        self.observed_peaks = [] # keeps track of all peaks over a certain intensities, identified and unidentified
        self.known_ions = []
        self.minimum_triggers = 300
        self.ppm_delta = 15
        self.t0 = timeit.default_timer()
        self.stats = { 'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0 }
        self.initial_tolerance = 20
        self.proton_mass = 1.007276
        self.isotope_mass = 1.00335
        self.crude_correction = 0
        self.has_correction_spline = False

        # for pdf graph output
        # rename to plot_output_rows/columns
        if params.rows is None or params.rows <= 0:
            self.rows = 5
        else:
            self.rows = params.rows
        # if params.columns is None or params.columns <= 0:
        self.columns = 3
        # else:
            # self.columns = params.columns
        if params.tolerance is None or params.tolerance <= 0:
            self.tolerance = 7
        else:
            self.tolerance = params.tolerance

        # verify valid file to read
        if params.mzml_file is None or params.mzml_file == "":
            print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
            return

        # tries to open the params file specified
        if not os.path.isfile(params.mzml_file):
            print(f"ERROR: File '{params.mzml_file}' not found or not a file")
            return

        self.file_name = params.mzml_file

        if self.file_name.endswith('.gz'):
            self.infile = gzip.open(self.file_name)
            self.peak_file = self.file_name[0:len(self.file_name) - 8]
        else:
            self.infile = open(self.file_name, 'rb')
            self.peak_file = self.file_name[0:len(self.file_name) - 5]

    def aggregate_spectra(self):   

        with mzml.read(self.infile) as reader:
            for spectrum in reader:

                #### Update counters and print progress
                self.stats['counter'] += 1

                # add most popular mz values to be printed out
                if spectrum['ms level'] == 1:
                    self.stats['ms1spectra'] += 1
                elif spectrum['ms level'] == 2:
                    self.stats['ms2spectra'] += 1
                    self.smallest_peak_intensity = sys.maxsize
                    for index in range(len(spectrum['m/z array'])):
                        peak_mz = spectrum['m/z array'][index]
                        if peak_mz > 400:
                            break
                        else:
                            intensity = spectrum['intensity array'][index]
                            self.all_peaks_intensities.append(intensity)
                            self.by_count[int(10000 * peak_mz + 0.5)] += 1
                            # self.by_intensity[int(10000 * peak_mz + 0.5)] += intensity
                            self.smallest_peak_intensity = min(self.smallest_peak_intensity, intensity)

                    # compare intensities to the smallest intensity
                    for index in range(len(spectrum['m/z array'])):
                        peak_mz = spectrum['m/z array'][index]
                        if peak_mz > 400:
                            break
                        else:
                            intensity = spectrum['intensity array'][index]
                            self.by_strength[int(10000 * peak_mz + 0.5)] += get_strength(intensity, self.smallest_peak_intensity)

                # updates terminal with the progress of reading peaks
                if self.stats['counter']/1000 == int(self.stats['counter']/1000):
                    print(f"  {self.stats['counter']}")

        print("finished aggregating spectra")

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
            delta_mzs = {}
            for variant in aa_immonium_losses[modification]:
                delta_mz = 0
                variant_add = variant.split('-')[0].split('+')[1:]
                variant_subtract = variant.split('-')[1:]
                for variant_formula in variant_add:
                    delta_mz += mass.calculate_mass(formula=variant_formula)
                for variant_formula in variant_subtract:
                    delta_mz -= mass.calculate_mass(formula=variant_formula)
                delta_mzs[variant] = delta_mz
            modification_deltas[modification] = delta_mzs
        
        # populate amino_acid_mass with theoretical ions and their mass
        for acid in amino_acids:
            for ion in ion_types:
                if len(acid) == 1:
                    base_acid_mz = mass.calculate_mass(sequence=acid, ion_type=ion, charge=1)
                else:
                    base_acid_mz = mass.calculate_mass(sequence=acid[0], ion_type=ion, charge=1) + amino_acid_modifications[acid]['mz']

                base_acid_mz = int(base_acid_mz * 100000 + 0.5) / 100000.0

                # add possible modifications, but only for a ions
                if ion == 'a':
                    # add base a ion
                    self.known_ions.append([base_acid_mz, 'I' + acid, False])
                    # add base a ion's isotope
                    self.known_ions.append([int(100000 * (base_acid_mz + self.isotope_mass)) / 100000, 'I' + acid + '+i', False])
                    # check if there are any possible modifications, if there are, add them too
                    for modification in modification_deltas[acid[0]]:
                        acid_mz = base_acid_mz
                        acid_mz += modification_deltas[acid[0]][modification]
                        acid_mz = int(acid_mz * 100000 + 0.5) / 100000.0
                        self.known_ions.append([acid_mz, 'I' + acid + modification, False])
                # add b and y ions
                else:
                    self.known_ions.append([base_acid_mz, ion + '-' + acid, False])
                    self.known_ions.append([int(100000 * (base_acid_mz + self.isotope_mass)) / 100000, ion + '-' + acid + '+i', False])

        # double nested for loops to identify pairs of amino acids (i.e. A and A)
        # first amino acid
        for identifier_index_1 in range(len(amino_acids)):
            for ion in ion_types:
                pair_acid_1 = amino_acids[identifier_index_1]
                # checks if it is Carbamidomethyl or Oxidation
                if len(pair_acid_1) == 1:
                    pair_mz_1 = mass.calculate_mass(sequence=pair_acid_1, ion_type=ion, charge=1)
                else:
                    pair_mz_1 = amino_acid_modifications[pair_acid_1]['mz'] + mass.calculate_mass(sequence=pair_acid_1[0], ion_type=ion, charge=1)
                
                # second amino acid
                for identifier_index_2 in range(len(amino_acids) - identifier_index_1):
                    identifier_index_2 += identifier_index_1
                    pair_acid_2 = amino_acids[identifier_index_2]
                    if len(pair_acid_2) == 1:
                        pair_mz_2 = mass.calculate_mass(sequence=pair_acid_2, ion_type='b', charge=0)
                    else:
                        pair_mz_2 = amino_acid_modifications[pair_acid_2]['mz'] + mass.calculate_mass(sequence=pair_acid_2[0], ion_type='b', charge=0)
                    
                    # add two amino acids together
                    pair_mz = int(100000 * (pair_mz_1 + pair_mz_2) + 0.5) / 100000
                    self.known_ions.append([pair_mz, ion + '-' + pair_acid_1 + pair_acid_2, False])

        # triple nested for loops to identify trios of amino acids (i.e. A and A)
        # first amino acid
        for identifier_index_1 in range(len(amino_acids)):
            for ion in ion_types:
                trio_acid_1 = amino_acids[identifier_index_1]
                # checks if it is Carbamidomethyl or Oxidation
                if len(trio_acid_1) == 1:
                    trio_mz_1 = mass.calculate_mass(sequence=trio_acid_1, ion_type=ion, charge=1)
                else:
                    trio_mz_1 = amino_acid_modifications[trio_acid_1]['mz'] + mass.calculate_mass(sequence=trio_acid_1[0], ion_type=ion, charge=1)
                
                # second amino acid
                for identifier_index_2 in range(len(amino_acids) - identifier_index_1):
                    identifier_index_2 += identifier_index_1
                    trio_acid_2 = amino_acids[identifier_index_2]
                    if len(trio_acid_2) == 1:
                        trio_mz_2 = mass.calculate_mass(sequence=trio_acid_2, ion_type='b', charge=0)
                    else:
                        trio_mz_2 = amino_acid_modifications[trio_acid_2]['mz'] + mass.calculate_mass(sequence=trio_acid_2[0], ion_type='b', charge=0)
                    
                    # third amino acid
                    for identifier_index_3 in range(len(amino_acids) - identifier_index_2):
                        identifier_index_3 += identifier_index_2
                        trio_acid_3 = amino_acids[identifier_index_3]
                        if len(trio_acid_3) == 1:
                            trio_mz_3 = mass.calculate_mass(sequence=trio_acid_3, ion_type='b', charge=0)
                        else:
                            trio_mz_3 = amino_acid_modifications[trio_acid_3]['mz'] + mass.calculate_mass(sequence=trio_acid_3[0], ion_type='b', charge=0)

                        # add all 3 together
                        trio_mz = int(100000 * (trio_mz_1 + trio_mz_2 + trio_mz_3) + 0.5) / 100000
                        if trio_mz <= 400:
                            self.known_ions.append([trio_mz, ion + '-' + trio_acid_1 + trio_acid_2 + trio_acid_3, False])
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
            if self.known_ions[index][1][0] == 'a' or self.known_ions[index][1][0] == 'b' or 'K' in self.known_ions[index][1]:
                self.known_ions.append([self.known_ions[index][0] + 229.162932, self.known_ions[index][1] + '+TMT', False])

        self.known_ions.append([mass.calculate_mass(formula="C3H8N2O2") + self.proton_mass, "C3H8N2O2", False]) # ppm = 1.18
        self.known_ions.append([mass.calculate_mass(formula="CH5N4O2") + self.proton_mass, "CH5N4O2", False]) # ppm = 2.38
        self.known_ions.append([mass.calculate_mass(formula="C5H12N2O") + self.proton_mass, "C5H12N2O", False]) # ppm = 1.18
        self.known_ions.append([mass.calculate_mass(formula="C4H10N2O2") + self.proton_mass, "C4H10N2O2", False]) # ppm = 1.46
        self.known_ions.append([126.127726, "TMT126", False])
        self.known_ions.append([127.124761, "TMT127N", False])
        self.known_ions.append([127.131081, "TMT127C", False])
        self.known_ions.append([128.128116, "TMT128N", False])
        self.known_ions.append([128.134436, "TMT128C", False])
        self.known_ions.append([129.131471, "TMT129N", False])
        self.known_ions.append([129.137790, "TMT129C", False])
        self.known_ions.append([130.134825, "TMT130N", False])
        self.known_ions.append([130.141145, "TMT130C", False])
        self.known_ions.append([131.138180, "TMT131N", False])
        self.known_ions.append([131.1445, "TMT131C", False])
        self.known_ions.append([132.141535, "TMT132N", False])
        self.known_ions.append([132.147855, "TMT132C", False])
        self.known_ions.append([133.14489, "TMT133N", False])
        self.known_ions.append([133.15121, "TMT133C", False])
        self.known_ions.append([134.148245, "TMT134N", False])
        self.known_ions.append([134.154565, "TMT134C", False])
        self.known_ions.append([135.151600, "TMT135N", False])
        self.known_ions.append([229.162932 + self.proton_mass, "TMT6Nterm", False])

        self.known_ions.sort()
        print("finished getting explanations")

    def determine_crude_correction(self):
        first_pass_peaks = {
            'IQ': 101.07094,
            'IE': 102.05495,
            'IM': 104.05285,
            'IH': 110.07127,
            'IR-NH3': 112.08692,
            'a-AA': 115.08659,
            'IR+H2O+H2O-N3H7': 116.0706,
            'IF': 120.08078,
            'b-AA-NH3': 126.05495,
            'IQ+CO': 129.06585,
            'IK+CO+H2O-NH3': 130.08625,
            'IC[Carbamidomethyl]': 133.04301,
            'IY': 136.07569,
            'a-AP': 141.10224,
            'b-AA': 143.0815,
            'a-AV': 143.11789,
            'b-GP': 155.0815,
            'IK+CO+H2ON2-NH3': 158.0924,
            'IW': 159.09167,
            'a-DP-H2O': 167.0815,
            'b-AP': 169.09715,
            'a-PV': 169.13354,
            'a-PT': 171.1128,
            'a-TV': 173.12845,
            'y-R': 175.11895,
            'a-LP': 183.14919,
            'b-PS': 185.09207,
            'b-PV': 197.12845
        }

        self.crude_xy_scatterplot = []
        first_pass_ppms = []
        for key in first_pass_peaks:
            peak_mz = first_pass_peaks[key]
            delta = self.initial_tolerance * first_pass_peaks[key] / 1e6
            lower_bound = peak_mz - delta
            upper_bound = peak_mz + delta
            best_match = [0, 0, 0]
            for i in range(len(self.by_count)):
                mz = i / 10000
                if mz > upper_bound:
                    break
                elif mz > lower_bound:
                    if self.by_count[i] > best_match[1]:
                        best_match = [mz, self.by_count[i], (peak_mz - mz) * 1e6 / peak_mz]
            
            if best_match[0] != 0:
                self.crude_xy_scatterplot.append(best_match)
                first_pass_ppms.append(best_match[2])

        self.crude_correction = np.median(first_pass_ppms)

    def find_initial_triggers(self):
        # scan through peaks to find all peaks with over 50
        previous_peak = [0, 0]
        counted = False
        for i in range(len(self.by_count)):
            if self.by_count[i] >= self.minimum_triggers:
                if self.by_count[i] < previous_peak[1] and not counted:
                    if self.by_count[i + 1] < self.by_count[i]:
                        self.triggered_peaks.append(previous_peak)
                        counted = True
                previous_peak = [i/10000, self.by_count[i]]
            else:
                if not counted and previous_peak[0] != 0 and self.by_count[i + 1] < self.by_count[i]:
                    self.triggered_peaks.append(previous_peak)
                counted = False
                previous_peak = [0, 0]
        
        print("finished finding initial triggers")

    def refine_triggered_peaks(self):
        for peak in self.triggered_peaks:
            peak_and_intensity = self.get_peak_fit_center(peak)
            if peak_and_intensity[0] == 0 and peak_and_intensity[1] == 0:
                continue
            else:
                peak_and_intensity = self.get_peak_fit_center(peak_and_intensity)
                if peak_and_intensity[0] == 0 and peak_and_intensity[1] == 0:
                    continue
                else:
                    if len(self.observed_peaks) > 0 and peak_and_intensity[0] == self.observed_peaks[len(self.observed_peaks) - 1][0]:
                        continue
                    else:
                        self.observed_peaks.append(peak_and_intensity)

        print("finished refining triggered peaks")
            
    def get_peak_fit_center(self, peak):
        mz_values = []
        intensity_values = []
        ppm_values = []
        mz_delta = int((self.ppm_delta / 1e6) * peak[0] * 10000)
        for index in range(mz_delta * 2 + 1):
            # change it to be 0
            add_index = int(peak[0] * 10000 + index - mz_delta - 1)
            mz_values.append(add_index / 10000 - peak[0])
            ppm_values.append((add_index / 10000 - peak[0]) * 1e6 / peak[0])
            intensity_values.append(self.by_strength[add_index])

        n = len(ppm_values)
        center = int(n/2)
        binsize = ppm_values[center]-ppm_values[center-1]

        try:
        #if 1:
            popt,pcov = curve_fit(gaussian_function,ppm_values,intensity_values,p0=[intensity_values[center],ppm_values[center],binsize])
        except:
            return [0, 0]

        peak_mz = round((peak[0] + (popt[1] * peak[0] / 1e6)), 5)
        return [peak_mz, round(popt[0], 1)]

    def identify_peaks(self):
        for index in range(len(self.known_ions)): # resets identifications
            if self.known_ions[index][2]:
                self.known_ions[index][2] = False

        if self.has_correction_spline:
            mz_values = []
            for peak in self.observed_peaks:
                mz_values.append(peak[0])
            x_values = np.array(mz_values)
            y_values = interpolate.BSpline(self.t, self.c, self.k)(x_values)
        
        for index in range(len(self.observed_peaks)):
            identifications = []
            peak = self.observed_peaks[index].copy()
            crude_correction_mz = self.crude_correction * peak[0] / 1e6
            peak[0] += crude_correction_mz

            if self.has_correction_spline:
                spline_correction_mz = y_values[index] * peak[0] / 1e6 # TODO: check if this should be before or after peak[0] += crude_correction_mz
                peak[0] += spline_correction_mz

            self.observed_peaks[index] = self.observed_peaks[index][0:1]
            for ion_index in range(len(self.known_ions)):
                amino_acid = self.known_ions[ion_index]
                # check tolerance based on ppm
                peak_tolerance = self.tolerance * peak[0] / 1e6
                if peak[0] > amino_acid[0] - peak_tolerance and peak[0] < amino_acid[0] + peak_tolerance and not amino_acid[2]:
                    identified_peak = [amino_acid[0], int(100000*(amino_acid[0] - peak[0]) + 0.5) / 100000, amino_acid[1]]
                    identifications.append(identified_peak)
                    self.known_ions[ion_index][2] = True
            
            identifications.sort(key = lambda x: abs(x[1]))
            for identification in identifications:
                self.observed_peaks[index].append(identification)

        print("finished identifying peaks")

    def plot_crude_calibrations(self):
        self.pdf = PdfPages(f'{self.peak_file}.calibration.pdf')
        self.fig, self.ax = plt.subplots(nrows=2, ncols=2, figsize=(8, 5))
        peak_mz = []
        calibration_ppm = []
        for closest_match in self.crude_xy_scatterplot:
            peak_mz.append(closest_match[0])
            calibration_ppm.append(closest_match[2])
        self.ax[0][0].plot(peak_mz, calibration_ppm, '.',c="g", markersize=2)
        self.ax[0][0].axhline(y = self.crude_correction, color = 'k', linestyle = 'dashed')
        self.ax[0][0].axhline(y = 0, color = 'k', linestyle = '-')
        self.ax[0][0].text(min(peak_mz), 0, f'tolerance: {self.initial_tolerance} ppm', fontsize='xx-small', ha='left', va='top')
        
        plt.tight_layout()

    def plot_peaks_strength(self):
        x = 0
        y = 0
        fig, ax = plt.subplots(self.rows, self.columns, figsize=(8.5,11))

        if self.has_correction_spline:
            mz_values = []
            for peak in self.observed_peaks:
                if len(peak) >= 3:
                    mz_values.append(peak[0])
            x_values = np.array(mz_values)
            y_values = interpolate.BSpline(self.t, self.c, self.k)(x_values)

        spline_index = 0
        for peak in self.observed_peaks:
            fig.subplots_adjust(top=0.98, bottom=0.05, left=0.12, right=0.98, hspace=0.2, wspace=0.38)
            # mz_values = []
            intensity_values = []
            ppm_values = []
            crude_correction_mz = self.crude_correction * peak[0] / 1e6
            mz_delta = int((self.ppm_delta / 1e6) * (peak[0] + crude_correction_mz) * 10000)

            for index in range(mz_delta * 2 + 1):
                # change it to be 0
                add_index = int((peak[0]) * 10000 + index - mz_delta - 1)
                # mz_values.append(add_index / 10000 - peak[0])
                ppm = (add_index / 10000 - peak[0]) * 1e6 / peak[0]
                intensity_values.append(self.by_strength[add_index])
                ppm_values.append(ppm)

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
                ax[x,y].plot(ppm_values,gaussian_function(ppm_values, *popt),'r:')
            except:
                # either skip it or show it, for now just skip it
                continue

            peak_fit_center = int(100000*(peak[0] + (popt[1] * peak[0] / 1e6))) / 100000

            if len(peak) >= 3:
                if self.has_correction_spline:
                    spline_correction_mz = y_values[spline_index] * peak[0] / 1e6
                    spline_index += 1
                else:
                    spline_correction_mz = 0
                line_x = (peak[2][0] - peak[0] - crude_correction_mz - spline_correction_mz) * 1e6 / peak[0]
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
                    self.pdf.savefig(fig)
                    plt.close('all')
                    fig, ax = plt.subplots(self.rows, self.columns,figsize=(8.5,11))
        
        if x != 0:
            self.pdf.savefig(fig) # maybe get rid of fig parameter in here
        self.pdf.close()
        plt.close()

        print("finished plotting peaks")

    def delta_scatterplots(self):
        for i in range(2):
            if i == 0: # the first one that i kept
                mz_values = []
                delta_values_ppm = []
                for peak in self.observed_peaks:
                    if len(peak) >= 3:
                        mz_values.append(peak[0])
                        delta_values_ppm.append((peak[2][1] * 1e6 / peak[2][0]) + self.crude_correction)
                    else:
                        continue
                self.ax[0,1].set_title(f"crude_correction: {self.crude_correction}", fontsize="xx-small")
                self.ax[0][1].scatter(mz_values, delta_values_ppm, 0.5)
                self.ax[0][1].axhline(y = self.crude_correction, color = 'k', linestyle = 'dashed')
                self.ax[0][1].text(min(mz_values), 0, f'tolerance: {self.tolerance} ppm', fontsize='xx-small', ha='left', va='top')
                self.ax[0][1].axhline(y = 0, color = 'k', linestyle = '-')
            elif i == 1:
                mz_values = []
                delta_values_ppm = []
                artificial_peak_range = []
                for peak in self.observed_peaks:
                    crude_correction_mz = self.crude_correction * peak[0] / 1e6
                    if len(peak) >= 3:
                        mz_values.append(peak[0] + crude_correction_mz)
                        delta_values_ppm.append((peak[2][1] * 1e6 / peak[2][0]))
                        if peak[0] >= 300 and peak[0] <= 400:
                            artificial_peak_range.append(peak[2][1] * 1e6 / peak[2][0])  
                    else:
                        continue
                
                # add an artificial peak at 500 with the median between 300 and 400 so the spline fit will flatten out
                index = mz_values[-1]
                while index < 500:
                    index += 10
                    mz_values.append(index)
                    delta_values_ppm.append(np.median(artificial_peak_range))
                self.ax[1][0].scatter(mz_values, delta_values_ppm, 0.5)
                self.ax[1][0].axhline(y = 0, color = 'k', linestyle = '-')
                
                mz_values.sort()

                if len(mz_values) <= 5:
                    print("There are too little data points to create a spline fit.")
                    self.ax[1,0].set_title("There are too little data points to create a spline fit", fontsize="xx-small")
                else:
                    x_new = np.linspace(0, 1, 5)[1:-1] # 5 = 3 knots + 2
                    q_knots = np.quantile(mz_values, x_new)
                    self.t, self.c, self.k = interpolate.splrep(mz_values, delta_values_ppm, t=q_knots, s=3)
                    self.has_correction_spline = True
                    x_regular = np.arange(mz_values[0], mz_values[-1])
                    y_values = interpolate.BSpline(self.t, self.c, self.k)(x_regular)
                    self.ax[1][0].plot(x_regular, y_values, 'b')
                    self.ax[1,0].set_title(f"t: {self.t}, c: {self.c}, k: {self.k}", fontsize="xx-small")
        print("finished creating delta scatterplot")
        plt.close()

    def plot_corrected_scatterplot(self):
        mz_values = []
        delta_values_ppm = []

        if self.has_correction_spline:
            mz_values = []
            for peak in self.observed_peaks:
                if len(peak) >= 3:
                    mz_values.append(peak[0])
            x_values = np.array(mz_values)
            y_values = interpolate.BSpline(self.t, self.c, self.k)(x_values)

        mz_values = []
        spline_index = 0
        for peak in self.observed_peaks:
            crude_correction_mz = self.crude_correction * peak[0] / 1e6
            
            if len(peak) >= 3:
                if self.has_correction_spline:
                    spline_correction_mz = y_values[spline_index]
                    spline_index += 1
                else:
                    spline_correction_mz = 0
                mz_values.append(peak[0] + crude_correction_mz + spline_correction_mz)
                delta_values_ppm.append((peak[2][1] * 1e6 / peak[2][0]))
            else:
                continue
        self.ax[1][1].scatter(mz_values, delta_values_ppm, 0.5)
        self.ax[1][1].axhline(y = 0, color = 'k', linestyle = '-')
        self.pdf.savefig(self.fig)

    def write_output(self):
        with open(f'{self.peak_file}.tsv', 'w') as file:
        # with open('common_peaks.tsv', 'w') as file:
            writer = csv.writer(file, delimiter='\t', lineterminator='\n')
            writer.writerows(self.observed_peaks)

        print("finished writing tsv file")

    def show_stats(self):
        t1 = timeit.default_timer()
        print(f"INFO: Read {self.stats['counter']} spectra from {self.file_name}")
        print(f"The number of ms1spectra is {self.stats['ms1spectra']}")
        print(f"The number of ms2spectra is {self.stats['ms2spectra']}")

        print(f"INFO: Elapsed time: {t1-self.t0}")
        print(f"INFO: Processed {self.stats['counter']/(t1-self.t0)} spectra per second")

    # Gaussian function used for curve fitting
def gaussian_function(x,a,x0,sigma):
    # a = amplitude
    # x0 = center point of the gaussian
    # sigma = measure of the width of the gaussian
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
    peak_finder = MSRunPeakFinder()
    # find the intensity of each m/z value
    peak_finder.aggregate_spectra()
    peak_finder.determine_crude_correction()
    peak_finder.plot_crude_calibrations()
    # create an array of all identifiable theoretical and their masses
    peak_finder.get_theoretical_ions() #
    # identify all the spectra with a certain intensity
    peak_finder.find_initial_triggers() #
    peak_finder.refine_triggered_peaks() # 
    peak_finder.identify_peaks() #
    # save the identified peaks to a file
    peak_finder.delta_scatterplots() #
    peak_finder.identify_peaks()
    peak_finder.plot_corrected_scatterplot()
    peak_finder.plot_peaks_strength()
    # print out data, including run time, number of peaks found, etc
    peak_finder.write_output() #
    peak_finder.show_stats()

if __name__ == "__main__": main()