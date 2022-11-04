#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python scan_ms2s.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML --rows 3 --columns 5

# close matches, but not close enough?
# W-NH3 - 142.0651
# line 59 in output - 142.1224, 67 appearances

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

def main():

    # sets up a module with a description
    argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
    argparser.add_argument('--mzml_file', action='store', help='Name of the mzML file to read')
    argparser.add_argument('--rows', action='store', type=int, help='Number of rows for graphs')
    argparser.add_argument('--columns', action='store', type=int, help='Number of columns for output')

    # runs the parse_args method to store all the results into the params variable
    params = argparser.parse_args()
    
    if params.mzml_file is None or params.mzml_file == "":
        print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
        return

    # tries to open the params file specified
    if not os.path.isfile(params.mzml_file):
        print(f"ERROR: File '{params.mzml_file}' not found or not a file")
        return

    if params.rows is None or params.rows == "":
        params.rows = 5
    
    if params.columns is None or params.columns == "":
        params.columns = 3

    peak_correction_factor = 0.0006

    # create three arrays of 4 million for bins
    by_count = np.zeros((4000000,), dtype=int)
    # by_intensity = np.zeros((4000000,), dtype=int)
    # by_strength = np.zeros((4000000,), dtype=int)

    check_spectra = True
    if check_spectra:
        all_peaks = []
        #### Read spectra from the file and isolate specific spectra
        t0 = timeit.default_timer()
        stats = { 'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0 }
        
        with mzml.read(params.mzml_file) as reader:
            for spectrum in reader:

                #### Update counters and print progress
                stats['counter'] += 1

                # start with a smaller sample size
                # if stats['counter'] == 101:
                    # break

                # add most popular mz values to be printed out
                satisfy_requirements = True
                if spectrum['ms level'] == 1:
                    stats['ms1spectra'] += 1
                elif spectrum['ms level'] == 2:
                    stats['ms2spectra'] += 1
                    smallest_peak_intensity = sys.maxsize
                    for index in range(len(spectrum['m/z array'])):
                        peak = spectrum['m/z array'][index]
                        peak -= peak_correction_factor
                        if peak > 400:
                            break
                        else:
                            intensity = spectrum['intensity array'][index]
                            all_peaks.append(intensity)
                            by_count[int(10000 * peak + 0.5)] += 1
                            # by_intensity[int(10000 * peak + 0.5)] += intensity
                            # smallest_peak_intensity = min(smallest_peak_intensity, intensity)

                    # for index in range(len(spectrum['m/z array'])):
                        # peak = spectrum['m/z array'][index]
                        # if peak > 400:
                            # break
                        # else:
                            # intensity = spectrum['intensity array'][index]
                            # by_strength[int(10000 * peak + 0.5)] += get_strength(intensity, smallest_peak_intensity)

                if stats['counter']/1000 == int(stats['counter']/1000):
                    print(f"  {stats['counter']}")

        #### Print final timing information
        t1 = timeit.default_timer()
        print(f"INFO: Read {stats['counter']} spectra from {params.mzml_file}")
        print(f"The number of ms1spectra is {stats['ms1spectra']}")
        print(f"The number of ms2spectra is {stats['ms2spectra']}")

    # create a histogram of all the intensities between 0 and 20000
    # plt.hist(all_peaks, bins=200, range=[0, 20000])
    # plt.show()

    # creates a dictionary of pyteomics mass for 60 amino acids, then sorts it
    # make this more complex to include the table in the second issue
    # while I loop through the amino acids, if the ion type is a, loop over the mass modifications of the issue, try to compute
    # the mass deltas based off the molecular compositions and add the new mass to the list
    # syntax to use:
    # for histodine: IH (instead of H type a)
    # for histodine b ions: IH+CO (instead of H type b)
    # for histodine y ions: IH+CO+H2O (instead of H type y)
    # for first pass, ignore duplicates to see it works properly (gets it through b & y roots, but also the special table)
    # in the future, can remove duplicates
    ion_types = ['a', 'b', 'y']
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
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
    isotopic_deltas = {}
    amino_acid_mass = []

    for isotope in aa_immonium_losses:
        delta_masses = {}
        for variant in aa_immonium_losses[isotope]:
            delta_mass = 0
            variant_add = variant.split('-')[0].split('+')[1:]
            variant_subtract = variant.split('-')[1:]
            for variant_formula in variant_add:
                delta_mass += mass.calculate_mass(formula=variant_formula)
            for variant_formula in variant_subtract:
                delta_mass -= mass.calculate_mass(formula=variant_formula)
            delta_masses[variant] = delta_mass
        isotopic_deltas[isotope] = delta_masses
    
    for acid in amino_acids:
        for ion in ion_types:
            acid_mass = mass.calculate_mass(sequence=acid, ion_type=ion, charge=1)
            acid_mass = int(acid_mass * 10000 + 0.5) / 10000.0
            amino_acid_mass.append([acid_mass, acid, ion])

            if ion == 'a':
                for isotope in isotopic_deltas[acid]:
                    print(acid + isotope)
                    acid_mass = mass.calculate_mass(sequence=acid, ion_type=ion, charge=1)
                    acid_mass += isotopic_deltas[acid][isotope]
                    acid_mass = int(acid_mass * 10000 + 0.5) / 10000.0
                    print(acid_mass)
                    amino_acid_mass.append([acid_mass, acid, ion + isotope])

    sorted(amino_acid_mass)
    # print(amino_acid_mass)

    # for key in amino_acid_mass:
        # for peak in tallest_peaks:
            # if peak[0] > key - 0.002 and peak[0] < key + 0.002:
                # print(str(peak[0]) + '\t' + str(peak[1]) + '\t' + str(peak[0] - key) + '\t' + str(amino_acid_mass[key]))
                # break

    # print out top mz values and number of appearances in a file
    if check_spectra:
        with open('popular_spectra.tsv', 'w') as file:
            writer = csv.writer(file, delimiter='\t', lineterminator='\n')
            tallest_peaks = []
            previous_peak = [0, 0]
            printed = False
            for i in range(len(by_count)):
                if by_count[i] >= 50:
                    if by_count[i] < previous_peak[1] and not printed:
                        for amino_acid in amino_acid_mass:
                            if i/10000 > amino_acid[0] - 0.002 and i/10000 < amino_acid[0] + 0.002:
                                if len(previous_peak) == 5 and previous_peak[4] != amino_acid[1:]:
                                    printed = False
                                    previous_peak = [i/10000, by_count[i]]
                                previous_peak.append(amino_acid[0])
                                previous_peak.append(int(10000 * (i/10000 - amino_acid[0])) / 10000)
                                previous_peak.append(amino_acid[1:])
                                if not printed:
                                    tallest_peaks.append(previous_peak)
                                    writer.writerow(previous_peak)
                                    printed = True
                        if not printed:
                            tallest_peaks.append(previous_peak)
                            writer.writerow(previous_peak)
                            printed = True
                    previous_peak = [i/10000, by_count[i]]
                else:
                    if not printed and previous_peak[0] != 0:
                        tallest_peaks.append(previous_peak)
                        writer.writerow(previous_peak)
                    printed = False
                    previous_peak = [0, 0]

        # extra step plotting for all tallest peaks
        want_plot = False
        if want_plot:
            for peak in tallest_peaks:
                mz_values = []
                intensity_values = []
                for index in range(31):
                    add_index = int(peak[0] * 10000 + index - 16)
                    mz_values.append(add_index / 10000)
                    intensity_values.append(by_count[add_index])
                plt.step(mz_values, intensity_values, where='mid')
                plt.show()

        save_plot = True
        if save_plot:
            x = 0
            y = 0
            pdf = PdfPages('common_peaks.pdf')
            fig, ax = plt.subplots(params.rows,params.columns,figsize=(8.5,11))
            for peak in tallest_peaks:
                fig.subplots_adjust(top=0.98, bottom=0.05, left=0.12, right=0.98, hspace=0.2, wspace=0.38)
                mz_values = []
                intensity_values = []
                for index in range(31):
                    # change it to be 0
                    add_index = int(peak[0] * 10000 + index - 16)
                    mz_values.append(add_index / 10000 - peak[0])
                    intensity_values.append(by_count[add_index])
                # change the axis labels to be a smaller text
                ax[x,y].step(mz_values, intensity_values, where='mid')
                ax[x,y].tick_params(axis='x', labelsize='xx-small')
                ax[x,y].locator_params(axis='x', nbins=5)
                # ax[x,y].set_xticks([-0.0015, -0.0007, 0, 0.0007, 0.0015])
                ax[x,y].tick_params(axis='y', labelsize='xx-small')
                if len(peak) == 5:
                    ax[x,y].axvline(x=peak[3], color='black', lw=1, linestyle='--')
                    ax[x,y].text(-0.0017, 45*max(intensity_values)/64, f'peak center: \n{peak[0]}\nion m/z: \n{peak[2]}\nion id: \n{peak[4][0]} type {peak[4][1]}', fontsize='xx-small')
                    # ax[x,y].set_title(f"Peak {peak[0]}, Amino Acid {peak[4][0]} type {peak[4][1]}", fontsize="xx-small")
                else:
                    ax[x,y].text(-0.0017, 15*max(intensity_values)/16, f'peak center: \n{peak[0]}', fontsize='xx-small')
                    # ax[x,y].set_title(f"Peak {peak[0]}", fontsize="xx-small")

                # creates a new figure if full
                y += 1
                y %= params.columns
                if y == 0:
                    x += 1
                    x %= params.rows
                    if x == 0:
                        pdf.savefig(fig)
                        fig, ax = plt.subplots(params.rows,params.columns,figsize=(8.5,11))
            
            pdf.savefig(fig)
            pdf.close()

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
