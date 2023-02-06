#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python spline_fit_test.py --tsv_file ..\data\HFX_9850_GVA_DLD1_2_180719_subset_delta_scatterplot_values.tsv
# to run the program (use terminal): python spline_fit_test.py --tsv_file ..\data\HFX_9850_GVA_DLD1_2_180719_delta_scatterplot_values.tsv
# to run the program (use terminal): python spline_fit_test.py --tsv_file ..\data\QEX03_210305_PTMscan_wt_IP_63_delta_scatterplot_values.tsv
# to run the program (use terminal): python spline_fit_test.py --tsv_file ..\data\06CPTAC_BCprospective_W_BI_20161116_BA_f17_delta_scatterplot_values.tsv
# to run the program (use terminal): python spline_fit_test.py --tsv_file ..\data\Q20181210_06_delta_scatterplot_values.tsv

import os
import argparse
import os.path
import matplotlib.pyplot as plt
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import operator

from scipy import interpolate
from collections import OrderedDict
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from numpy import exp

class SplineFitter:

    def __init__(self):

        argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
        argparser.add_argument('--tsv_file', action='store', help='Name of the mzML file to read')

        params = argparser.parse_args()

        # variables declaration
        self.mz_values = []
        self.delta_ppm_values = []

        # verify valid file to read
        if params.tsv_file is None or params.tsv_file == "":
            print('ERROR: Parameter --tsv_file must be provided. See --help for more information')
            return

        # tries to open the params file specified
        if not os.path.isfile(params.tsv_file):
            print(f"ERROR: File '{params.tsv_file}' not found or not a file")
            return

        self.file_name = params.tsv_file
        self.infile = open(self.file_name, 'rb')
        self.peak_file = self.file_name[0:len(self.file_name) - 29]

    def update_arrays(self):
        with open(self.file_name) as file:
            lines = [line.rstrip() for line in file]
            for pairs in lines:
                line_split = [i for i in pairs.split("\t")]
                self.mz_values.append(float(line_split[0]))
                self.delta_ppm_values.append(float(line_split[1]))

    def fit_spline(self):
        knots = [2, 3, 4, 5]
        i = 0
        pdf = PdfPages(f'{self.peak_file}_spline_fits.pdf')

        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8, 5))
        
        for row in range(2):
            for col in range(2):
                self.refine_scatterplot()
                ax[row][col].plot(self.refined_mz_values, self.refined_delta_ppm_values, '.',c="g", markersize=2)
                x = self.refined_mz_values
                x_new = np.linspace(0, 1, knots[i]+2)[1:-1]
                q_knots = np.quantile(x, x_new)
                t, c, k = interpolate.splrep(x, self.refined_delta_ppm_values, t=q_knots, s=3)
                # instead, there should be a yfit value for each x value in self.mz_values which is then graphed
                yfit = interpolate.BSpline(t, c, k)(x)
                ax[row][col].plot(x, yfit, 'r')
                x_regular = np.arange(x[0], x[-1])
                yfit = interpolate.BSpline(t, c, k)(x_regular)
                ax[row][col].plot(x_regular, yfit, 'b')
                ax[row][col].set_title("Knots number - "+str(knots[i]))
                ax[row][col].grid()
                i=i+1
                
        plt.tight_layout()

        pdf.savefig(fig)
        pdf.close()

    def refine_scatterplot(self):
        # xy refers to the x and y axis on the scatterplot, x being mz values and y being delta ppm
        # values
        bin_size = 25 # (self.mz_values[len(self.mz_values) - 1] - self.mz_values[0])/12
        refined_xy = {}
        xy = {}
        removable_index = []
        upper_bound = self.mz_values[0] + bin_size
        for index in range(len(self.mz_values)):
            xy[self.mz_values[index]] = self.delta_ppm_values[index]
            if self.mz_values[index] > upper_bound:
                sorted_xy = sorted(refined_xy.items(), key=operator.itemgetter(1))
                length = len(sorted_xy)
                i = 1 / length
                if length >= 10:
                    for tuple in sorted_xy:
                        if i <= 0.1 or i >= 0.9:
                            removable_index.append(tuple[0])
                        i += 1 / length
                elif length >= 5:
                    for tuple in sorted_xy:
                        if i <= 0.2 or i >= 0.8:
                            removable_index.append(tuple[0])
                        i += 1 / length
                else:
                    for tuple in sorted_xy:
                        removable_index.append(tuple[0])
                upper_bound += bin_size
                refined_xy = {}
            refined_xy[self.mz_values[index]] = self.delta_ppm_values[index]
        
        for key in removable_index:
            xy.pop(key)

        self.refined_mz_values = []
        self.refined_delta_ppm_values = []
        for key in xy:
            self.refined_mz_values.append(int(key))
            self.refined_delta_ppm_values.append(xy[key])

def main():
    splinefitter = SplineFitter()
    splinefitter.update_arrays()
    splinefitter.fit_spline()

if __name__ == "__main__": main()
