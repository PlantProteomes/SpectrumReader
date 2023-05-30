#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python determine_isotope.py --chemical_formula C9H9NO

import argparse
import csv
import matplotlib.pyplot as plt

from pyteomics import mzml, auxiliary, mass

class isotopeFinder:

    def __init__(self, chemical_formula):
        self.isotopic_elements = {
            "C": [[12, 0.9884], [13.003354835336, 0.0096]],
            "H": [[1.007825031898, 0.99972], [2.014101777844, 0.00001]],
            "O": [[15.994914619257, 0.99738], [16.999131755953, 0.000367], [17.999159612136, 0.00187]],
            "N": [[14.003074004251, 0.99578], [15.000108898266, 0.00337]]
        }
        self.chemical_formula = chemical_formula
        self.proton_mass = 1.007276

    def find_isotopes(self):
        self.monoisotope_mass = round(mass.calculate_mass(formula=self.chemical_formula) - mass.calculate_mass(formula="CO") + self.proton_mass, 5)
        self.isotope_masses = {}
        for isotope in self.isotopic_elements:
            self.isotope_masses[f"{isotope}{int(round(self.isotopic_elements[isotope][1][0], 0))}"] = round(self.isotopic_elements[isotope][1][0] - self.isotopic_elements[isotope][0][0], 5)

    def create_table(self):
        self.table = [[self.monoisotope_mass, 0, "monoisotope", 1.0]]
        for isotope in self.isotope_masses:
            isotope_element = isotope[0]
            if isotope_element == "C" or isotope_element == "O":
                lost_elements = 1
            else:
                lost_elements = 0
            isotope_proportion = self.isotopic_elements[isotope_element][1][1]
            index = self.chemical_formula.index(isotope_element)
            if index == len(self.chemical_formula) - 1:
                num_of_element = 1 - lost_elements
            else:
                if self.chemical_formula[index + 1].isnumeric():
                    num_of_element = int(self.chemical_formula[index + 1]) - lost_elements
                else:
                    num_of_element = 1 - lost_elements
            if num_of_element == 0:
                continue
            else:
                self.table.append([round(self.monoisotope_mass + self.isotope_masses[isotope], 5), self.isotope_masses[isotope], isotope, isotope_proportion * num_of_element])

    def write_table(self):
        table = self.table.copy()
        table.sort(key=lambda x:x[0])
        table.insert(0, ["m/z", "delta from monoisotope", "isotope", "intensity"])
        with open('isotopes.tsv', 'w') as file:
            writer = csv.writer(file, delimiter='\t', lineterminator='\n')
            writer.writerows(table)

    def print_graph(self):
        isotopes = self.table.copy()
        isotopes.pop(0)
        ppm_x = self.to_ppm([i[0] for i in isotopes])
        # mz_x = [i[0] for i in isotopes]
        intensity_y = [j[3] for j in isotopes]
        plt.bar(ppm_x, intensity_y, width = 0.5)
        # plt.axvline(x=ppm_x[0]+1, color='black', lw=1, linestyle='--')
        plt.savefig("isotopes.pdf")
        plt.close

    def to_ppm(self, mz_array):
        ppm_array = []
        reference_mz = mz_array[0]
        for mz in mz_array:
            ppm = (mz - reference_mz) / reference_mz * 1e6
            ppm_array.append(ppm)

        return ppm_array

def main():
    argparser = argparse.ArgumentParser(description='An example program that calculate mass of isotopes of a compound')
    argparser.add_argument('--chemical_formula', action='store', default=5, help='Number of rows for graphs')

    params = argparser.parse_args()
    isotope_finder = isotopeFinder(params.chemical_formula)
    isotope_finder.find_isotopes()
    isotope_finder.create_table()
    isotope_finder.write_table()
    isotope_finder.print_graph()

if __name__ == "__main__": main()