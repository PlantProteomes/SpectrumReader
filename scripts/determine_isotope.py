#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal):

import argparse

class isotopeFinder:

    def __init__(self, chemical_formula):
        self.isotopic_elements = {
            "C": [[12, 0.9884], [13.003354835336, 0.0096]],
            "H": [[1.007825031898, 0.99972], [2.014101777844, 0.00001]],
            "O": [[15.994914619257, 0.99738], [16.999131755953, 0.000367], [17.999159612136, 0.00187]],
            "N": [[14.003074004251, 0.99578], [15.000108898266, 0.00337]]
        }
        self.chemical_formula = chemical_formula

def main():
    argparser = argparse.ArgumentParser(description='An example program that calculate mass of isotopes of a compound')
    argparser.add_argument('--chemical_formula', action='store', default=5, help='Number of rows for graphs')

    params = argparser.parse_args()
    isotope_finder = isotopeFinder(params.chemical_formula)
    isotope_finder.find_

if __name__ == "__main__": main()