#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python make_heatmap.py --tsv_file ..\data\combined_list.tsv

import matplotlib.pyplot as plt
import argparse
import numpy

from ast import literal_eval

class HeatMapCreator:

    def __init__(self):
        # Takes in all the tsv file names in a string format separated by commas
        argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
        argparser.add_argument('--tsv_file', action='store', help='Name of the mzML file to read')
        
        params = argparser.parse_args()
        self.tsv_file = params.tsv_file
        self.heat_map_matrix = numpy.zeros((20,25))
        self.residues = ['+CO', '+H2O', '+H2ON2', '+NH3', '-CHN', '-CH2', '-CH2N', '-CH2N2', '-CH3N', '-C2H2', '-C2H3N3', '-C2H4', '-C2H4N2', '-C3H6', '-C3H6N2', '-C4H6N2', '-C4H7N', '-CH5N3', '-CH6N2', '-CO', '-H2O', '-NH', '-NH2', '-NH3', '-N3H7']
        self.amino_acids = ['A', 'C', 'D' , 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    def read_file(self):
        with open(self.tsv_file) as file:
            lines = [line.rstrip() for line in file]
            lines = lines[1:]
            for peak in lines:
                line_split = [i for i in peak.split("\t")]
                if line_split[6] != '?':
                    identification = line_split[6]
                    if identification[0] == 'I':
                        amino_acid_index = self.amino_acids.index(identification[1])
                        for residue_index in range(len(self.residues)):
                            if self.residues[residue_index] in identification:
                                self.heat_map_matrix[amino_acid_index][residue_index] = 1
                if line_split[7] != '?':
                    identifications = line_split[7]
                    identifications = literal_eval(identifications)
                    for identification in identifications:
                        if identification[0] == 'I':
                            amino_acid_index = self.amino_acids.index(identification[1])
                            for residue_index in range(len(self.residues)):
                                if self.residues[residue_index] in identification:
                                    self.heat_map_matrix[amino_acid_index][residue_index] = 1

    def generate_heatmap(self):
        # plot heatmap
        plt.imshow(self.heat_map_matrix, cmap = 'autumn')
        plt.colorbar()

        # Set tick labels
        plt.xticks(range(len(self.residues)), self.residues, rotation = 90)
        plt.yticks(range(len(self.amino_acids)), self.amino_acids)
        
        # save the figure
        plt.savefig("heatmap.pdf")
        plt.close

def main():
    heat_map = HeatMapCreator()
    heat_map.read_file()
    heat_map.generate_heatmap()

if __name__ == "__main__": main()