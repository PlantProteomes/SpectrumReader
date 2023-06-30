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
        # self.heat_map_matrix = numpy.zeros((20,25))
        self.residues = {} # make a nested dictionary
        # 
        self.scanned_residues = []
        self.amino_acids = ['A', 'C', 'C[Carbamidomethyl]', 'D' , 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'M[Oxidation]', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    def read_file(self):
        with open(self.tsv_file) as file:
            lines = [line.rstrip() for line in file]
            lines = lines[1:] # removes the header line
            for peak in lines:
                line_split = [i for i in peak.split("\t")]
                if line_split[6] != '?':
                    identification = line_split[6]
                    if identification[0] == 'I' and "+i" not in identification:
                        residue = identification[2:]
                        if not residue in self.scanned_residues:
                            self.residues[residue] = {}
                            self.scanned_residues.append(residue)
                        self.residues[residue][identification[1]] = float(line_split[1]) # adds the intensity
                if line_split[7] != '?':
                    identifications = line_split[7]
                    identifications = literal_eval(identifications)
                    for identification in identifications:
                        if identification[0] == 'I' and "+i" not in identification:
                            residue = identification[2:]
                            if not residue in self.scanned_residues:
                                self.residues[residue] = {}
                                self.scanned_residues.append(residue)
                            self.residues[residue][identification[1]] = float(line_split[1]) # adds the intensity

    def generate_heatmap(self):
        self.heat_map_matrix = numpy.zeros((20,len(self.residues)))

        residue_index = 0
        for residue in self.residues:
            for amino_acid in self.residues[residue]:
                amino_acid_index = self.amino_acids.index(amino_acid)
                self.heat_map_matrix[amino_acid_index][residue_index] = self.residues[residue][amino_acid]
            residue_index += 1

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