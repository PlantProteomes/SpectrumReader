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
        argparser = argparse.ArgumentParser(description='An example program that creates heatmaps using tsv files')
        argparser.add_argument('--tsv_file', action='store', help='Name of the tsv file to read')
        
        params = argparser.parse_args()
        self.tsv_file = params.tsv_file
        # self.heat_map_matrix = numpy.zeros((20,25))
        self.residues = {'no loss': {}} # make a nested dictionary
        # 
        self.scanned_residues = []
        self.amino_acids = ['IA', 'IC', 'IC[Carbamidomethyl]', 'ID' , 'IE', 'IF', 'IG', 'IH', 'II', 'IK', 'IL', 'IM', 'IM[Oxidation]', 'IN', 'IP', 'IQ', 'IR', 'IS', 'IT', 'IV', 'IW', 'IY']

    def read_file(self):
        with open(self.tsv_file) as file:
            lines = [line.rstrip() for line in file]
            lines = lines[1:] # removes the header line
            for peak in lines:
                line_split = [i for i in peak.split("\t")]
                if line_split[6] != '?':
                    identification = line_split[6]
                    if identification[0] == 'I' and "+i" not in identification:
                        if "[Acetyl]" in identification:
                            continue
                        elif "[Oxidation]" in identification:
                            amino_acid = identification[:13]
                            residue = identification[13:]
                        elif "[Carbamidomethyl]" in identification:
                            amino_acid = identification[:19]
                            residue = identification[19:]
                        else:
                            amino_acid = identification[:2]
                            residue = identification[2:]
                        if not residue in self.scanned_residues:
                            self.residues[residue] = {}
                            self.scanned_residues.append(residue)
                        self.residues[residue][amino_acid] = float(line_split[1]) # adds the intensity
                if len(line_split) >= 8 and line_split[7] != '?':
                    identifications = line_split[7].split(", ")
                    for identification in identifications:
                        if len(identification) > 0 and identification[0] == 'I' and "+i" not in identification:
                            if "[Acetyl]" in identification:
                                continue
                            elif "[Oxidation]" in identification:
                                amino_acid = identification[:13]
                                residue = identification[13:]
                            elif "[Carbamidomethyl]" in identification:
                                amino_acid = identification[:19]
                                residue = identification[19:]
                            else:
                                amino_acid = identification[:2]
                                residue = identification[2:]
                            if not residue in self.scanned_residues:
                                self.residues[residue] = {}
                                self.scanned_residues.append(residue)
                            self.residues[residue][amino_acid] = float(line_split[1]) # adds the intensity

    def to_strength(self, intensity):
        if intensity == 0:
            return 0
        elif intensity <= 2.5:
            return 0.2
        elif intensity <= 20:
            return 0.5
        else:
            return 0.8

    def generate_heatmap(self):
        self.residues["no loss"] = self.residues.pop('')
        self.heat_map_matrix = numpy.zeros((len(self.amino_acids),len(self.residues)))

        residue_list = list(self.residues.items())
        residue_list.sort(key = lambda x: -1 * sum(x[1].values()))
        residue_list.sort(key = lambda x: -1 * len(x[1]))
        self.residues = dict(residue_list)

        residue_index = 0
        for residue in self.residues:
            for amino_acid in self.residues[residue]:
                amino_acid_index = self.amino_acids.index(amino_acid)
                self.heat_map_matrix[amino_acid_index][residue_index] = self.to_strength(self.residues[residue][amino_acid])
            residue_index += 1

        # plot heatmap
        fig, ax = plt.subplots(nrows=1, ncols=2, gridspec_kw={'width_ratios': [1, 10]})
        plt.xticks(range(len(self.residues)), self.residues, rotation = 90, fontsize = 5)
        plt.yticks(range(len(self.amino_acids)), self.amino_acids, fontsize = 5)
        ax[1].imshow(self.heat_map_matrix, cmap = 'BuPu')
        # im_ratio = self.heat_map_matrix.shape[0]/self.heat_map_matrix.shape[1]
        # plt.colorbar(fraction=0.047*im_ratio)

        colorBar = numpy.array([[0, 2], [15, 25]])
        ax[0].imshow(colorBar, cmap = 'BuPu')
        ax[0].set_title("Color Scale", fontsize = 7)
        ax[0].text(0, 0, "0",
                    ha="center", va="center", color="r", fontsize = 5)
        ax[0].text(0, 1, "<2.5",
                    ha="center", va="center", color="r", fontsize = 5)
        ax[0].text(1, 0, "<=20",
                    ha="center", va="center", color="r", fontsize = 5)
        ax[0].text(1, 1, ">20",
                    ha="center", va="center", color="r", fontsize = 5)
        plt.setp(ax[0].get_xticklabels(), visible=False)
        plt.setp(ax[0].get_yticklabels(), visible=False)
        plt.xticks()
        ax[0].tick_params(
            axis='both',
            length=0,
            width=0,
            labelbottom=False)
        
        # save the figure
        plt.tight_layout()
        plt.savefig("heatmap.pdf")
        plt.close

def main():
    heat_map = HeatMapCreator()
    heat_map.read_file()
    heat_map.generate_heatmap()

if __name__ == "__main__": main()