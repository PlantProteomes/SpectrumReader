
#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python find_nist_identification.py --file_name ..\data\IARPA3_best_tissue_add_info.msp --unknown_ions_file_name ..\data\HFX_9850_GVA_DLD1_2_180719_subset_unidentified.tsv

import argparse
import os
import re
import csv

class IdentificationFinder:
    def __init__(self):
        # Takes in all the tsv file names in a string format separated by commas
        argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
        argparser.add_argument('--file_name', type=str, help='Filenames of one or more files to read')
        argparser.add_argument('--unknown_ions_file_name', type=str, help='Filenames of tsv files containing unknown ions to read')
        
        params = argparser.parse_args()
        self.filename = params.file_name
        self.unknown_ions_file_name = params.unknown_ions_file_name

        self.unknown_mzs = []
        self.possible_identifications = []

    def get_unidentified(self):
        if not os.path.isfile(self.unknown_ions_file_name):
            print(f"ERROR: File '{self.unknown_ions_file_name}' not founds or not a file")
            return
        if self.unknown_ions_file_name[-3:] != "tsv":
            print(f"ERROR: File '{self.unknown_ions_file_name}' is not a tsv file")
            return
        
        print("starting to find all unidentified/chemical formula identifications from provided list")
        with open(self.unknown_ions_file_name) as file:
            lines = [line.rstrip() for line in file]
            lines = lines[1:] # Removes header from tsv file
            for unknown_ions in lines:
                line_split = [i for i in unknown_ions.split("\t")]
                self.unknown_mzs.append(float(line_split[0]))
    
    def read_file(self):
        print("starting to search for possible identifications")
        unknown_mz_index = 0
        with open(self.filename) as infile:
            for line in infile:
                match = re.match(r'\d', line)
                unknown_mz = self.unknown_mzs[unknown_mz_index]
                mz = 0
                if match:
                    columns = line.strip().split("\t") # use RegEx to check if the line starts with a number
                    mz = float(columns[0])
                    while unknown_mz + 0.001 < mz and unknown_mz_index != len(self.unknown_mzs) - 1:
                        unknown_mz_index = min(unknown_mz_index + 1, len(self.unknown_mzs) - 1)
                        unknown_mz = self.unknown_mzs[unknown_mz_index]
                    # if mz < unknown_mz + 0.001 and mz > unknown_mz - 0.001: #is close enough to what i'm looking for, print the line
                    if mz < unknown_mz + 0.001 and mz > unknown_mz - 0.001 and columns[2][1] != "?":
                        possible_identification = ""
                        for i in line:
                            possible_identification += i
                        self.possible_identifications.append([possible_identification[0:-2]])
                        # print(line, end = '')
                else:
                    unknown_mz_index = 0

    def write_identifications(self):
        print("starting to write possible identifications")
        self.possible_identifications.sort()
        output_file = self.unknown_ions_file_name[0:-17]
        with open(f'{output_file}_possible_identifications.tsv', 'w') as file:
            writer = csv.writer(file, delimiter='\t', lineterminator='\n')
            writer.writerows(self.possible_identifications)

def main():
    find_identifications = IdentificationFinder()
    # Populates the aggregated_observed_peaks dictionary
    find_identifications.get_unidentified()
    find_identifications.read_file()
    find_identifications.write_identifications()

if __name__ == "__main__": main()