#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python find_unidentified.py ..\data\HFX_9850_GVA_DLD1_2_180719_subset.tsv
# to run the program (use terminal): python find_unidentified.py ..\data\HFX_9850_GVA_DLD1_2_180719.tsv
# to run the program (use terminal): python find_unidentified.py ..\data\QEX03_210305_PTMscan_wt_IP_63.tsv
# to run the program (use terminal): python find_unidentified.py ..\data\06CPTAC_BCprospective_W_BI_20161116_BA_f17.tsv
# to run the program (use terminal): python find_unidentified.py ..\data\Q20181210_06.tsv
# to run the program (use terminal): python find_unidentified.py ..\data\OR13_A2_20161014_CMP_Ecoli_glycerol_exp_2mg_IMAC_R2.tsv

import csv
import argparse
import os
import os.path
import statistics

class CombineList:
    
    def __init__(self):
        # Takes in all the tsv file names in a string format separated by commas
        argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
        argparser.add_argument('tsv_files', type=str, nargs='+', help='Filenames of one or more tsv files to read')
        
        params = argparser.parse_args()

        # Creates an array of all the file names for file in params.files:
        for file in params.tsv_files:
            if not os.path.isfile(file):
                print(f"ERROR: File '{file}' not founds or not a file")
                return
            if file[-3:] != "tsv":
                print(f"ERROR: File '{file}' is not a tsv")
                return
        
        self.file_names = params.tsv_files
            
        # Creates a dictionary with the key being the file name and the value being a list of all the
        # unidentified/formulaic peaks from the tsv file
        self.aggregated_unidentified_peaks = {}

    def read_files(self):
        # Goes through each file to add all observed peaks
        for file_name in self.file_names:
            # Creates a list to store all the observed peaks
            unidentified_peaks = [["mz", "intensity", "", "closest identification"]]
            # Splits the file into separate lines. In this case, each line is a new observed peak
            with open(file_name) as file:
                lines = [line.rstrip() for line in file]
                lines = lines[1:] # Removes header from tsv file
                # Splits each observed peak by tabs, which includes mz and intensity and if applicable,
                # the identifications
                for observed_peak in lines:
                    line_split = [i for i in observed_peak.split("\t")]
                    if file_name != "combined_list.tsv":
                        line_split = line_split[1:] # Removes the uncorrected m/z value
                    # Adds the relevant information (mz, intensity, identifications) to the observed peak
                    # list. The first two elements are the mz and intensity for unknown ions
                    # For known ions, the intensity is the second element whereas the mz is from the first
                    # identification list. Since the identification list is a string, it will be split by
                    # commas then the mz will be used and the identification is stored. If there are other
                    # identifications, those are stored as well 
                    if len(line_split) > 3:
                        if file_name == "combined_list.tsv":
                            primary_identification = line_split[4]
                            mz = round(float(line_split[0]), 5)
                        else:
                            primary_identification = line_split[3].split(", ")[2][1:-2] # -1 accounts for the closing bracket and single quotation mark
                            mz = round(float(line_split[3].split(", ")[0][1:]), 5) # 1 accounts for the opening bracket
                        # If it is a formula, then add it to the list. Otherwise, no need to add it to the list
                        if primary_identification[0] != 'I' and primary_identification[0] != "a" and primary_identification[0] != "b"  and primary_identification[0] != "y":
                            intensity = float(line_split[1])
                            unidentified_peaks.append([mz, intensity, "", primary_identification])
                    else:
                        mz = float(line_split[0])
                        intensity = float(line_split[1])
                        unidentified_peaks.append([mz, intensity, "", "?"])

                self.aggregated_unidentified_peaks[file_name] = unidentified_peaks

        print("finished reading files")

    def rewrite_tsvs(self):
        for key in self.aggregated_unidentified_peaks:
            with open(f'{key}_unidentified.tsv', 'w') as file:
                writer = csv.writer(file, delimiter='\t', lineterminator='\n')
                writer.writerows(self.aggregated_unidentified_peaks[key])

            print(f"finished finding unidentified/formulaic peaks in {key}")

def main():
    combine_list = CombineList()
    # Populates the aggregated_observed_peaks dictionary
    combine_list.read_files()
    combine_list.rewrite_tsvs()

if __name__ == "__main__": main()