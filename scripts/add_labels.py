
#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python add_labels.py --tsv_file ..\data\combined_list.tsv

import csv
import argparse
import os
import os.path
import statistics

class LabelAdder:
    def __init__(self):
        # Takes in all the tsv file names in a string format separated by commas
        argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
        argparser.add_argument('--tsv_file', action='store', help='Name of the tsv file to read')
        
        params = argparser.parse_args()

        # Creates an array of all the file names for file in params.files:
        if not os.path.isfile(params.tsv_file):
            print(f"ERROR: File '{params.tsv_file}' not founds or not a file")
            return
        if params.tsv_file[-3:] != "tsv":
            print(f"ERROR: File '{params.tsv_file}' is not a tsv")
            return
        
        self.file_name = params.tsv_file
        self.all_peaks = []

    def read_files(self):
        first_line = True
        with open(self.file_name) as file:
            lines = [line.rstrip() for line in file]
            for observed_peak in lines:
                line_split = [i for i in observed_peak.split("\t")]
                row = []
                for column in line_split:
                    row.append(column)

                if first_line:
                    row.append("labels")
                    first_line = False
                else:
                    is_common = False
                    is_immonium = False
                    fraction = row[2].split("/")
                    if float(fraction[0])/float(fraction[1]) == 1.0 or float(row[1]) >= 5000.0:
                        is_common = True
                    if row[6][0] == "I":
                        is_immonium = True
                    if is_common and is_immonium:
                        row.append("COM, IMM")
                    elif is_common:
                        row.append("COM")
                    elif is_immonium:
                        row.append("IMM")
                    else:
                        row.append("")

                self.all_peaks.append(row)

    def rewrite_file(self):
        with open(f'{self.file_name}', 'w') as file:
        # with open('common_peaks.tsv', 'w') as file:
            writer = csv.writer(file, delimiter='\t', lineterminator='\n')
            writer.writerows(self.all_peaks)

def main():
    add_label = LabelAdder()
    # Populates the aggregated_observed_peaks dictionary
    add_label.read_files()
    add_label.rewrite_file()

if __name__ == "__main__": main()