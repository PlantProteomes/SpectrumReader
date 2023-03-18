#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python combine_list.py --tsv_files "..\data\HFX_9850_GVA_DLD1_2_180719_subset.tsv, ..\data\HFX_9850_GVA_DLD1_2_180719.tsv, ..\data\QEX03_210305_PTMscan_wt_IP_63.tsv, ..\data\Q20181210_06.tsv, ..\data\OR13_A2_20161014_CMP_Ecoli_glycerol_exp_2mg_IMAC_R2.tsv"

import csv
import argparse
import os
import os.path
import statistics

class CombineList:
    
    def __init__(self):
        # Takes in all the tsv file names in a string format separated by commas
        argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
        argparser.add_argument('--tsv_files', action='store', type=str, help='Name of the mzML file to read')
        
        params = argparser.parse_args()

        # Creates an array of all the file names
        self.file_names = params.tsv_files.split(", ")

        # Confirms all file names will work
        for file in self.file_names:
            if file is None or file == "":
                print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
                return

            # Tries to open the params file specified. If it fails, an error is printed
            if not os.path.isfile(file):
                print(f"ERROR: File '{file}' not found or not a file")
                return
            
        # Creates a dictionary with the key being the file name and the value being a list of all the
        # observed peak from the tsv file
        self.aggregated_observed_peaks = {}
        self.all_peaks = [["mz", "intensity", "ms runs", "primary identification", "other identifications"]]

    def read_files(self):
        # Goes through each file to add all observed peaks
        removed = 0
        for index in range(len(self.file_names)):
            file_name = self.file_names[index - removed]
            # Creates a list to store all the observed peaks
            observed_peaks = []
            has_histidine = False
            # Splits the file into separate lines. In this case, each line is a new observed peak
            with open(file_name) as file:
                lines = [line.rstrip() for line in file]
                # Splits each observed peak by tabs, which includes mz and intensity and if applicable,
                # the identifications
                for observed_peak in lines:
                    line_split = [i for i in observed_peak.split("\t")]
                    # Adds the relevant information (mz, intensity, identifications) to the observed peak
                    # list. The first two elements are the mz and intensity for unknown ions
                    # For known ions, the intensity is the second element whereas the mz is from the first
                    # identification list. Since the identification list is a string, it will be split by
                    # commas then the mz will be used and the identification is stored. If there are other
                    # identifications, those are stored as well 
                    if len(line_split) > 2:
                        intensity = float(line_split[1])
                        mz = round(float(line_split[2].split(", ")[0][1:]), 5) # 1 accounts for the opening bracket
                        primary_identification = line_split[2].split(", ")[2][1:-2] # -1 accounts for the closing bracket and single quotation mark

                        # If there is no histidine, I end this process early and don't consider it
                        # TODO: Remove it in the future
                        if not has_histidine and primary_identification == "IH":
                            has_histidine = True
                        elif not has_histidine and mz > 110.1:
                            break

                        other_identification = []
                        for index in range(len(line_split) - 3):
                            other_identification.append(line_split[index + 3].split(", ")[2][1:-2])
                        observed_peaks.append([mz, intensity, primary_identification, other_identification])
                    else:
                        mz = float(line_split[0])
                        if not has_histidine and mz > 110.1:
                            break
                        intensity = float(line_split[1])
                        observed_peaks.append([mz, intensity, "?", "?"])

                if has_histidine:
                    self.aggregated_observed_peaks[file_name] = observed_peaks
                else:
                    print(f"removed {file_name} since it did not find a histidine")
                    removed += 1
                    self.file_names.remove(file_name)
                    # TODO: change track the largest intensity/other peak, can use that as the normalized intensity

        print("finished reading files")

    def find_histidines(self): # can eventually rename to find_relative_intensity or something
        # Goes through the aggregated_observed_peaks to find the histidine and stores the intensity to be
        # used as the normalized intensity
        self.histidine_intensity = []
        for file_name in self.file_names:
            for peak in self.aggregated_observed_peaks[file_name]:
                if len(peak) > 2 and peak[2] == "IH":
                    self.histidine_intensity.append(peak[1])

        print("finished finding histidines")

    def merge_lists(self):
        peaks = [self.aggregated_observed_peaks[file_name][0] for file_name in self.file_names]

        for index in range(len(peaks)):
            peaks[index][1] = round(peaks[index][1] * 10000 / self.histidine_intensity[index], 1)

        peak_mzs = [item[0] for item in peaks]
        added_all = False
        while not added_all:
            smallest_mz = min(peak_mzs)
            upper_bound = smallest_mz + 5 * smallest_mz / 1e6

            same_peaks = []
            for index in range(len(self.file_names)):
                file_name = self.file_names[index]
                if peaks[index] == -1:
                    continue
                elif peak_mzs[index] <= upper_bound:
                    same_peaks.append(peaks[index])
                    self.aggregated_observed_peaks[file_name].pop(0)
                    if len(self.aggregated_observed_peaks[file_name]) == 0:
                        peaks[index] = -1
                        peak_mzs[index] = float('inf')
                        del self.aggregated_observed_peaks[file_name]
                    else:
                        peaks[index] = self.aggregated_observed_peaks[file_name][0]
                        peaks[index][1] = round(peaks[index][1] * 10000 / self.histidine_intensity[index], 1)
                        peak_mzs[index] = peaks[index][0]

            peak = []
            if len(same_peaks) == 1:
                peak = same_peaks[0]
                peak.insert(2, 1)
            else:
                peak.append(len(same_peaks)) # adds number of ms runs it appears in
                same_peaks.sort(key = lambda x: -x[1])
                peak.append(same_peaks[0][2]) # adds primary explanation
                peak.append(same_peaks[0][3]) # adds other explanations
                mzs = [item[0] for item in same_peaks]
                mzs.sort() # sort by mz value
                peak.insert(0, round(statistics.median(mzs), 5)) # add the median mz value
                intensities = [item[1] for item in same_peaks]
                intensities.sort()
                peak.insert(1, round(statistics.median(intensities), 1)) # adds the intensity value
            self.all_peaks.append(peak)

            if len(self.aggregated_observed_peaks) == 0:
                added_all = True

        print("finished merging lists")

    def write_combined_list(self):
        with open('combined_list.tsv', 'w') as file:
            writer = csv.writer(file, delimiter='\t', lineterminator='\n')
            writer.writerows(self.all_peaks)

        print("finished writing tsv file")

def main():
    combine_list = CombineList()
    # Populates the aggregated_observed_peaks dictionary
    combine_list.read_files()
    combine_list.find_histidines()
    combine_list.merge_lists()
    combine_list.write_combined_list()

if __name__ == "__main__": main()