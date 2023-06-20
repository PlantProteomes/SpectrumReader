#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python combine_list.py ..\data\HFX_9850_GVA_DLD1_2_180719_subset.tsv, ..\data\HFX_9850_GVA_DLD1_2_180719.tsv, ..\data\QEX03_210305_PTMscan_wt_IP_63.tsv, ..\data\Q20181210_06.tsv, ..\data\OR13_A2_20161014_CMP_Ecoli_glycerol_exp_2mg_IMAC_R2.tsv

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
        # observed peak from the tsv file
        self.aggregated_observed_peaks = {}
        self.all_peaks = [["mz", "intensity", "ms runs", "percentage", "theoretical m/z", "delta PPM", "primary identification", "other identifications"]]

    def read_files(self):
        # Goes through each file to add all observed peaks
        removed = 0
        self.smallest_mz = []
        for index in range(len(self.file_names)):
            file_name = self.file_names[index - removed]
            # Creates a list to store all the observed peaks
            observed_peaks = []
            has_histidine = False
            # Splits the file into separate lines. In this case, each line is a new observed peak
            with open(file_name) as file:
                lines = [line.rstrip() for line in file]
                lines = lines[1:] # Removes header from tsv file
                # Splits each observed peak by tabs, which includes mz and intensity and if applicable,
                # the identifications
                for observed_peak in lines:
                    line_split = [i for i in observed_peak.split("\t")]
                    line_split = line_split[1:] # Removes the uncorrected m/z value
                    # Adds the relevant information (mz, intensity, identifications) to the observed peak
                    # list. The first two elements are the mz and intensity for unknown ions
                    # For known ions, the intensity is the second element whereas the mz is from the first
                    # identification list. Since the identification list is a string, it will be split by
                    # commas then the mz will be used and the identification is stored. If there are other
                    # identifications, those are stored as well 
                    if len(line_split) > 3:
                        intensity = float(line_split[1])
                        mz = round(float(line_split[0]), 5)
                        percent = float(line_split[2][0:-1])
                        primary_identification = [line_split[3].split(", ")[2][1:-2], line_split[3].split(", ")[0][1:]] # -1 accounts for the closing bracket and single quotation mark

                        # If there is no histidine, I end this process early and don't consider it
                        # TODO: Remove it in the future
                        if not has_histidine and primary_identification[0] == "IH":
                            has_histidine = True
                        elif not has_histidine and mz > 110.1:
                            break

                        other_identification = []
                        for index in range(len(line_split) - 4):
                            other_identification.append(line_split[index + 4].split(", ")[2][1:-2])
                        observed_peaks.append([mz, intensity, percent, primary_identification, other_identification])
                    else:
                        mz = float(line_split[0])
                        if not has_histidine and mz > 110.1:
                            break
                        intensity = float(line_split[1])
                        percent = float(line_split[2][0:-1])
                        observed_peaks.append([mz, intensity, percent, "?", "?"])

                if has_histidine:
                    self.aggregated_observed_peaks[file_name] = observed_peaks
                    lines = lines[0].split("\t")
                    self.smallest_mz.append(float(lines[1]) - 0.01) # if there are multiple peaks that are
                    # the same, the median will be taken for the combined_list plot. if there are two peaks,
                    # the average is taken. if this is the case for the very first peak, comparing the
                    # peak to the minimum of each file will not work, as the median value will now be
                    # less than one but greater than the other. subtracting 0.01 will ensure taking the median
                    # will not affect the accuracy of how many ms runs the peak may exist in, but is small enough
                    # that it will not affect any other calculations
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
                if len(peak) > 3 and peak[3][0] == "IH":
                    self.histidine_intensity.append(peak[1])

        print("finished finding histidines")

    def merge_lists(self):
        # accesses the first peak in all files
        peaks = [self.aggregated_observed_peaks[file_name][0] for file_name in self.file_names]

        # round the intensity to 1 decimal point after normalizing it to the histidine intensity
        for index in range(len(peaks)):
            peaks[index][1] = round(peaks[index][1] * 10000 / self.histidine_intensity[index], 1)

        # take the m/z and find those with the same m/z
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
                mz = peak[0]
                if len(peak[3]) > 1: # accounts for if peak is not "?" as it also has the m/z value
                    theoretical_mz = float(peak[3][1])
                    peak[3] = peak[3][0]
                    peak.insert(3, theoretical_mz) # adds the theoretical m/z value
                    peak.insert(4, round((theoretical_mz - mz) / mz * 1e6, 2)) # adds delta PPM
                else: # accounts for if peak is "?"
                    peak.insert(3, "?") # adds the theoretical m/z value
                    peak.insert(4, "?") # adds delta PPM
                peak.insert(2, f"1/{self.find_denominator(mz)}")
                peak[3] = f"{peak[3]}%"
            else:
                mzs = [item[0] for item in same_peaks]
                mzs.sort() # sort by mz value
                median_mz = round(statistics.median(mzs), 5)
                peak.append(f"{len(same_peaks)}/{self.find_denominator(median_mz)}") # adds number of ms runs it appears in
                same_peaks.sort(key = lambda x: -x[1]) # sorts the array of all peaks by intensity
                if len(same_peaks[0][3]) > 1: # accounts for if peak is not "?"
                    theoretical_mz = float(same_peaks[0][3][1]) # takes the most intense peak to use for data
                    peak.append(theoretical_mz) # adds the theoretical m/z value
                    peak.append(round((theoretical_mz - median_mz) / median_mz * 1e6, 2)) # adds delta PPM
                    peak.append(same_peaks[0][3][0]) # adds primary explanation
                else: # accounts for if peak is "?"
                    peak.append("?")
                    peak.append("?")
                    peak.append("?")
                peak.append(same_peaks[0][4]) # adds other explanations
                peak.insert(0, median_mz) # add the median mz value
                intensities = [item[1] for item in same_peaks]
                intensities.sort()
                peak.insert(1, round(statistics.median(intensities), 1)) # adds the intensity value
                percentages = [item[2] for item in same_peaks]
                peak.insert(3, f"{round(statistics.mean(percentages), 1)}%")
            self.all_peaks.append(peak)

            if len(self.aggregated_observed_peaks) == 0:
                added_all = True

        print("finished merging lists")

    def find_denominator(self, mz):
        # returns the number of ms runs the peak could have appeared in. for example, if the peak is 90 m/z
        # and appears in 5/6 ms runs but the 6th ms run starts data collection at 100 m/z, the 6th ms run
        # is not counted
        possible_ms_runs = 0
        for minimum_peak in self.smallest_mz:
            if mz > minimum_peak:
                possible_ms_runs += 1
        return possible_ms_runs

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