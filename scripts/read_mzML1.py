#!/usr/bin/env python3

# must first change to scripts folder through "cd .\scripts\"
# to run the program (use terminal): python read_mzML1.py --mzml_file ..\data\HFX_9850_GVA_DLD1_2_180719.mzML
# can do python read_mzML1.py --help and look at the usage statements for what I can do

import os
import argparse
import os.path
import timeit
import matplotlib.pyplot as plt
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup

from pyteomics import mzml, auxiliary

def main():

    # sets up a module with a description
    argparser = argparse.ArgumentParser(description='An example program that reads an mzML file sequentially')
    # tells there is a mzml file to access and the user can input data that params will store through the '--mzml_file ..\data\HFX...' in the terminal
    argparser.add_argument('--mzml_file', action='store', help='Name of the mzML file to read')
    # argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    # can use params.precursor_mz to access what the precursor_mz value is
    argparser.add_argument('--precursor_mz', action='store', type=float, help='Precursor m/z to select')
    # runs the parse_args method to store all the results into the params variable
    params = argparser.parse_args()

    #### Ensure that mzml_file was passed, look in the params variable for the specified mzml_file. if you want to open precursor, it would be params.precursor_mz
    if params.mzml_file is None or params.mzml_file == "":
        print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
        return

    # tries to open the params file specified
    if not os.path.isfile(params.mzml_file):
        print(f"ERROR: File '{params.mzml_file}' not found or not a file")
        return


    #### Read spectra from the file
    t0 = timeit.default_timer()
    stats = { 'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0 }
    with mzml.read(params.mzml_file) as reader:
        for spectrum in reader:
            if stats['counter'] == 0:
                auxiliary.print_tree(spectrum)

            #### Update counters and print progress
            stats['counter'] += 1

            # check to see how i can access the ms level


            if stats['counter'] == 2:
                print(spectrum)
                print(len(spectrum['intensity array']))
                fig, ax = plt.subplots(figsize=(12, 6))
                # creating variables for the required arguments
                precursor_mz = (((((((spectrum['precursorList'])['precursor']))[0])['selectedIonList'])['selectedIon'])[0])['selected ion m/z']
                # precursor_mz = 0.0
                precursor_charge = (((((((spectrum['precursorList'])['precursor']))[0])['selectedIonList'])['selectedIon'])[0])['charge state']
                # precursor_charge = 0
                plotSpectrum = sus.MsmsSpectrum(f"Spectrum {stats['counter']}", precursor_mz, precursor_charge, spectrum['m/z array'], spectrum['intensity array'])
                plotSpectrum = plotSpectrum.set_mz_range(min_mz=min(spectrum['m/z array']), max_mz=max(spectrum['m/z array']))
                sup.spectrum(plotSpectrum)
                plt.title(f"Spectrum {stats['counter']}, precursor mz = {precursor_mz}, precursor charge = {precursor_charge}")
                plt.show()
                plt.close()


            if spectrum['ms level'] == 1:
                stats['ms1spectra'] += 1
            elif spectrum['ms level'] == 2:
                stats['ms2spectra'] += 1

            if stats['counter']/1000 == int(stats['counter']/1000):
                print(f"  {stats['counter']}")


    #### Print final timing information
    t1 = timeit.default_timer()
    print(f"INFO: Read {stats['counter']} spectra from {params.mzml_file}")
    print(f"The number of ms1spectra is {stats['ms1spectra']}")
    print(f"The number of ms2spectra is {stats['ms2spectra']}")
    print(f"INFO: Elapsed time: {t1-t0}")
    print(f"INFO: Processed {stats['counter']/(t1-t0)} spectra per second")


if __name__ == "__main__": main()
