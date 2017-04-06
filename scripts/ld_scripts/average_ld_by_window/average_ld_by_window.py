"""The purpose of this script is to take pairwise LD output from plink and calculate average LD in set windows
	Updated: 3/30/17
	Daniel Cotter
"""

# import required packages
from sys import argv
import csv
import numpy as np
import time
import ld_analysis

############################
## Set argument variables ##
############################
""" input file: pairwise LD file generated in space-delimited table format 
				by plink using the command 'with-freqs' (essential so the columns are correctly indexed)
	windows: bed file with desired windows for the analysis
	ld_bin_size: length in either direction from the focal position in which we should consider LD (in kilobases)
	analysis type: type of summary file to be generated using the desired options
	output_file: path to an output file
"""
script, input_LD, windows_input, ld_bin_length, output_file = argv
ld_bin_length = int(ld_bin_length)
#script, input_file, windows, ld_bin_size, analysis_type, output_file = argv

######################################################
## Read in windows file and set starting conditions ##
######################################################
# Open the window file as bed coordinates so that every window_coordinates[i][1] corresponds to the start position and [2] to the end position.
with open(windows_input, 'rU') as f:
	window_file = list(csv.reader(f, delimiter = '\t'))
	for window in window_file:
			window[1] = float(window[1])
			window[2] = float(window[2]) 
			del window[3]

############################################################
## Open the LD file and perform the analysis line by line ##
############################################################
ld_calculations = ld_analysis.LD_loop(input_LD, window_file, ld_bin_length)

with open(output_file, 'w') as OUT:
	csvwriter = csv.writer(OUT, delimiter='\t')
	for line in ld_calculations:
		csvwriter.writerow(line)