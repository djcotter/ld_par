"""The purpose of this script is to take pairwise LD output from plink and calculate average LD in set windows
	Updated: 3/20/17
	Daniel Cotter
"""

# import required packages
from sys import argv
import csv
import numpy as np

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
script, input_file, windows, ld_bin_size = argv
#script, input_file, windows, ld_bin_size, analysis_type, output_file = argv

######################################################
## Read in windows file and set starting conditions ##
######################################################
# Open the window file as bed coordinates so that every window_coordinates[i][1] corresponds to the start position and [2] to the end position.
with open(windows, 'rU') as f:
	window_coordinates = list(csv.reader(f, delimiter = '\t'))
	for window in window_coordinates:
			window[1] = float(window[1])
			window[2] = float(window[2]) 
			del window[3]

win_num = 0 
windows = window_coordinates[win_num]
R2_values = {}
reverse_R2 = {} #The information that is stored as pairwise is not duplicated. We have to create a matrix that contains values for BP_B > window[i][2] but BP_B - BP_A < ld_bin
ld_bin = float(ld_bin_size) * 1000
fp = 0 #initalize the focal position as 0
results = [] # initialize an empty results list


###########################################
## Functions for compiling summary stats ##
###########################################
def mean_LD(window_LD_dictionary):
	"""Calculate mean LD in a window by taking all focal positions and their associated distance and R2 values
		The structure of the data is {fp: [[distance, R2], [distance, R2] ...]}
	"""
	values = []
	for focal_position in window_LD_dictionary:
		for pair in window_LD_dictionary[focal_position]:
			values.append(pair[0])
	mean_LD = np.mean(values)
	return mean_LD
	
def median_LD(window_LD_dictionary):
	"""Calculate median LD in a window by taking all focal positions and their associated distance and R2 values
		The structure of the data is {fp: [[distance, R2], [distance, R2] ...]}
	"""
	values = []
	for focal_position in window_LD_dictionary:
		for pair in window_LD_dictionary[focal_position]:
			values.append(pair[0])
	median_LD = np.median(values)
	return median_LD

def LD_percentiles(window_LD_dictionary, lower, upper):
	"""Calculate percintiles of LD values in a window by taking all focal positions and their associated distance and R2 values
		The structure of the data is {fp: [[distance, R2], [distance, R2] ...]}
	"""
	values = []
	for focal_position in window_LD_dictionary:
		for pair in window_LD_dictionary[focal_position]:
			values.append(pair[0])
	quartiles = [np.percentile(values, lower), np.percentile(values, upper)]
	return quartiles

def bootstrap_CI(window_LD_dictionary, replicates):
	pass

def decay_value(window_LD_dictionary):
	pass

############################################################
## Open the LD file and perform the analysis line by line ##
############################################################
with open(input_file, 'rU') as f:
	for line in f:
		row = line.split()
		
		# skip the header line of the file
		if row[0] == 'CHR_A':
			continue

		# initializes the focal position for the first line of the file	
		if fp == 0: 
			fp = float(row[1])
			R2_values[fp] = []
			LD_bin = [fp - ld_bin, fp + ld_bin]

		# change script paramaters and initialize a new list when reaching a new focal position
		if float(row[1]) != fp: 
			fp = float(row[1])

			#if the focal position lies outside of the current window, change the window index, calculate a summary statistic, and clear R2_values from memory
			if (fp > windows[2]):
				results.append([windows[1], windows[2], mean_LD(R2_values)])
				print windows #temp
				
				R2_values = {}
				win_num += 1
				windows = window_coordinates[win_num]

			R2_values[fp] = []


		# if the site being compared is within ld_bin of the focal position, store the distance and R2 value in the dictionary as a list
		if (abs(float(row[5]) - float(row[1])) < ld_bin):
			R2_values[fp].append([float(row[5])-float(row[1]), float(row[8])])
		
		# if the site in position 2 lies outside the current window, but (BP_B - BP_A) lies within ld_bin, we need to save the information for analyses of future windows
		if (float(row[5]) > windows[2]) and (abs(float(row[5]) - float(row[1])) < ld_bin):
			if float(row[5]) in reverse_R2:
				reverse_R2[float(row[5])].append([float(row[5])-float(row[1]), float(row[8])])	
			else:
				reverse_R2[float(row[5])] = [[float(row[5])-float(row[1]), float(row[8])]]

		# check the reverse list to see if there are any saved data for the current focal position, then combine the two lists
		if fp in reverse_R2:
			R2_values[fp] = R2_values[fp] + reverse_R2[fp]
			del reverse_R2[fp]

# once the last line has been reached, R2_values will have all the information correspoding to the last window
# add the summary of these values to the results file
results.append([windows[1], windows[2], mean_LD(R2_values)])

for line in results:
	print line
