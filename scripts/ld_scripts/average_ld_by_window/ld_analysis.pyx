
import numpy as np
import time

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
			values.append(pair[1])
	mean_LD = np.mean(values)
	return mean_LD


def median_LD(window_LD_dictionary):
	"""Calculate median LD in a window by taking all focal positions and their associated distance and R2 values
		The structure of the data is {fp: [[distance, R2], [distance, R2] ...]}
	"""
	values = []
	for focal_position in window_LD_dictionary:
		for pair in window_LD_dictionary[focal_position]:
			values.append(pair[1])
	median_LD = np.median(values)
	return median_LD


def LD_percentiles(window_LD_dictionary, lower, upper):
	"""Calculate percintiles of LD values in a window by taking all focal positions and their associated distance and R2 values
		The structure of the data is {fp: [[distance, R2], [distance, R2] ...]}
	"""
	values = []
	for focal_position in window_LD_dictionary:
		for pair in window_LD_dictionary[focal_position]:
			values.append(pair[1])
	quartiles = [np.percentile(values, lower), np.percentile(values, upper)]
	return quartiles


def R2_values_to_array(window_LD_dictionary):
	"""Function to take the R2_Values dictionary for a specific window and return a numpy 
		array of all R2 values within the window
	"""
	values = []
	for fp in window_LD_dictionary:
		for pair in window_LD_dictionary[fp]:
			values.append(pair[1])
	return np.asarray(values)


def rand_samples(data_array):
	"""
	Create a randomly dsitributed set of data based on a given array and the length of the given data array
	"""
	array_len = len(data_array)
	indices = np.random.random_integers(0, array_len - 1, array_len)
	return np.take(data_array, indices)


def bootstrap_CI_mean(data_array, replicates):
	resamples = []
	for i in range(replicates):
		resamples.append(np.mean(rand_samples(data_array)))
	return (np.nanpercentile(resamples, 2.5), np.nanpercentile(resamples, 97.5))


def decay_value(window_LD_dictionary):
	pass


#Loop for passing through the LD summary file

def LD_loop(input_file, window_coordinates, int ld_bin_size):
	cdef int win_num, fp, ld_bin, row_1, row_5
	cdef float row_8

	win_num = 0 
	windows = window_coordinates[win_num]
	R2_values = {}
	reverse_R2 = {} #The information that is stored as pairwise is not duplicated. We have to create a matrix that contains values for BP_B > window[i][2] but BP_B - BP_A < ld_bin
	ld_bin = ld_bin_size * 1000
	fp = 0 #initalize the focal position as 0
	results = [] # initialize an empty results list
	with open(input_file, 'rU') as f:
		for line in f:
			row = line.split()

			# skip the header line of the file
			if row[0] == 'CHR_A':
				continue

			#declare all relevant values
			a1 = int(row[1]) #BP_A
			a5 = int(row[5]) #BP_B
			a8 = float(row[8]) #R2
			row_1 = a1
			row_5 = a5
			row_8 = a8

			# initializes the focal position for the first line of the file	
			if fp == 0: 
				fp = row_1
				R2_values[fp] = []

			# change script paramaters and initialize a new list when reaching a new focal position
			if row_1 != fp: 
				fp = row_1

				#if the focal position lies outside of the current window, change the window index, calculate a summary statistic, and clear R2_values from memory
				if (fp > windows[2]):
					bootstrap = bootstrap_CI_mean(R2_values_to_array(R2_values), 1000)
					results.append([windows[1], windows[2], mean_LD(R2_values), bootstrap[0], bootstrap[1]])
					print windows[1]
					R2_values = {}
					win_num += 1
					windows = window_coordinates[win_num]

				R2_values[fp] = []

			# if the site being compared is within ld_bin of the focal position, store the distance and R2 value in the dictionary as a list
			if (abs(row_5 - row_1) < ld_bin):
				R2_values[fp].append([row_5 - row_1, row_8])
			
			# if the site in position 2 lies outside the current window, but (BP_B - BP_A) lies within ld_bin, we need to save the information for analyses of future windows
			if (row_5 > windows[2]) and (abs(row_5 - row_1) < ld_bin):
				if row_5 in reverse_R2:
					reverse_R2[row_5].append([row_5 - row_1, row_8])	
				else:
					reverse_R2[row_5] = [[row_5 - row_1, row_8]]

			# check the reverse list to see if there are any saved data for the current focal position, then combine the two lists
			if fp in reverse_R2:
				R2_values[fp] = R2_values[fp] + reverse_R2[fp]
				del reverse_R2[fp]

	# once the last line has been reached, R2_values will have all the information correspoding to the last window
	# add the summary of these values to the results file
	bootstrap = bootstrap_CI_mean(R2_values_to_array(R2_values), 1000)
	results.append([windows[1], windows[2], mean_LD(R2_values), bootstrap[0], bootstrap[1]])
	return results

