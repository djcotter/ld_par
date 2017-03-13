from sys import argv
import csv

script, input_file = argv
#script, input_file, windows, ld_window_size, analysis_type, output_file = argv

# Open the window file as bed coordinates so that every window_coordinates[i][1] corresponds to the start position and [2] is the end position.
with open(windows, 'rU') as f:
	window_coordinates = list(csv.reader(f, delimiter = '\t'))

#open the input file to start the analysis and leave it open to save space
with open(input_file, 'rU') as f:
	for line in f:
		row = line.split()
