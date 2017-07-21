from sys import argv
import csv

script, input_file, output_file = argv

with open(input_file, 'rU') as f:
	positions = list(csv.reader(f, delimiter = '\t'))
	
with open(output_file, 'wb') as csvfile:
	writer = csv.writer(csvfile, delimiter = '\t')
	for info in positions:
		writer.writerow(['chr'+info[0], str(int(info[1])-1), info[1], info[2]])