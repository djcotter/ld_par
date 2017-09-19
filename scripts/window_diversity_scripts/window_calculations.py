from sys import argv
import csv

script, pop_diversity, called_sites, windows, output_file = argv

# Open the window file as bed coordinates so that every window_coordinates[i][1] corresponds to the start position and [2] is the end position.
with open(windows, 'rU') as f:
	window_coordinates = list(csv.reader(f, delimiter = '\t'))
	
with open(pop_diversity, 'rU') as f:
	diversity = list(csv.reader(f, delimiter = '\t'))
	
with open(called_sites, 'rU') as f:
	callable = list(csv.reader(f, delimiter = '\t'))
	
data = []

for c in window_coordinates:
	sum_called = 0
	count = 0
	for r in callable:
		called = 0
		if int(r[1]) >= int(c[1]) and int(r[2]) <= int(c[2]):
			called = int(r[2]) - int(r[1])
			sum_called += called
		if int(r[1]) >= int(c[1]) and int(r[1]) < int(c[2]) and int(r[2]) > int(c[2]):
			called = int(c[2]) - int(r[1])
			sum_called += called
		if int(r[1]) < int(c[1]) and int(r[2]) <= int(c[2]) and int(r[2]) > int(c[1]):
			called = int(r[2]) - int(c[1])
			sum_called += called
		else:
			sum_called += 0
	pi_sum = float(0)
	for d in diversity:
		if int(d[2]) > int(c[1]) and int(d[2]) <= int(c[2]):
			pi_sum += float(d[3])
			count += 1
	if sum_called > 0:
		data.append(['chrX', c[1], c[2], float(pi_sum/sum_called), sum_called, count])
	else:
		data.append(['chrX', c[1], c[2], pi_sum, "NA"])
	
		
with open(output_file, 'wb') as csvfile:
	writer = csv.writer(csvfile, delimiter = '\t')
	for info in data:
		writer.writerow(info)