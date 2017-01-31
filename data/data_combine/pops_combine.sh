#!/bin/bash

MYTMPDIR=$(mktemp -d)

function clean_up {
	# Perform program exit housekeeping
	rm -rf "$MYTMPDIR"
}
trap clean_up EXIT

pops_file=subpopulation_list.txt
declare -a pops

# Load file into array.
let i=0
while IFS=$'\n' read -r line_data; do
    pops[i]="${line_data}"
    ((++i))
done < $pops_file

for i in ${pops[@]}
do
	echo $i > $MYTMPDIR/${i}_temp_pi.txt
	awk '{print (($5~/NA/)?"NA":$4)}' ../${i}_females_100kb_no_diversity.txt >> $MYTMPDIR/${i}_temp_pi.txt
done

paste $(for i in ${pops[@]}; do echo -n $MYTMPDIR/${i}_temp_pi.txt ""; done) > output_file.txt