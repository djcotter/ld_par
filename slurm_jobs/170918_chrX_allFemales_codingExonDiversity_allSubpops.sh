#!/bin/bash

#SBATCH --job-name=diversityArrayJob
#SBATCH --output=170118_diversityArrayJob_%A_%a.slurm.out
#SBATCH --error=170118_diversityArrayJob_%A_%a.slurm.err
#SBATCH --array=0-25
#SBATCH -n 4                       # number of cores
#SBATCH -t 3-00:00                  # wall time (D-HH:MM)
#SBATCH --mem-per-cpu=4000
#SBATCH -A mwilsons					# Account to charge job to
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=dcotter1@asu.edu # send-to address

date
echo -e "\n"

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

#Print corresponding population code
# This file is a newline seperated file with all of the three letter population codes to be analyzed. 
# Can be set to a list of all populations and subpopulations for an all-in-one analysis.
pops_file=subpopulation_list.txt
readarray -t pops < $pops_file
i=${pops[$SLURM_ARRAY_TASK_ID]}
echo "My POP code: ${i}"

echo -e "\n"

echo "Starting Analysis..."

cd ~/projects/1000_genomes/window_diversity_analysis/

#set window coordinates
window_file=merged_refSeq_codingExons_hg19.bed
window_path=~/projects/1000_genomes/window_diversity_analysis/inputs/windows/
windows=${window_path}${window_file}
win=refseq_Exons

# This will convert all of the provided population data into bed format
if [ -e temps/raw_diversity_by_site/${i}_females_pi_per_site.bed ]; then
	echo "${i} data already converted to BED format..."
else
	python scripts/bedConvert.py inputs/diversity_by_sites/${i}_females_pi_output_per_site.txt temps/raw_diversity_by_site/${i}_females_pi_per_site.bed 
	echo "Converted ${i} data to bed format..."
fi

# The next section will perform the window anlysis using the provided window file, the filtered_called_sites file, and the filtered pop_diversity files.
echo "Ready to perform window analysis using: ${window_file}..."

python scripts/window_calculations.py temps/raw_diversity_by_site/${i}_females_filtered_pi_per_site.bed inputs/callable_sites/20141020.strict_mask.chrX.bed $windows outputs/${i}_females_${win}_diversity.txt
echo "Performed window calculations on ${i} data..."

echo -e "\n"

########################################
echo "Finished Window Analysis..."
echo "Job Complete"

echo -e "\n"
date 

