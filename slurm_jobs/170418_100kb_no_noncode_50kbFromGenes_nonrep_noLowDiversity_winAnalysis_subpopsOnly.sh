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


# This file is the filter to be used for this analysis.
filter_file=noncode_nonrep_50kb_flanked_hg19_filter.bed
filter_path=~/projects/1000_genomes/window_diversity_analysis/inputs/filters/
filter=${filter_path}${filter_file}

#This file will be the windows that the analysis will use.
window_file=chrX_100kb_windows.bed
window_path=~/projects/1000_genomes/window_diversity_analysis/inputs/windows/
windows=${window_path}${window_file}
win=100kb_no

# This will filter the callable sites file and save the filtered version to a temp directory.
~/tools/bedtools/bin/bedtools subtract -a inputs/callable_sites/20141020.strict_mask.chrX.bed -b $filter > temps/filtered_callable_sites/filtered_called_sites.bed
echo "Filtered the callable_sites mask using ${filter_file}..."

# This will convert all of the provided population data into bed format
if [ -e temps/raw_diversity_by_site/${i}_females_pi_per_site.bed ]; then
	echo "${i} data already converted to BED format..."
else
	python scripts/bedConvert.py inputs/diversity_by_sites/${i}_females_pi_output_per_site.txt temps/raw_diversity_by_site/${i}_females_pi_per_site.bed 
	echo "Converted ${i} data to bed format..."
fi


# This will filter the converted population bed files by removing sites that overlap the filter.
~/tools/bedtools/bin/bedtools intersect -a temps/raw_diversity_by_site/${i}_females_pi_per_site.bed -b temps/filtered_callable_sites/filtered_called_sites.bed > temps/filtered_diversity_by_site/${i}_females_filtered_pi_per_site.bed
echo "Filtered ${i} data using filtered_called_sites file..."


# The next section will perform the window anlysis using the provided window file, the filtered_called_sites file, and the filtered pop_diversity files.
echo "Ready to perform window analysis using: ${window_file}..."

python scripts/window_calculations.py temps/filtered_diversity_by_site/${i}_females_filtered_pi_per_site.bed temps/filtered_callable_sites/filtered_called_sites.bed $windows outputs/${i}_females_${win}_50kb_nearGenes_diversity.txt
echo "Performed window calculations on ${i} data..."

echo -e "\n"

########################################
echo "Finished Window Analysis..."
echo "Job Complete"

echo -e "\n"
date 
