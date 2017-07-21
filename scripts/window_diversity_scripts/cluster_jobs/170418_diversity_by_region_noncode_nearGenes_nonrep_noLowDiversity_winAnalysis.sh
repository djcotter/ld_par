#!/bin/bash

#SBATCH --job-name=diversityArrayJob
#SBATCH --output=170118_diversityArrayJob_%A_%a.slurm.out
#SBATCH --error=170118_diversityArrayJob_%A_%a.slurm.err
#SBATCH --array=0-2
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

echo -e "\n"

echo "Starting Analysis..."

cd ~/projects/1000_genomes/window_diversity_analysis/

# This file is the filter to be used for this analysis.
nearGenes=(10 25 50)
flank=${nearGenes[$SLURM_ARRAY_TASK_ID]}
filter_file=noncode_nonrep_${flank}kb_flanked_noLowDiverstiyAmp_hg19_filter.bed
filter_path=~/projects/1000_genomes/window_diversity_analysis/inputs/filters/
filter=${filter_path}${filter_file}

#label pop code + flank 
i=ALL
echo "My POP code: ${i}"

#This file will be the windows that the analysis will use.
window_file=PAR_Regions
# window_path=~/projects/1000_genomes/window_diversity_analysis/inputs/windows/
# windows=${window_path}${window_file}
win=by_region

# This will filter the callable sites file and save the filtered version to a temp directory.
~/tools/bedtools/bin/bedtools subtract -a inputs/callable_sites/20141020.strict_mask.chrX.bed -b $filter > temps/filtered_callable_sites/${flank}_filtered_called_sites.bed
echo "Filtered the callable_sites mask using ${filter_file}..."

# This will convert all of the provided population data into bed format
if [ -e temps/raw_diversity_by_site/${i}_females_pi_per_site.bed ]; then
	echo "${i} data already converted to BED format..."
else
	python scripts/bedConvert.py inputs/diversity_by_sites/${i}_females_pi_output_per_site.txt temps/raw_diversity_by_site/${i}_females_pi_per_site.bed 
	echo "Converted ${i} data to bed format..."
fi

# This will filter the converted population bed files by removing sites that overlap the filter.
~/tools/bedtools/bin/bedtools intersect -a temps/raw_diversity_by_site/${i}_females_pi_per_site.bed -b temps/filtered_callable_sites/${flank}_filtered_called_sites.bed > temps/filtered_diversity_by_site/${i}_${flank}_females_filtered_pi_per_site.bed
echo "Filtered ${i} data using filtered_called_sites file..."

# The next section will perform the window anlysis using the provided window file, the filtered_called_sites file, and the filtered pop_diversity files.
echo "Ready to perform window analysis using: ${window_file}..."

# The difference is that this script has built in windows and addresses the overlap of the nonPAR with the XTR; the windows are built in to the script.
python scripts/window_calculations_PAR_regions.py temps/filtered_diversity_by_site/${i}_${flank}_females_filtered_pi_per_site.bed temps/filtered_callable_sites/${flank}_filtered_called_sites.bed outputs/ALL_females_diversity_${flank}kb_nearGenes_${win}.txt
echo "Performed window calculations on ${i} data..."

echo -e "\n"

########################################
echo "Finished Window Analysis..."
echo "Job Complete"

echo -e "\n"
date 
