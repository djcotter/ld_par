#!/bin/bash

#SBATCH --job-name=ld_calculations
#SBATCH --array=0-25
#SBATCH -n 4                       # number of cores
#SBATCH -t 3-00:00                  # wall time (D-HH:MM)
#SBATCH --mem-per-cpu=4000
#SBATCH -A mwilsons					# Account to charge job to
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=dcotter1@asu.edu # send-to address

date

cd ~/projects/1000_genomes/outputs/all_individuals_subpop_diversity_by_site/chrX
source activate exome_diversity

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

#Print corresponding population code
# This file is a newline seperated file with all of the three letter population codes to be analyzed. 
# Can be set to a list of all populations and subpopulations.
pops_file=~/projects/1000_genomes/populations/subpopulation_list.txt
readarray -t pops < $pops_file
i=${pops[$SLURM_ARRAY_TASK_ID]}

echo "My POP code: ${i}"

echo -e "\n"

echo "Starting Analysis..."

python ~/projects/1000_genomes/scripts/170721_Diversity_from_VCF_cyvcf_Output_pi_per_site_per_population_ignoreINDELs_Hohenlohe_check_for_haploid_sites.py --vcf ~/projects/1000_genomes/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf --population_lists ~/projects/1000_genomes/populations/subpops/${i}_individuals --chrom_inc 8

date