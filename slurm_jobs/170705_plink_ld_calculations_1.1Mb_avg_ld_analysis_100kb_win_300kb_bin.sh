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

cd ~/projects/1000_genomes/ld_par_boundary/
export PATH=$PATH:~/tools/bcftools/bin/
export PATH=$PATH:~/tools/plink2/bin/

MYTMPDIR=$(mktemp -p /scratch/dcotter1 -d)
function clean_up {
	# Perform program exit housekeeping
	ls -lh $MYTMPDIR
	rm -rf "$MYTMPDIR"
}
trap clean_up EXIT

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

bcftools view --samples-file ~/projects/1000_genomes/populations/subpops/${i}_females -Ov -o $MYTMPDIR/${i}_females_unfiltered_variants.vcf ~/projects/1000_genomes/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf
bcftools view -m2 -M2 -v snps -Ov -o $MYTMPDIR/${i}_females_biallelic_variants.vcf $MYTMPDIR/${i}_females_unfiltered_variants.vcf
vcftools --vcf $MYTMPDIR/${i}_females_biallelic_variants.vcf --mac 1 --recode --out $MYTMPDIR/${i}_females_biallelic_variants_mac_filtered
cp $MYTMPDIR/${i}_females_biallelic_variants_mac_filtered.log log_files/${i}_females_biallelic_variants_mac_filtered.log

plink --vcf $MYTMPDIR/${i}_females_biallelic_variants_mac_filtered.recode.vcf --out plink_files/${i}_females_biallelic_mac_filtered
plink --bfile plink_files/${i}_females_biallelic_mac_filtered --ld_window 700000 --ld_window-kb 1100 --memory 16000 --r2 with-freqs --out $MYTMPDIR/${i}_females_biallelic_mac_filtered_ld_R2
cp $MYTMPDIR/${i}_females_biallelic_mac_filtered_ld_R2.log log_files/${i}_females_biallelic_mac_filtered_ld_R2.log

python ld_scripts/170622_average_ld_by_window/average_ld_by_window/average_ld_by_window.py $MYTMPDIR/${i}_females_biallelic_mac_filtered_ld_R2.ld ~/projects/1000_genomes/window_diversity_analysis/inputs/windowschrX_100kb_windows.bed 300 results/170705_avg_LD_100kb_no_300kb_LDbins/${i}_females_avg_R2_100kb_windows_300kb_ldBins_95bootstrapCI.txt

date

