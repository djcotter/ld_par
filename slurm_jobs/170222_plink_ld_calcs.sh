#!/bin/bash

#SBATCH --job-name=ld_calculations
#SBATCH -n 8                        # number of cores
#SBATCH -t 3-00:00                  # wall time (D-HH:MM)
#SBATCH -A mwilsons					# Account to charge job to
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=dcotter1@asu.edu # send-to address

date

cd ~/projects/1000_genomes/ld_par_boundary
export PATH=$PATH:~/tools/bcftools/bin/
export PATH=$PATH:~/tools/plink2/bin/

MYTMPDIR=$(mktemp -d)
function clean_up {
	# Perform program exit housekeeping
	rm -rf "$MYTMPDIR"
}
trap clean_up EXIT

pops_file=population_list.txt
readarray -t pops < $pops_file

for i in ${pops[@]}
do
	bcftools view --samples-file ~/projects/1000_genomes/populations/subpops/${i}_females -Ov -o $MYTMPDIR/${i}_females_unfiltered_variants.vcf ~/projects/1000_genomes/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf
	bcftools view -m2 -M2 -v snps -Ov -o $MYTMPDIR/${i}_females_biallelic_variants.vcf $MYTMPDIR/${i}_females_unfiltered_variants.vcf
	plink --vcf $MYTMPDIR/${i}_females_biallelic_variants.vcf --out plink_files/${i}_females_biallelic_variants.vcf
	#plink --bfile plink_files/${i}_females_biallelic_variants.vcf --out 
done

date 