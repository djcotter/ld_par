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

~/tools/plink2/bin/plink --bfile 

date 