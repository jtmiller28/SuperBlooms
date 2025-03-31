#!/bin/bash

#SBATCH --job-name=phenovision-desert-filtering     # Job name
#SBATCH --mail-type=ALL                  # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jtmiller@ucsb.edu    # Where to send mail
#SBATCH --output=%j.log                  # Standard output and error log
#SBATCH --nodes=1                        # Run all processes on a single node
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem=100gb                       # Job memory request
#SBATCH --time=00-48:00:00               # Time limit days-hrs:min:sec
#SBATCH --account=guralnick
#SBATCH --qos=guralnick-b
pwd; hostname; date

#load modules

module load R/4.3

#do some (or alot) of coding
Rscript --vanilla /blue/guralnick/millerjared/SuperBlooms/finished-code/r/001-data-subsetting.R

## example for running: sbatch /blue/guralnick/millerjared/SuperBlooms/finished-code/bash/001-data-subsetting-scheduler.sh
