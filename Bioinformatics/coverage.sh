#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=cov
#################
#a file for job output, you can check job progress
#SBATCH --output=COV.out
#################
# a file for errors from the job
#SBATCH --error=COV.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 12:00:00
#################
#quality of service; think of it as job priority
#SBATCH -p RM-shared
#################
#number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks=2
#################
#SBATCH --mem=2G
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=ALL
#################
#who to send email to; please change to your email
#################
#now run normal batch commands
##################
#echo commands to stdout

set -x

mapfile -t FILENAMES < bamlist.txt
b=${FILENAMES[$SLURM_ARRAY_TASK_ID]}

echo ${b}
BAMDIR=/ocean/projects/deb200006p/rachbay/Zostera/dups

../../programs/mosdepth --fast-mode --mapq 20 --by 5000 ../coverage/$b $BAMDIR/${b}_markdup.bam
