#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=DUP
#################
#a file for job output, you can check job progress
#SBATCH --output=DUP.out
#################
# a file for errors from the job
#SBATCH --error=DUP.err
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
#SBATCH -A bio150014p
#################
#SBATCH --mem=3G
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=ALL
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=rachaelbay@gmail.com
#################
#now run normal batch commands
##################
#echo commands to stdout

set -x

module load java
module load picard

BAMDIR=/ocean/projects/deb200006p/rachbay/Zostera/bam
DUPDIR=/ocean/projects/deb200006p/rachbay/Zostera/dups

mapfile -t FILENAMES < samplist.txt
b=${FILENAMES[$SLURM_ARRAY_TASK_ID]}

bam=$BAMDIR/${b}.bam

picard MarkDuplicates \
INPUT=$bam \
OUTPUT=$DUPDIR/dup/${b}_markdup.bam \
METRICS_FILE=$DUPDIR/metrics/${b}.txt \
VALIDATION_STRINGENCY=SILENT \
REMOVE_DUPLICATES=false \
ASSUME_SORTED=true
