#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=map
#################
#a file for job output, you can check job progress
#SBATCH --output=MAP.out
#################
# a file for errors from the job
#SBATCH --error=MAP.err
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
#SBATCH --ntasks=8
#################
#SBATCH --mem=15G
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

module load BWA/0.7.3a
module load samtools/1.11.0

REF=/ocean/projects/deb200006p/lmschieb/transferred/TestSlurm/Zostera_marina.mainGenome.fasta

mapfile -t FILENAMES < samps.txt
b=${FILENAMES[$SLURM_ARRAY_TASK_ID]}

f1=$b\_1.fq.gz
f2=$b\_2.fq.gz
ID=$(echo "$b" | cut -f8 -d'/')

bwa mem -M -t 8 -R "@RG\tID:$ID\tSM:$ID\tPL:Illumina" $REF ${f1} ${f2} | samtools view -bhS - | \
	samtools sort -o ./bam/$ID.bam

samtools index ./bam/$ID.bam
