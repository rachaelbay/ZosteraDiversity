#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=GVCF
#################
#a file for job output, you can check job progress
#SBATCH --output=GVCF.out
#################
# a file for errors from the job
#SBATCH --error=GVCF.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 24:00:00
#################
#quality of service; think of it as job priority
#SBATCH -p RM-shared
#################
#number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks=4
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

module load python
module load java
module load samtools
module load GATK

DUPDIR=/ocean/projects/deb200006p/rachbay/Zostera/dups
REFDIR=/ocean/projects/deb200006p/lmschieb/transferred/TestSlurm/Zostera_marina.mainGenome.fasta
GVCFDIR=/ocean/projects/deb200006p/rachbay/Zostera/gvcf

chr=$1
bam=${b}.markeddups.rg.bam

$GATKDIR/gatk --java-options "-Xmx7g" HaplotypeCaller \
-R $REFDIR \
-I $DUPDIR/$bam \
-L $chr \
-O $GVCFDIR/${b}_markeddups.g.vcf.gz \
-ERC GVCF \
