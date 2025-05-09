#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=GATK
#################
#a file for job output, you can check job progress
#SBATCH --output=GATK2.out
#################
# a file for errors from the job
#SBATCH --error=GATK2.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 48:00:00
#################
#quality of service; think of it as job priority
#SBATCH -p RM-shared
#################
#number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks=10
#################
#SBATCH --mem=16G
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=ALL
#################
#now run normal batch commands
##################
#echo commands to stdout

set -x

module load python
module load java
module load GATK

REFDIR=/ocean/projects/deb200006p/lmschieb/transferred/TestSlurm/Zostera_marina.mainGenome.fasta
SNPDIR=/ocean/projects/deb200006p/rachbay/Zostera/snp


chr=$1

gatk --java-options "-Xmx12g" \
GenotypeGVCFs \
-R $REFDIR \
-V gendb://$SNPDIR/chr$chr \
-O $SNPDIR/../genotypes/$chr.vcf.gz
