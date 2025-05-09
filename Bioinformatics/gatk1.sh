#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=GenomicsDB
#################
#a file for job output, you can check job progress
#SBATCH --output=GATK1.out
#################
# a file for errors from the job
#SBATCH --error=GATK1.err
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

module load python
module load java
module load samtools
module load GATK

GVCFDIR=/ocean/projects/deb200006p/rachbay/Zostera/gvcf
REFDIR=/ocean/projects/deb200006p/lmschieb/transferred/TestSlurm/Zostera_marina.mainGenome.fasta
SNPDIR=/ocean/projects/deb200006p/rachbay/Zostera/snp

chr=$1

ls $GVCFDIR/*.vcf.gz > gvcf.list


gatk --java-options "-Xmx16g" \
GenomicsDBImport \
-R $REFDIR \
-V gvcf.list \
-L $chr \
--genomicsdb-workspace-path ../snp/chr$chr

##To add to existing database
#$GATKDIR/gatk --java-options "-Xmx16g" \
#GenomicsDBImport \
#-V chr$chr.list \
#--genomicsdb-update-workspace-path ../bqsr/snp.round2/chr$chr
