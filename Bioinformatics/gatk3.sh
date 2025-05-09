#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=GATK
#################
#a file for job output, you can check job progress
#SBATCH --output=GATK3.out
#################
# a file for errors from the job
#SBATCH --error=GATK3.err
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
module load picard
module load vcftools
module load GATK



REFDIR=/ocean/projects/deb200006p/lmschieb/transferred/TestSlurm/Zostera_marina.mainGenome.fasta
SNPDIR=ocean/projects/deb200006p/rachbay/Zostera/genotypes

ls $VCFDIR/*.vcf.gz | sort -V > vcf.list

gatk GatherVcfs \
-I vcf.list \
-O $SNPDIR/combined_raw.vcf

gatk --java-options "-Xmx8g" \
SelectVariants \
-R $REFDIR \
--select-type-to-include SNP \
-V $SNPDIR/combined_raw.vcf \
-O $SNPDIR/combined_snps.vcf

gatk VariantFiltration \
-R $REFDIR \
-V $SNPDIRcombined_snps.vcf \
--filter-name "FS" \
--filter "FS > 30.0" \
--filter-name "QD" \
--filter "QD < 2.0" \
-O snp/variants/combined_filttags.vcf

gatk SelectVariants \
-R $REFDIR \
-V $SNPDIR/combined_filttags.vcf \
-select 'vc.isNotFiltered()' \
-O $SNPDIR/combined_filtered.vcf

vcftools --vcf $SNPDIRcombined_filtered.vcf \
	--min-alleles 2 \
	--max-alleles 2 \
	--remove-indels --max-missing 0.8 \
	--minGQ 30 \
	--min-meanDP 10 \
	--max-meanDP 60 \
	--maf 0.01 \
	--recode --out snp/variants/combined_vcffilter
