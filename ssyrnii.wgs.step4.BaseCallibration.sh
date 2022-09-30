#!/bin/sh
#
#SBATCH --job-name=Ssyr.wgs              # Job Name
#SBATCH --nodes=1             # 40 nodes
#SBATCH --ntasks-per-node=25               # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/panfs/pfs.local/home/e001o910/../../scratch/bi/e001o910/genomes/DevonPipelineresults/newGATKss	# Set working d$
#SBATCH --mem-per-cpu=5gb            # memory requested
#SBATCH --time=168:00:00
#SBATCH --mail-user=emily.ostrow@ku.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=ssyr-step4_%j.log

#load modules that you will need for executing the pipeline
module load java
module load R
source activate gatk

#save list of input file names
samples="ENO297FSs_S1
ENO298FSs_S2
ENO299FSs_S3
ENO299MSs_S4
ENO301FSs_S5
ENO302FSs_S6
ENO303FSs_S7
ENO304FSs_S8
ENO305FSs_S9
ENO306FSs_S10
ENO307FSs_S11
ENO308FSs_S12
ENO310FSs_S13
ENO311FSs_S14
ENO312FSs_S15
ENO313FSs_S25
ENO320FSs_S19
ENO321FSs_S20
ENO322MSs_S24
ENO327FSs_S29
ENO330FSs_S32
ENO331FSs_S33
ENO332FSs_S35
ENO333FSs_S36
LHD1849FSs_S22
LHD1850MSs_S23"

##path to your reference
ref="../syrnii_eno30710kb_spades_scaffolds.fasta"
##name of combined vcf files for output
comb="ssyrnii_10kb"


###############################################
###  Base quality recalibration ###
###############################################

for i in $samples
do

/panfs/pfs.local/work/bi/e001o910/programs/gatk-4.2.6.1/gatk --java-options "-Xmx8g" BaseRecalibrator \
-I ${i}.merged.sorted.dups.bam \
-R $ref \
--known-sites $comb.bqsr_indels.vcf.gz \
--known-sites $comb.bqsr_snps.vcf.gz \
-O ${i}.table

## ApplyBQSR
/panfs/pfs.local/work/bi/e001o910/programs/gatk-4.2.6.1/gatk --java-options "-Xmx8g" ApplyBQSR \
-I ${i}.merged.sorted.dups.bam \
-R $ref  \
--bqsr-recal-file ${i}.table \
-O ${i}.recal.bam

## Index the recalibrated bam files
##don't think this is necessary
#java -Xmx10g -jar /home/d669d153/work/picard/build/libs/picard.jar BuildBamIndex \
#I=${i}.recal.bam VALIDATION_STRINGENCY=LENIENT

done










