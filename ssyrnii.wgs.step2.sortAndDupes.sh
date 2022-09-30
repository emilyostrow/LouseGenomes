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
#SBATCH --output=ssyr-step2_%j.log

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

######################
#### Variant Refinement ###
######################
##may need to add -Xmx flag (e.g java -Xmx48g -jar) if you run out of memory on the cluster during these steps
##this code points to the executable for picard v.2.26.11
#
## prep reference genome for GATK (get .fai and .dict files)
java -jar /home/d669d153/work/picard/build/libs/picard.jar CreateSequenceDictionary R= $ref O= ../syrnii_eno30710kb_spades_scaffolds.ref.dict
/panfs/pfs.local/work/bi/bin/samtools-1.3.1/bin/samtools faidx $ref
#
## add read group information
for i in $samples
do 
java -jar /home/d669d153/work/picard/build/libs/picard.jar AddOrReplaceReadGroups INPUT=${i}_sorted.bam OUTPUT=${i}_sortedRG.bam RGID=1 RGLB= library1 RGPL=illumina RGPU=R1 RGSM=${i}_sorted
done
#
## mark duplicates
for i in $samples
do 
java -jar /home/d669d153/work/picard/build/libs/picard.jar MarkDuplicates INPUT=${i}_sortedRG.bam OUTPUT=${i}.sorted.marked.bam METRICS_FILE=${i}_sorted.metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
done

