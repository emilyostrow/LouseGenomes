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
#SBATCH --output=ssyr-step1_%j.log

#load modules that you will need for executing the pipeline
module load java
module load R

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

################################
### TRIM AND REMOVE ADAPTORS ###
################################
#AdapterRemoval v.2.3.2 is already installed via conda, so can be called without specifying the full path
#run in a loop over all samples
for i in $samples
do 
AdapterRemoval --file1 ../../soraw.fastq/${i}_R1_001.fastq.gz --file2 ../../soraw.fastq/${i}_R2_001.fastq.gz --basename ${i} --trimns --trimqualities --collapse --threads 25
done
##used fastp to remove adapters for atram previously

##################################################
### ALIGN TO REFERENCE GENOME USING Bowtie2 ###
##################################################
#Bowtie 2 v.2.3.3.1 is already installed and in path
#bowtie settings info:
#-p = number of threads
#--very-sensitive-local = following presets: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
#-x = indexed reference genome input
#--rg-id = set a read group ID
#--rg = add <text> ("lab:value") to @RG line of SAM header.
#-1 = mate pair 1
#-2 = mate pair 2
#-U = input overlapping and unpaired sequences for the given sample

# index the reference genome
bowtie2-build -f ../syrnii_eno30710kb_spades_scaffolds.fasta syrnii_ref_bowtie

#run in a loop over all samples
for i in $samples
do 
bowtie2 -p 25 --very-sensitive-local -x oculatus_ref_bowtie --rg-id ${i} --rg SM:${i} -1 ${i}.pair1.truncated -2 ${i}.pair2.truncated -U ${i}.collapsed,${i}.collapsed.truncated,${i}.singleton.truncated -S ${i}.sam
/panfs/pfs.local/work/bi/bin/samtools-1.3.1/bin/samtools view --threads 25 -bS ${i}.sam > ${i}.bam
/panfs/pfs.local/work/bi/bin/samtools-1.3.1/bin/samtools sort --threads 25 ${i}.bam -o ${i}_sorted.bam
rm ${i}.sam
rm ${i}.bam
done

###############################################
### ASSESS ALIGNMENT TO REFERENCE: QUALIMAP ###
###############################################
#run qualimap on each sample to ensure decent mapping
for i in $samples
do 
/home/d669d153/work/qualimap_v2.2.1/qualimap bamqc -bam ${i}_sorted.bam -outfile ${i}_sorted.pdf 
done


### This is the end of the mapping pipeline. Samples will randomly throw errors during the sam > bam conversion, resulting in nearly the entire alignment getting dumped.
### At this point make sure that each _sorted.bam file is of reasonable size, if a sample has a tiny _sorted.bam file, it's usually a good indication that an error occurred.
### Manually inspect file sizes by running 'ls -lh *_sorted.bam' in the target directory before moving on to variant calling.
### If samples failed, simply create a copy of this pipeline in the same directory, remove all the samples that worked from the $samples variable created above, and rerun, which should solve the issue.
