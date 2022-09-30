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
#SBATCH --output=ssyr-step5_%j.log

#load modules that you will need for executing the pipeline
module load java
module load R
module load vcftools
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

# following: https://www.melbournebioinformatics.org.au/tutorials/tutorials/variant_calling_gatk1/variant_calling_gatk1/

##you need to create a gatk environment, so if you are trying to use gatk and it doesn't work, start there.
#cd path/to/gatk 
#conda env create -f gatkcondaenv.yml
##path to your reference
ref="../syrnii_eno30710kb_spades_scaffolds.fasta"
##name of combined vcf files for output
comb="ssyrnii_10kb"

###############################################
### VARIANT CALLING: GATK HaplotypeCaller ###
###############################################

for i in $samples
do
/panfs/pfs.local/work/bi/e001o910/programs/gatk-4.2.6.1/gatk --java-options "-Xmx7g" HaplotypeCaller \
    -I ${i}.recal.bam \
    -R $ref \
    -O ${i}.recal.g.vcf.gz
done

###############################################
### CombineGVCFs ###
###############################################

/panfs/pfs.local/work/bi/e001o910/programs/gatk-4.2.6.1/gatk --java-options "-Xmx7g" CombineGVCFs \
    -R $ref \
	-V ENO297FSs_S1.recal.g.vcf.gz \
	-V ENO298FSs_S2.recal.g.vcf.gz \
	-V ENO299FSs_S3.recal.g.vcf.gz \
	-V ENO299MSs_S4.recal.g.vcf.gz \
	-V ENO301FSs_S5.recal.g.vcf.gz \
	-V ENO302FSs_S6.recal.g.vcf.gz \
	-V ENO303FSs_S7.recal.g.vcf.gz \
	-V ENO304FSs_S8.recal.g.vcf.gz \
	-V ENO305FSs_S9.recal.g.vcf.gz \
	-V ENO306FSs_S10.recal.g.vcf.gz \
	-V ENO307FSs_S11.recal.g.vcf.gz \
	-V ENO308FSs_S12.recal.g.vcf.gz \
	-V ENO310FSs_S13.recal.g.vcf.gz \
	-V ENO311FSs_S14.recal.g.vcf.gz \
	-V ENO312FSs_S15.recal.g.vcf.gz \
	-V ENO313FSs_S25.recal.g.vcf.gz \
	-V ENO320FSs_S19.recal.g.vcf.gz \
	-V ENO321FSs_S20.recal.g.vcf.gz \
	-V ENO322MSs_S24.recal.g.vcf.gz \
	-V ENO327FSs_S29.recal.g.vcf.gz \
	-V ENO330FSs_S32.recal.g.vcf.gz \
	-V ENO331FSs_S33.recal.g.vcf.gz \
	-V ENO332FSs_S35.recal.g.vcf.gz \
	-V ENO333FSs_S36.recal.g.vcf.gz \
	-V LHD1849FSs_S22.recal.g.vcf.gz \
	-V LHD1850MSs_S23.recal.g.vcf.gz \
    -O $comb.recal.g.vcf.gz


###############################################
### GenotypeGVCFs ###
###############################################

/panfs/pfs.local/work/bi/e001o910/programs/gatk-4.2.6.1/gatk --java-options "-Xmx7g" GenotypeGVCFs \
    -R $ref \
    -V $comb.recal.g.vcf.gz \
    -O $comb.recal.vcf.gz


###############################################
### Export SNPs ###
###############################################

    
/panfs/pfs.local/work/bi/e001o910/programs/gatk-4.2.6.1/gatk SelectVariants \
    -V $comb.recal.vcf.gz \
    -select-type SNP \
    -O $comb.recal.snps.vcf.gz

    ##options are BOTH, INDEL, and SNP
###############################################
### FinalFiltering ###
###############################################

/panfs/pfs.local/work/bi/e001o910/programs/gatk-4.2.6.1/gatk VariantFiltration \
    -V $comb.recal.snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "SOR > 4.0" --filter-name "SOR4" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O $comb.recal.snps_filtered.vcf.gz


# remove filtered sites and add additional filters with VCFtools
#VCFtools v.0.1.15 already in path, no need to specify full location
# maximum missing = 20%, MAF = 5%, minimum depth = 3X, maximum mean depth = 100X, biallelic
vcftools --vcf $comb.snps.filtered.vcf --remove-filtered-all --max-missing 0.80 --maf 0.05 --min-meanDP 3 --max-meanDP 100 --min-alleles 2 --max-alleles 2 --recode --out $comb.vcftools.filtered


###############################################
### VariantsToTable - not sure why you would need this ###
###############################################

#/panfs/pfs.local/work/bi/e001o910/programs/gatk-4.2.6.1/gatk VariantsToTable \
#    -R $ref \
#    -V $comb.vqsr.varfilter.pass.vcf.gz \
#    -F CHROM -F POS -F FILTER -F TYPE -GF AD -GF DP \
#    --show-filtered \
#    -O $comb.vqsr.varfilter.pass.tsv


