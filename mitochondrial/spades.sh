#!/bin/bash
#SBATCH --job-name=batchspades
#SBATCH --partition=bi
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=168:00:00
#SBATCH --mem=60G
#SBATCH --mail-user=emily.ostrow@ku.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o
#SBATCH -chdir=/panfs/pfs.local/home/e001o910/scratch/
#SBATCH --output=batchplasmidspades_%j.log


module load spades

##working with trimmed files (adapters removed) above 35x coverage

#spades.py -1 ENO321FSs_S20_R1.fastq.gz -2 ENO321FSs_S20_R2.fastq.gz -o ENO321 -t 12 --plasmid
#spades.py -1 LHD1850MSs_S23_R1.fastq.gz -2 LHD1850MSs_S23_R2.fastq.gz -o LHD1850 -t 12 --plasmid
#spades.py -1 ENO307FSs_S11_R1.fastq.gz -2 ENO307FSs_S11_R2.fastq.gz -o ENO307 -t 12 --plasmid
#spades.py -1 ENO331FSs_S33_R1.fastq.gz -2 ENO331FSs_S33_R2.fastq.gz -o ENO331 -t 12 --plasmid
#spades.py -1 ENO327FSs_S29_R1.fastq.gz -2 ENO327FSs_S29_R2.fastq.gz -o ENO327 -t 12 --plasmid
#spades.py -1 ENO330FSs_S32_R1.fastq.gz -2 ENO330FSs_S32_R2.fastq.gz -o ENO330 -t 12 --plasmid
#spades.py -1 ENO308FSs_S12_R1.fastq.gz -2 ENO308FSs_S12_R2.fastq.gz -o ENO308 -t 12 --plasmid
#spades.py -1 ENO313FSs_S25_R1.fastq.gz -2 ENO313FSs_S25_R2.fastq.gz -o ENO313 -t 12 --plasmid
#spades.py -1 ENO320FSs_S19_R1.fastq.gz -2 ENO320FSs_S19_R2.fastq.gz -o ENO320 -t 12 --plasmid


sosamples="ENO316FSo_S17
ENO317MSo_S18
ENO321FSo_S21
ENO322FSo_S16
ENO325MSo_S26
ENO326FSo_S27
ENO327FSo_S28
ENO328FSo_S30
ENO329FSo_S31
ENO331FSo_S34"
sssamples="ENO297FSs_S1
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



for i in $sosamples
do
seqtk sample ../soTrimmed/${i}_R1.fastp.fastq.gz 0.1 > ${i}_10p_R1.fq
seqtk sample ../soTrimmed/${i}_R2.fastp.fastq.gz 0.1 > ${i}_10p_R2.fq
done
for i in $sssamples
do
seqtk sample ../ssTrimmed/${i}_R1.fastp.fastq.gz 0.1 > ${i}_10p_R1.fq
seqtk sample ../ssTrimmed/${i}_R2.fastp.fastq.gz 0.1 > ${i}_10p_R2.fq
done

dirPath="/panfs/pfs.local/home/e001o910/../../scratch/bi/e001o910/genomes/spades"

for i in $sosamples
do 
echo '#!/bin/sh' | tee ${i}_plasmid_spades.sh
echo "#" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --job-name=${i}_plasmid_spades              # Job Name" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --nodes=1             # 40 nodes" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --ntasks-per-node=12               #1 CPU allocation" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --partition=bi            # Name of the Slurm partition used" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --chdir=$dirPath	# Set working d$" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --mem-per-cpu=5gb            # memory requested" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --time=168:00:00" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --mail-user=emily.ostrow@ku.edu" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --mail-type=END,FAIL" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --output=plasmid_spades${i}_%j.log" | tee -a ${i}_plasmid_spades.sh
echo "module load spades" | tee -a ${i}_plasmid_spades.sh
echo "spades.py -1 ${i}_10p_R1.fq -2 ${i}_10p_R2.fq -o ${i}_10p -t 12 --plasmid" | tee -a ${i}_plasmid_spades.sh
sbatch ${i}_plasmid_spades.sh
done

for i in $sssamples
do 
echo '#!/bin/sh' | tee ${i}_plasmid_spades.sh
echo "#" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --job-name=${i}_plasmid_spades              # Job Name" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --nodes=1             # 40 nodes" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --ntasks-per-node=12               #1 CPU allocation" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --partition=bi            # Name of the Slurm partition used" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --chdir=$dirPath	# Set working d$" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --mem-per-cpu=5gb            # memory requested" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --time=168:00:00" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --mail-user=emily.ostrow@ku.edu" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --mail-type=END,FAIL" | tee -a ${i}_plasmid_spades.sh
echo "#SBATCH --output=plasmid_spades${i}_%j.log" | tee -a ${i}_plasmid_spades.sh
echo "module load spades" | tee -a ${i}_plasmid_spades.sh
echo "spades.py -1 ${i}_10p_R1.fq -2 ${i}_10p_R2.fq -o ${i}_10p -t 12 --plasmid" | tee -a ${i}_plasmid_spades.sh
sbatch ${i}_plasmid_spades.sh
done