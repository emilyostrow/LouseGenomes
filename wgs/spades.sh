#!/bin/bash
#SBATCH --job-name=LHD1850spades
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
#SBATCH --output=LHD1850spades_%j.log


module load spades

##working with trimmed files (adapters removed)

#spades.py -1 ENO321FSs_S20_R1.fastq.gz -2 ENO321FSs_S20_R2.fastq.gz -o ENO321 -t 12
spades.py -1 LHD1850MSs_S23_R1.fastq.gz -2 LHD1850MSs_S23_R2.fastq.gz -o LHD1850 -t 12