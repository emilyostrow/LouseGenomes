#!/bin/sh
#
#SBATCH --job-name=ENO297FSs_S1_plasmid_spades              # Job Name
#SBATCH --nodes=1             # 40 nodes
#SBATCH --ntasks-per-node=12               #1 CPU allocation
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/panfs/pfs.local/home/e001o910/../../scratch/bi/e001o910/genomes/spades	# Set working d$
#SBATCH --mem-per-cpu=5gb            # memory requested
#SBATCH --time=168:00:00
#SBATCH --mail-user=emily.ostrow@ku.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=plasmid_spadesENO297FSs_S1_%j.log
module load spades
spades.py -1 ENO297FSs_S1_10p_R1.fq -2 ENO297FSs_S1_10p_R2.fq -o ENO297FSs_S1_10p -t 12 --plasmid
