#!/bin/bash
#SBATCH -p smallcpu
#SBATCH -A nhr-tpchange
#SBATCH --time=10:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=60G
#SBATCH --job-name=write_WCB_files
#SBATCH --output=slurmfiles/%x_%j.out
#SBATCH --error=slurmfiles/%x_%j.err

#load CDO module
module load data/CDO/1.9.8-intel-2019b

#change mod to executable and execute (this takes hours)
chmod +x MERGE_TRAJ_FILES.sh
./MERGE_TRAJ_FILES.sh

#THIS ONLY TAKES 30 MINUTES
cd /path/to/julia/code/ #ADD PATH HERE
srun julia --project=. run_WCB_selection.jl
