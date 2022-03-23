#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --output=julia_test_1.%A_%a.log
#SBATCH -p konings,owners
#SBATCH --array=0-2

argnames=("o1AMPM" "o6AMPM" "oAll")
outnames=("o1AMPM_c1/" "o6AMPM_c1/" "oAll_c1/")
in_names=("obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv")
# load the module
ml julia/1.5.1

# run the Julia application
srun julia loop_mode.jl ${argnames[$SLURM_ARRAY_TASK_ID]} ${outnames[$SLURM_ARRAY_TASK_ID]} ${in_names[$SLURM_ARRAY_TASK_ID]}