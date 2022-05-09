#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --output=julia_test_1.%A_%a.log
#SBATCH -p konings,owners
#SBATCH --array=0-2

argnames=("o6AMPM" "o6AMPM" "o6AMPM")
outnames=("o6AMPM_c1/" "o6AMPM_c2/" "o6AMPM_c3/")
in_names=("obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv")
init_names=("o6AMPM_c1/post_par.csv" "o6AMPM_c2/post_par.csv" "o6AMPM_c3/post_par.csv")

# load the module
ml julia/1.7.2

# run the Julia application
srun julia loop_continue.jl ${argnames[$SLURM_ARRAY_TASK_ID]} ${outnames[$SLURM_ARRAY_TASK_ID]} ${in_names[$SLURM_ARRAY_TASK_ID]} ${init_names[$SLURM_ARRAY_TASK_ID]}
