#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --output=julia_test_1.log
#SBATCH -p konings,owners
#SBATCH --array=0-3

argnames=("o1AMPM" "o6AMPM" "o1AM" "oAll")
outnames=("o1AMPM/" "o6AMPM/" "o1AM/" "oAll/")
# load the module
ml julia/1.5.1

# run the Julia application
srun julia loop_coupledET.jl ${argnames[$SLURM_ARRAY_TASK_ID]} ${outnames[$SLURM_ARRAY_TASK_ID]}
