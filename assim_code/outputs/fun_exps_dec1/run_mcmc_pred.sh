#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --output=julia_test_pred.log
#SBATCH -p konings,owners
#SBATCH --array=0-8

#argnames=("o1AMPM" "o6AMPM" "o1AM" "oAll" "o1AMPM" "o6AMPM" "o1AM" "oAll")
outnames=("o1AMPM_c3/" "o6AMPM_c3/" "oAll_c3/" "o1AMPM_c2/" "o6AMPM_c2/" "oAll_c2/" "o1AMPM_c1/" "o6AMPM_c1/" "oAll_c1/")
# load the module
ml julia/1.5.1

# run the Julia application
srun julia predict_coupledET_noK.jl ${outnames[$SLURM_ARRAY_TASK_ID]}
