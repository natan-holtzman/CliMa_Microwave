#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --output=julia_test_pred2.%A_%a..log
#SBATCH -p konings,owners
#SBATCH --array=9-14

#argnames=("o1AMPM" "o6AMPM" "o1AM" "oAll" "o1AMPM" "o6AMPM" "o1AM" "oAll")
outnames=("o1AMPM_c3/" "o6AMPM_c3/" "oAll_c3/" "o1AMPM_c2/" "o6AMPM_c2/" "oAll_c2/" "o1AMPM_c1/" "o6AMPM_c1/" "oAll_c1/" "o1and6_c1/" "o1and6_c2/" "o1and6_c3/" "o16offset_c1/" "o16offset_c2/" "o16offset_c3/")

# load the module
ml julia/1.7.2

# run the Julia application
srun julia predict_5000_v2.jl ${outnames[$SLURM_ARRAY_TASK_ID]}
