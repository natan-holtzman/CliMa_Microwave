#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --output=julia_test_1.%A_%a.log
#SBATCH -p konings,owners
#SBATCH --array=9-14

argnames=("o1AMPM" "o6AMPM" "oAll" "o1AMPM" "o6AMPM" "oAll" "o1AMPM" "o6AMPM" "oAll" "o1and6" "o1and6" "o1and6" "o16offset" "o16offset" "o16offset")
outnames=("o1AMPM_c1/" "o6AMPM_c1/" "oAll_c1/" "o1AMPM_c2/" "o6AMPM_c2/" "oAll_c2/" "o1AMPM_c3/" "o6AMPM_c3/" "oAll_c3/" "o1and6_c1/" "o1and6_c2/" "o1and6_c3/" "o16offset_c1/" "o16offset_c2/" "o16offset_c3/")
in_names=("obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv" "obsTB_witherr_1.csv")
init_names=("o1AMPM_c1/opt_par.csv" "o6AMPM_c1/opt_par.csv" "oAll_c1/opt_par.csv" "o1AMPM_c1/opt_par.csv" "o6AMPM_c1/opt_par.csv" "oAll_c1/opt_par.csv" "o1AMPM_c1/opt_par.csv" "o6AMPM_c1/opt_par.csv" "oAll_c1/opt_par.csv" "o1and6_c1/opt_par.csv" "o1and6_c1/opt_par.csv" "o1and6_c1/opt_par.csv" "o16offset_c1/opt_par.csv" "o16offset_c1/opt_par.csv" "o16offset_c1/opt_par.csv")


# load the module
ml julia/1.7.2

# run the Julia application
srun julia loop_simple.jl ${argnames[$SLURM_ARRAY_TASK_ID]} ${outnames[$SLURM_ARRAY_TASK_ID]} ${in_names[$SLURM_ARRAY_TASK_ID]} ${init_names[$SLURM_ARRAY_TASK_ID]}
