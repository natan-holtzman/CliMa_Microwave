#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --output=julia_test_1.log
#SBATCH -p konings,owners

#argnames=("o1AMPM" "o6AMPM" "o1AM" "oAll")
#outnames=("o1AMPM/" "o6AMPM/" "o1AM/" "oAll/")
# load the module
ml julia/1.5.1

# run the Julia application
srun julia make_noisyTB.jl 
