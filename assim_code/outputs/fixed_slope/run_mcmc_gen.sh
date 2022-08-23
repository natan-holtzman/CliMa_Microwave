#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --output=julia_test_1.log
#SBATCH -p konings,owners

# load the module
ml julia/1.7.2

# run the Julia application
srun julia make_noisyTB.jl 
