#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --output=julia_test_3b.log
#SBATCH -p konings,owners

# load the module
ml julia/1.5.1

# run the Julia application
srun julia assim_realdata_8par.jl 
