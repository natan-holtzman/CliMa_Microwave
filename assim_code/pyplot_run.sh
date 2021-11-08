#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --output=julia_test_1.log
#SBATCH -p konings,owners

# load the module
ml viz
ml py-matplotlib/3.2.1_py36
ml devel
ml py-pandas/1.0.3_py36 
# run the Julia application
python3 pyplot_post.py
