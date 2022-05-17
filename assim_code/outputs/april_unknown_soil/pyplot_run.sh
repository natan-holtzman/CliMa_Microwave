#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --output=grstat_may15.log
#SBATCH -p konings,owners

# load the module
ml viz
ml py-matplotlib/3.2.1_py36
ml devel
ml py-pandas/1.0.3_py36 
# run the Julia application
#python3 distrib_stats2.py
#python3 violin_plot10k.py
python3 GMstat.py
#python3 violin_12par.py
#python3 plot_vod.py
