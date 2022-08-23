#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --output=stats_hourly_full.log
#SBATCH -p konings,owners

# load the module
ml viz
ml py-matplotlib/3.2.1_py36
ml devel
ml py-pandas/1.0.3_py36 
ml py-scipy/1.4.1_py36

# run the Julia application
#python3 distrib_stats3.py
#python3 GMstat.py
#python3 box_may3.py
#python3 box_hourly.py
#python3 diurnal_avg.py
#python3 box_prior_3run.py
#python3 violin_5run.py
#python3 violin_3x5.py
#python3 plot_vod.py
python3 u_test_hourly.py
#python3 violin_withVOD_3run.py
