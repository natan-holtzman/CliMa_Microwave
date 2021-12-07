#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --output=julia_test_rz.log
#SBATCH -p konings,owners

# load the module
ml viz
ml py-matplotlib/3.2.1_py36
ml devel
ml py-pandas/1.0.3_py36 
# run the Julia application
#python3 range_diurnal.py
#python3 stats_range.py
#python3 soilpot_range.py
#python3 meanstats_drought2.py
#python3 meanstats_all.py
#python3 summerstats.py
#python3 meanstats_2012.py 
python3 plot3_10000v2.py
#python3 plot_vod.py
