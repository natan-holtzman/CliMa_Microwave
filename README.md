# CliMa_Microwave
This repository contains code related to the paper "Constraining plant hydraulics with microwave radiometry in a land surface model: Impacts of temporal resolution" by Holtzman et al., submitted to Water Resources Research in 2023.

The "simulation_code" folder contains code for running the modified version of the CliMA Land model that is described in the paper.

The "assim_code" folder contains code for model-data-fusion, as well as analysis and plotting results.

The "data" folder contains input data for the model and scripts that were used to process that data.

For an example of the modified CliMA Land model that was presented in the paper, you can run assim_code/plot_fig2.jl or assim_code/plot_fig3_4.jl

To generate the full results, carry out the following command-line workflow within the assim_code/outputs/run_nov_2022 folder:

#1. Generate synthetic brightness temperature observations from the "true" model

  julia make_noisyTB.jl

#2. Find a point estimate of model parameters that approximately matches the HOURLY observations. 
#This will be used as an initial guess in the MCMC algorithm for all observation scenarios.

  julia loop_mode.jl "oALL" "opt_par.csv" "obsTB_witherr_1.csv"

#3. Run the model-data fusion algorithm for the various observing scenarios
#Here this is shown as a loop for clarity, but in practice it should be implemented as a batch job in parallel,
#as in the example Slurm script run_mcmc_shareTB.sh

  argnames=("o1AMPM" "o6AMPM" "oAll" "o1and6" "o1AMPM" "o6AMPM" "oAll" "o1and6" "o1AMPM" "o6AMPM" "oAll" "o1and6")
  dirnames=("o1AMPM_c1/" "o6AMPM_c1/" "oAll_c1/" "o1and6_c1/" "o1AMPM_c2/" "o6AMPM_c2/" "oAll_c2/" "o1and6_c2/" "o1AMPM_c3/" "o6AMPM_c3/" "oAll_c3/" "o1and6_c3")
  mkdir dirnames
  for i in {1..12}; do
    julia loop_simple.jl argnames[$i] dirnames[$i] "obsTB_witherr_1.csv" "opt_par.csv"
  done

#4. Make predictions over 13 years from each observation scenario, using samples from the posterior parameter distributions

  for i in {1..12}; do
    julia predict_5000_v2.jl dirnames[$i]
  done

#4. Run the analysis and plotting code

  julia vod_quantiles_save.jl
  julia plot_fig5.jl
  python3 save_pars.py
  python3 do_stats_all.jl
  python3 plot_fig6.py
  python3 plot_fig7_8.py
  python3 plot_fig9.py

