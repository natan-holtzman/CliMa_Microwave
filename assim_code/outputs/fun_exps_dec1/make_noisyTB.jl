include("../../../get_home_dir.jl")

using DataFrames
using CSV
#using Plots
using Random
using StatsBase



include(string(PROJECT_HOME,"/assim_code/mironov.jl"));
include(string(PROJECT_HOME,"/assim_code/tau_omega_funs.jl"));
include(string(PROJECT_HOME,"/assim_code/time_averaging.jl"))

df_raw = CSV.read(string(PROJECT_HOME,"/data/moflux_land_data_skipyear_hourly2.csv"), DataFrame);
df_raw[!,"RAIN"] *= 2;

N = 24*365
istart = 24*365*2 + 1; 
soil0 = 0.30;

include(string(PROJECT_HOME,"/simulation_code/sim_vary_new_stomata_med.jl"));

deltaT = FT(60*60);
alpha = FT(50);
nsoil = FT(1.2);

function run_sim_2(vcmax_par::FT, k_frac::FT, k_plant::FT,  z_soil::FT, weibB::FT, vol_factor::FT, g1::FT)
	return run_sim_vary(vcmax_par, k_frac, weibB, FT(2), k_plant, FT(1e-4), z_soil, istart, N, soil0, vol_factor, 1e-5, 3, df_raw, g1, deltaT, alpha, nsoil);
end

sim_res1 = run_sim_2(FT(31),FT(0.5), FT(16), FT(2000),FT(5),FT(0.33),FT(506));

#vcmax_par, k_frac, k_plant, k_soil, z_soil, istart, N, smc0, slope_index,use_flux, weibB, weibC

noise_std = 1.3;
noise_var = noise_std^2

tsoil = sim_res1[1].T_SOIL .+ 273.15;
tair = sim_res1[1].T_AIR .+ 273.15;
tcan = (sim_res1[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
laiM = sim_res1[1].LAI_modis;
canpot_true = mean(sim_res1[3],dims=2);

using Distributions
using LabelledArrays
using LinearAlgebra

vodA = 0.067;
vodB = 0.82;
vodC = 0.051;

trueTB = get_TB_2(sim_res1[2][:,1], tsoil, tcan, canpot_true, laiM, vodA, vodB, vodC);
obsH = trueTB[1] + noise_std*randn(Float64, length(trueTB[1]));
obsV = trueTB[2] + noise_std*randn(Float64, length(trueTB[1]));

obsDF = DataFrame(hpol=obsH, vpol=obsV);
CSV.write("obsTB_witherr_1.csv",obsDF);

obsH2 = trueTB[1] + noise_std*randn(Float64, length(trueTB[1]));
obsV2 = trueTB[2] + noise_std*randn(Float64, length(trueTB[1]));

obsDF = DataFrame(hpol=obsH2, vpol=obsV2);
CSV.write("obsTB_witherr_2.csv",obsDF);

obsH3 = trueTB[1] + noise_std*randn(Float64, length(trueTB[1]));
obsV3 = trueTB[2] + noise_std*randn(Float64, length(trueTB[1]));

obsDF = DataFrame(hpol=obsH3, vpol=obsV3);
CSV.write("obsTB_witherr_3.csv",obsDF);




