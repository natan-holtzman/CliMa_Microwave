using DataFrames
using CSV
using PyPlot
using Random
using StatsBase
using Statistics

using Distributions
using LabelledArrays
using LinearAlgebra

include("../get_home_dir.jl")

pygui(true)

df_raw = CSV.read(string(PROJECT_HOME,"/data/moflux_land_data_skipyear_hourly2.csv"), DataFrame);
df_raw[!,"RAIN"] *= 2; #this is because rain was originally in mm/half hour time step

include(string(PROJECT_HOME,"/simulation_code/full_model.jl"));
include("time_averaging.jl")
include("mironov.jl")
include("tau_omega_funs.jl")

N = 24*365*1
istart = 24*365*2 + 1; 
deltaT = FT(60*60);
soil0 = 0.375;
alpha = FT(10);
nsoil = FT(1.5);

function run_sim_3layer(vcmax_par::FT, k_frac::FT, k_plant::FT,
      z_soil::FT, weibB::FT, vol_factor::FT, g1::FT, k_soil::FT,
      smc_runoff::FT, exp_root_dist::FT,
      canopy_pvslope::FT, trunk_pvslope::FT)
	return run_sim_varyB(vcmax_par, k_frac, weibB, FT(2), k_plant, 
    k_soil, z_soil, istart, N, soil0, vol_factor, 1e-5, 3, df_raw, g1, deltaT, alpha, nsoil,0,0,
    smc_runoff, exp_root_dist,canopy_pvslope, trunk_pvslope);
end

pars0 = convert(Array{FT}, [40, 0.5, 12, 2000, 8, 1.0, 250,0.4e-6,
         0.35,2,1/20,1/20])
sim_res1 = run_sim_3layer(pars0...);
cm1 = mean(sim_res1[6][:,1:3],dims=2);
pdLWP = Float64.(replace(sim_res1[1].LWP_predawn, missing => NaN));


tsoil = sim_res1[1].T_SOIL .+ 273.15;
tair = sim_res1[1].T_AIR .+ 273.15;
tcan = (sim_res1[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
laiM = sim_res1[1].LAI_modis;

vodA = 0.067;
vodB = 0.82;
vodC = 0.051;

tb1 = get_TB_2(sim_res1[2][:,1], tsoil, tcan, cm1, laiM, vodA, vodB, vodC);