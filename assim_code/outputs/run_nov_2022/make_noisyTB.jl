include("../../../get_home_dir.jl")
using Pkg
Pkg.activate(string(PROJECT_HOME,"/feb_j171"));

using DataFrames
using CSV
#using Plots
using Random
using StatsBase



include(string(PROJECT_HOME,"/assim_code/mironov.jl"));
include(string(PROJECT_HOME,"/assim_code/tau_omega_funs.jl"));
include(string(PROJECT_HOME,"/assim_code/time_averaging.jl"))

df_raw = CSV.read(string(PROJECT_HOME,"/data/moflux_fluxnet_data.csv"), DataFrame);
df_raw[!,"RAIN"] *= 2;


include(string(PROJECT_HOME,"/simulation_code/full_model_newP63.jl"));


N = 24*365*1
istart = 24*365*2 + 1;
soil0 = FT(0.4);

deltaT = FT(60*60);

soil_curve_1 = VanGenuchten{FT}(stype = "Ozark",α = FT(10.9),   n = 1.5, Θs = 0.55,   Θr = 0.067);
smc1 = soil_swc(soil_curve_1,FT(-1));

function run_sim_0(vcmax_par::FT, weibB_stom::FT, k_plant::FT,
    z_soil::FT, xylem_extraB::FT, vol_factor::FT, g1::FT, k_soil::FT,
    slope_runoff::FT, exp_root_dist::FT,nsoil::FT,soilfac::FT
    )
weibC_plant = FT(2);
soil_curve_2 =  VanGenuchten{FT}(stype = "Ozark",α = FT(10.9),   n = nsoil, Θs = 0.55,   Θr = 0.067);
pot1b = soil_p_25_swc(soil_curve_2,smc1);
alpha2 = FT(10.9*soilfac)*abs(pot1b);
return run_sim_varyB(vcmax_par, weibB_stom/soilfac, (weibB_stom+xylem_extraB)/soilfac, weibC_plant, k_plant*soilfac,
  k_soil*soilfac, z_soil, istart, N, soil0, vol_factor*soilfac, 1e-5, 3, df_raw, g1, deltaT, alpha2, nsoil,0,0,
  slope_runoff, exp_root_dist,FT(1/20),FT(1/20));
end



pars0 = convert(Array{FT}, [90, 3, 10, 2000, 1, 1.0, 300, 0.4e-6,
       0.4,2,1.5,1]);





sim_res1 = run_sim_0(pars0...);

#vcmax_par, k_frac, k_plant, k_soil, z_soil, istart, N, smc0, slope_index,use_flux, weibB, weibC

noise_std = 1.3;
noise_var = noise_std^2

tsoil = sim_res1[1].T_SOIL .+ 273.15;
tair = sim_res1[1].T_AIR .+ 273.15;
tcan = (sim_res1[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
laiM = sim_res1[1].LAI_modis;
canpot_true = mean(sim_res1[6][:,1:3],dims=2);

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




