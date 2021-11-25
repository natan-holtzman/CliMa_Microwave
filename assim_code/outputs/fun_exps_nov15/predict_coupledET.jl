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

N = 24*365*6
istart = 24*365*0 + 1;
soil0 = 0.40;

include(string(PROJECT_HOME,"/simulation_code/sim_vary_new_stomata_med.jl"));

deltaT = FT(60*60);
alpha = FT(1.368);
nsoil = FT(2.6257);

function run_sim_2(vcmax_par::FT, k_frac::FT, k_plant::FT, k_soil::FT, z_soil::FT, weibB::FT, weibC::FT, vol_factor::FT, g1::FT)
        return run_sim_vary(vcmax_par, k_frac, weibB, weibC, k_plant, k_soil, z_soil, istart, N, soil0, vol_factor, 1e-5, 3, df_raw, g1, deltaT, alpha, nsoil);
end

sim_res1 = run_sim_2(FT(31),FT(0.5), FT(2),FT(0.4e-5), FT(800),FT(5),FT(2),FT(1),FT((16-0.0*9.3)*sqrt(1000)));


noise_std = 1.3;
noise_var = noise_std^2

tsoil = sim_res1[1].T_SOIL .+ 273.15;
tair = sim_res1[1].T_AIR .+ 273.15;
tcan = (sim_res1[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
laiM = sim_res1[1].LAI_modis;
canpot_true = mean(sim_res1[3],dims=2);
trueTB = get_TB_2(sim_res1[2][:,1], tsoil, tcan, canpot_true, laiM, 0.1, 1, 0.05);
obsH = trueTB[1] + noise_std*randn(Float64, length(trueTB[1]));
obsV = trueTB[2] + noise_std*randn(Float64, length(trueTB[1]));

function write_post(mydir)

par_post = convert(Array,CSV.read(string(mydir,"post_par.csv"),DataFrame))[2001:25:end,1:9];
outdir = mydir

nday = convert(Int, N/24);
npost = size(par_post)[1];

post_RZ = zeros(nday,npost+1);
#post_Surf = zeros(nday,npost+1);
post_ET = zeros(nday,npost+1);
post_LWP = zeros(nday*24,npost+1);
post_branch = zeros(nday*24,npost+1);
post_trunk = zeros(nday*24,npost+1);

for j in 1:npost
	println(j)
	a = par_post[j,:]
try
	sim_res0 = run_sim_2(convert(Array{FT}, exp.(a))...); 

	post_RZ[:,j] = get_daily(sim_res0[2][:,end],24);
	post_ET[:,j] = get_daily(sim_res0[1].ETmod,24);
	post_LWP[:,j] = mean(sim_res0[3],dims=2)[:,1];
	post_branch[:,j] = mean(sim_res0[4],dims=2)[:,1];
	post_trunk[:,j] = sim_res0[5];
catch err
	println("Error")
end
end

	#place the results of the original true run in the last column
	post_RZ[:,end] = get_daily(sim_res1[2][:,end],24);
	post_ET[:,end] = get_daily(sim_res1[1].ETmod,24);
	post_LWP[:,end] = mean(sim_res1[3],dims=2)[:,1];
	post_branch[:,end] = mean(sim_res1[4],dims=2)[:,1];
	post_trunk[:,end] = sim_res1[5];
	
#outdir = "y2/"



#post_RZ, post_ET, post_LWP, post_branch, post_trunk;
CSV.write(string(outdir,"postRZ.csv"), DataFrame(post_RZ) );
CSV.write(string(outdir,"postET.csv"), DataFrame(post_ET));
CSV.write(string(outdir,"postLeaf.csv"), DataFrame(post_LWP));
CSV.write(string(outdir,"postBranch.csv"), DataFrame(post_branch));
CSV.write(string(outdir,"postTrunk.csv"), DataFrame(post_trunk));

end

write_post(ARGS[1])

