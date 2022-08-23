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

N = 24*365*12
istart = 24*365*0 + 1;
soil0 = 0.4;

include(string(PROJECT_HOME,"/simulation_code/full_model_newP63.jl"));
include(string(PROJECT_HOME,"/simulation_code/full_model_newP63_ballberry.jl"));

deltaT = FT(60*60);

soil_curve_1 = VanGenuchten{FT}(stype = "Ozark",α = FT(10.9),   n = 1.5, Θs = 0.55,   Θr = 0.067);
smc1 = soil_swc(soil_curve_1,FT(-1));


function run_sim_2(vcmax_par::FT, weibB_stom::FT, k_plant::FT,
    z_soil::FT, xylem_extraB::FT, vol_factor::FT, g1::FT, k_soil::FT,
    slope_runoff::FT, exp_root_dist::FT,nsoil::FT,soilfac::FT
    )
weibC_plant = FT(2);
soil_curve_2 =  VanGenuchten{FT}(stype = "Ozark",α = FT(10.9),   n = nsoil, Θs = 0.55,   Θr = 0.067);
pot1b = soil_p_25_swc(soil_curve_2,smc1);
alpha2 = FT(10.9*soilfac)*abs(pot1b);
return run_sim_BallBerry(vcmax_par, weibB_stom/soilfac, (weibB_stom+xylem_extraB)/soilfac, weibC_plant, k_plant*soilfac,
  k_soil*soilfac, z_soil, istart, N, soil0, vol_factor*soilfac, 1e-5, 3, df_raw, g1, deltaT, alpha2, nsoil,0,0,
  slope_runoff, exp_root_dist,FT(1/20),FT(1/20));
end

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

noise_std = 1.3;
noise_var = noise_std^2

tsoil = sim_res1[1].T_SOIL .+ 273.15;
tair = sim_res1[1].T_AIR .+ 273.15;
tcan = (sim_res1[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
laiM = sim_res1[1].LAI_modis;

function write_post(mydir)

par_post = convert(Array,CSV.read(string(mydir,"post_par.csv"),DataFrame))[(end-4000+1):100:end,1:12];
outdir = mydir

nday = convert(Int, N/24);
npost = size(par_post)[1];

#post_RZ = zeros(nday,npost+1);
#post_Surf = zeros(nday,npost+1);
post_ET = zeros(nday*24,npost+1);
post_LWP = zeros(nday*24,npost+1);
#post_branch = zeros(nday*24,npost+1);
#post_trunk = zeros(nday*24,npost+1);
post_gsw = zeros(nday*24,npost+1);
post_TWS = zeros(nday,npost+1);
#post_TWP = zeros(nday,npost+1);
post_GPP = zeros(nday*24,npost+1);

for j in 1:npost
	println(j)
	a = par_post[j,:]
try
	sim_res0 = run_sim_2(convert(Array{FT}, exp.(a))...); 

	#post_RZ[:,j] = get_daily(sim_res0[2][:,end],24);
        post_TWS[:,j] = get_daily(sim_res0[1].ColumnSMC,24);
        #post_TWP[:,j] = get_daily(sim_res0[1].ColumnSWP,24);
        post_GPP[:,j] = sim_res0[1].GPP;
	post_ET[:,j] = sim_res0[1].ETmod;
	post_gsw[:,j] = sim_res0[1].glw;
	post_LWP[:,j] =  mean(sim_res0[6][:,1:3],dims=2); 
	#post_branch[:,j] = mean(sim_res0[4],dims=2)[:,1];
	#post_trunk[:,j] = sim_res0[5];
catch err
	println("Error", a)
end
end

	#place the results of the original true run in the last column
	#post_RZ[:,end] = get_daily(sim_res1[2][:,end],24);
        post_TWS[:,end] = get_daily(sim_res1[1].ColumnSMC,24);
        #post_TWP[:,end] = get_daily(sim_res1[1].ColumnSWP,24);
	post_GPP[:,end] = sim_res1[1].GPP;
	post_ET[:,end] = sim_res1[1].ETmod;
	post_gsw[:,end] = sim_res1[1].glw;
	post_LWP[:,end] = mean(sim_res1[6][:,1:3],dims=2);
	#post_branch[:,end] = mean(sim_res1[4],dims=2)[:,1];
	#post_trunk[:,end] = sim_res1[5];
	
#outdir = "y2/"



#post_RZ, post_ET, post_LWP, post_branch, post_trunk;
#CSV.write(string(outdir,"postRZ.csv"), DataFrame(post_RZ) );
CSV.write(string(outdir,"postSWS.csv"), DataFrame(post_TWS) );
#CSV.write(string(outdir,"postSWP.csv"), DataFrame(post_TWP) );
CSV.write(string(outdir,"postET.csv"), DataFrame(post_ET));
CSV.write(string(outdir,"postLeaf.csv"), DataFrame(post_LWP));
#CSV.write(string(outdir,"postBranch.csv"), DataFrame(post_branch));
#CSV.write(string(outdir,"postTrunk.csv"), DataFrame(post_trunk));
CSV.write(string(outdir,"postGSW.csv"), DataFrame(post_gsw));
CSV.write(string(outdir,"postGPP.csv"), DataFrame(post_GPP));
end

write_post(ARGS[1])

