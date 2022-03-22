include("../../../get_home_dir.jl")

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

N = 24*365*1
istart = 24*365*2 + 1;
soil0 = 0.4;

include(string(PROJECT_HOME,"/simulation_code/full_model_newP63.jl"));

deltaT = FT(60*60);
alpha = FT(10.9);
nsoil = FT(1.5);

function run_sim_2(vcmax_par::FT, weibB_stom::FT, k_plant::FT,
      z_soil::FT, xylem_extraB::FT, vol_factor::FT, g1::FT, k_soil::FT,
      slope_runoff::FT, exp_root_dist::FT,weibC_plant::FT
      )
return run_sim_varyB(vcmax_par, weibB_stom, weibB_stom+xylem_extraB, weibC_plant, k_plant,
    k_soil, z_soil, istart, N, soil0, vol_factor, 1e-5, 3, df_raw, g1, deltaT, alpha, nsoil,0,0,
    slope_runoff, exp_root_dist,FT(1/20),FT(1/20));
end

pars0 = convert(Array{FT}, [90, 3, 10, 2000, 1, 1.0, 300, 0.4e-6,
         0.4,2,4]);

sim_res1 = run_sim_2(pars0...);

#vcmax_par, k_frac, k_plant, k_soil, z_soil, istart, N, smc0, slope_index,use_flux, weibB, weibC

noise_std = 1.3;
noise_var = noise_std^2

tsoil = sim_res1[1].T_SOIL .+ 273.15;
tair = sim_res1[1].T_AIR .+ 273.15;
tcan = (sim_res1[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
laiM = sim_res1[1].LAI_modis;

N = 24*100;
istart = 24*365*2 + 24*125 + 1; 

tsoil = tsoil[(24*125+1):(24*125+N)];
tair = tair[(24*125+1):(24*125+N)];
tcan = tcan[(24*125+1):(24*125+N)];
laiM = laiM[(24*125+1):(24*125+N)];


using Distributions
using LabelledArrays
using LinearAlgebra
using Optim


prior_min = [10, 0.75, 0.1, 500, 0.01,  0.01,100,0.1e-6,0.3,0.01,1];
prior_max = [120,15,  50, 3000, 10, 10,1000,    4e-6,0.5,5,8];
prior_dist = Product(Uniform.(log.(prior_min), log.(prior_max)));


struct MySimplexer <: Optim.Simplexer end
Optim.simplexer(S::MySimplexer, initial_x) = [rand(prior_dist) for i = 1:length(initial_x)+1]

function do_assim(obs_mask_full, norm_factor, outdir, vodA, vodB, vodC,obs_file)

obsdata = CSV.read(obs_file,DataFrame);

obsH = convert(Array,obsdata.hpol)[(24*125+1):(24*125+N)];
obsV = convert(Array,obsdata.vpol)[(24*125+1):(24*125+N)];

obs_mask = obs_mask_full[(24*125+1):(24*125+N)];

function log_p_nolabel(a, k0, alpha0, beta0, maxiter, stop_crit, error_var)

p = -50000*(sum(a .< log.(prior_min)) + sum(a .> log.(prior_max)));

sim_res = 1;
try
	sim_res = run_sim_2(convert(Array{FT}, exp.(a))...);
catch err_occurred
	return -50000*11*2;
end

	canpot = mean(sim_res[6][:,1:3],dims=2); 
	soilsurf = sim_res[2][:,1];
        n_obs = length(obsH[obs_mask]);
        n_to_skip = 5;#trunc(Int, n_obs / 240);
        optim_res = optim_vod_stages2(obsH[obs_mask], obsV[obs_mask], canpot[obs_mask], soilsurf[obs_mask],
                                      tsoil[obs_mask], tcan[obs_mask], laiM[obs_mask], k0, alpha0, beta0, maxiter, stop_crit, 1e-5, n_to_skip)
        sim_TB = get_TB_2(soilsurf, tsoil, tcan, canpot, laiM,optim_res[1], optim_res[2],optim_res[3]);

	p += norm_factor*logpdf(MvNormal(sim_TB[1][obs_mask], sqrt(error_var)), obsH[obs_mask]);
	p += norm_factor*logpdf(MvNormal(sim_TB[2][obs_mask], sqrt(error_var)), obsV[obs_mask]);		

	return p
end

function mywrapper(a)
    -log_p_nolabel(a,vodA,vodB,vodC,10000,10^-6, 1.3^2)
end

#=
if isfile(init_name)
        ipar_LL = Array(CSV.read(init_name, DataFrame))
else
        ipar_LL = zeros(100,NPAR+1);
        for j in 1:100
                println(j)
                xi, LLi = init_works();
                ipar_LL[j,1:NPAR] = xi;
                ipar_LL[j,end] = LLi;
        end
        CSV.write(init_name, DataFrame(ipar_LL));
end

a01 = ipar_LL[argmax(ipar_LL[:,end]),1:NPAR];
=#
a01 = mean(prior_dist);


myopt = optimize(mywrapper, a01,NelderMead(initial_simplex = MySimplexer()),   
 Optim.Options(show_trace=true,iterations=500,store_trace=true))

dfmin = DataFrame(mymin=Optim.minimizer(myopt));
CSV.write(string(outdir,"opt_par.csv"),dfmin);

dftrace = DataFrame(mytrace=Optim.f_trace(myopt));
CSV.write(string(outdir,"opt_err.csv"),dftrace);

end
