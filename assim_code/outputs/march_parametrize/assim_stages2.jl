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


using Distributions
using LabelledArrays
using LinearAlgebra

function do_assim(obs_mask, norm_factor, outdir, vodA, vodB, vodC,obs_file,init_file)

obsdata = CSV.read(obs_file,DataFrame);

obsH = convert(Array,obsdata.hpol);
obsV = convert(Array,obsdata.vpol);

prior_min = [10, 0.75, 0.1, 500, 0.01,  0.01,100,0.1e-6,0.3,0.01,1];
prior_max = [120,15,  50, 3000, 10, 10,1000,    4e-6,0.5,5,8];
prior_dist = Product(Uniform.(log.(prior_min), log.(prior_max)));


function log_p_nolabel(a, k0, alpha0, beta0, maxiter, stop_crit, error_var, old_TB, old_par)
	p = logpdf(prior_dist, a);
	p0 = logpdf(prior_dist, old_par);
	p0 += norm_factor*logpdf(MvNormal(old_TB[1][obs_mask], sqrt(error_var)), obsH[obs_mask]);
	p0 += norm_factor*logpdf(MvNormal(old_TB[2][obs_mask], sqrt(error_var)), obsV[obs_mask]);		

	if isinf(p)
		return -Inf, 0,0,0,zeros(3);
	else
		try
			sim_res = run_sim_2(convert(Array{FT}, exp.(a))...);

			if isnan(mean(sim_res[1].leafpot))
				return -Inf, 0, 0, 0, zeros(3);

                        else
				canpot = mean(sim_res[6][:,1:3],dims=2); 
				soilsurf = sim_res[2][:,1];
                                n_obs = length(obsH[obs_mask]);
                                n_to_skip = 5;#trunc(Int, n_obs / 240);
                                optim_res = optim_vod_stages2(obsH[obs_mask], obsV[obs_mask], canpot[obs_mask], soilsurf[obs_mask],
                                                                        tsoil[obs_mask], tcan[obs_mask], laiM[obs_mask], k0, alpha0, beta0, maxiter, stop_crit, 1e-5, n_to_skip)
                                sim_TB = get_TB_2(soilsurf, tsoil, tcan, canpot, laiM,optim_res[1], optim_res[2],optim_res[3]);

				p += norm_factor*logpdf(MvNormal(sim_TB[1][obs_mask], sqrt(error_var)), obsH[obs_mask]);
				p += norm_factor*logpdf(MvNormal(sim_TB[2][obs_mask], sqrt(error_var)), obsV[obs_mask]);		

				return p,p0,sim_res,sim_TB,optim_res;
			end
		catch err_occurred
			return -Inf, 0,0,0,zeros(3);
		end
	end
end


function runAMH(x_init, niter,start_full,burn_len)

	k0 = exp(rand(Uniform(log(0.01),log(0.5))));
	alpha0 = exp(rand(Uniform(log(0.25),log(2))));
	beta0 = exp(rand(Uniform(log(0.01),log(0.5))));

	pars = Array(x_init);

	err0 = 1.3^2;

	nPar= length(pars);

	p,p0,sim_res0,tb0,optim_res0 = log_p_nolabel(pars,k0,alpha0,beta0,10000,10^-6, err0, (obsH,obsV),pars);
	ll0 = p;
	mycov = I*0.01/length(pars);
	chain = zeros(niter, nPar*2+2+3);
	
	nday = Int(length(sim_res1[1].leafpot)/24);
	
	for j in 1:niter


		if j == 2500
			CSV.write(string(outdir,"first2500.csv"),DataFrame(chain[1:2500,:]));
		end
		if j % 2 == 0
			println(j)
		end

		if j % 25 == 0
			println(chain[j-1,:])
		end

		#each row of chain is [pars old, pars proposed, old log likelihood, old err estimate, old VOD par]
		
		chain[j,1:nPar] = pars;
		chain[j,nPar*2+1] = ll0;
		chain[j,nPar*2+2] = err0;
		chain[j,(nPar*2+3):end] = optim_res0[1:3];


		if j > start_full+burn_len
			mycov = (1-10^-6)*cov(chain[start_full:j,1:nPar])  + 10^-6*I 
		end

		parsP = rand(MvNormal(pars, mycov),1)[:,1];
	    chain[j,(nPar+1):(nPar*2)] = parsP;
		p,p0,sim_resP,tbP,optim_resP = log_p_nolabel(parsP,optim_res0[1],optim_res0[2],optim_res0[3],10000,10^-6, err0, tb0,pars);
		mhratio = exp(p-p0);
		randI = rand();
		if randI < mhratio
			ll0 = p;
			pars = parsP;
			optim_res0 = optim_resP;
			tb0 = tbP;
			sim_res0 = sim_resP;
		end
		
	end
	
	#place the results of the original true run in the last column
	
	return chain #, post_RZ, post_ET, post_LWP, post_branch, post_trunk;
end

function init_works()

	k0 = exp(rand(Uniform(log(0.01),log(0.5))));
	alpha0 = exp(rand(Uniform(log(0.25),log(2))));
	beta0 = exp(rand(Uniform(log(0.01),log(0.5))));

	x0 = rand(prior_dist, 1)[:,1];
	ll_init,p0,sim_resP,tbP,optim_resP = log_p_nolabel(Array(x0), k0,alpha0,beta0,10000,10^-6,1.3^2,(obsH,obsV),x0);
	while isinf(ll_init) | isnan(ll_init)
		x0 = rand(prior_dist, 1)[:,1];
		ll_init,p0,sim_resP,tbP,optim_resP = log_p_nolabel(Array(x0), k0,alpha0,beta0,10000,10^-6,1.3^2,(obsH,obsV),x0);
	end
	return x0, ll_init
end

NPAR = length(prior_min);

opt_sol = Array(CSV.read(init_file,DataFrame).mymin);
opt_sol_rel = (opt_sol .- log.(prior_min)) ./ (log.(prior_max) - log.(prior_min));
perturb_start_rel = opt_sol_rel + 0.1*(rand(NPAR) .- 0.5);
perturb_start_rel[perturb_start_rel .<= 0] .= 0.05;
perturb_start_rel[perturb_start_rel .>= 1] .= 0.95;
 
a01 = perturb_start_rel .* (log.(prior_max) - log.(prior_min)) .+ log.(prior_min); 

c1 = runAMH(a01,10000,1000,500);

dfc = DataFrame(c1);
CSV.write(string(outdir,"post_par.csv"),dfc);

end
