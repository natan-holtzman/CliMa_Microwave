
using DataFrames
using CSV
#using Plots
using Random
using StatsBase

include("mironov.jl");
include("tau_omega_funs.jl");
include("time_averaging.jl")

df_raw = CSV.read("../data/moflux_land_data_newnames_7.csv", DataFrame);

N = 48*365
istart = 48*365*2 - 52*48 #+ 230*48
soil0 = 0.39;

#include("../simulation_code/rebuild_sim_Scrit_drain_plugin_ET_setsoil.jl");
#function run_sim_2(vcmax_par::FT, k_frac::FT, k_plant::FT, k_soil::FT, z_soil::FT, weibB::FT, weibC::FT)
#	return convert_sim(run_sim(vcmax_par, k_frac, k_plant, k_soil, z_soil, istart, N, soil0, FT(0), 0, weibB, weibC));
#end

include("../simulation_code/sim_vary_new_stomata.jl");
function run_sim_2(vcmax_par::FT, k_frac::FT, k_plant::FT, k_soil::FT, z_soil::FT, weibB::FT, weibC::FT)
	return convert_sim(run_sim_vary(vcmax_par, k_frac, weibB, weibC, k_plant, k_soil, z_soil, istart, N, soil0, 1, 1e-5, 1, df_raw));
end

sim_res1 = run_sim_2(FT(22),FT(0.33), FT(15.0),FT(1e-5), FT(800),FT(1.5),FT(1));
#vcmax_par, k_frac, k_plant, k_soil, z_soil, istart, N, smc0, slope_index,use_flux, weibB, weibC

noise_std = 1.3;
noise_var = noise_std^2

trueTB = get_TB(sim_res1, mean(sim_res1[3],dims=2), 0.25, 1, 1/50);
obsH = trueTB[1] + noise_std*randn(Float64, length(trueTB[1]));
obsV = trueTB[2] + noise_std*randn(Float64, length(trueTB[1]));

obsH2 = trueTB[1] + noise_std*randn(Float64, length(trueTB[1]));
obsV2 = trueTB[2] + noise_std*randn(Float64, length(trueTB[1]));

using Distributions
using LabelledArrays
using LinearAlgebra

obs_mask = zeros(length(obsH));
obs_mask[2:(24*3):end] .= 1;
obs_mask[(2+12):(24*3):end] .= 1;

#obs_mask = ones(length(obsH));

obs_mask = obs_mask .== 1;

norm_factor = FT(1);

prior_min = [10, 0.01, 0.1,  1e-6, 500, 0.75, 0.75];
prior_max = [120,0.75,  100, 2e-5, 3000, 10, 8];
prior_dist = Product(Uniform.(log.(prior_min), log.(prior_max)));


function log_p_nolabel(a, k0, alpha0, beta0, maxiter, stop_crit, error_var, old_TB, old_par)
	p = logpdf(prior_dist, a);
	p0 = logpdf(prior_dist, old_par);
	p0 += norm_factor*logpdf(MvNormal(old_TB[1][obs_mask], sqrt(error_var)), obsH[obs_mask]);
	p0 += norm_factor*logpdf(MvNormal(old_TB[2][obs_mask], sqrt(error_var)), obsV[obs_mask]);		

	if isinf(p)
		return -Inf, 0,0,0,0;
	else
		try
			sim_res = run_sim_2(convert(Array{FT}, exp.(a))...);

			if isnan(mean(sim_res[1].leafpot))
				return -Inf, 0, 0, 0, 0;
			else
				optim_res = optim_vod_stages(obsH[obs_mask], obsV[obs_mask],
							mean(sim_res[3],dims=2)[obs_mask], sim_res[2][:,1][obs_mask],
							sim_res[1][obs_mask,:], k0, alpha0, beta0, maxiter, stop_crit,10^-5, 1);
				sim_TB = get_TB(sim_res, mean(sim_res[3],dims=2), optim_res[1], optim_res[2],optim_res[3]);

				p += norm_factor*logpdf(MvNormal(sim_TB[1][obs_mask], sqrt(error_var)), obsH[obs_mask]);
				p += norm_factor*logpdf(MvNormal(sim_TB[2][obs_mask], sqrt(error_var)), obsV[obs_mask]);		

				return p,p0,sim_res,sim_TB,optim_res;
			end
		catch err_occurred
			return -Inf, 0,0,0,0;
		end
	end
end


function log_p_err(error_var, new_error_var, pred_TB)
	p0 = 0;
	p1 = 0;
	
	p0 += logpdf(truncated(Cauchy(), 0, Inf), error_var);
	p1 += logpdf(truncated(Cauchy(), 0, Inf), new_error_var);

	p0 += norm_factor*-0.5*(sum( (pred_TB[1]-obsH) .^2 ./ error_var .+ log(error_var)) + sum( (pred_TB[2]-obsV) .^2 ./ error_var .+ log(error_var)))
	p1 += norm_factor*-0.5*(sum( (pred_TB[1]-obsH) .^2 ./ new_error_var .+ log(new_error_var)) + sum( (pred_TB[2]-obsV) .^2 ./ new_error_var .+ log(new_error_var)))

	return p0, p1
end



function runAMH(x_init, niter, burnlen)

	k0 = exp(rand(Uniform(log(0.01),log(0.75))));
	alpha0 = exp(rand(Uniform(log(0.25),log(2))));
	beta0 = exp(rand(Uniform(log(1/100),log(1/25))));

	pars = Array(x_init);

	err0 = 1.7;

	nPar= length(pars);

	p,p0,sim_res0,tb0,optim_res0 = log_p_nolabel(pars,k0,alpha0,beta0,10000,10^-6, err0, (obsH,obsV),pars);
	ll0 = p;
	mycov = I*0.01/length(pars);
	chain = zeros(niter, nPar*2+2+3);
	
	nday = Int(length(sim_res1[1].leafpot)/24);
	
	post_RZ = zeros(nday,Int(niter/25)+1);
	#post_Surf = zeros(nday,Int(niter/25)+1);
	post_ET = zeros(nday,Int(niter/25)+1);
	post_LWP = zeros(nday*24,Int(niter/25)+1);
	post_branch = zeros(nday*24,Int(niter/25)+1);
	post_trunk = zeros(nday*24,Int(niter/25)+1);

	
	for j in 1:niter
		if j % 2 == 0
			println(j)
		end

		if j % 25 == 0
			println(chain[j-1,:])
			post_RZ[:,Int(j/25)] = get_daily(sim_res0[2][:,end],24);
			#post_Surf[:,Int(j/25)] = get_daily(sim_res0[2][:,1],24);
			post_ET[:,Int(j/25)] = get_daily(sim_res0[1].ETmod,24);
			post_LWP[:,Int(j/25)] = mean(sim_res0[3],dims=2)[:,1];
			post_branch[:,Int(j/25)] = mean(sim_res0[4],dims=2)[:,1];
			post_trunk[:,Int(j/25)] = sim_res0[5];
		end

		#each row of chain is [pars old, pars proposed, old log likelihood, old err estimate, old VOD par]
		
		chain[j,1:nPar] = pars;
		chain[j,nPar*2+1] = ll0;
		chain[j,nPar*2+2] = err0;
		chain[j,(nPar*2+3):end] = optim_res0[1:3];

		if j > burnlen
			mycov = (1-10^-6)*cov(chain[(j-burnlen):j,1:nPar])  + 10^-6*I;
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
		
		err_prop = rand(LogNormal(log(err0), 0.01));

		p0_new, p1_new = log_p_err(err0, err_prop, tb0);
		mhratio_err = exp(p1_new - p0_new) * err_prop/err0;
		randIe = rand();
		if randIe < mhratio_err
			err0 = err_prop;
		end
	end
	
	#place the results of the original true run in the last column
	post_RZ[:,end] = get_daily(sim_res1[2][:,end],24);
	#post_Surf[:,end] = get_daily(sim_res1[2][:,1],24);
	post_ET[:,end] = get_daily(sim_res1[1].ETmod,24);
	post_LWP[:,end] = mean(sim_res1[3],dims=2)[:,1];
	post_branch[:,end] = mean(sim_res1[4],dims=2)[:,1];
	post_trunk[:,end] = sim_res1[5];
	
	return chain, post_RZ, post_ET, post_LWP, post_branch, post_trunk;
end

function init_works()

	k0 = exp(rand(Uniform(log(0.01),log(0.75))));
	alpha0 = exp(rand(Uniform(log(0.25),log(2))));
	beta0 = exp(rand(Uniform(log(1/100),log(1/25))));

	x0 = rand(prior_dist, 1)[:,1];
	ll_init,p0,sim_resP,tbP,optim_resP = log_p_nolabel(Array(x0), k0,alpha0,beta0,10000,10^-6,1.7,(obsH,obsV),x0);
	while isinf(ll_init) | isnan(ll_init)
		x0 = rand(prior_dist, 1)[:,1];
		ll_init,p0,sim_resP,tbP,optim_resP = log_p_nolabel(Array(x0), k0,alpha0,beta0,10000,10^-6,1.7,(obsH,obsV),x0);
	end
	return x0, ll_init
end

NPAR = length(prior_min);

init_name = "obs2x3/init_par.csv"

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
c1 = runAMH(a01, 10000, 500);

dfc = DataFrame(c1[1]);

CSV.write("obs2x3/post_par.csv",dfc);

#post_RZ, post_ET, post_LWP, post_branch, post_trunk;
CSV.write("obs2x3/postRZ.csv", DataFrame(c1[2]));
CSV.write("obs2x3/postET.csv", DataFrame(c1[3]));
CSV.write("obs2x3/postLeaf.csv", DataFrame(c1[4]));
CSV.write("obs2x3/postBranch.csv", DataFrame(c1[5]));
CSV.write("obs2x3/postTrunk.csv", DataFrame(c1[6]));







