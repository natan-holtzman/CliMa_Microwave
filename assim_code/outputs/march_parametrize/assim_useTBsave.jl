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

N = 24*365*2
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

function do_assim(obs_mask, norm_factor, outdir, vodA, vodB, vodC,obs_file)

#trueTB = get_TB_2(sim_res1[2][:,1], tsoil, tcan, canpot_true, laiM, vodA, vodB, vodC);
#obsH = trueTB[1] + noise_std*randn(Float64, length(trueTB[1]));
#obsV = trueTB[2] + noise_std*randn(Float64, length(trueTB[1]));

#obsDF = DataFrame(hpol=obsH, vpol=obsV);
obsdata = CSV.read(obs_file,DataFrame);

obsH = convert(Array,obsdata.hpol);
obsV = convert(Array,obsdata.vpol);

#obs_mask = zeros(length(obsH));
#obs_mask[7:(24*3):end] .= 1;
#obs_mask[(7+12):(24*3):end] .= 1;

#obs_mask = ones(length(obsH));

#obs_mask = obs_mask .== 1;

#norm_factor = FT(1);

prior_min = [10, 0.01, 0.1, 500, 0.75,  0.01,100];
prior_max = [120,0.9,  50, 3000, 10, 10,1000];
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
                                canpot = mean(sim_res[3],dims=2);
				soilsurf = sim_res[2][:,1];
                                n_obs = length(obsH[obs_mask]);
                                n_to_skip = trunc(Int, n_obs / 240);
                                optim_res = optim_vod_stages2(obsH[obs_mask], obsV[obs_mask], canpot[obs_mask], soilsurf[obs_mask],
                                                                        tsoil[obs_mask], tcan[obs_mask], laiM[obs_mask], k0, alpha0, beta0, maxiter, stop_crit, 1e-5, n_to_skip)
                                sim_TB = get_TB_2(soilsurf, tsoil, tcan, canpot, laiM,optim_res[1], optim_res[2],optim_res[3]);

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



function runAMH(x_init, niter, burnlen,full_cov_len)

	k0 = exp(rand(Uniform(log(0.01),log(0.5))));
	alpha0 = exp(rand(Uniform(log(0.25),log(2))));
	beta0 = exp(rand(Uniform(log(0.01),log(0.5))));

	pars = Array(x_init);

	err0 = 1.7;

	nPar= length(pars);

	p,p0,sim_res0,tb0,optim_res0 = log_p_nolabel(pars,k0,alpha0,beta0,10000,10^-6, err0, (obsH,obsV),pars);
	ll0 = p;
	mycov = I*0.01/length(pars);
	chain = zeros(niter, nPar*2+2+3);
	
	nday = Int(length(sim_res1[1].leafpot)/24);
	
	for j in 1:niter


		#if j == 5000
		#	CSV.write(string(outdir,"first5000.csv"),DataFrame(chain[1:5000,:]));
		#end
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

		if j > burnlen
			mycov = (1-10^-6)*cov(chain[(j-burnlen):j,1:nPar])  + 10^-6*I;
		end

		if j > burnlen+full_cov_len
			mycov = (1-10^-6)*cov(chain[full_cov_len:j,1:nPar])  + 10^-6*I 
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
	
	return chain #, post_RZ, post_ET, post_LWP, post_branch, post_trunk;
end

function init_works()

	k0 = exp(rand(Uniform(log(0.01),log(0.5))));
	alpha0 = exp(rand(Uniform(log(0.25),log(2))));
	beta0 = exp(rand(Uniform(log(0.01),log(0.5))));

	x0 = rand(prior_dist, 1)[:,1];
	ll_init,p0,sim_resP,tbP,optim_resP = log_p_nolabel(Array(x0), k0,alpha0,beta0,10000,10^-6,1.7,(obsH,obsV),x0);
	while isinf(ll_init) | isnan(ll_init)
		x0 = rand(prior_dist, 1)[:,1];
		ll_init,p0,sim_resP,tbP,optim_resP = log_p_nolabel(Array(x0), k0,alpha0,beta0,10000,10^-6,1.7,(obsH,obsV),x0);
	end
	return x0, ll_init
end

NPAR = length(prior_min);

#outdir = "./"
init_name = string(outdir,"init_par.csv");


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
c1 = runAMH(a01, 5000, 500, 5000);

dfc = DataFrame(c1);
CSV.write(string(outdir,"post_par.csv"),dfc);

end
