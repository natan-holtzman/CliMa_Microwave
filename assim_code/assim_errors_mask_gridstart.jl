
using DataFrames
using CSV
#using Plots
using Random
using StatsBase

include("mironov.jl");
include("tau_omega_funs.jl");
include("../simulation_code/rebuild_sim_N.jl");
include("time_averaging.jl")


#N = 48*3
#istart = 48*365*2 - 52*48 + 230*48


N = 48*365
istart = 48*365*2 - 52*48 #+ 230*48
soil0 = 0.39;

sim_res1 = convert_sim(run_sim(FT(22),FT(1.5), FT(4.0), FT(1e-4), FT(2.5), FT(2.0), FT(700),FT(0.45),istart,N, soil0));


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

function log_p_nolabel(a, k0, alpha0, beta0, maxiter, stop_crit, error_var, old_TB)

	v = LVector(logVcmax=a[1], logPcrit=a[2], logKmaxPlant=a[3], logKmaxSoil=a[4], logZsoil=a[5], logNsoil=a[6]);

	p = 0.0;
	p += logpdf(Uniform(log(10), log(120)), v.logVcmax);
	p += logpdf(Uniform(log(0.1), log(100)), v.logKmaxPlant);
	p += logpdf(Uniform(log(0.75), log(5)), v.logPcrit);
        p += logpdf(Uniform(log(1e-5), log(2e-3)), v.logKmaxSoil);
        p += logpdf(Uniform(log(0.4), log(0.6)), v.logNsoil);
        p += logpdf(Uniform(log(500), log(3000)), v.logZsoil);

		
sim_res = convert_sim(run_sim(FT(exp(v.logVcmax)),FT(exp(v.logPcrit)), FT(exp(v.logKmaxPlant)), FT(exp(v.logKmaxSoil)), FT(2.5), FT(2.0), FT(exp(v.logZsoil)), 
			  FT(exp(v.logNsoil)), istart, N, soil0));

	if isnan(mean(sim_res[1].leafpot))
		return -Inf;

	else
	
	p0 = p+0;
	p0 += logpdf(MvNormal(old_TB[1][obs_mask], sqrt(error_var)), obsH[obs_mask]);
	p0 += logpdf(MvNormal(old_TB[2][obs_mask], sqrt(error_var)), obsV[obs_mask]);		


	optim_res = optim_vod_stages(obsH[obs_mask], obsV[obs_mask],
			    mean(sim_res[3],dims=2)[obs_mask], sim_res[2][:,1][obs_mask],
			    sim_res[1][obs_mask,:], k0, alpha0, beta0, maxiter, stop_crit,10^-5, 1);
		sim_TB = get_TB(sim_res, mean(sim_res[3],dims=2), optim_res[1], optim_res[2],optim_res[3]);

		p += logpdf(MvNormal(sim_TB[1][obs_mask], sqrt(error_var)), obsH[obs_mask]);
		p += logpdf(MvNormal(sim_TB[2][obs_mask], sqrt(error_var)), obsV[obs_mask]);		


		return p, optim_res, sim_TB, p0
	end
end


function log_p_err(error_var, new_error_var, pred_TB)
	p0 = 0;
	p1 = 0;
	
	p0 += logpdf(truncated(Cauchy(), 0, Inf), error_var);
	p1 += logpdf(truncated(Cauchy(), 0, Inf), new_error_var);

	p0 -= 0.5*(sum( (pred_TB[1]-obsH) .^2 ./ error_var .+ log(error_var)) + sum( (pred_TB[2]-obsV) .^2 ./ error_var .+ log(error_var)))
	p1 -= 0.5*(sum( (pred_TB[1]-obsH) .^2 ./ new_error_var .+ log(new_error_var)) + sum( (pred_TB[2]-obsV) .^2 ./ new_error_var .+ log(new_error_var)))

	return p0, p1
end



function runAMH(x_init, niter, burnlen)

k0 = exp(rand(Uniform(log(0.01),log(0.75))));
alpha0 = exp(rand(Uniform(log(0.25),log(2))));
beta0 = exp(rand(Uniform(log(1/100),log(1/25))));

#pars = log.([FT(22),FT(1.5), FT(4.0), FT(1e-4), FT(700),FT(0.45)]);
pars = Array(x_init);

err0 = 1.7;

nPar= length(pars);

ll0,vodpar0,tb1,ll0prev = log_p_nolabel(pars,k0,alpha0,beta0,10000,10^-6, err0, (obsH,obsV));
mycov = I*0.01/length(pars);
chain = zeros(niter, nPar*2+11+2);

for j in 1:niter
	if j % 2 == 0
		println(j)
	end

        if j % 25 == 0
                println(chain[j-1,:])
        end



	
	chain[j,1:nPar] = pars;
	chain[j,(nPar+1):(nPar+3)] = vodpar0[1:3];
	chain[j,nPar+4] = ll0

	if j > burnlen
		mycov = (1-10^-6)*cov(chain[(j-burnlen):j,1:nPar])  + 10^-6*I;
	end


    parsP = rand(MvNormal(pars, mycov),1);
	#llP = log_p_nolabel(parsP,0,0,0,10000,10^-5);
    llP,vodparP,tbP,ll0prev = log_p_nolabel(parsP,vodpar0[1],vodpar0[2],vodpar0[3],10000,10^-6, err0, tb1);
   chain[j,nPar+5] = llP;
   chain[j,(nPar+6):(nPar+10)] = vodparP;
    mhratio = exp(llP - ll0prev);
    randI = rand();
   chain[j,nPar+11] = randI;
   chain[j,(nPar+12):(nPar+12+nPar-1)] = parsP;
    if randI < mhratio
		ll0 = llP*1;
		pars = parsP[:,1]*1;
		vodpar0 = vodparP*1;
		tb1 = tbP;
	end
	
	err_prop = rand(LogNormal(log(err0), 0.01));
	chain[j,(nPar+12+nPar):end] = [err0, err_prop];

	p0_new, p1_new = log_p_err(err0, err_prop, tb1);
	mhratio_err = exp(p1_new - p0_new) * err_prop/err0;
	randIe = rand();
	if randIe < mhratio_err
		err0 = err_prop;
	end
	
	
end

return chain
end

function init_works()

k0 = exp(rand(Uniform(log(0.01),log(0.75))));
alpha0 = exp(rand(Uniform(log(0.25),log(2))));
beta0 = exp(rand(Uniform(log(1/100),log(1/25))));


x0 = LVector(logVcmax=rand(Uniform(log(10), log(120))),
			 logPcrit=rand(Uniform(log(0.75), log(5))),
			 logKmaxPlant=rand(Uniform(log(0.1), log(100))),
			 logKmaxSoil=rand(Uniform(log(1e-5), log(2e-3))), 
 logZsoil=rand(Uniform(log(500), log(3000))),
 logNsoil=rand(Uniform(log(0.4), log(0.6)))
			 );


ll_init, vodpar_init,a1,a2 = log_p_nolabel(Array(x0), k0,alpha0,beta0,10000,10^-6,1.7,(obsH,obsV));

if isnan(ll_init)

while isnan(ll_init)



x0 = LVector(logVcmax=rand(Uniform(log(10), log(120))),
                         logPcrit=rand(Uniform(log(0.75), log(5))),
                         logKmaxPlant=rand(Uniform(log(0.1), log(100))),
                         logKmaxSoil=rand(Uniform(log(1e-5), log(2e-3))),
 logZsoil=rand(Uniform(log(500), log(3000))),
 logNsoil=rand(Uniform(log(0.4), log(0.6)))
                         );




ll_init, vodpar_init,a1,a2 = log_p_nolabel(Array(x0), k0,alpha0,beta0,10000,10^-6,1.7,(obsH,obsV));
end
end

return x0, ll_init

end


init_name = "init_par_LL_2x3_6par.csv"

if isfile(init_name)
        ipar_LL = Array(CSV.read(init_name, DataFrame))
else
	ipar_LL = zeros(100,7);

	for j in 1:100
		println(j)
		xi, LLi = init_works();
		ipar_LL[j,1:6] = Array(xi);
		ipar_LL[j,7] = LLi;
	end
end

a01 = ipar_LL[argmax(ipar_LL[:,7]),1:7];
c1 = runAMH(a01, 10000, 500);

dfc = DataFrame(c1);

CSV.write("post_2x3obs_newvar_gridstart.csv",dfc);



