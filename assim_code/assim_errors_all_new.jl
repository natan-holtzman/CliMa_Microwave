
using DataFrames
using CSV
#using Plots
using Random
using StatsBase

include("mironov.jl");

include("../simulation_code/rebuild_sim_N.jl");




function avg2(x)
if (typeof(x[1])==Float64) | (typeof(x[1])==Float32)
        return (x[1:2:end] + x[2:2:end])/2;
else
        return x[1:2:end]
end

end


function avg2_2d(x)

return (x[1:2:end,:] + x[2:2:end,:])/2;

end


function avg2_df(x)

ncol = Integer(size(x)[2]);
nrow = Integer(size(x)[1]/2);
ans = deepcopy(x[1:nrow,:]);

for j in 1:ncol
#println(x[1:10,j])
#println(avg2(x[1:10,j]))
        ans[:,j] = avg2(x[:,j]);
end

return ans;

end

function convert_sim(x)
a = avg2_df(x[1]);
b = avg2_2d(x[2]);
c = avg2_2d(x[3]);
d = avg2_2d(x[4]);

return a,b,c,d
end




#N = 48*3
#istart = 48*365*2 - 52*48 + 230*48


N = 48*365
istart = 48*365*2 - 52*48 #+ 230*48
soil0 = 0.39;

sim_res1 = convert_sim(run_sim(FT(22),FT(1.5), FT(4.0), FT(1e-4), FT(2.5), FT(2.0), FT(700),FT(0.45),istart,N, soil0));


function get_TB(sim_res, leafpot, pvs, alpha, beta)

smc_mod = sim_res[2][:,1];

eps_real = 0*smc_mod; eps_im = 0*smc_mod;
for i in 1:length(eps_real)
        eps_i = mironov(1.4e9, smc_mod[i], 0.15);
        eps_real[i] = eps_i[1];
        eps_im[i] = eps_i[2];
end

eps = eps_real + eps_im*1im;

theta = 40*pi/180;
mycostheta = cos(theta);
mysin2theta = sin(theta)^2;

fH = abs2.((mycostheta .- sqrt.(eps .- mysin2theta)) ./ (mycostheta .+ sqrt.(eps .- mysin2theta)));
fV = abs2.((eps*mycostheta .- sqrt.(eps .- mysin2theta)) ./ (eps*mycostheta .+ sqrt.(eps .- mysin2theta)));

#leafpot = mean(sim_res1[3],dims=2);
leafvol_rel = 1 .+ pvs*leafpot;
laiM = sim_res[1].LAI_modis;

vod = leafvol_rel .* (alpha .+ laiM*beta);

gamma = exp.(-1*vod / mycostheta);
omega = 0.05
rhfac =  exp(-0.16*mycostheta^2)

Tsoil = sim_res[1].T_SOIL .+ 273.15;
Tcan = (sim_res[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
predH = gamma .* (1 .- rhfac*fH) .* Tsoil + (1-omega) * (1 .- gamma) .* (1 .+ gamma .* (rhfac*fH)) .* Tcan;
predV = gamma .* (1 .- rhfac*fV) .* Tsoil + (1-omega) * (1 .- gamma) .* (1 .+ gamma .* (rhfac*fV)) .* Tcan;

return predH[:,1], predV[:,1], vod

end


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



function optim_vod_while(obsHi, obsVi, leafpot, smc_mod, simtab, k, alpha, beta, maxiter, stop_crit, step_size)

#smc_mod = sim_res[2][:,1];

eps_real = zeros(length(smc_mod)); eps_im = zeros(length(smc_mod));
for i in 1:length(eps_real)
        eps_i = mironov(1.4e9, smc_mod[i], 0.15);
        eps_real[i] = eps_i[1];
        eps_im[i] = eps_i[2];
end

#eps = eps_real + eps_im*1im;

eps = complex.(eps_real, eps_im);


theta = 40*pi/180;
mycostheta = cos(theta);
mysin2theta = sin(theta)^2;

fH = abs2.((mycostheta .- sqrt.(eps .- mysin2theta)) ./ (mycostheta .+ sqrt.(eps .- mysin2theta)));
fV = abs2.((eps*mycostheta .- sqrt.(eps .- mysin2theta)) ./ (eps*mycostheta .+ sqrt.(eps .- mysin2theta)));

omega = 0.05;
rhfac =  exp(-0.16*mycostheta^2);

Tsoil = simtab.T_SOIL .+ 273.15;
Tcan = (simtab.LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
laiM = simtab.LAI_modis;


#k = exp(rand(Uniform(log(0.01),log(0.25))));
#alpha = exp(rand(Uniform(log(0.25),log(2))));
#beta = exp(rand(Uniform(log(1/100),log(1/25))));

#k = 0.051; alpha = 1.01; beta = 1/55;

#step_size = 10^-4;

function overall_grad(diff_overall, r_rough, dgamma, gamma_star)

	t1 = (1 .- r_rough) .* Tsoil .* dgamma;
	t2 = (1-omega) * (1 .+ gamma_star .* r_rough) .* Tcan .* (-1* dgamma);
	t3 = (1-omega) * (1 .- gamma_star) .* Tcan .* r_rough .* dgamma;
	ans = mean(diff_overall .* (t1+t2+t3));
	return ans;
end



#niterO = 500
#results = zeros(maxiter,5)

total_err_old = 1000;
err_diff = 1000;
itercount = 0;

while (itercount < maxiter) & (err_diff > stop_crit)
	
	itercount = itercount + 1;
	leafvol_rel = 1 .+ k*leafpot;
	vod = leafvol_rel .* (alpha .+ laiM * beta);
	gamma = exp.(-1*vod / mycostheta);

	predH_oi = gamma .* (1 .- rhfac*fH) .* Tsoil + (1-omega) * (1 .- gamma) .* (1 .+ gamma .* (rhfac*fH)) .* Tcan;
	predV_oi = gamma .* (1 .- rhfac*fV) .* Tsoil + (1-omega) * (1 .- gamma) .* (1 .+ gamma .* (rhfac*fV)) .* Tcan;

	dgamma_dk = gamma .* (-1 / mycostheta) .* leafpot;
	dgamma_dalpha = gamma .* (-1 / mycostheta) .* leafvol_rel;
	dgamma_dbeta = gamma .* (-1 / mycostheta) .* leafvol_rel .* laiM;

	diffH = (predH_oi .- obsHi); diffV = (predV_oi .- obsVi);

	grad_k = overall_grad(diffH, rhfac*fH, dgamma_dk, gamma) + overall_grad(diffV, rhfac*fV, dgamma_dk, gamma);
	grad_alpha = overall_grad(diffH, rhfac*fH, dgamma_dalpha, gamma) + overall_grad(diffV, rhfac*fV, dgamma_dalpha, gamma);
	grad_beta = overall_grad(diffH, rhfac*fH, dgamma_dbeta, gamma) + overall_grad(diffV, rhfac*fV, dgamma_dbeta, gamma);


	k = k - step_size*grad_k;
	alpha = alpha - step_size*grad_alpha;
	beta = beta - step_size*grad_beta;

	#k = min(max(0.01,k - step_size*grad_k),0.25);
	#alpha = min(max(0.25,alpha - step_size*grad_alpha),2);
	#beta = min(max(1/100,beta - step_size*grad_beta),1/25);
	
	total_err_new = (sqrt(mean(diffH .^ 2)) + sqrt(mean(diffV .^ 2)))/2;
	err_diff = total_err_old-total_err_new;
	total_err_old = total_err_new*1;

	#results[itercount,:] = [k, alpha, beta, sqrt(mean(diffH .^ 2)), err_diff]
end

        k = min(max(0.00,k),0.75);
        alpha = min(max(0.0,alpha),2);
        beta = min(max(0,beta ),1/25);

	
return [k, alpha, beta, total_err_old, itercount];
end


function optim_vod_stages(obsHi, obsVi, leafpot, smc_mod, simtab, k, alpha, beta, maxiter, stop_crit,step_size, nskip)
	#stage1 = optim_vod_while(obsHi[1:nskip:end], obsVi[1:nskip:end], leafpot[1:nskip:end], smc_mod[1:nskip:end], simtab[1:nskip:end,:], k, alpha, beta, maxiter, stop_crit, step_size);
	stage2 = optim_vod_while(obsHi, obsVi, leafpot, smc_mod, simtab, k,alpha, beta, maxiter, stop_crit, step_size);
	return stage2;
end

####################

#obs_mask = zeros(length(obsH));
#obs_mask[2:(24*3):end] .= 1;
#obs_mask[(2+12):(24*3):end] .= 1;

obs_mask = ones(length(obsH));

obs_mask = obs_mask .== 1;

function log_p_nolabel(a, k0, alpha0, beta0, maxiter, stop_crit, error_var, old_TB)

	v = LVector(logVcmax=a[1], logPcrit=a[2], logKmaxPlant=a[3], logKmaxSoil=a[4], logBsoil=a[5], logP20=a[6], logZsoil=a[7], logNsoil=a[8]);

	p = 0.0;
	p += logpdf(Uniform(log(10), log(120)), v.logVcmax);
	p += logpdf(Uniform(log(0.1), log(100)), v.logKmaxPlant);
	p += logpdf(Uniform(log(0.75), log(5)), v.logPcrit);
        p += logpdf(Uniform(log(1e-5), log(2e-3)), v.logKmaxSoil);
        #p += logpdf(Uniform(log(1), log(8)), v.logBsoil);
        #p += logpdf(Uniform(log(0.25), log(5)), v.logP20);
        p += logpdf(Uniform(log(0.4), log(0.6)), v.logNsoil);
        p += logpdf(Uniform(log(500), log(3000)), v.logZsoil);

		
sim_res = convert_sim(run_sim(FT(exp(v.logVcmax)),FT(exp(v.logPcrit)), FT(exp(v.logKmaxPlant)), FT(exp(v.logKmaxSoil)), FT(2.5), FT(2.0), FT(exp(v.logZsoil)), 
			  FT(exp(v.logNsoil)), istart, N, soil0));

	if isnan(mean(sim_res[1].leafpot))
		return -Inf;

	else
	
	p0 = p+0;
	p0 += 1/36.0*logpdf(MvNormal(old_TB[1][obs_mask], sqrt(error_var)), obsH[obs_mask]);
	p0 += 1/36.0*logpdf(MvNormal(old_TB[2][obs_mask], sqrt(error_var)), obsV[obs_mask]);		


	optim_res = optim_vod_stages(obsH[obs_mask], obsV[obs_mask],
			    mean(sim_res[3],dims=2)[obs_mask], sim_res[2][:,1][obs_mask],
			    sim_res[1][obs_mask,:], k0, alpha0, beta0, maxiter, stop_crit,10^-5, 24);
		sim_TB = get_TB(sim_res, mean(sim_res[3],dims=2), optim_res[1], optim_res[2],optim_res[3]);

		p += 1/36.0*logpdf(MvNormal(sim_TB[1][obs_mask], sqrt(error_var)), obsH[obs_mask]);
		p += 1/36.0*logpdf(MvNormal(sim_TB[2][obs_mask], sqrt(error_var)), obsV[obs_mask]);		


		return p, optim_res, sim_TB, p0
	end
end


function log_p_err(error_var, new_error_var, pred_TB)
	p0 = 0;
	p1 = 0;
	
	p0 += logpdf(truncated(Cauchy(), 0, Inf), error_var);
	p1 += logpdf(truncated(Cauchy(), 0, Inf), new_error_var);

	p0 -= 1/36.0*0.5*(sum( (pred_TB[1]-obsH) .^2 ./ error_var .+ log(error_var)) + sum( (pred_TB[2]-obsV) .^2 ./ error_var .+ log(error_var)))
	p1 -= 1/36.0*0.5*(sum( (pred_TB[1]-obsH) .^2 ./ new_error_var .+ log(new_error_var)) + sum( (pred_TB[2]-obsV) .^2 ./ new_error_var .+ log(new_error_var)))

	return p0, p1
end



function runAMH(x_init, niter, burnlen)

k0 = exp(rand(Uniform(log(0.01),log(0.75))));
alpha0 = exp(rand(Uniform(log(0.25),log(2))));
beta0 = exp(rand(Uniform(log(1/100),log(1/25))));

pars = log.([FT(22),FT(1.5), FT(4.0), FT(1e-4), FT(2.5), FT(2.0), FT(700),FT(0.45)]);
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
 logBsoil=rand(Uniform(log(1), log(8))),
 logP20=rand(Uniform(log(0.25), log(5))),
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
 logBsoil=rand(Uniform(log(1), log(8))),
 logP20=rand(Uniform(log(0.25), log(5))),
 logZsoil=rand(Uniform(log(500), log(3000))),
 logNsoil=rand(Uniform(log(0.4), log(0.6)))
                         );




ll_init, vodpar_init,a1,a2 = log_p_nolabel(Array(x0), k0,alpha0,beta0,10000,10^-6,1.7,(obsH,obsV));
end
end

return x0

end


x01 = init_works();

c1 = runAMH(x01, 10000, 500);

dfc = DataFrame(c1);

CSV.write("post_allnorm_newvar2.csv",dfc);


