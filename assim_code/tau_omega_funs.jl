
using DataFrames
using CSV
#using Plots
using Random
using StatsBase

include("mironov.jl");


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
leafvol_rel[leafvol_rel .< 0] .= 0;
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



function get_TB_2(smc_mod, Tsoil, Tcan, leafpot, laiM, pvs, alpha, beta)

#smc_mod = sim_res[2][:,1];

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
leafvol_rel[leafvol_rel .< 0] .= 0;
#laiM = sim_res[1].LAI_modis;

vod = leafvol_rel .* (alpha .+ laiM*beta);

gamma = exp.(-1*vod / mycostheta);
omega = 0.05
rhfac =  exp(-0.16*mycostheta^2)

#Tsoil = sim_res[1].T_SOIL .+ 273.15;
#Tcan = (sim_res[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
predH = gamma .* (1 .- rhfac*fH) .* Tsoil + (1-omega) * (1 .- gamma) .* (1 .+ gamma .* (rhfac*fH)) .* Tcan;
predV = gamma .* (1 .- rhfac*fV) .* Tsoil + (1-omega) * (1 .- gamma) .* (1 .+ gamma .* (rhfac*fV)) .* Tcan;

return predH[:,1], predV[:,1], vod

end



using LinearAlgebra


function optim_vod_while(obsHi, obsVi, leafpot, smc_mod, simtab, k, alpha, beta, maxiter, stop_crit, step_size)

eps_real = zeros(length(smc_mod)); eps_im = zeros(length(smc_mod));
for i in 1:length(eps_real)
        eps_i = mironov(1.4e9, smc_mod[i], 0.15);
        eps_real[i] = eps_i[1];
        eps_im[i] = eps_i[2];
end

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


function overall_grad(diff_overall, r_rough, dgamma, gamma_star)

	t1 = (1 .- r_rough) .* Tsoil .* dgamma;
	t2 = (1-omega) * (1 .+ gamma_star .* r_rough) .* Tcan .* (-1* dgamma);
	t3 = (1-omega) * (1 .- gamma_star) .* Tcan .* r_rough .* dgamma;
	ans = mean(diff_overall .* (t1+t2+t3));
	return ans;
end

total_err_old = 1000;
err_diff = 1000;
itercount = 0;

while (itercount < maxiter) & (err_diff > stop_crit)
	
	itercount = itercount + 1;
	leafvol_rel = 1 .+ k*leafpot;
	leafvol_rel[leafvol_rel .< 0] .= 0;
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

end

        k = min(max(0.00,k),0.75);
        alpha = min(max(0.0,alpha),2);
        beta = min(max(0,beta ),0.2);

	
return [k, alpha, beta, total_err_old, itercount];
end


function optim_vod_while2(obsHi, obsVi, leafpot, smc_mod, Tsoil, Tcan, laiM, k, alpha, beta, maxiter, stop_crit, step_size)

eps_real = zeros(length(smc_mod)); eps_im = zeros(length(smc_mod));
for i in 1:length(eps_real)
        eps_i = mironov(1.4e9, smc_mod[i], 0.15);
        eps_real[i] = eps_i[1];
        eps_im[i] = eps_i[2];
end

eps = complex.(eps_real, eps_im);


theta = 40*pi/180;
mycostheta = cos(theta);
mysin2theta = sin(theta)^2;

fH = abs2.((mycostheta .- sqrt.(eps .- mysin2theta)) ./ (mycostheta .+ sqrt.(eps .- mysin2theta)));
fV = abs2.((eps*mycostheta .- sqrt.(eps .- mysin2theta)) ./ (eps*mycostheta .+ sqrt.(eps .- mysin2theta)));

omega = 0.05;
rhfac =  exp(-0.16*mycostheta^2);

#Tsoil = simtab.T_SOIL .+ 273.15;
#Tcan = (simtab.LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
#laiM = simtab.LAI_modis;


function overall_grad(diff_overall, r_rough, dgamma, gamma_star)

	t1 = (1 .- r_rough) .* Tsoil .* dgamma;
	t2 = (1-omega) * (1 .+ gamma_star .* r_rough) .* Tcan .* (-1* dgamma);
	t3 = (1-omega) * (1 .- gamma_star) .* Tcan .* r_rough .* dgamma;
	ans = mean(diff_overall .* (t1+t2+t3));
	return ans;
end

total_err_old = 1000;
err_diff = 1000;
itercount = 0;

while (itercount < maxiter) & (err_diff > stop_crit)
	
	itercount = itercount + 1;
	leafvol_rel = 1 .+ k*leafpot;
	leafvol_rel[leafvol_rel .< 0] .= 0;
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

end

        k = min(max(0.00,k),0.75);
        alpha = min(max(0.0,alpha),2);
        beta = min(max(0,beta ),0.2);

	
return [k, alpha, beta, total_err_old, itercount];
end

function optim_vod_stages(obsHi, obsVi, leafpot, smc_mod, simtab, k, alpha, beta, maxiter, stop_crit,step_size, nskip)
	stage1 = optim_vod_while(obsHi[1:nskip:end], obsVi[1:nskip:end], leafpot[1:nskip:end], smc_mod[1:nskip:end], simtab[1:nskip:end,:], k, alpha, beta, maxiter, stop_crit, step_size);
	k1, alpha1, beta1, e1, i1 = stage1;
	stage2 = optim_vod_while(obsHi, obsVi, leafpot, smc_mod, simtab, k1, alpha1, beta1, maxiter, stop_crit, step_size);
	return stage2;
end



function optim_vod_stages2(obsHi, obsVi, leafpot, smc_mod, Tsoil, Tcan, laiM, k, alpha, beta, maxiter, stop_crit,step_size, nskip)
	stage1 = optim_vod_while2(obsHi[1:nskip:end], obsVi[1:nskip:end], leafpot[1:nskip:end], smc_mod[1:nskip:end], Tsoil[1:nskip:end], Tcan[1:nskip:end], laiM[1:nskip:end], k, alpha, beta, maxiter, stop_crit, step_size);
	k1, alpha1, beta1, e1, i1 = stage1;
	stage2 = optim_vod_while2(obsHi, obsVi, leafpot, smc_mod, Tsoil, Tcan, laiM, k1, alpha1, beta1, maxiter, stop_crit, step_size);
	return stage2;
end