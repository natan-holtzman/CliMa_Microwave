using DataFrames
using CSV
#using Plots
using Random
using StatsBase
using Zygote
using Optim

include("mironov.jl");

theta = 40*pi/180;
mycostheta = cos(theta);
mysin2theta = sin(theta)^2;

roughness = 0.16
rhfac =  exp(-roughness*mycostheta^2)

slopes0 = ones(24);
intercepts0 = randn(24);

#k0, alpha0, beta0, omega0 = (0.067, 0.82, 0.051, 0.05);


function get_soil_reflect(smc_mod)

eps_real = 0*smc_mod; eps_im = 0*smc_mod;
for i in 1:length(eps_real)
        eps_i = mironov(1.4e9, smc_mod[i], 0.15);
        eps_real[i] = eps_i[1];
        eps_im[i] = eps_i[2];
end

eps = eps_real + eps_im*1im;

fH = abs2.((mycostheta .- sqrt.(eps .- mysin2theta)) ./ (mycostheta .+ sqrt.(eps .- mysin2theta)));
fV = abs2.((eps*mycostheta .- sqrt.(eps .- mysin2theta)) ./ (eps*mycostheta .+ sqrt.(eps .- mysin2theta)));

return rhfac*fH, rhfac*fV;
end

function transmissivity_forward(leafpot,laiM,pvs,alpha,beta)
#leafvol_rel = max.(1 .+ pvs*leafpot,0);
leafvol_rel = 1 .+ pvs*leafpot;
vod = leafvol_rel .* (alpha .+ laiM*beta);
gamma = exp.(-1*vod / mycostheta);
return gamma;
end

function tau_omega_forward(gamma,rH,rV,omega,Tsoil,Tcan)
predH = gamma .* (1 .- rH) .* Tsoil + (1-omega) * (1 .- gamma) .* (1 .+ gamma .* rH) .* Tcan;
predV = gamma .* (1 .- rV) .* Tsoil + (1-omega) * (1 .- gamma) .* (1 .+ gamma .* rV) .* Tcan;
return predH, predV;
end

hourmat = zeros(365*24,24);
for i in 1:24
    hourmat[i:24:end,i] .= 1;
end

function tcan_model_forwardK(Tsoil,slopes,intercepts)
    return hourmat*intercepts + Tsoil .* (hourmat*slopes);
end

function tcan_model_forward_constmean(Tsoil,slopes,intercepts,mask)
    t1 = hourmat[mask,:]*intercepts + Tsoil[mask] .* (hourmat[mask,:]*slopes);
    return t1;# .- mean(t1) .+ mean(Tsoil[mask]);
end



function overall_pred_mask(pvs,alpha,beta,omega,slopes,intercepts,
                            mask,tsoil,canpot1,laiM,rH,rV)
    tcan_pred = tcan_model_forward_constmean(tsoil,slopes,intercepts,mask);
    gamma = transmissivity_forward(canpot1[mask],laiM[mask],pvs,alpha,beta);
    predTB = tau_omega_forward(gamma,rH[mask],rV[mask],omega,tsoil[mask],tcan_pred);
    return predTB;
end


function obs_opt(soilsurf1,canpot1,tsoil,laiM,obs_mask,k0, alpha0, beta0, omega0)
    rH, rV = get_soil_reflect(soilsurf1);


    function get_loss_flat(pars_in)
        pvs,alpha,beta,omega = pars_in[1:4]
        slopes = pars_in[5:(5+23)];
        intercepts = pars_in[(5+23+1):end];
        forward_pred = overall_pred_mask(pvs,alpha,beta,omega,slopes,intercepts,
        obs_mask,tsoil,canpot1,laiM,rH,rV);
        errsH = mean((obsH[obs_mask] - forward_pred[1]) .^ 2);
        errsV = mean((obsV[obs_mask] - forward_pred[2]) .^ 2);
        return errsH+errsV;
    end

    function overall_grad_flat!(G,mypars)
        G[1:end] = gradient(get_loss_flat,mypars)[1];
    end

    obs_pars0 = [k0, alpha0, beta0, omega0,slopes0,intercepts0];
    initpar_flat = vcat(obs_pars0...);

    min2 = Optim.minimizer(optimize(get_loss_flat, overall_grad_flat!, initpar_flat, BFGS()));

    #slopes2 = min2[5:(5+23)];
    #intercepts2 = min2[(5+23+1):end];

    #pred1 = overall_pred(min2[1:4]...,slopes2,intercepts2,obs_mask);
    loss2 = get_loss_flat(min2);
    return min2, loss2;
end

#=
init_radpar = [0.067, 0.82, 0.051, 0.05];
opt6 = obs_opt(soilsurf1,canpot1,tsoil,laiM,mask6,init_radpar...);
opt1 = obs_opt(soilsurf1,canpot1,tsoil,laiM,mask1,init_radpar...);
opt1and6 = obs_opt(soilsurf1,canpot1,tsoil,laiM,mask1and6,init_radpar...);
opt9 = obs_opt(soilsurf1,canpot1,tsoil,laiM,mask9,init_radpar...);
opt0 = obs_opt(soilsurf1,canpot1,tsoil,laiM,mask0,init_radpar...);


tcan_pred6 = tcan_model_forward_constmean(tsoil,opt6[1][5:(5+23)],opt6[1][(5+23+1):end],mask6);
tcan_pred1 = tcan_model_forward_constmean(tsoil,opt1[1][5:(5+23)],opt1[1][(5+23+1):end],mask1);
tcan_pred1and6 = tcan_model_forward_constmean(tsoil,opt1and6[1][5:(5+23)],opt1and6[1][(5+23+1):end],mask1and6);
tcan_pred9 = tcan_model_forward_constmean(tsoil,opt9[1][5:(5+23)],opt9[1][(5+23+1):end],mask9);
tcan_pred0 = tcan_model_forward_constmean(tsoil,opt0[1][5:(5+23)],opt0[1][(5+23+1):end],mask0);



figure()
plot(get_diurnal(tcan,24).-mean(tcan),"o-");
plot([6,6+12],get_diurnal(tcan_pred6,2) .- mean(tcan_pred6),"o-");
plot([1,1+12],get_diurnal(tcan_pred1,2) .- mean(tcan_pred1),"o-");
plot([9,9+12],get_diurnal(tcan_pred9,2) .- mean(tcan_pred9),"o-");
plot([0,0+12],get_diurnal(tcan_pred0,2) .- mean(tcan_pred0),"o-");


figure()
plot(tcan[mask6] .- mean(tcan[mask6]));
plot(tcan_pred6 .- mean(tcan_pred6));

rHall, rVall = get_soil_reflect(soilsurf1);

pars_in = opt1[1];
pvs,alpha,beta,omega = pars_in[1:4];
slopes = pars_in[5:(5+23)];
intercepts = pars_in[(5+23+1):end];

tb1 = overall_pred_mask(pvs,alpha,beta,omega,slopes,intercepts,
                            mask1,tsoil,canpot1,laiM,rHall,rVall);
pars_in = opt6[1];
pvs,alpha,beta,omega = pars_in[1:4];
slopes = pars_in[5:(5+23)];
intercepts = pars_in[(5+23+1):end];

tb6 = overall_pred_mask(pvs,alpha,beta,omega,slopes,intercepts,
                            mask6,tsoil,canpot1,laiM,rHall,rVall);
=#

