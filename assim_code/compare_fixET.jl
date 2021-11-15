include("../get_home_dir.jl")

using DataFrames
using CSV
#using Plots
using Random
using StatsBase

include("mironov.jl");
include("tau_omega_funs.jl");
include("time_averaging.jl")

df_raw = avg2_df(CSV.read("../data/moflux_land_data_newnames_7.csv", DataFrame));
df_raw[!,"RAIN"] *= 2;

N = 24*365
istart = 24*365*2 - 52*24 #+ 230*48
soil0 = 0.39;

include("../simulation_code/sim_vary_new_stomata_med_fixET.jl");

deltaT = FT(60*60);
alpha = FT(1.368);
nsoil = FT(2.6257);
#m = 1 - 1/nsoil;
#alpha = (FT(0.1) ^ (-1/m) - 1) ^(1/nsoil)/pot1;

function run_sim_2(vcmax_par::FT, k_plant::FT, k_soil::FT, z_soil::FT, weibB::FT, weibC::FT, vol_factor::FT, alpha_mult::FT=FT(1))
#vcmax_par, k_weibB, k_weibC, k_plant, k_soil, z_soil, istart, N, smc0, storage_mult, buffrate,scheme_number,df_raw,deltaT,alpha,n)
	return run_sim_vary(vcmax_par, weibB, weibC, k_plant, k_soil, z_soil, istart, N, soil0, vol_factor, 1e-5, 3, df_raw, deltaT, alpha*alpha_mult, nsoil);
end

baseK = 2
baseC = 1

sim_res1 = run_sim_2(FT(22), FT(baseK),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC));
#sim_resCbig = run_sim_2(FT(22), FT(baseK),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC*2));
#sim_resCsmall = run_sim_2(FT(22), FT(baseK),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC/2));

#sim_resKbig = run_sim_2(FT(22), FT(baseK*2),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC));
#sim_resKsmall = run_sim_2(FT(22), FT(baseK/2),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC));
#sim_resKsmall2 = run_sim_2(FT(22), FT(baseK/2),FT(1e-5), FT(2000),FT(4.0),FT(2),FT(baseC));

#sim_resCsmall_Kbig = run_sim_2(FT(22), FT(baseK*1.25),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC/2));
#sim_resCbig_Ksmall = run_sim_2(FT(22), FT(baseK/1.25),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC*2));

#sim_resCsmall_Kbig2 = run_sim_2(FT(22), FT(baseK*1.25),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC/2),FT(0.9));
#sim_resCbig_Ksmall2 = run_sim_2(FT(22), FT(baseK/1.25),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC*2),FT(1.1));

#sim_resKbig2 = run_sim_2(FT(22), FT(baseK*2),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC),FT(0.85));
#sim_resKbig3 = run_sim_2(FT(22), FT(baseK*2),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC/8),FT(0.85));

#=
plot(get_diurnal(sim_res1[1].leafpot,24),"k", label="Original")
plot(get_diurnal(sim_resKbig[1].leafpot,24),color="orange",label="Adjust K")
plot(get_diurnal(sim_resKbig2[1].leafpot,24),"--",color="orange",label="Adjust K, alpha")
plot(get_diurnal(sim_resKbig3[1].leafpot,24),"--",color="blue",label="Adjust K, alpha, C")
xlabel("Time of day")
ylabel("LWP (MPa)")
legend()

leafpot_ratio = mean(sim_resKbig[1].leafpot[1:1:end]) / mean(sim_res1[1].leafpot[1:1:end]);
mydiffs = sim_resKbig[1].leafpot/leafpot_ratio - sim_res1[1].leafpot;
rmse_list = zeros(24);
for i in 1:24
	rmse_list[i] = sqrt(mean(mydiffs[i:24:end] .^ 2)); 
end
=#

sim_resKbig = run_sim_2(FT(22), FT(baseK*1.5),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC));

leafpot_ratio = mean(sim_resKbig[3][1:end]) / mean(sim_res1[3][1:end]);

sim_TB1 = get_TB(sim_res1, mean(sim_res1[3],dims=2), 0.25, 1, 0.05);
sim_TB_Kbig = get_TB(sim_resKbig, mean(sim_resKbig[3],dims=2), 0.25/leafpot_ratio, 1, 0.05);

mydiffs = sim_TB1[1] - sim_TB_Kbig[1];
rmse_list = zeros(24);
for i in 1:24
	rmse_list[i] = sqrt(mean(mydiffs[i:24:end] .^ 2)); 
end

tsoil = sim_res1[1].T_SOIL .+ 273.15;
tair = sim_res1[1].T_AIR .+ 273.15;
tcan = (sim_res1[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;

tdiff = tcan - tair;
tdiff_reshape = reshape(tdiff,(24,:));
tdiff_day = mean(tdiff_reshape[6:18,:],dims=1)[1,:];

et_flux = sim_res1[1].LE/44200;
et_norm = et_flux ./ sim_res1[1].LAI_modis;

et_reshape = reshape(et_flux,(24,:));
et_day = mean(et_reshape[6:18,:],dims=1)[1,:];

et_norm_reshape = reshape(et_norm,(24,:));
et_norm_day = mean(et_norm_reshape[6:18,:],dims=1)[1,:];

soilsurf = sim_res1[2][:,1];
canpot = mean(sim_res1[3],dims=2);
laiM = sim_res1[1].LAI_modis;
#smc_mod, Tsoil, Tcan, leafpot, laiM, pvs, alpha, beta)
trueTB = get_TB_2(soilsurf, tsoil, tcan, canpot, laiM, 0.1, 1, 0.05);
sim_TB_air = get_TB_2(soilsurf, tsoil, tair, canpot, laiM, 0.1, 1, 0.05);
#sim_TBc = get_TB_2(soilsurf, tsoil, tcan, canpot, laiM, 0.1, 1, 0.05);

noise_std = 1.3;
noise_var = noise_std^2

obsH = trueTB[1] + noise_std*randn(Float64, length(trueTB[1]));
obsV = trueTB[2] + noise_std*randn(Float64, length(trueTB[1]));

optim_res_good = optim_vod_stages2(obsH, obsV, canpot, soilsurf, tsoil, tcan, laiM, 0.1, 1, 0.05, 10000, 1e-6, 1e-5, 1);

obs_mask = zeros(length(obsH));
obs_mask[2:(24*3):end] .= 1;
obs_mask[(2+12):(24*3):end] .= 1;
obs_mask = obs_mask .== 1;

optim_res_good_1 = optim_vod_stages2(obsH[obs_mask], obsV[obs_mask], canpot[obs_mask], soilsurf[obs_mask],
									tsoil[obs_mask], tcan[obs_mask], laiM[obs_mask], 0.1, 1, 0.05, 10000, 1e-6, 1e-5, 1);
optim_res_air = optim_vod_stages2(obsH, obsV, canpot, soilsurf, tsoil, tair, laiM, 0.1, 1, 0.05, 10000, 1e-6, 1e-5, 1);

optim_res_air_1 = optim_vod_stages2(obsH[obs_mask], obsV[obs_mask], canpot[obs_mask], soilsurf[obs_mask],
									tsoil[obs_mask], tair[obs_mask], laiM[obs_mask], 0.1, 1, 0.05, 10000, 1e-6, 1e-5, 1);

#=
plot(get_diurnal(sim_resCbig_Ksmall2[1].leafpot,24),color="purple")

plot(get_diurnal(sim_resKbig[1].leafpot,24),"r")
plot(get_diurnal(sim_resKsmall[1].leafpot,24),"b")

plot(get_diurnal(sim_resCbig[1].leafpot,24),"r--")
plot(get_diurnal(sim_resCsmall[1].leafpot,24),"b--")
=#

#now do without adjusting soil but with adjusting Tb

#=
sim_res1 = run_sim_2(FT(22), FT(baseK),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC));
sim_resKbig = run_sim_2(FT(22), FT(baseK*2),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC));
sim_resKbig_Csmall = run_sim_2(FT(22), FT(baseK*2),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC/2));


#sim_resKbig2 = run_sim_2(FT(22), FT(baseK*2),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC),FT(0.85));
#sim_resKbig3 = run_sim_2(FT(22), FT(baseK*2),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(baseC/8),FT(0.85));

kmean_ratio = mean(sim_resKbig[3]) / mean(sim_res1[3])
kmean_ratio2 = mean(sim_resKbig_Csmall[3]) / mean(sim_res1[3])

sim_TB1 = get_TB(sim_res1, mean(sim_res1[3],dims=2), 0.25, 1, 0.05);
sim_TB_Kbig = get_TB(sim_resKbig, mean(sim_resKbig[3],dims=2), 0.25, 1, 0.05);
sim_TB_Kbig2 = get_TB(sim_resKbig, mean(sim_resKbig[3],dims=2), 0.25/kmean_ratio, 1, 0.05);
sim_TB_Kbig3 = get_TB(sim_resKbig_Csmall, mean(sim_resKbig_Csmall[3],dims=2), 0.25/kmean_ratio2, 1, 0.05);

plot(get_diurnal(sim_TB1[3],24),"k", label="Original")
plot(get_diurnal(sim_TB_Kbig[3],24),color="orange",label="Adjust K")
plot(get_diurnal(sim_TB_Kbig2[3],24),"--",color="orange",label="Adjust K, alpha")
plot(get_diurnal(sim_TB_Kbig3[3],24),"--",color="blue",label="Adjust K, alpha, C")
xlabel("Time of day")
ylabel("VOD")
legend()
=#
