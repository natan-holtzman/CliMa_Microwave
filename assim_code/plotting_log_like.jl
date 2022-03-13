
using DataFrames
using CSV
using PyPlot
using Random
using StatsBase
using Statistics

using Distributions
using LabelledArrays
using LinearAlgebra

include("../get_home_dir.jl")

pygui(true)


df_raw = CSV.read(string(PROJECT_HOME,"/data/moflux_land_data_skipyear_hourly2.csv"), DataFrame);
df_raw[!,"RAIN"] *= 2;

avg_et_daily = mean(df_raw.LE/44200)*18.02/1000 * 60*60*24;
avg_rain_daily = mean(df_raw.RAIN)*24;

include(string(PROJECT_HOME,"/simulation_code/sim_vary_new_stomata_med_varyLayer.jl"));
include("time_averaging.jl")
include("mironov.jl")
include("tau_omega_funs.jl")



N = 24*365*1
istart = 24*365*2 + 1; 
soil0 = 0.375;


deltaT = FT(60*60);
#alpha = FT(50);
#alpha = FT(20);

#alpha = FT(5);
#nsoil = FT(1.2);

#alpha = FT(50)#FT(50);
#nsoil = FT(1.6);


alpha = FT(10)#FT(50);
nsoil = FT(1.5);


parnames = ["Kplant","Zsoil","Capacitance","Ksoil",
"Terrain_slope","root_dist_factor",
"Vcmax","P63","g1"];
par_indices = [3,4,6,8,9,10,
				1,5,7 ];

function run_sim_3layer(vcmax_par::FT, k_frac::FT, k_plant::FT,
      z_soil::FT, weibB::FT, vol_factor::FT, g1::FT, k_soil::FT,
      smc_runoff::FT, exp_root_dist::FT,
      canopy_pvslope::FT, trunk_pvslope::FT)
	return run_sim_varyB(vcmax_par, k_frac, weibB, FT(4), k_plant, 
    k_soil, z_soil, istart, N, soil0, vol_factor, 1e-5, 3, df_raw, g1, deltaT, alpha, nsoil,0,0,
    smc_runoff, exp_root_dist,canopy_pvslope, trunk_pvslope);
end

c1 = 1

pars0 = convert(Array{FT}, [40, 0.5, 12, 2000, 8, 1.0, 250,0.4e-6,
         0.35,2,1/20,1/20])
sim_res1 = run_sim_3layer(pars0...);
cm1 = mean(sim_res1[6][:,1:3],dims=2);
pdLWP = Float64.(replace(sim_res1[1].LWP_predawn, missing => NaN));


tsoil = sim_res1[1].T_SOIL .+ 273.15;
tair = sim_res1[1].T_AIR .+ 273.15;
tcan = (sim_res1[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
laiM = sim_res1[1].LAI_modis;

vodA = 0.067;
vodB = 0.82;
vodC = 0.051;

tb1 = get_TB_2(sim_res1[2][:,1], tsoil, tcan, cm1, laiM, vodA, vodB, vodC);


mynormfac = 1/36;
synthetic_ll = zeros(1000)
for ri in 1:1000
	noisyH = tb1[1] + 1.3*randn(length(tb1[1]));
	noisyV = tb1[2] + 1.3*randn(length(tb1[2]));
	obsLL_H = logpdf(MvNormal(tb1[1],1.3),noisyH);
	obsLL_V = logpdf(MvNormal(tb1[2],1.3),noisyV);
	obsLL_total = mynormfac * (obsLL_H + obsLL_V);
	synthetic_ll[ri] = obsLL_total;
end

prior_min = [10, 0.01, 0.1, 500, 0.75,  0.01,100,1e-6,0.01,0.01];
prior_max = [120,0.9,  50, 3000, 10, 10,1000,    5e-4,50.0,5];
prior_dist = Product(Uniform.(log.(prior_min), log.(prior_max)));
constLL = logpdf(prior_dist,log.(pars0[1:10]));

#=
chaindata = CSV.read("assim_code/sample_likelihood/post_par.csv",DataFrame);
chaindata2 = CSV.read("assim_code/sample_likelihood/allobs_c2/post_par.csv",DataFrame);
chaindata3 = CSV.read("assim_code/sample_likelihood/allobs_c3/post_par.csv",DataFrame);


chaindata = CSV.read("assim_code/sample_likelihood/only6/post_par.csv",DataFrame);
chaindata2 = CSV.read("assim_code/sample_likelihood/only6_c2/post_par.csv",DataFrame);
chaindata3 = CSV.read("assim_code/sample_likelihood/only6_c3/post_par.csv",DataFrame);
=#


chaindata1 = CSV.read("assim_code/sample_like_2/all1/post_par.csv",DataFrame);
chaindata2 = CSV.read("assim_code/sample_like_2/all2/post_par.csv",DataFrame);
chaindata3 = CSV.read("assim_code/sample_like_2/all3/post_par.csv",DataFrame);

post_LL = vcat([x[2000:25:end,21] for x in [chaindata1,chaindata2,chaindata3]]...);
post_errvar = vcat([x[2000:25:end,22] for x in [chaindata1,chaindata2,chaindata3]]...);

#=
figure()
hist(synthetic_ll .+ constLL,density=1,label="Theoretical")
hist(post_LL,density=1,alpha=0.5,label="From MCMC")
xlabel("Log likelihood")
ylabel("Density")
legend()
=#

sse_mcmc = (-(post_LL .- constLL) - 243*log.(sqrt.(post_errvar)) .- 243/2*log(2*pi)) .* (2*post_errvar);
sse_theory = (-synthetic_ll .- 243*log.(1.3) .- 243/2*log(2*pi)) .* (2*1.3^2);
sse_mcmc2 = (-(post_LL .- constLL) .- 243*log.(1.3) .- 243/2*log(2*pi)) .* (2*1.3^2);


post_tab = vcat([x[2000:50:end,:] for x in [chaindata1,chaindata2,chaindata3]]...);
#=
figure()
parI = 7
hist(post_tab[:,parI])
xlim(log(prior_min[parI]),log(prior_max[parI]))
=#
prior_mean = (Array(chaindata1[1,1:10]) + Array(chaindata2[1,1:10]) + Array(chaindata3[1,1:10]))/3;
post_mean = mean(Array(post_tab[:,1:10]),dims=1);

post_samp = Array(post_tab[end,1:10]);


sim_res_prior = run_sim_3layer(convert(Array{FT},exp.(prior_mean))...,FT(1/20),FT(1/20));
sim_res_post = run_sim_3layer(convert(Array{FT},exp.(post_mean))...,FT(1/20),FT(1/20));
sim_res_samp = run_sim_3layer(convert(Array{FT},exp.(post_samp))...,FT(1/20),FT(1/20));



cm_prior = mean(sim_res_prior[6][:,1:3],dims=2);
cm_post = mean(sim_res_post[6][:,1:3],dims=2);
cm_samp = mean(sim_res_samp[6][:,1:3],dims=2);

N = 24*365*12
istart = 24*365*0 + 1; 

#=
println("long 1")
sim_res1long = run_sim_3layer(pars0...);
cm1_long = mean(sim_res1long[6][:,1:3],dims=2);
#println("long 2")
#sim_res_post_long = run_sim_3layer(convert(Array{FT},exp.(post_mean))...,FT(1/20),FT(1/20));
println("long 3")
par_samp = convert(Array{FT},exp.(post_samp));
sim_res_samp_long = run_sim_3layer(convert(Array{FT},exp.(post_samp))...,FT(1/20),FT(1/20));

cm2_long = mean(sim_res_samp_long[ 6][:,1:3],dims=2);
=#
#parametrize as stomatal p63 and quantity abs(xylem p63)-abs(stomatal p63), which must be more than 0
psi_list = collect(-8:0.01:0);
wS = exp.(-(psi_list/2) .^ 2);
wX = exp.(-(psi_list/2.5) .^ 2);
#amount of water leaving via ET in this time step is proportional to wS*(psoil - psi_list)
#dV = k dP
#amount of water incoming is proportional to wX
#change in potential is inversely proportional to wX
#runaway will occur when there is a positive feedback between ET and dP
#i.e. if the change in 
function plot_tree(x,y)
	leaf_levels = [12,16,18]
	trunk_level = 9
	root_levels = collect(-8:1:-1)
	figure()
	for i in 1:3
		plot([x[i],x[i+3],x[7]],[leaf_levels[i]+1,leaf_levels[i],trunk_level],
			 "ko-")
	end
	plot([1,1]*x[7],[0,trunk_level],"ko-")
	for i in 1:8
		plot([y[i],x[7+i],x[7]],[root_levels[i],root_levels[i],0],
			 "ko-")
	end
end

#=
figure()
plot(sim_res1[6][:,1:3],color="green",alpha=0.5)
plot(pdLWP,"ko")
plot(sim_res1[6][:,4:6],color="orange",alpha=0.5)
plot(sim_res1[6][:,7],color="blue",alpha=0.5)
plot(sim_res1[6][:,8:15],color="brown",alpha=0.5)

daylist = collect(1:length(sim_res1long[1].VPD))/24;
figure()
plot(daylist, sim_res1long[6][:,1:3],color="green",alpha=0.5)
#plot(pdLWP2,"ko")
plot(daylist, sim_res1long[6][:,4:6],color="orange",alpha=0.5)
plot(daylist, sim_res1long[6][:,7],color="blue",alpha=0.5)
plot(daylist,sim_res1long[6][:,8:15],color="brown",alpha=0.5)

figure()
plot(daylist[1:(24*5):end],get_daily(cm1_long,24*5))
xlabel("Time (days since 1 Jan 2005)")
ylabel("5-day average LWP (MPa)")

ylabel("LWP (MPa)")

=#
#=
pars1 = [3.1927648046462975,-2.6467226953772807,2.099458400464035,7.6636431524168485,2.200646613730597,-1.1787635223907147,5.74270506922163,-14.480768091978758,-1.1286712435746784,0.673368997239951,3.3204521717948343,-3.764499807981834,1.9399104710601602,7.329693071448318,2.3019307796709465,1.2814998731299516,5.5803917979572315,-14.278778878664784,-1.0954914081021205,0.7479192727655555,-826.554645047097,1.6900000000000002,0.05087957867351999,0.802651566774747,0.04858293737916987];
pars1b = convert(Array{FT}, exp.(pars1[1:10]))

post1 = convert(Array,CSV.read("assim_code/sample_output_mar9/postLeaf.csv",DataFrame));
plot(daylist,post1,"r",alpha=0.33)
plot(daylist,cm1_long,"k")
xlim(365,365*2)
ylabel("LWP (MPa)")
xlabel("Time (days since 1 Jan 2005)")
=#
figure()
subplot(2,2,1)
plot(sim_res1[1].ETmod[(244*24):(246*24)])
plot(sim_res1[1].ETmod[(253*24):(255*24)])
ylim(-0.001,0.008)
#xlabel("Time (hours)")
title("Model ET (mol/m2/s)")

subplot(2,2,2)
plot(sim_res1[1].LE[(244*24):(246*24)]/42400)
plot(sim_res1[1].LE[(253*24):(255*24)]/42400)
ylim(-0.001,0.008)
#xlabel("Time (hours)")
title("Flux tower ET (mol/m2/s)")

subplot(2,2,3)
plot(mean(sim_res1[10][(244*24):(246*24),:],dims=2))
plot(mean(sim_res1[10][(253*24):(255*24),:],dims=2))
#xlabel("Time (hours)")
title("Model stomatal\nconductance (mol/m2/s)")

subplot(2,2,4)
plot(cm1[(244*24):(246*24)])
plot(cm1[(253*24):(255*24)])
#xlabel("Time (hours)")
title("Model LWP (MPa)")


figure()
subplot(2,2,1)
plot(mean(sim_res1[8][(244*24):(246*24),:],dims=2))
plot(mean(sim_res1[8][(253*24):(255*24),:],dims=2))
#xlabel("Time (hours)")
title("APAR (umol/m2/s)")

subplot(2,2,2)
plot(sim_res1[1].RelHum[(244*24):(246*24)])
plot(sim_res1[1].RelHum[(253*24):(255*24)])
title("Relative humidity")

subplot(2,2,3)
plot(sim_res1[1].T_AIR[(244*24):(246*24)])
plot(sim_res1[1].T_AIR[(253*24):(255*24)])
title("Temperature (deg C)")

subplot(2,2,4)
plot(sim_res1[1].VPD[(244*24):(246*24)])
plot(sim_res1[1].VPD[(253*24):(255*24)])
#xlabel("Time (hours)")
title("VPD (Pa)")

sval = collect(FT,0.1:0.001:0.55);
scurve = VanGenuchten{FT}(stype = "Ozark",       
                                           α = FT(10),  
                                           n = FT(1.5), 
                                          Θs = 0.55,    
                                          Θr = 0.067);
spotval = [soil_p_25_swc(scurve,x) for x in sval];

figure()
plot(sim_res1[1].SMC,pdLWP,"o")
plot(sval,spotval)
