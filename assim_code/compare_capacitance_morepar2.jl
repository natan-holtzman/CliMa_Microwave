include("../get_home_dir.jl")

using DataFrames
using CSV
using PyPlot
using Random
using StatsBase
using Statistics

df_raw = CSV.read(string(PROJECT_HOME,"/data/moflux_land_data_skipyear_hourly2.csv"), DataFrame);
df_raw[!,"RAIN"] *= 2;

include(string(PROJECT_HOME,"/simulation_code/sim_vary_new_stomata_med_varyLayer.jl"));
include("time_averaging.jl")
include("mironov.jl")
include("tau_omega_funs.jl")



N = 24*365*1
istart = 24*365*2 + 1; 
soil0 = 0.30;


deltaT = FT(60*60);
alpha = FT(50);
#alpha = FT(20);

alpha = FT(5);
nsoil = FT(1.2);
#nsoil = FT(1.4);


parnames = ["vcmax","Kplant","Zsoil","p63","Capacitance","g1","Ksoil",
"Terrain_slope","root_dist_factor"];
par_indices = [1,3,4,5,6,7,8,9,10 ];

function run_sim_3layer(vcmax_par::FT, k_frac::FT, k_plant::FT,
      z_soil::FT, weibB::FT, vol_factor::FT, g1::FT, k_soil::FT,
      slope_runoff::FT, exp_root_dist::FT,
      canopy_pvslope::FT, trunk_pvslope::FT)
	return run_sim_varyB(vcmax_par, k_frac, weibB, FT(2), k_plant, 
    k_soil, z_soil, istart, N, soil0, vol_factor, 1e-5, 3, df_raw, g1, deltaT, alpha, nsoil,0,0,
    slope_runoff, exp_root_dist,canopy_pvslope, trunk_pvslope);
end

c1 = 1

pars0 = convert(Array{FT}, [60, 0.5, 12, 2000, 5, 0.4, 250,1e-4,
         0.05,2,1/20,1/5])
sim_res1 = run_sim_3layer(pars0...);

#factor_list = [1/3,1/2,1/1.5,1.5,2,3];

factor_list = [1/2,2]

sim_list = [];

#push all of these to sim_list


#effect of doubling pv curve slope for trunk is same as effect of doubling storage volume

for parI in 1:length(par_indices)
    println(parnames[parI])
    for factorI in 1:length(factor_list)
        newpars = pars0*1;
        newpars[par_indices[parI]] *= factor_list[factorI]
        simIJ = run_sim_3layer(newpars...);
        push!(sim_list,simIJ);
    end
end

cm1 = mean(sim_res1[6][:,1:3],dims=2);
cmx = [mean(x[6][:,1:3],dims=2) for x in sim_list];

pdLWP = Float64.(replace(sim_res1[1].LWP_predawn, missing => NaN));


tsoil = sim_res1[1].T_SOIL .+ 273.15;
tair = sim_res1[1].T_AIR .+ 273.15;
tcan = (sim_res1[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
laiM = sim_res1[1].LAI_modis;

using Distributions
using LabelledArrays
using LinearAlgebra

vodA = 0.067;
vodB = 0.82;
vodC = 0.051;

tb1 = get_TB_2(sim_res1[2][:,1], tsoil, tcan, cm1, laiM, vodA, vodB, vodC);


function optim_vod_mask(obsHi, obsVi, leafpot, smc_mod, Tsoil, Tcan, laiM, k, alpha, beta, maxiter, stop_crit,step_size, nskip, mask)
	return optim_vod_stages2(obsHi[mask], obsVi[mask], leafpot[mask], smc_mod[mask], Tsoil[mask], Tcan[mask], laiM[mask], k, alpha, beta, maxiter, stop_crit,step_size, nskip)
end


hodlist = (1:length(tb1[1])) .% 24;

mask6 = ((hodlist .== 7) + (hodlist .== (7+12))) .== 1;
mask1 = ((hodlist .== 2) + (hodlist .== (2+12))) .== 1;

tbx_all = []
tbx_6 = []
tbx_1 = []

for i in 1:length(sim_list)
	println(i)
	sim_res2 = sim_list[i];
	cm2 = cmx[i]
	ores2_all = optim_vod_stages2(tb1[1], tb1[2], cm2, sim_res2[2][:,1], tsoil, tcan, laiM, vodA, vodB, vodC, 10^5, 1e-6, 1e-5, 24);
	ores2_6 = optim_vod_mask(tb1[1], tb1[2], cm2, sim_res2[2][:,1], tsoil, tcan, laiM, vodA, vodB, vodC, 10^5, 1e-6, 1e-5, 2, mask6);
	ores2_1 = optim_vod_mask(tb1[1], tb1[2], cm2, sim_res2[2][:,1], tsoil, tcan, laiM, vodA, vodB, vodC, 10^5, 1e-6, 1e-5, 2, mask1);

	tb2_all = get_TB_2(sim_res2[2][:,1], tsoil, tcan, cm2, laiM, ores2_all[1:3]...);
	tb2_1 = get_TB_2(sim_res2[2][:,1], tsoil, tcan, cm2, laiM, ores2_1[1:3]...);
	tb2_6 = get_TB_2(sim_res2[2][:,1], tsoil, tcan, cm2, laiM, ores2_6[1:3]...);
	
	push!(tbx_all,tb2_all)
	push!(tbx_6,tb2_1)
	push!(tbx_1,tb2_6)
end

function cat_temp(x)
	return vcat(x[1],x[2])
end

tdiff_all = [sqrt.(mean((cat_temp(tb1)-cat_temp(x)) .^ 2)) for x in tbx_all];
tdiff_6 = [sqrt.(mean((cat_temp(tb1)-cat_temp(x))[vcat(mask6,mask6)] .^ 2)) for x in tbx_6];
tdiff_1 = [sqrt.(mean((cat_temp(tb1)-cat_temp(x))[vcat(mask1,mask1)] .^ 2)) for x in tbx_1];


significant_rmse_1yr = 1.3*sqrt(quantile.(Chisq.(365), 0.95)/(365) - 1);
significant_rmse_10yr = 1.3*sqrt(quantile.(Chisq.(365*10), 0.95)/(365*10) - 1);

npar = length(par_indices)
nfac = length(factor_list)

factor_list4 = repeat(factor_list, outer=npar);

pygui(true)


ETdiffD = [sqrt.(mean(get_daily(sim_res1[1].ETmod - x[1].ETmod, 24) .^ 2)) for x in sim_list]*18/1000*60*60*24;
GPPdiffD = [sqrt.(mean(get_daily(sim_res1[1].GPP - x[1].GPP, 24) .^ 2)) for x in sim_list];
LWPdiffD = [sqrt.(mean(get_daily(cm1 - x, 24) .^ 2)) for x in cmx];

DstdET = std(get_daily(sim_res1[1].ETmod,24))*18/1000*60*60*24;
DstdGPP = std(get_daily(sim_res1[1].GPP,24));
DstdLWP = std(get_daily(cm1,24));
DstdSMC = std(get_daily(sim_res1[1].ColumnSMC,24));


SMCdiffD = [sqrt.(mean(get_daily(sim_res1[1].ColumnSMC - x[1].ColumnSMC, 24) .^ 2)) for x in sim_list];

sSMCdiffD = [sqrt.(mean(get_daily(sim_res1[2][:,1] - x[2][:,1], 24) .^ 2)) for x in sim_list];

figure()
plot(tdiff_all,tdiff_1,"o",label="1 AM/PM")
plot(tdiff_all,tdiff_6,"o",label="6 AM/PM")
plot([0,1],[0,1],"k--")
xlabel("Full diurnal RMSE (K)")
ylabel("Twice-daily RMSE (K)")
legend()


function diurnal_pattern(x)
	z1 = get_diurnal(x,24);
	return (z1 .- mean(z1)) / std(z1);
end
shapeLWPdiff = [sqrt.(mean((diurnal_pattern(cm1) - diurnal_pattern(x)) .^ 2)) for x in cmx];

all_diurnal = [diurnal_pattern(x) for x in cmx];

figure()
plot(hcat(all_diurnal...),alpha=0.25)
plot(diurnal_pattern(cm1),"k")
xlabel("Time of day")
ylabel("Normalized LWP")


figure()
plot(tdiff_all, 1 .- LWPdiffD / DstdLWP,"o",label="LWP")
plot(tdiff_all, 1 .- ETdiffD / DstdET,"o",label="ET")
plot(tdiff_all, 1 .- SMCdiffD / DstdSMC,"o",label="SMC")
plot([0],[1],"ko")
xlabel("RMSE of Tb (K)")
ylabel("R2 of model state")
plot([1,1]*significant_rmse_10yr,[-2,1],"k--")
plot([1,1]*significant_rmse_1yr,[-2,1],"k--")
plot([0,maximum(tdiff_all)], [0,0],"k")
text(significant_rmse_10yr+0.01,-1,"Significant\nover 10 years",rotation="vertical")
text(significant_rmse_1yr+0.01,-1,"Significant\nover 1 year",rotation="vertical")
legend()

figure()
bar(1:npar, 0.5*(tdiff_1[1:2:end]+tdiff_1[2:2:end]))
bar(1:npar, 0.5*(tdiff_all[1:2:end]+tdiff_all[2:2:end]),alpha=0.5)
xticks(1:npar,parnames)

figure()
bar(1:npar, 0.5*(LWPdiffD[1:2:end]+LWPdiffD[2:2:end]))
xticks(1:npar,parnames)

figure()
bar(1:npar, 0.5*(ETdiffD[1:2:end]+ETdiffD[2:2:end]))
xticks(1:npar,parnames)