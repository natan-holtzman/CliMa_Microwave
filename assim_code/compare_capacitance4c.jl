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
istart = 24*365*3 + 1; 
soil0 = 0.30;


deltaT = FT(60*60);
alpha = FT(50);
#alpha = FT(20);

alpha = FT(5);
nsoil = FT(1.2);
#nsoil = FT(1.4);

function run_sim_3layer(vcmax_par::FT, k_frac::FT, k_plant::FT,  z_soil::FT, weibB::FT, vol_factor::FT, g1::FT)
	return run_sim_varyB(vcmax_par, k_frac, weibB, FT(2), k_plant, FT(1e-4), z_soil, istart, N, soil0, vol_factor, 1e-5, 3, df_raw, g1, deltaT, alpha, nsoil,0,0);
end

c1 = 1


sim_res1 = run_sim_3layer(FT(60),FT(0.5), FT(12), FT(2000),FT(5),FT(0.4),FT(250));

factor_list = [1/3,1/2,1/1.5,1.5,2,3];

sim_list = [];

#push all of these to sim_list

parnames = ["g1","p63","Kplant","Capacitance","Zsoil"];

println("Set 0 of 4")
#change g1
partial_list = [run_sim_3layer(FT(60),FT(0.5), FT(12), FT(2000),FT(5),FT(0.4),FT(250*factor_list[i])) for i in 1:length(factor_list)];
append!(sim_list, partial_list);

println("Set 1 of 4")
#change stomatal and xylem p50
partial_list = [run_sim_3layer(FT(60),FT(0.5), FT(12), FT(2000),FT(5*factor_list[i]),FT(0.4),FT(250)) for i in 1:length(factor_list)];
append!(sim_list, partial_list);

println("Set 2 of 4")
#change conductance
partial_list = [run_sim_3layer(FT(60),FT(0.5), FT(12*factor_list[i]), FT(2000),FT(5),FT(0.4),FT(250)) for i in 1:length(factor_list)];
append!(sim_list, partial_list);

println("Set 3 of 4")
#change cap
partial_list = [run_sim_3layer(FT(60),FT(0.5), FT(12), FT(2000),FT(5),FT(0.4*factor_list[i]),FT(250)) for i in 1:length(factor_list)];
append!(sim_list, partial_list);


#change cap and cond
#sim_list = [run_sim_3layer(FT(60),FT(0.5), FT(12*factor_list[i]), FT(2000),FT(5),FT(0.4*factor_list[i]),FT(250)) for i in 1:length(factor_list)];

println("Set 4 of 4")
#change soil depth
partial_list = [run_sim_3layer(FT(60),FT(0.5), FT(12), FT(2000*factor_list[i]),FT(5),FT(0.4),FT(250)) for i in 1:length(factor_list)];
append!(sim_list, partial_list);

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

factor_list4 = repeat(factor_list, outer=5);

significant_rmse_1yr = 1.3*sqrt(quantile.(Chisq.(365), 0.95)/(365) - 1);
significant_rmse_10yr = 1.3*sqrt(quantile.(Chisq.(365*10), 0.95)/(365*10) - 1);


figure()

subplot(1,3,1)
plot(reshape(factor_list4,(6,5)),reshape(tdiff_all,(6,5)),"o-",label=parnames)
#xlabel("Parameter multiplier")
#ylabel("Tb RMSE (K)")
xscale("log")
ylim(0,0.9)
title("Full diurnal")
legend()
xticks(factor_list,round.(factor_list,digits=2))
plot([1/3,3],significant_rmse_1yr*[1,1],"k--")
plot([1/3,3],significant_rmse_10yr*[1,1],"k--")


subplot(1,3,2)
plot(reshape(factor_list4,(6,5)),reshape(tdiff_1,(6,5)),"o-",label=parnames)
xscale("log")
ylim(0,0.9)
title("1 AM/PM")
xticks(factor_list,round.(factor_list,digits=2))
plot([1/3,3],significant_rmse_1yr*[1,1],"k--")
plot([1/3,3],significant_rmse_10yr*[1,1],"k--")

subplot(1,3,3)
plot(reshape(factor_list4,(6,5)),reshape(tdiff_6,(6,5)),"o-",label=parnames)
xscale("log")
ylim(0,0.9)
title("6 AM/PM")
xticks(factor_list,round.(factor_list,digits=2))
plot([1/3,3],significant_rmse_1yr*[1,1],"k--")
plot([1/3,3],significant_rmse_10yr*[1,1],"k--")



figure()
plot(factor_list4, tdiff_all,"o",label="All")
plot(factor_list4, tdiff_6,"o",label="6 AM/PM")
plot(factor_list4, tdiff_1,"o",label="1 AM/PM")
plot([1],[0],"ko")

significant_rmse_1yr = 1.3*sqrt(quantile.(Chisq.(365), 0.95)/(365) - 1);
plot([factor_list[1],factor_list[end]],[1,1]*significant_rmse_1yr,"k--")


significant_rmse_10yr = 1.3*sqrt(quantile.(Chisq.(365*10), 0.95)/(365*10) - 1);
plot([factor_list[1],factor_list[end]],[1,1]*significant_rmse_10yr,"k--")

xlabel("Parameter change")
ylabel("Tb RMSE")
xscale("log")


#npoint_list = 120:10:(360*10);
#chilist = [1.3*sqrt(quantile.(Chisq.(i), 0.95)/i - 1) for i in npoint_list];


ETdiff = [sqrt.(mean((sim_res1[1].ETmod - x[1].ETmod) .^ 2)) for x in sim_list]*18/1000*60*60*24;
GPPdiff = [sqrt.(mean((sim_res1[1].GPP - x[1].GPP) .^ 2)) for x in sim_list];
LWPdiff = [sqrt.(mean((cm1 - x) .^ 2)) for x in cmx];
LWPcor = [cor(cm1,x)[1,1] for x in cmx];



stdET = std(sim_res1[1].ETmod)*18/1000*60*60*24;
stdGPP = std(sim_res1[1].GPP);
stdLWP = std(cm1);

figure()
plot(factor_list4, 1 .- ETdiff ./ stdET, "o", label="ET")
plot(factor_list4, 1 .- GPPdiff ./ stdGPP, "o", label="GPP")
plot(factor_list4, 1 .- LWPdiff ./ stdLWP, "o", label="LWP")
plot([1],[1],"ko")
xlabel("Parameter change")
ylabel("R2")
legend()
title("Hourly")
xscale("log")



ETdiffD = [sqrt.(mean(get_daily(sim_res1[1].ETmod - x[1].ETmod, 24) .^ 2)) for x in sim_list]*18/1000*60*60*24;
GPPdiffD = [sqrt.(mean(get_daily(sim_res1[1].GPP - x[1].GPP, 24) .^ 2)) for x in sim_list];
LWPdiffD = [sqrt.(mean(get_daily(cm1 - x, 24) .^ 2)) for x in cmx];

DstdET = std(get_daily(sim_res1[1].ETmod,24))*18/1000*60*60*24;
DstdGPP = std(get_daily(sim_res1[1].GPP,24));
DstdLWP = std(get_daily(cm1,24));
DstdSMC = std(get_daily(sim_res1[1].ColumnSMC,24));


figure()
plot(factor_list4, 1 .- ETdiffD ./ DstdET, "o", label="ET")
plot(factor_list4, 1 .- GPPdiffD ./ DstdGPP, "o", label="GPP")
plot(factor_list4, 1 .- LWPdiffD ./ DstdLWP, "o", label="LWP")
plot([1],[1],"ko")
xlabel("Parameter change")
ylabel("R2")
legend()
title("Daily")
xscale("log")

figure()
plot(tdiff_all, LWPdiffD,"o",label="Full diurnal")
#plot(tdiff_1, LWPdiffD,"o",label="1 AM/PM")
#plot(tdiff_6, LWPdiffD,"o",label="6 AM/PM")
plot([0],[0],"ko")
xlabel("RMSE of Tb (K)")
ylabel("RMSE of LWP (MPa)")
plot([1,1]*significant_rmse_10yr,[0,maximum(LWPdiffD)],"k--")
plot([1,1]*significant_rmse_1yr,[0,maximum(LWPdiffD)],"k--")
text(significant_rmse_10yr+0.01,maximum(LWPdiffD)*0.75,"Significant\nover 10 years",rotation="vertical")
text(significant_rmse_1yr+0.01,maximum(LWPdiffD)*0.75,"Significant\nover 1 year",rotation="vertical")
#legend()

figure()
plot(tdiff_all, ETdiffD,"o",label="Full diurnal")
#plot(tdiff_1, ETdiffD,"o",label="1 AM/PM")
#plot(tdiff_6, ETdiffD,"o",label="6 AM/PM")
plot([0],[0],"ko")
xlabel("RMSE of Tb (K)")
ylabel("RMSE of ET (mm/day)")
plot([1,1]*significant_rmse_10yr,[0,maximum(ETdiffD)],"k--")
plot([1,1]*significant_rmse_1yr,[0,maximum(ETdiffD)],"k--")
text(significant_rmse_10yr+0.01,maximum(ETdiffD)*0.75,"Significant\nover 10 years",rotation="vertical")
text(significant_rmse_1yr+0.01,maximum(ETdiffD)*0.75,"Significant\nover 1 year",rotation="vertical")
#legend()


SMCdiffD = [sqrt.(mean(get_daily(sim_res1[1].ColumnSMC - x[1].ColumnSMC, 24) .^ 2)) for x in sim_list];

sSMCdiffD = [sqrt.(mean(get_daily(sim_res1[2][:,1] - x[2][:,1], 24) .^ 2)) for x in sim_list];

myY = SMCdiffD;

figure()
plot(tdiff_all, myY,"o",label="Full diurnal")
#plot(tdiff_1, myY,"o",label="1 AM/PM")
#plot(tdiff_6, myY,"o",label="6 AM/PM")
plot([0],[0],"ko")
xlabel("RMSE of Tb (K)")
ylabel("RMSE of Column SMC")
plot([1,1]*significant_rmse_10yr,[0,maximum(myY)],"k--")
plot([1,1]*significant_rmse_1yr,[0,maximum(myY)],"k--")
text(significant_rmse_10yr+0.01,maximum(myY)*0.75,"Significant\nover 10 years",rotation="vertical")
text(significant_rmse_1yr+0.01,maximum(myY)*0.75,"Significant\nover 1 year",rotation="vertical")
#legend()

figure()
plot(tdiff_all,tdiff_1,"o",label="1 AM/PM")
plot(tdiff_all,tdiff_6,"o",label="6 AM/PM")
plot([0,1],[0,1],"k--")
xlabel("Full diurnal RMSE (K)")
ylabel("Twice-daily RMSE (K)")
legend()



#ETdiff_obs = [sqrt.(mean((sim_res1[1].ETmod - sim_res1[1].LE/44200) .^ 2)) for x in sim_list]*18/1000*60*60;

figure()

subplot(1,3,1)
plot(reshape(factor_list4,(6,5)),reshape(LWPdiffD,(6,5)),"o-",label=parnames)
#xlabel("Parameter multiplier")
#ylabel("Tb RMSE (K)")
xscale("log")
#ylim(0,0.9)
title("LWP (MPa)")
legend()
xticks(factor_list,round.(factor_list,digits=2))

subplot(1,3,2)
plot(reshape(factor_list4,(6,5)),reshape(ETdiffD,(6,5)),"o-",label=parnames)
xscale("log")
title("ET (mm/day)")
xticks(factor_list,round.(factor_list,digits=2))

subplot(1,3,3)
plot(reshape(factor_list4,(6,5)),reshape(SMCdiffD,(6,5)),"o-",label=parnames)
xscale("log")
title("Column SMC")
xticks(factor_list,round.(factor_list,digits=2))

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

subplot(1,2,1)
plot(reshape(factor_list4,(6,5)),reshape(LWPdiffD,(6,5)),"o-",label=parnames)
#xlabel("Parameter multiplier")
#ylabel("Tb RMSE (K)")
xscale("log")
#ylim(0,0.9)
title("Daily LWP (MPa)")
legend()
xticks(factor_list,round.(factor_list,digits=2))

subplot(1,2,2)
plot(reshape(factor_list4,(6,5)),reshape(shapeLWPdiff,(6,5)),"o-",label=parnames)
#xlabel("Parameter multiplier")
#ylabel("Tb RMSE (K)")
xscale("log")
#ylim(0,0.9)
title("LWP diurnal shape")
xticks(factor_list,round.(factor_list,digits=2))




tbx_H = []

for i in 1:length(sim_list)
	println(i)
	sim_res2 = sim_list[i];
	cm2 = cmx[i]

	tb2_all = get_TB_2(sim_res2[2][:,1], tsoil, tcan, cm2, laiM, vodA, vodB, vodC);
	
	push!(tbx_H,tb2_all)
end

tdiff_allH = [sqrt.(mean((cat_temp(tb1)-cat_temp(x)) .^ 2)) for x in tbx_H];
tdiff_6H = [sqrt.(mean((cat_temp(tb1)-cat_temp(x))[vcat(mask6,mask6)] .^ 2)) for x in tbx_H];
tdiff_1H = [sqrt.(mean((cat_temp(tb1)-cat_temp(x))[vcat(mask1,mask1)] .^ 2)) for x in tbx_H];

figure()
plot(tdiff_all, shapeLWPdiff,"o",label="Full diurnal")
#plot(tdiff_1, shapeLWPdiff, "o",label="1 AM/PM")
#plot(tdiff_6, shapeLWPdiff, "o",label="6 AM/PM")
plot([0],[0],"ko")
xlabel("RMSE of Tb (K)")
ylabel("RMSE of LWP diurnal shape")
plot([1,1]*significant_rmse_10yr,[0,maximum(myY)],"k--")
plot([1,1]*significant_rmse_1yr,[0,maximum(myY)],"k--")
text(significant_rmse_10yr+0.01,maximum(myY)*0.75,"Significant\nover 10 years",rotation="vertical")
text(significant_rmse_1yr+0.01,maximum(myY)*0.75,"Significant\nover 1 year",rotation="vertical")
#legend()



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

