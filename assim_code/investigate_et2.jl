include("../get_home_dir.jl")

using DataFrames
using CSV
using PyPlot
using Random
using StatsBase



include(string(PROJECT_HOME,"/assim_code/mironov.jl"));
include(string(PROJECT_HOME,"/assim_code/tau_omega_funs.jl"));
include(string(PROJECT_HOME,"/assim_code/time_averaging.jl"))

df_raw = CSV.read(string(PROJECT_HOME,"/data/moflux_land_data_skipyear_hourly2.csv"), DataFrame);
df_raw[!,"RAIN"] *= 2; #rain is mm/half hour, need to convert it to mm/time step

N = 24*365*1
istart = 24*365*8+1
soil0 = 0.42;

include(string(PROJECT_HOME,"/simulation_code/sim_vary_new_stomata_med.jl"));

deltaT = FT(60*60);
#alpha = FT(1.368);
#nsoil = FT(2.6257);

#alpha = FT(2);
#nsoil = FT(1.8);


#nsoil = FT(1.7);
#alpha = FT(3.2);

#alpha = FT(1.9);
#nsoil = FT(3.2);

#alpha = FT(204);
#nsoil = FT(1.4);


alpha = FT(2.5);
nsoil = FT(1.2);


#alpha = FT(1);
#nsoil = FT(1.8);

#alpha = FT(204/4);
#nsoil = FT(1.6);


function run_sim_2(vcmax_par::FT, k_frac::FT, k_plant::FT, z_soil::FT, weibB::FT, weibC::FT, vol_factor::FT, g1::FT)
        return run_sim_vary(vcmax_par, k_frac, weibB, weibC, k_plant, FT(2e-5), z_soil, istart, N, soil0, vol_factor, 1e-5, 3, df_raw, g1, deltaT, alpha, nsoil);
end

#sim_res1 = run_sim_2(FT(60),FT(0.25), FT(2),FT(1e-5), FT(600),FT(3),FT(2),FT(1),FT(506));
#sim_res1 = run_sim_2(FT(31),FT(0.25), FT(2),FT(0.4e-6), FT(800),FT(4),FT(2),FT(1),FT((16-0.0*9.3)*sqrt(1000)));

sim_res1 = run_sim_2(FT(31),FT(0.5), FT(2), FT(2000),FT(5),FT(2),FT(1),FT((16-0.0*9.3)*sqrt(1000)));

logpar2 = [3.420552042,-3.005631752,0.140492204,7.77502719,1.384741914,1.146245499,0.559185057,6.697470399]
par2 = convert(Array{FT}, exp.(logpar2));
sim_res2 = run_sim_2(par2...);
#sim_res2 = run_sim_2(FT(31),FT(0.5), FT(2), FT(500),FT(5),FT(2),FT(1),FT(506));


gravity_factor = (18.5+9)/2 * 1000/mpa2mm;

leafpot = mean(sim_res1[3],dims=2);
simTB = get_TB(sim_res1, leafpot, 0.041, 0.82, 0.051);
simTB0 = get_TB(sim_res1, leafpot, 0, 0.82, 0.051);
simTB2 = get_TB(sim_res1, leafpot, 0.082, 0.82, 0.051);


#%%
noon_ET = sim_res1[1].LE[12:24:end]/44200;
noon_VPD = sim_res1[1].VPD[12:24:end];
noon_PATM = sim_res1[1].P_ATM[12:24:end]*1000;
noon_LAI = sim_res1[1].LAI_modis[12:24:end];
noon_GLW = noon_ET ./ (noon_VPD ./ noon_PATM) ./ noon_LAI;

noon_ET_mod = sim_res1[1].ETmod[12:24:end];
noon_GLW_mod = noon_ET_mod ./ (noon_VPD ./ noon_PATM) ./ noon_LAI;


noon_GLW_select = noon_GLW[noon_LAI .> 3];
noon_GLW_mod_select = noon_GLW_mod[noon_LAI .> 3];

noon_GLW_select2 = noon_GLW_select[noon_GLW_select .> 0];
noon_GLW_mod_select2 = noon_GLW_mod_select[noon_GLW_select .> 0];


pdLWP = Float64.(replace(sim_res1[1].LWP_predawn, missing => NaN));
mdLWP = Float64.(replace(sim_res1[1].LWP_midday, missing => NaN));

#nsoil = 1.4
#alpha = 200
msoil = 1-1/nsoil;
swc_eq = ((-pdLWP * alpha) .^ nsoil .+ 1) .^ (-msoil) * (0.45 - 0.067) .+ 0.067;


figure()
plot(leafpot[7:24:end])
#plot(sim_res1[5][7:24:end])

#plot(sim_res1[5][7:24:end])
plot(pdLWP[7:24:end] .- gravity_factor, "o")

#=
s_arr = collect(0.14:0.005:0.45);
n = 1.4; m = 1-1/n; a = 40;
psi_arr = (((s_arr .- 0.067) / (0.45 - 0.067)) .^ (-1/m) .-1) .^ (1/n) / a;
psi_arr *= 1.0/psi_arr[s_arr .== 0.25];

figure()
scatter(mean(sim_res1[2],dims=2), pdLWP, c = sim_res1[1].YEAR)
#scatter(sim_res1[2][:,1], pdLWP)
#scatter(sim_res1[2][:,4], pdLWP)
plot(s_arr, -psi_arr,"k")

n = 2.0; m = 1-1/n; a = 40;
psi_arr = (((s_arr .- 0.067) / (0.45 - 0.067)) .^ (-1/m) .-1) .^ (1/n) / a;
psi_arr *= 1.75/psi_arr[s_arr .== 0.2];
plot(s_arr, -psi_arr,color="orange")
=#

figure()
plot((mean(sim_res1[2][:,3:4],dims=2) .- 0.067) / (0.45-0.067), pdLWP, "o")
srel_arr = collect(0.6:0.01:0.95);
n = 1.2; m = 1-1/n; a = 200;
psi_arr = (srel_arr .^ (-1/m) .-1) .^ (1/n);
psi_fac = 1.5/psi_arr[srel_arr .== 0.75];
plot(srel_arr, -psi_arr*psi_fac,"k")


#=
figure()
plot(-get_daily(sim_res1[1].LE/44200, 24),color="grey")
twinx()
plot(pdLWP[7:24:end],"ro")
twinx()
plot(sim_res1[1].SMC[1:24:end],"k")
plot(sim_res1[2][1:24:end,1],color="orange")
plot(sim_res1[2][1:24:end,4],color="blue")
=#

#%%

#plot(log.(vpdobs[12:24:end]),    log.(etratio)[12:24:end], "o")
figure()
plot(get_daily(sim_res1[1].LE/44200, 24*5))
plot(get_daily(sim_res1[1].ETmod, 24*5))
#%%
#=
figure()
plot(cumsum(get_daily(sim_res1[1].LE/44200, 24)))
plot(cumsum(get_daily(sim_res1[1].ETmod, 24)))
=#
#%%
figure()
plot(get_diurnal(sim_res1[1].LE/44200, 24))
plot(get_diurnal(sim_res1[1].ETmod, 24))
#=
function run_sim_1day(vcmax_par::FT, k_frac::FT, k_plant::FT, k_soil::FT, z_soil::FT, weibB::FT, weibC::FT, vol_factor::FT, g1::FT, hnew)
	df2 = deepcopy(df_raw[(150*24):(151*24),:]);
	df2.RelHum = 0*df2.RelHum .+ hnew;
    return run_sim_vary(vcmax_par, k_frac, weibB, weibC, k_plant, k_soil, z_soil, 1, 24, soil0, vol_factor, 1e-5, 3, df2, g1, deltaT, alpha, nsoil);
end

h_list = collect(0.2:0.1:1);
sim_list = [run_sim_1day(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(24), h) for h in h_list];
et_list = [mean(x[1].ETmod) for x in sim_list];
noon_glw = [x[1].glw[12] for x in sim_list];

=#

#=
pars_post = CSV.read("C:/Users/natan/OneDrive - Leland Stanford Junior University/Documents/moflux_docs/pars_to_experiment/obs6/post_par.csv", DataFrame);
logpars_ret = convert(Array,pars_post[end,1:9]);
#logpars_ret = [3.2180858900299945,-0.886292105647866,0.7870081261568632,-11.047899629216229,6.545214727186734,1.3461444453928548,1.195044099092091,0.17049235312958116,6.7906980609243295];
pars_ret = convert(Array{FT}, exp.(logpars_ret));
sim_res2 = run_sim_2(pars_ret...);

etdaily1 = get_daily(sim_res1[1].ETmod, 24);
etdaily2 = get_daily(sim_res2[1].ETmod, 24);
figure()
plot(etdaily1,label="True");
plot(etdaily2,label="Retrieved");

figure()
plot(get_diurnal(sim_res1[1].ETmod,24),label="True");
plot(get_diurnal(sim_res2[1].ETmod,24),label="Retrieved");

figure()
plot(get_diurnal(sim_res1[1].leafpot,24),label="True");
plot(get_diurnal(sim_res2[1].leafpot,24),label="Retrieved");

figure()
plot(sim_res1[1].leafpot[1:24:end],label="True");
plot(sim_res2[1].leafpot[1:24:end],label="Retrieved");

include(string(PROJECT_HOME,"/simulation_code/sim_sub_LWPvcmax.jl"));

function run_sim_3(vcmax_par::FT, k_frac::FT, k_plant::FT, k_soil::FT, z_soil::FT, weibB::FT, weibC::FT, vol_factor::FT, g1::FT)
        return run_sim_plugin(vcmax_par, k_frac, weibB, weibC, k_plant, k_soil, z_soil, istart, N, soil0, vol_factor, 1e-5, 3, df_raw, g1, deltaT, alpha, nsoil, sim_res1[3]);
end

#sim_res1b = run_sim_3(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));
sim_res2b = run_sim_3(pars_ret...);


etdaily2b = get_daily(sim_res2b[1].ETmod, 24);

figure()
plot(etdaily1,label="True");
#plot(etdaily2,label="Retrieved");
plot(etdaily2b,label="Retrieved 2");

figure()
plot(get_diurnal(sim_res1[1].ETmod,24),label="True");
plot(get_diurnal(sim_res2b[1].ETmod,24),label="Retrieved");

figure()
plot(get_diurnal(sim_res1[1].leafpot,24),label="True");
plot(get_diurnal(sim_res2b[1].leafpot,24),label="Retrieved");

figure()
plot(sim_res1[1].leafpot[1:24:end],label="True");
plot(sim_res2b[1].leafpot[1:24:end],label="Retrieved");
=#
