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

N = 24*365*12
istart = 1 #+ 365*24 #24*365*2 - 52*24 #+ 230*48
soil0 = 0.42;

include(string(PROJECT_HOME,"/simulation_code/sim_vary_new_stomata_med.jl"));

deltaT = FT(60*60);
alpha = FT(1.368);
nsoil = FT(2.6257);

function run_sim_2(vcmax_par::FT, k_frac::FT, k_plant::FT, k_soil::FT, z_soil::FT, weibB::FT, weibC::FT, vol_factor::FT, g1::FT)
        return run_sim_vary(vcmax_par, k_frac, weibB, weibC, k_plant, k_soil, z_soil, istart, N, soil0, vol_factor, 1e-5, 3, df_raw, g1, deltaT, alpha, nsoil);
end

#sim_res1 = run_sim_2(FT(60),FT(0.25), FT(2),FT(1e-5), FT(600),FT(3),FT(2),FT(1),FT(506));

sim_res1 = run_sim_2(FT(60),FT(0.25), FT(2),FT(1e-6), FT(800),FT(3),FT(2),FT(1),FT(506));


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




#%%

#plot(log.(vpdobs[12:24:end]),    log.(etratio)[12:24:end], "o")
figure()
plot(get_daily(sim_res1[1].LE/44200, 24*5))
plot(get_daily(sim_res1[1].ETmod, 24*5))
#%%
figure()
plot(cumsum(get_daily(sim_res1[1].LE/44200, 24)))
plot(cumsum(get_daily(sim_res1[1].ETmod, 24)))
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
