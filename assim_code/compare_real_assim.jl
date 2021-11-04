
using DataFrames
using CSV
#using Plots
using Random
using StatsBase

#include("../simulation_code/rebuild_sim_medlyn_0g0_plugin_ET.jl");


include("time_averaging.jl")

#N = 48*3
#istart = 48*365*2 - 52*48 + 230*48


N = 48*365*1
istart = 48*365*2 - 52*48 #+ 230*48
soil0 = 0.39;

post_medlyn = Array(CSV.read("samples_emp0_onlyET.csv", DataFrame));
post_wang = Array(CSV.read("samples_opt_onlyET.csv", DataFrame));

#include("../simulation_code/rebuild_sim_medlyn_0g0_plugin_ET.jl");

include("../simulation_code/rebuild_sim_Scrit_drain_plugin_ET.jl");


#fixdepth = FT(1650);
fixdepth = FT(800);

soil_cut = FT(0.2);
drain_rate = FT(0);
porosity = FT(0.55);
cutoff = FT(1);
plantK = FT(5);
soilB = FT(1.5);
soilP = FT(1.5);
soil_K = FT(0.0001);

#a_med = convert(Array{FT},exp.(post_medlyn[end,1:11]))
#sim_res_med = convert_sim(run_sim(FT(22), FT(0.15), plantK, soil_K, soilB, soilP, fixdepth, porosity,istart,N,soil0,drain_rate, a_med[10], cutoff,0));

a_wang = convert(Array{FT},exp.(post_wang[end,1:9])) 
sim_res_wang = convert_sim(run_sim(FT(22), FT(0.1), plantK, soil_K, soilB, soilP, fixdepth, porosity,istart,N,soil0,drain_rate,0));
sim_res_wang_t = @timed convert_sim(run_sim(FT(22), FT(0.1), plantK, soil_K, soilB, soilP, fixdepth, porosity,istart,N,soil0,drain_rate,0));


#sim_res_wang_t = @timed convert_sim(run_sim(FT(22), FT(0.1), plantK, soil_K, soilB, soilP, fixdepth, porosity,istart,N,soil0,drain_rate,0,0));
#sim_res_wang = sim_res_wang_t.value;

#sim_res_med_t = @timed convert_sim(run_sim(FT(22), FT(0.1), plantK, soil_K, soilB, soilP, fixdepth, porosity,istart,N,soil0,drain_rate,0,1));
#sim_res_med = sim_res_med_t.value;

include("../simulation_code/rebuild_sim_ij_par_vol.jl");
sim_res_med = convert_sim(run_sim(FT(22), FT(0.1), plantK, soil_K, soilB, soilP, fixdepth, porosity,istart,N,soil0,drain_rate,1));
sim_res_med_t = @timed convert_sim(run_sim(FT(22), FT(0.1), plantK, soil_K, soilB, soilP, fixdepth, porosity,istart,N,soil0,drain_rate,0.1));



obsET = sum(reshape(sim_res_med[1].LE/44200, (24,:)),dims=1)[1,:]/24;
etWang = sum(reshape(sim_res_wang[1].ETmod, (24,:)),dims=1)[1,:]/24;
etMed = sum(reshape(sim_res_med[1].ETmod, (24,:)),dims=1)[1,:]/24;

obsSMC = sum(reshape(sim_res_med[1].SMC, (24,:)),dims=1)[1,:]/24;
smcWang = sum(reshape(sim_res_wang[2][:,1], (24,:)),dims=1)[1,:]/24;
smcMed = sum(reshape(sim_res_med[2][:,1], (24,:)),dims=1)[1,:]/24;


using PyPlot

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 1;
rcParams["font.size"] = 24;

figure()
plot(obsET,"k",label="Obs"); plot(etWang,color="blue",label="Opt"); plot(etMed,color="red", label="Emp")
legend()
xlabel("Day of year")
ylabel("ET (mol/m2/s)")

figure()
plot(obsSMC,"k",label="Obs"); plot(smcWang,color="blue",label="Opt"); plot(smcMed,color="red", label="Emp")
legend()
xlabel("Day of year")
ylabel("Surface SMC")

day_arr = collect(1:(24*365))/24;
figure()
plot(day_arr, sim_res_wang[1].leafpot,color="blue",label="Opt"); plot(day_arr, sim_res_med[1].leafpot, color="red", label="Emp")
legend()
xlabel("Day of year")
ylabel("LWP (MPa)")



figure()
plot(day_arr, mean(sim_res_wang[3],dims=2),label="Leaf")
plot(day_arr, mean(sim_res_wang[4],dims=2),label="Branch")
plot(day_arr,sim_res_wang[5],label="Trunk")
legend()
xlabel("Day of year")
ylabel("LWP (MPa)")
title("No storage")

figure()
plot(day_arr, mean(sim_res_med[3],dims=2),label="Leaf")
plot(day_arr, mean(sim_res_med[4],dims=2),label="Branch")
plot(day_arr, sim_res_med[5],label="Trunk")
legend()
xlabel("Day of year")
ylabel("LWP (MPa)")
title("With storage")


