
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

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 2;
rcParams["font.size"] = 18;
rcParams["mathtext.default"] = "regular";

df_raw = CSV.read(string(PROJECT_HOME,"/data/moflux_fluxnet_data.csv"), DataFrame);

df_raw[!,"RAIN"] *= 2; #convert per half hour to per hour

include(string(PROJECT_HOME,"/simulation_code/full_model.jl"));
include(string(PROJECT_HOME,"/assim_code/time_averaging.jl"))



N = 24*365*1
istart = 24*365*6 + 1; 
soil0 = 0.35;


deltaT = FT(60*60);
alpha = FT(10.9)
nsoil = FT(1.5);

function run_sim_3layer(vcmax_par::FT, k_frac::FT, k_plant::FT,
      z_soil::FT, weibB::FT, vol_factor::FT, g1::FT, k_soil::FT,
      smc_runoff::FT, exp_root_dist::FT,
      canopy_pvslope::FT, trunk_pvslope::FT)
	return run_sim_varyB(vcmax_par, k_frac, weibB, FT(4), k_plant, 
    k_soil, z_soil, istart, N, soil0, vol_factor, 1e-5, 3, df_raw, g1, deltaT, alpha, nsoil,0,0,
    smc_runoff, exp_root_dist,canopy_pvslope, trunk_pvslope);
end


pars0 = convert(Array{FT}, [90, 0.25, 10, 2000, 4, 1.0, 400, 1e-6,
         0.4,2,1/20,1/20]);
sim_res1 = run_sim_3layer(pars0...);
cm1 = mean(sim_res1[6][:,1:3],dims=2);
pdLWP = Float64.(replace(sim_res1[1].LWP_predawn, missing => NaN));


function plot_tree(x,y,node_example)
	leaf_levels = [12,16,18.5]
	trunk_level = 9
	soil_bounds = node_example.soil_bounds
	root_levels = (soil_bounds[2:end] + soil_bounds[1:(end-1)])/2
	figure()
	for i in 1:3
		plot([x[i],x[i+3],x[7]],[leaf_levels[i]+1,leaf_levels[i],trunk_level],
			 "gX-")
	end
	for i in 1:8
		plot([y[i],x[7+i],x[7]],[root_levels[i],root_levels[i],0],
			 "s-",color="brown")
	end
	plot([1,1]*x[7],[0,trunk_level],"ko-")
	xlabel("Water potential (MPa)")
	ylabel("Height (m)")
end

inight = 180*24 + 4;
plot_tree(sim_res1[6][inight,:], sim_res1[7][inight,:],sim_res1[end])
xlim(-3.5,0)

hourlist = collect(1:(365*24))/24;

figure()
plot(hourlist,sim_res1[6][:,1:3],color="green",alpha=0.5)
plot(hourlist,sim_res1[6][:,4:6],color="blue",alpha=0.5)
plot(hourlist,sim_res1[6][:,7],color="k",alpha=0.5)
plot(hourlist,sim_res1[6][:,8:15],color="brown",alpha=0.5)
xlabel("Time (day of year)")
ylabel("Water potential (MPa)")
plot([],[],color="green",alpha=0.5,label="Leaves")
plot([],[],color="blue",alpha=0.5,label="Branches")
plot([],[],color="k",alpha=0.5,label="Trunk")
plot([],[],color="brown",alpha=0.5,label="Roots")
legend()


daylist = collect(1:365);

figure()
plot(daylist, get_daily(sim_res1[1].LE/40650,24)*18/1000*60*60*24,label="Eddy covariance")
plot(daylist, get_daily(sim_res1[1].ETmod,24)*18/1000*60*60*24,label="Model")
xlabel("Time (day of year)")
ylabel("ET (mm/day)")
legend()

figure()
plot(daylist, get_daily(sim_res1[1].GPP_night,24),label="Eddy covariance")
plot(daylist, get_daily(sim_res1[1].GPP,24),label="Model")
xlabel("Time (day of year)")
ylabel("GPP (umol/m2/s)")
legend()


figure()
plot(cm1,pdLWP,"o",fillstyle="none")
plot([0,-4],[0,-4],color="grey")
xlabel("Modeled (MPa)")
ylabel("Observed (MPa)")
title("Pre-dawn leaf water potential")



figure()
plot(daylist,cm1[7:24:end],label="Model")
plot(daylist,pdLWP[7:24:end],"o",label="Observations")
xlabel("Time (day of year)")
ylabel("Water potential (MPa)")
title("Pre-dawn leaf water potential")
legend()