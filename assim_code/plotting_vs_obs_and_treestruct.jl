
using DataFrames
using CSV
using PyPlot
using Dates

include("../get_home_dir.jl")

pygui(true)

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 2;
rcParams["font.size"] = 18;
rcParams["mathtext.default"] = "regular";


#include(string(PROJECT_HOME,"/assim_code/mironov.jl"));
#include(string(PROJECT_HOME,"/assim_code/tau_omega_new2_only_omega.jl"));
include(string(PROJECT_HOME,"/assim_code/time_averaging.jl"))

df_raw = CSV.read(string(PROJECT_HOME,"/data/moflux_fluxnet_data.csv"), DataFrame);
df_raw[!,"RAIN"] *= 2;

#N = 24*365*12
#istart = 24*365*0 + 1; 

N = 24*365*1;
istart = 24*365*2 + 1; 

soil0 = 0.4;

include(string(PROJECT_HOME,"/simulation_code/full_model_newP63.jl"));
#include(string(PROJECT_HOME,"/simulation_code/full_model_newP63_ballberry.jl"));

deltaT = FT(60*60);

soil_curve_1 = VanGenuchten{FT}(stype = "Ozark",α = FT(10.9),   n = 1.5, Θs = 0.55,   Θr = 0.067);
smc1 = soil_swc(soil_curve_1,FT(-1));


function run_sim_0(vcmax_par::FT, weibB_stom::FT, k_plant::FT,
    z_soil::FT, xylem_extraB::FT, vol_factor::FT, g1::FT, k_soil::FT,
    slope_runoff::FT, exp_root_dist::FT,nsoil::FT,soilfac::FT
    )
weibC_plant = FT(2);
soil_curve_2 =  VanGenuchten{FT}(stype = "Ozark",α = FT(10.9),   n = nsoil, Θs = 0.55,   Θr = 0.067);
pot1b = soil_p_25_swc(soil_curve_2,smc1);
alpha2 = FT(10.9*soilfac)*abs(pot1b);
return run_sim_varyB(vcmax_par, weibB_stom/soilfac, (weibB_stom+xylem_extraB)/soilfac, weibC_plant, k_plant*soilfac,
  k_soil*soilfac, z_soil, istart, N, soil0, vol_factor*soilfac, 1e-5, 3, df_raw, g1, deltaT, alpha2, nsoil,0,0,
  slope_runoff, exp_root_dist,FT(1/20),FT(1/20));
end



pars0 = convert(Array{FT}, [90, 3, 10, 2000, 1, 1.0, 300, 0.4e-6,
       0.4,2,1.5,1]);


sim_res1 = run_sim_0(pars0...);


cm1 = mean(sim_res1[6][:,1:3],dims=2);
pdLWP = Float64.(replace(sim_res1[1].LWP_predawn, missing => NaN));

mydates = DateTime.(sim_res1[1].DateTime, "yyyy-mm-dd HH:MM:SS");
daylist = Date.(mydates)[1:24:end];


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

inight = 237*24 + 4;
plot_tree(sim_res1[6][inight,:], sim_res1[7][inight,:],sim_res1[end])
xlim(-3.5,0)
title("5 AM")

inight = 237*24 + 12;
plot_tree(sim_res1[6][inight,:], sim_res1[7][inight,:],sim_res1[end])
xlim(-3.5,0)
title("Noon")


figure()
plot(daylist[1:5:end],get_daily(sim_res1[1].LE/40650,24*5)*18/1000*60*60*24,"k",label="Eddy covariance")
plot(daylist[1:5:end],get_daily(sim_res1[1].ETmod,24*5)*18/1000*60*60*24,"r",label="Model",alpha=0.75)
xlabel("Time")
ylabel("ET (mm/day)")
legend()

figure()
plot(daylist[1:5:end],get_daily(sim_res1[1].GPP_night,24*5),"k",label="Eddy covariance")
plot(daylist[1:5:end],get_daily(sim_res1[1].GPP,24*5),"r",label="Model",alpha=0.75)
xlabel("Time")
ylabel("GPP (umol/m^2/s)")
legend()


figure()
plot(daylist[1:5:end],get_daily(sim_res1[1].SMC,24*5),"k",label="Observed")
plot(daylist[1:5:end],get_daily(sim_res1[2][:,1],24*5),"r",label="Model",alpha=0.75)
#plot(daylist[1:5:end]/365,get_daily(sim_res1long[1].ColumnSMC,24*5),"b",label="Model surface",alpha=0.75)
xlabel("Time")
ylabel("Surface soil moisture")
legend()

pdLWP = Float64.(replace(sim_res1[1].LWP_predawn, missing => NaN));

figure()
plot(daylist,cm1[7:24:end],"r",alpha=0.75,label="Model")
plot(daylist,pdLWP[7:24:end],"ko",label="Observations")
#plot(daylist[1:5:end]/365,get_daily(sim_res1long[1].ColumnSMC,24*5),"b",label="Model surface",alpha=0.75)
xlabel("Time")
ylabel("Water potential (MPa)")
legend()
