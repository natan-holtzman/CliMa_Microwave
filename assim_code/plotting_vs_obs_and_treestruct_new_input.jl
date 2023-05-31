
using DataFrames
using CSV
using PyPlot
using Dates

include("../get_home_dir.jl")

pygui(true)

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 1;
rcParams["font.size"] = 10;
rcParams["mathtext.default"] = "regular";


#include(string(PROJECT_HOME,"/assim_code/mironov.jl"));
#include(string(PROJECT_HOME,"/assim_code/tau_omega_new2_only_omega.jl"));
include(string(PROJECT_HOME,"/assim_code/time_averaging.jl"))

df_raw = CSV.read(string(PROJECT_HOME,"/data/moflux_fluxnet_data_nov2022_lef.csv"), DataFrame);
df_raw[!,"RAIN"] *= 2;

#df_raw.LAI_modis = 0.9 .+ 5/3*(df_raw.LAI_modis .- 0.9);

#N = 24*365*12
#istart = 24*365*0 + 1; 

N = 24*365*13;
istart = 24*365*0 + 1; 

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

#=
function run_sim_1(vcmax_par::FT, weibB_stom::FT, k_plant::FT,
    z_soil::FT, xylem_extraB::FT, vol_factor::FT, g1::FT, k_soil::FT,
    slope_runoff::FT, exp_root_dist::FT,nsoil::FT,soilfac::FT
    )
weibC_plant = FT(2);
soil_curve_2 =  VanGenuchten{FT}(stype = "Ozark",α = FT(10.9),   n = nsoil, Θs = 0.55,   Θr = 0.067);
pot1b = soil_p_25_swc(soil_curve_2,smc1);
alpha2 = FT(10.9*soilfac)*abs(pot1b);
return run_sim_BallBerry(vcmax_par, weibB_stom/soilfac, (weibB_stom+xylem_extraB)/soilfac, weibC_plant, k_plant*soilfac,
  k_soil*soilfac, z_soil, istart, N, soil0, vol_factor*soilfac, 1e-5, 3, df_raw, g1, deltaT, alpha2, nsoil,0,0,
  slope_runoff, exp_root_dist,FT(1/20),FT(1/20));
end
=#

pars0 = convert(Array{FT}, [90, 3, 10, 2000, 1, 1.0, 300, 0.4e-6,
       0.4,2,1.5,1]);

pars1 = convert(Array{FT}, [60, 3, 10, 2000, 1, 1.0, 500, 0.4e-6,
       0.4,2,1.5,1]);
pars1b = convert(Array{FT}, [60, 3, 10, 2000, 1, 1.0, 400, 0.4e-6,
       0.4,2,1.5,1]);

sim_res1 = run_sim_0(pars0...);


cm1 = mean(sim_res1[6][:,1:3],dims=2);
pdLWP = Float64.(replace(sim_res1[1].LWP_predawn, missing => NaN));

mydates = DateTime.(sim_res1[1].DateTime, "yyyy-mm-dd HH:MM:SS");
daylist = Date.(mydates)[1:24:end];

dry_year = [2005,2012,2013,2014];
drystart = [Date(x,6,1) for x in dry_year];
dryend = [Date(x,9,30) for x in dry_year];

#150 - 275
#year 1, 3, 7


function plot_tree(plant_psi,soil_psi,node_example)
	leaf_levels = [12,16,18.5]
	trunk_level = 9
	soil_bounds = node_example.soil_bounds
	root_levels = (soil_bounds[2:end] + soil_bounds[1:(end-1)])/2
	figure()
	for i in 1:3
		plot([plant_psi[i],plant_psi[i+3],plant_psi[7]],[leaf_levels[i]+1,leaf_levels[i],trunk_level],
			 "gX-")
	end
	for i in 1:8
		plot([soil_psi[i],plant_psi[7+i],plant_psi[7]],[root_levels[i],root_levels[i],0],
			 "s-",color="brown")
	end
	plot([1,1]*plant_psi[7],[0,trunk_level],"ko-")
	xlabel("Water potential (MPa)")
	ylabel("Height (m)")
end

function plot_tree_relative(plant_psi,soil_psi,node_example)
	leaf_levels = [12,16,18.5]
	trunk_level = 9
	soil_bounds = node_example.soil_bounds
	root_levels = (soil_bounds[2:end] + soil_bounds[1:(end-1)])*2
	figure()
	for i in 1:3
		plot([plant_psi[i],plant_psi[i+3],plant_psi[7]],[leaf_levels[i]+1,leaf_levels[i],trunk_level],
			 "gX-")
	end
	for i in 1:8
		plot([soil_psi[i],plant_psi[7+i],plant_psi[7]],[root_levels[i],root_levels[i],0],
			 "s-",color="brown")
	end
	plot([1,1]*plant_psi[7],[0,trunk_level],"ko-",markersize=10)
	xlabel("Water potential (MPa)")
	ylabel("Relative height (not to scale)")
	yticks([],[])
end


function plot_tree_1trunk(plant_psi,soil_psi,node_example)
	leaf_levels = [8,12,16]
	trunk_level = 9/2;
	soil_bounds = node_example.soil_bounds
	root_levels = (soil_bounds[2:end] + soil_bounds[1:(end-1)])*2
	figure()
	for i in 1:3
		plot([plant_psi[i],plant_psi[i+3],plant_psi[7]],[leaf_levels[i]+1,leaf_levels[i],trunk_level],
			 "gX-")
	end
	for i in 1:8
		plot([soil_psi[i],plant_psi[7+i],plant_psi[7]],[root_levels[i],root_levels[i],trunk_level],
			 "s-",color="brown")
	end
	plot([1]*plant_psi[7],[trunk_level],"ko",markersize=10)
	xlabel("Water potential (MPa)")
	ylabel("Relative height (not to scale)")
	yticks([],[])
end



function plot_tree_1trunk_fig(ax,plant_psi,soil_psi,node_example)
	leaf_levels = [8,12,16]
	trunk_level = 9/2;
	soil_bounds = node_example.soil_bounds
	root_levels = (soil_bounds[2:end] + soil_bounds[1:(end-1)])*3
	for i in 1:3
		ax.plot([plant_psi[i],plant_psi[i+3],plant_psi[7]],[leaf_levels[i]+1,leaf_levels[i],trunk_level],
			 "gX-",markersize=4)
	end
	for i in 1:8
		ax.plot([soil_psi[i],plant_psi[7+i],plant_psi[7]],[root_levels[i],root_levels[i],trunk_level],
			 "s-",color="brown",markersize=4)
	end
	ax.plot([1]*plant_psi[7],[trunk_level],"ko",markersize=8)
	ax.set_yticks([],[])
end

fig, axes = subplots(2,1)
inight = 275*24 + 6;
plot_tree_1trunk_fig(axes[1,1],sim_res1[6][inight,:], sim_res1[7][inight,:],sim_res1[end]);
axes[1,1].set_xlim(-2.5,0)
axes[1,1].set_title("(a) 5 AM",fontsize=20,loc="left")


inight = 275*24 + 13;
plot_tree_1trunk_fig(axes[2,1],sim_res1[6][inight,:], sim_res1[7][inight,:],sim_res1[end]);
axes[2,1].set_xlim(-2.5,0)
axes[2,1].set_title("(b) Noon",fontsize=20,loc="left")

#xlabel("Water potential (MPa)")
#ylabel("Relative height (not to scale)")



function get_daily_missing(x,navg)
	if (typeof(x[1])==Float64) | (typeof(x[1])==Float32)
            xres = reshape(x,(navg,:));
            return [mean(skipmissing(xres[:,j])) for j in 1:size(xres)[2]];
	else
			return x[1:navg:end]
	end
end


etmiss = get_daily_missing(sim_res1[1].LE/44200,24*5)*18/1000*60*60*24;
etnan = Float64.(replace(etmiss, missing => NaN));



figure(figsize=(10,8))
subplot(2,1,1)
plot(daylist[1:5:end],etnan,"k",label="Eddy covariance")
plot(daylist[1:5:end],get_daily(sim_res1[1].ETmod,24*5)*18/1000*60*60*24,"r",label="CliMA Land",alpha=0.75)
xlabel("Time")
#ylabel("ET (mm/day)")
ylim(-0.1,8.5)
xlim(daylist[1],daylist[end])
legend()
title("(a)",loc="left",fontsize=18)
title("ET (mm/day)",loc="center",fontsize=18)


pdLWP = Float64.(replace(sim_res1[1].LWP_predawn, missing => NaN));

subplot(2,1,2)
plot(daylist,cm1[7:24:end],"r",alpha=0.75,label="CliMA Land")
plot(daylist,pdLWP[7:24:end],"ko",label="Observations",markersize=5,fillstyle="none",markeredgewidth=2)
#plot(daylist[1:5:end]/365,get_daily(sim_res1long[1].ColumnSMC,24*5),"b",label="Model surface",alpha=0.75)
xlabel("Time")
title("Predawn leaf water potential (MPa)",loc="center",fontsize=18)
axvspan(Date(2007,1,1),Date(2007,12,31),color="green",alpha=0.33,label="Model-data fusion period")
#=
axvspan(drystart[1],dryend[1],color="grey",alpha=0.33,label="Dry summers")
for j in 2:length(drystart)
	axvspan(drystart[j],dryend[j],color="grey",alpha=0.33)
end
=#
xlim(daylist[1],daylist[end])
legend()
title("(b)",loc="left",fontsize=18)

tight_layout()

savefig("data_comparison_may31.png")



lwpdaily = get_daily(cm1,24);
smdaily = get_daily(sim_res1[1].ColumnSMC,24);
etdaily = get_daily(sim_res1[1].ETmod,24)*18/1000*60*60*24;
gppdaily = get_daily(sim_res1[1].GPP,24);
vars_daily = hcat(lwpdaily, smdaily, etdaily, gppdaily);
var_means = mean(vars_daily,dims=1);
var_std = std(vars_daily,dims=1);




#=
figure()
plot(sim_res1[1].glw[13:24:end],(sim_res1[1].GPP ./ sim_res1[1].LAI_modis)[13:24:end],"o")
plot(get_daily(sim_res1[1].glw,24),get_daily(sim_res1[1].GPP ./ sim_res1[1].LAI_modis,24),"o")
plot([0,0.4],[0,0.4*50],"k")
#plot([0,0.4],[0,0.4*100],"k")
=#
#=
figure()
plot(daylist[1:5:end],get_daily(sim_res1[1].GPP_night,24*5),"k",label="Eddy covariance")
plot(daylist[1:5:end],get_daily(sim_res1[1].GPP,24*5),"r",label="Model",alpha=0.75)
xlabel("Time")
ylabel("GPP (umol/m^2/s)")
ylim(-1,18)
legend()


smiss = get_daily_missing(sim_res1[1].SMC,24*5);
snan = Float64.(replace(smiss, missing => NaN));
snan[450:510] .= NaN;


figure()
plot(daylist[1:5:end],snan,"k",label="Observed")
plot(daylist[1:5:end],get_daily(sim_res1[2][:,1],24*5),"r",label="Model",alpha=0.75)
#plot(daylist[1:5:end]/365,get_daily(sim_res1long[1].ColumnSMC,24*5),"b",label="Model surface",alpha=0.75)
xlabel("Time")
ylabel("Surface soil moisture")
legend()
=#
