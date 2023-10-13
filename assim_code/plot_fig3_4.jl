
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

include(string(PROJECT_HOME,"/assim_code/time_averaging.jl"))

df_raw = CSV.read(string(PROJECT_HOME,"/data/moflux_fluxnet_data_nov2022_lef.csv"), DataFrame);
df_raw[!,"RAIN"] *= 2;


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

pars0 = convert(Array{FT}, [90, 3, 10, 2000, 1, 1.0, 300, 0.4e-6,
       0.4,2,1.5,1]);

sim_res1 = run_sim_0(pars0...);

cm1 = mean(sim_res1[6][:,1:3],dims=2);
pdLWP = Float64.(replace(sim_res1[1].LWP_predawn, missing => NaN));

mydates = DateTime.(sim_res1[1].DateTime, "yyyy-mm-dd HH:MM:SS");
daylist = Date.(mydates)[1:24:end];

dry_year = [2005,2012,2013,2014];
drystart = [Date(x,6,1) for x in dry_year];
dryend = [Date(x,9,30) for x in dry_year];

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
axes[1,1].set_title("(a) 5 AM",fontsize=16,loc="left")


inight = 275*24 + 13;
plot_tree_1trunk_fig(axes[2,1],sim_res1[6][inight,:], sim_res1[7][inight,:],sim_res1[end]);
axes[2,1].set_xlim(-2.5,0)
axes[2,1].set_title("(b) Noon",fontsize=16,loc="left")

axes[2,1].set_xlabel("Water potential (MPa)",fontsize=16)

axes[2,1].text(-2.45,12,"Leaves",color="g",fontweight="bold")
axes[2,1].text(-1.65,12,"Branches",color="g",fontweight="bold")
axes[2,1].text(-1.8,-0.5,"Trunk",color="k",fontweight="bold")
axes[2,1].text(-1.15,2,"Roots",color="brown",fontweight="bold")
axes[2,1].text(-0.45,2,"Soil layers",color="brown",fontweight="bold")


fig.supylabel("Relative height (not to scale)",fontsize=16)

tight_layout()
savedir = "C:/Users/natan/OneDrive - Stanford/Documents/clima paper writing/coauthor edits/wrr reviews/fig_pdfs/";
savefig(string(savedir,"fig3.pdf"))

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
plot(daylist[1:5:end],get_daily(sim_res1[1].ETmod,24*5)*18/1000*60*60*24,"r",label="CliMA Land",alpha=0.75)
plot(daylist[1:5:end],etnan,"k",label="Eddy covariance")

#xlabel("Time",fontsize=24)
#ylabel("ET (mm/day)")
ylim(-0.1,8)
xlim(daylist[1],daylist[end])
legend(fontsize=15,loc="upper left",ncol=2)
title("(a)",loc="left",fontsize=20)
title("Evapotranspiration",loc="center",fontsize=20)
ylabel("ET (mm/day)",fontsize=16)


pdLWP = Float64.(replace(sim_res1[1].LWP_predawn, missing => NaN));

subplot(2,1,2)
plot(daylist,cm1[7:24:end],"r",alpha=0.75,label="CliMA Land")
plot(daylist,pdLWP[7:24:end],"ko",label="Observations",markersize=5,fillstyle="none",markeredgewidth=2)
xlabel("Time",fontsize=20)
title("Predawn leaf water potential",loc="center",fontsize=20)
ylabel(raw"$\Psi$ (MPa)",fontsize=20)

axvspan(Date(2007,1,1),Date(2007,12,31),color="green",alpha=0.33,label="Model-data fusion period")

ylim(-5,0)
xlim(daylist[1],daylist[end])
legend(fontsize=15,ncol=3)
title("(b)",loc="left",fontsize=16)

tight_layout()

savefig(string(savedir,"fig4.pdf"))

#savefig("data_comparison_fig4c.png")

