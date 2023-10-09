
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


N = 24*365*1;
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

mydates = DateTime.(sim_res1[1].DateTime, "yyyy-mm-dd HH:MM:SS");
daylist = Date.(mydates)[1:24:end];

pars2 = convert(Array{FT}, [90, 3, 10, 2000, 1, 2.0, 300, 0.4e-6,
       0.4,2,1.5,1]);

sim_res2 = run_sim_0(pars2...);
cm2 = mean(sim_res2[6][:,1:3],dims=2);

figure()
plot(mydates,cm1,color="black",label="Baseline",linewidth=2)
plot(mydates,cm2,"--",color="red",label="2x capacitance",linewidth=2,alpha=0.75)
xlabel("Time",fontsize=20)
ylabel(raw"Leaf water potential (MPa)",fontsize=20)

#ylim(-5,0)
xlim(mydates[180*24],mydates[(180+5)*24])
legend(fontsize=20,ncol=2)
xticks(fontsize=16)
yticks(fontsize=16)
ylim(-3.5,0)

