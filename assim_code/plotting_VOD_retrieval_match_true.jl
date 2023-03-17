
using DataFrames
using CSV
using PyPlot
using Dates

include("../get_home_dir.jl")

pygui(true)

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 1;
rcParams["lines.markersize"] = 2;

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

N = 24*365*1;
istart = 24*365*2 + 1; 

soil0 = 0.4;

include(string(PROJECT_HOME,"/simulation_code/full_model_newP63.jl"));
include(string(PROJECT_HOME,"/simulation_code/full_model_newP63_ballberry.jl"));

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


pars0 = convert(Array{FT}, [90, 3, 10, 2000, 1, 1.0, 300, 0.4e-6,
       0.4,2,1.5,1]);
lplist = [[  4.99405767,   0.80343857,   2.45331971,   7.53797763,
0.47409801,   0.2804809 ,   2.92590326, -14.88419025,
-0.89269341,   0.15087817,   0.41616954,   0.21764421],
[  5.30183155,   0.71639358,   2.60355008,   7.57660079,
0.78765758,   0.60152583,   3.01529433, -16.06559991,
-0.80184993,  -1.50657843,   0.47462958,   0.78405518],
[  4.89091545,   0.88239074,   2.60353234,   7.65274827,
0.76215424,   0.91681717,   2.86015231, -15.40996107,
-0.94303818,   0.24113653,   0.50159069,   0.76009207],
[  5.21425855,   0.67373583,   2.48586249,   7.46280712,
          1.03727476,   0.46119343,   2.85151803, -15.82363334,
         -0.87158386,   0.24401006,   0.48198944,   0.50062545]];

sim_res0 = run_sim_0(pars0...);

vodpars = [[0.08044635, 0.82301716, 0.05181002],
[0.17225362, 0.83603052, 0.05339761],
[0.16877381, 0.83558454, 0.05491833],
[0.14339696, 0.83141694, 0.05564335]];
#vodB, vodA, vodC


include(string(PROJECT_HOME,"/assim_code/mironov.jl"));
include(string(PROJECT_HOME,"/assim_code/tau_omega_new2_only_omega.jl"));

function getVOD(sim_resA,vodA,vodB,vodC,omega)
    canpot_true = mean(sim_resA[6][:,1:3],dims=2);

    tsoil = sim_resA[1].T_SOIL .+ 273.15;
    tcan = (sim_resA[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
    laiM = sim_resA[1].LAI_modis;
    canpot_true = mean(sim_resA[6][:,1:3],dims=2);
    surfsoil = sim_resA[2][:,1];
    
    #trueTB = get_TB_2(surfsoil, tsoil, tcan, canpot_true, laiM, vodA, vodB, vodC);
    trueTB = get_TB_2(surfsoil, tsoil, tcan, canpot_true, laiM, vodA, vodB, vodC,omega);
    return trueTB;
end

mydates = DateTime.(sim_res0[1].DateTime, "yyyy-mm-dd HH:MM:SS");
daylist = Date.(mydates)[1:24:end];

vod0 = getVOD(sim_res0,0.067,0.82,0.051,0.05);

#=
parslast = [  5.69853422,   0.64423276,   2.5652525 ,   7.75276212,
-0.58018373,   0.25647054,   3.11366709, -14.49991223,
-0.87973907,   0.40173571,   0.38305563,  -0.87020028];
=#
parslast = [  5.04098848,   0.66518166,   2.53994516,   7.7394941 ,
1.47762442,   0.60178479,   2.74880588, -14.87316503,
-0.87175587,   0.16035905,   0.36830784,  -0.81136998];

simresH = run_sim_1(convert(Array{FT},exp.(parslast))...);
#vodparslast = [0.02467557, 0.8141346 , 0.0479153 ,0.04905261];
vodparslast = [0.0397103 , 0.81123016, 0.05394451, 0.04797286];

vodH = getVOD(simresH,vodparslast...);




#sI = 1;
#simresH = run_sim_1(convert(Array{FT},exp.(lplist[sI]))...);
#vodH = getVOD(simresH,vodpars[sI]...);

#=
optim_res = optim_vod_stages2(obsH, obsV, canpot, soilsurf,
    tsoil, tcan, laiM, k0, alpha0, beta0,omega0, maxiter, stop_crit, 1e-5, 5);
=#


#=
sI = 2;
simres1 = run_sim_1(convert(Array{FT},exp.(lplist[sI]))...);
vod1 = getVOD(simres1,vodpars[sI]...);

sI = 3;
simres6 = run_sim_1(convert(Array{FT},exp.(lplist[sI]))...);
vod6 = getVOD(simres6,vodpars[sI]...);

sI = 4;
simres16 = run_sim_1(convert(Array{FT},exp.(lplist[sI]))...);
vod16 = getVOD(simres16,vodpars[sI]...);

figure();
plot(mydates,vod0[3][:,1],"k");
plot(mydates,vodH[3][:,1],"r");
plot(mydates,vod1[3][:,1],"b");

figure();
plot(mydates,vod0[1],"k");
plot(mydates,vodH[1],"r");
#%%
=#

tcan0 = (sim_res0[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;

#=
function myr2(ret,obs)
  return 1 - mean((ret-obs) .^ 2) / var(obs)
end
=#

function myr2(ret,obs)
  return sqrt(mean((ret-obs) .^ 2));
end

retrieved = deepcopy(vodH);

a0 = 0.2

fig,axes = subplots(2,2);
ax = axes[1,1];
x = retrieved[3][:,1];
y = vod0[3][:,1];
ax.plot(x,y,".",alpha=a0)
println(myr2(x,y))
ax.set_title("VOD")
ax.plot([minimum(x),maximum(x)],[minimum(x),maximum(x)],"r--")
#ax.set_ylabel("True model")
#ax.text(0.1,0.9,raw"R^2")

ax = axes[1,2];
x = sim_res0[2][:,1];
y = simresH[2][:,1];
ax.plot(x,y,".",alpha=a0)
println(myr2(x,y))
ax.set_title(raw"Surface soil moisture $(m^3/m^3$)")
ax.plot([minimum(x),maximum(x)],[minimum(x),maximum(x)],"r--")
#ax.set_ylabel("True model")

ax = axes[2,1];
x = retrieved[1];
y = vod0[1];
ax.plot(x,y,".",alpha=a0)
println(myr2(x,y))
ax.set_title(raw"$T_{BV}$ (K)")
ax.plot([minimum(x),maximum(x)],[minimum(x),maximum(x)],"r--")
#ax.set_xlabel("Retrieved")
#ax.set_ylabel("True model")

ax = axes[2,2];
x = retrieved[2];
y = vod0[2];
ax.plot(x,y,".",alpha=a0)
println(myr2(x,y))
ax.set_title(raw"$T_{BH}$ (K)")
ax.plot([minimum(x),maximum(x)],[minimum(x),maximum(x)],"r--")
#ax.set_xlabel("Retrieved")
#ax.set_ylabel("True model")
tight_layout()




