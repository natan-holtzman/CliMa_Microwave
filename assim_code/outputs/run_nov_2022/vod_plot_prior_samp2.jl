include("../../../get_home_dir.jl")
using Pkg
Pkg.activate(string(PROJECT_HOME,"/feb_j171"));

using DataFrames
using CSV
#using PyPlot
using Dates


#pygui(true)

#=
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 1;
rcParams["lines.markersize"] = 2;

rcParams["font.size"] = 15;
rcParams["mathtext.default"] = "regular";
=#

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

prior_min = [10, 0.75, 0.1, 500, 0.01,  0.01,1,0.5e-7,0.3,0.01,1.1,0.125];
prior_max = [300,15,  50, 3000, 10, 10,120,    20e-6,0.5,5,1.9,8];

prior_mean = 0.5*(log.(prior_min)+log.(prior_max));

println("true run")
sim_res0 = run_sim_0(pars0...);

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
obsH0 = 1*vod0[1];
obsV0 = 1*vod0[2];
obsTB = CSV.read("obsTB_witherr_1.csv",DataFrame);
obsH = Array(obsTB.hpol);
obsV = Array(obsTB.vpol);

postpar_tab = Array(CSV.read("postpar_hourly.csv",DataFrame,header=false));
priorpar_tab = transpose(Array(CSV.read("prior_sample_inlog.csv",DataFrame,header=true)));

function solveVOD(sim_resA)
    canpot_true = mean(sim_resA[6][:,1:3],dims=2);

    tsoil = sim_resA[1].T_SOIL .+ 273.15;
    tcan = (sim_resA[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
    laiM = sim_resA[1].LAI_modis;
    canpot_true = mean(sim_resA[6][:,1:3],dims=2);
    surfsoil = sim_resA[2][:,1];
    vodopt = optim_vod_stages2(obsH, obsV, canpot_true, surfsoil,
                    tsoil, tcan, laiM,
                0.067, 0.82, 0.051, 0.05, 10000,10^-6, 1e-5, 5);
    #trueTB = get_TB_2(surfsoil, tsoil, tcan, canpot_true, laiM, vodA, vodB, vodC);
    trueTB = get_TB_2(surfsoil, tsoil, tcan, canpot_true, laiM, vodopt[1:4]...);
    return trueTB;
end

println("prior run")
#simresPrior = [run_sim_1(convert(Array{FT},exp.(px))...) for px in priorpars];
soilPrior = [];
radPrior = [];

for i in 1:120
    println(i)
    simI = run_sim_1(convert(Array{FT},exp.(priorpar_tab[i,:]))...);
    if sum(isnan.(simI[2][:,1])) == 0
        push!(soilPrior,simI[2][:,1]);
        radI = solveVOD(simI);
        push!(radPrior,radI);
    end
end

println("posterior run")
soilPost = [];
radPost = [];

for i in 1:120
    println(i)
    simI = run_sim_1(convert(Array{FT},exp.(postpar_tab[i,:]))...);
    push!(soilPost,simI[2][:,1]);
    radI = solveVOD(simI);
    push!(radPost,radI);
end

#tcan0 = (sim_res0[1].LW_OUT / 0.97 / K_STEFAN ) .^ 0.25;
#=
function myr2(ret,obs)
  return sqrt(mean((ret-obs) .^ 2));
end
=#
#function that takes in an ax, an averaging fun, and 3 input datasets

#=
function myplot3_old(ax,afun,data0,dataP,dataH,extract_fun,xvals)
    vP = hcat([afun(extract_fun(x),24) for x in dataP]...);
    vH = hcat([afun(extract_fun(x),24) for x in dataH]...);
    v0 = afun(extract_fun(data0),24)
    ax.plot(xvals,vP,"b",alpha=0.67,label="Prior")
    ax.plot(xvals,vH,"r",alpha=0.67,label="HOURLY retrieval")
    ax.plot(xvals,v0,"k",label="True model")
    ax.set_xlim(xvals[1],xvals[end])
end

function row_quant(x,q)
    return [quantile(x[i,:],q) for i in 1:size(x)[1]]
end


function myplot3(ax,afun,data0,dataP,dataH,extract_fun,xvals)
    vP = hcat([afun(extract_fun(x),24) for x in dataP]...);
    vH = hcat([afun(extract_fun(x),24) for x in dataH]...);
    v0 = afun(extract_fun(data0),24)
    ax.plot(xvals,v0,"k--",label="True model")
    ax.fill_between(xvals,row_quant(vP,0.25),row_quant(vP,0.75),color="b",alpha=0.5,label="Prior",linewidth=1)
    ax.fill_between(xvals,row_quant(vH,0.25),row_quant(vH,0.75),color="r",alpha=0.5,label="HOURLY retrieval",linewidth=1)
    ax.set_xlim(xvals[1],xvals[end])
end
=#

function pickVOD(x)
    return x[3][:,1]
end

function pickTbh(x)
    return x[1]
end

function pickTbv(x)
    return x[2]
end

function pickSSM(x)
    return x
end

function diurnal_nomean(x,d)
    y = get_diurnal(x,d)
    return y; # .- mean(y)
  end
#=
vodP = deepcopy(radPrior);
vodH = deepcopy(radPost);


dayser = 1:365;
hser = 0:23;
println("plotting")
fig, axes = subplots(4,2,figsize=(11,15))
axes[1,1].plot([],[],"k--",label="True model")
axes[1,1].plot([],[],"b",alpha=0.67,label="Prior distribution")
axes[1,1].plot([],[],"r",alpha=0.67,label="HOURLY retrieval")
fig.legend(loc="upper center",bbox_to_anchor=(0.5,0.97),ncol=3)

myplot3(axes[1,1],get_daily,vod0,vodP,vodH,pickVOD,dayser)
myplot3(axes[1,2],diurnal_nomean,vod0,vodP,vodH,pickVOD,hser)
axes[1,1].set_ylabel("VOD")
axes[1,1].set_title("Daily means")
axes[1,2].set_title("Mean diurnal cycle")


ssm_true = sim_res0[2][:,1]
myplot3(axes[2,1],get_daily,ssm_true,soilPrior,soilPost,pickSSM,dayser)
myplot3(axes[2,2],diurnal_nomean,ssm_true,soilPrior,soilPost,pickSSM,hser)

smclab = string("Surface soil\nmoisture ",raw"$(m^3/m^3)$");

axes[2,1].set_ylabel(smclab)


myplot3(axes[3,1],get_daily,vod0,vodP,vodH,pickTbh,dayser)
myplot3(axes[3,2],diurnal_nomean,vod0,vodP,vodH,pickTbh,hser)
axes[3,1].set_ylabel(raw"$T_{BH}$ (K)")

myplot3(axes[4,1],get_daily,vod0,vodP,vodH,pickTbv,dayser)
myplot3(axes[4,2],diurnal_nomean,vod0,vodP,vodH,pickTbv,hser)
axes[4,1].set_ylabel(raw"$T_{BV}$ (K)")

axes[4,1].set_xlabel("Day of year 2007")
axes[4,2].set_xlabel("Hour of day")
tight_layout()

fig.subplots_adjust(top=0.9)

fig.savefig("testfig_mar22e.png")
=#

dataP = deepcopy(radPrior);
dataH = deepcopy(radPost);

all_vod_P =  hcat([pickVOD(x) for x in dataP]...);
all_vod_H =  hcat([pickVOD(x) for x in dataH]...);
CSV.write("vod_prior.csv",DataFrame(all_vod_P, :auto));
CSV.write("vod_hourly.csv",DataFrame(all_vod_H, :auto));


all_vod_P =  hcat([pickTbh(x) for x in dataP]...);
all_vod_H =  hcat([pickTbh(x) for x in dataH]...);
CSV.write("tbh_prior.csv",DataFrame(all_vod_P, :auto));
CSV.write("tbh_hourly.csv",DataFrame(all_vod_H, :auto));


all_vod_P =  hcat([pickTbv(x) for x in dataP]...);
all_vod_H =  hcat([pickTbv(x) for x in dataH]...);
CSV.write("tbv_prior.csv",DataFrame(all_vod_P, :auto));
CSV.write("tbv_hourly.csv",DataFrame(all_vod_H, :auto));


CSV.write("ssm_prior.csv",DataFrame(soilPrior, :auto));
CSV.write("ssm_hourly.csv",DataFrame(soilPost, :auto));







