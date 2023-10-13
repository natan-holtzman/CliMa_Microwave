
using DataFrames
using CSV
using PyPlot
using Dates

include("../get_home_dir.jl")

pygui(true)

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 1;
rcParams["lines.markersize"] = 2;

rcParams["font.size"] = 12;
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

prior_min = [10, 0.75, 0.1, 500, 0.01,  0.01,1,0.5e-7,0.3,0.01,1.1,0.125];
prior_max = [300,15,  50, 3000, 10, 10,120,    20e-6,0.5,5,1.9,8];

prior_mean = 0.5*(log.(prior_min)+log.(prior_max));

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
#obsH0 = 1*vod0[1];
#obsV0 = 1*vod0[2];
trueVOD = 1*vod0[3][:,1];

figure();
plot(mydates,trueVOD);
xlim(daylist[200],daylist[210]);
ylim(0.775,1);


etdaily = get_daily(sim_res0[1].ETmod,24)*18/1000*24*60*60;
hiday = 160;
loday = 161;

colors_list = ["#E69F00","#56B4E9","#009E73","#CC79A7","#0072B2"];
obs_names = ["HOURLY", "1 AM/PM", "6 AM/PM","1+6","1+6 offset","True model"];

#=
figure()
hivod = trueVOD[(hiday-1)*24 .+ (1:24)];
lovod = trueVOD[(loday-1)*24 .+ (1:24)];
xtime = Array(0:23);

plot(xtime,hivod,color=colors_list[1],label="Hourly")
plot(xtime,lovod,"--",color=colors_list[1])

hsel = 2:12:24;
plot(xtime[hsel], hivod[hsel],"o-",color=colors_list[2],label="1 AM/PM")
plot(xtime[hsel], lovod[hsel],"o--",color=colors_list[2])


hsel = 7:12:24;
plot(xtime[hsel], hivod[hsel],"o-",color=colors_list[3],label="6 AM/PM")
plot(xtime[hsel], lovod[hsel],"o--",color=colors_list[3])


hsel = [2,7,14,19];
plot(xtime[hsel], hivod[hsel],"o-",color=colors_list[4],label="1+6 AM/PM")
plot(xtime[hsel], lovod[hsel],"o--",color=colors_list[4])


legend()
xlabel("Hour of day")
ylabel("VOD")
#plot the hourly in color of hourly
#then plot "-o" with other colors same as bar plots




allvod = trueVOD[(155-1)*24 .+ (1:(24*18))];
xtime = mydates[(155-1)*24 .+ (1:(24*18))];

figure()
plot(xtime,allvod,color=colors_list[1],label="Hourly")

hsel = zeros(size(xtime));
hsel[2:(24*3):end] .= 1;
hsel[(2+12):(24*3):end] .= 1;
hsel = hsel .== 1;
plot(xtime[hsel], allvod[hsel],"o--",color=colors_list[2],label="1 AM/PM")


hsel = zeros(size(xtime));
hsel[7:(24*3):end] .= 1;
hsel[(7+12):(24*3):end] .= 1;
hsel = hsel .== 1;
plot(xtime[hsel], allvod[hsel],"o--",color=colors_list[3],label="6 AM/PM")


hsel = zeros(size(xtime));
hsel[2:(24*3):end] .= 1;
hsel[(2+12):(24*3):end] .= 1;
hsel[7:(24*3):end] .= 1;
hsel[(7+12):(24*3):end] .= 1;
hsel = hsel .== 1;
plot(xtime[hsel], allvod[hsel],"o--",color=colors_list[4],label="1+6")


legend()
xlabel("Time")
ylabel("VOD")


=#

allvod = trueVOD[(155-1)*24 .+ (1:(24*18))];
xtime = mydates[(155-1)*24 .+ (1:(24*18))];


obsTB = CSV.read("obsTB_witherr_1.csv",DataFrame);
obsH = Array(obsTB.hpol);
obsV = Array(obsTB.vpol);

figure(figsize=(8,10))
subplot(6,1,1)
plot(xtime,allvod,color=colors_list[1],label="Hourly",markersize=4)
xlim(xtime[24*2],xtime[end-(24*7)])
ylim(0.77,1)
xticks([],[])
#title("Hourly")
title("(a) HOURLY VOD",loc="left")

subplot(6,1,2)
hsel = zeros(size(xtime));
hsel[2:(24*3):end] .= 1;
hsel[(2+12):(24*3):end] .= 1;
hsel = hsel .== 1;
plot(xtime[hsel], allvod[hsel],"o--",color=colors_list[2],label="1 AM/PM",markersize=4)
xlim(xtime[24*2],xtime[end-(24*7)])
ylim(0.77,1)
xticks([],[])
#title("1 AM/PM")
title("(b) 1 AM/PM VOD",loc="left")


subplot(6,1,3)
hsel = zeros(size(xtime));
hsel[7:(24*3):end] .= 1;
hsel[(7+12):(24*3):end] .= 1;
hsel = hsel .== 1;
plot(xtime[hsel], allvod[hsel],"o--",color=colors_list[3],label="6 AM/PM",markersize=4)
xlim(xtime[24*2],xtime[end-(24*7)])
ylim(0.77,1)
xticks([],[])
#title("6 AM/PM")
title("(c) 6 AM/PM VOD",loc="left")


subplot(6,1,4)
hsel = zeros(size(xtime));
hsel[2:(24*3):end] .= 1;
hsel[(2+12):(24*3):end] .= 1;
hsel[7:(24*3):end] .= 1;
hsel[(7+12):(24*3):end] .= 1;
hsel = hsel .== 1;
plot(xtime[hsel], allvod[hsel],"o--",color=colors_list[4],label="1+6",markersize=4)
xlim(xtime[24*2],xtime[end-(24*7)])
ylim(0.77,1)
xticks([],[])

title("(d) 1+6 VOD",loc="left")


etsamp  = sim_res0[1].ETmod[(155-1)*24 .+ (1:(24*18))] * (18/1000 * 60*60);

subplot(6,1,5)
plot(xtime,etsamp,"blue")
#ylim(0.77,1)
#title("Hourly")
title("(e) ET (mm/hour)",loc="left")
dtlocs = xtime[1:48:end];
DateTick = Dates.format.(dtlocs, "U d");
#xticks(dtlocs,DateTick);
xticks([],[])

xlim(xtime[24*2],xtime[end-(24*7)]);


subplot(6,1,6)
plot([],[],"k",alpha=0.75,label="Without noise")
plot([],[],"r",alpha=0.75,label="After adding noise")


plot(xtime,obsH[(155-1)*24 .+ (1:(24*18))],"r",alpha=0.75)
plot(xtime,vod0[1][(155-1)*24 .+ (1:(24*18))],"k",alpha=0.75)
#ylim(0.77,1)
#title("Hourly")
title(raw"(f) $T_{BH} (K)$",loc="left")
dtlocs = xtime[1:48:end];
DateTick = Dates.format.(dtlocs, "U d");
xticks(dtlocs,DateTick);
xlim(xtime[24*2],xtime[end-(24*7)]);
legend(loc=(0.3,1.1),ncol=2,fontsize=12)
xlabel("Time in 2007")

tight_layout()

savedir = "C:/Users/natan/OneDrive - Stanford/Documents/clima paper writing/coauthor edits/wrr reviews/fig_pdfs/";
savefig(string(savedir,"fig2.pdf"))
#supylabel("VOD")



