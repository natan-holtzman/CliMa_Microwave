using DataFrames
using CSV
using PyPlot
using Random
using StatsBase
using Statistics
using Dates

include("../get_home_dir.jl")
include(string(PROJECT_HOME,"/assim_code/time_averaging.jl"))


pygui(true)

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 2;
rcParams["font.size"] = 18;
rcParams["mathtext.default"] = "regular";


pall = convert(Array,CSV.read("../../files_april17/runs_may30/post_par_2slope_ac1.csv",DataFrame));
p1 = convert(Array,CSV.read("../../files_april17/runs_may30/post_par_2slope_1c1.csv",DataFrame));

pall_samp = exp.(pall[6001:100:end,1:12]);
p1_samp = exp.(p1[6001:100:end,1:12]);

refpot_all = transpose(pall_samp[:,12]);
refpot_1 = transpose(p1_samp[:,12]);


leaf_6_all = convert(Array,CSV.read("../../files_april17/runs_may30/postLeaf_2slope_ac1.csv",DataFrame));
leaf_6_c1 = convert(Array,CSV.read("../../files_april17/runs_may30/postLeaf_2slope_1c1.csv",DataFrame));


istart = (365*2+180)*24;
iend = (365*2+190)*24;

mpa2mm = 10^6/9.8;
grav_pot = 13.5*1000/mpa2mm

df_raw = CSV.read(string(PROJECT_HOME,"/data/moflux_land_data_skipyear_hourly2.csv"), DataFrame);
df_raw[!,"RAIN"] *= 2; #this is because rain was originally in mm/half hour time step

mydates = DateTime.(df_raw.DateTime, "yyyy-mm-dd HH:MM:SS");
mydates_onlydate = Date.(mydates);

leafnorm_all = (leaf_6_all[:,1:40] .+ grav_pot) .* refpot_all .- grav_pot;
leafnorm_1 = (leaf_6_c1[:,1:40] .+ grav_pot) .* refpot_1 .- grav_pot;

figure()
plot(istart:iend, leafnorm_all[istart:iend,1:40] ,color="tab:blue",alpha=0.25);
#plot(istart:iend, leafnorm_1[istart:iend,1:40] ,color="tab:orange",alpha=0.25);
plot([],[],color="tab:blue",alpha=0.67,label="Retrieval");
#plot([],[],color="tab:orange",alpha=0.67,label="1 AM/PM retrieval");
plot(istart:iend,leaf_6_all[istart:iend,end],"k",label="True model");
xticks(istart:(24*2):iend .+ 1,mydates_onlydate[istart:(24*2):iend .+ 1]);#, mydates)
xlim(istart+18,istart+18+24*5);
ylabel("Leaf water potential (MPa)");
xlabel("Time");
legend();
