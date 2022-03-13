#include("../get_home_dir.jl")

using DataFrames
using CSV
using PyPlot
using Random
using StatsBase

#include(string(PROJECT_HOME,"/assim_code/time_averaging.jl"))

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 1;
rcParams["font.size"] = 24;


mpa2mm = 10^6/9.8;
gravity_fac = (18.5+9)/2 * 1000/mpa2mm;

pmean = CSV.read("chain3_range/lwp_means.csv",DataFrame, header=0);
p25 = CSV.read("chain3_range/lwp_lower.csv",DataFrame,header=0) ;
p75 = CSV.read("chain3_range/lwp_upper.csv",DataFrame,header=0) ;

ptrue = CSV.read("chain3_range/lwp_true.csv",DataFrame,header=0) ;


colors_list = ["tab:blue","tab:green","tab:orange","tab:red"];
label_list = ["All obs.", "1 AM/PM obs.", "6 AM/PM obs."];

year_x = collect(1:(size(pmean)[1]))/365 .+ 2006;

figure()

plot(year_x, ptrue[:,1],color="k",linewidth=3, label="Truth")
for i in 1:3
	plot(year_x, pmean[:,i],color=colors_list[i],linewidth=3, label=label_list[i])
end

for i in 1:3
	fill_between(year_x, p25[:,i], p75[:,i] , color=colors_list[i],alpha=0.33)
end

legend()
xlabel("Time (years)")
ylabel("Leaf water potential (MPa)")
#xlim(2011, 2013.5)

#=

spmean = CSV.read("swp_means.csv",DataFrame, header=0);
sp25 = CSV.read("swp_lower.csv",DataFrame,header=0) ;
sp75 = CSV.read("swp_upper.csv",DataFrame,header=0) ;
sptrue = CSV.read("swp_true.csv",DataFrame,header=0) ;
#=
plot(year_x,ptrue[:,1],"k", label="True leaf")
#plot(year_x[1:(end-1)],(pmean[1:(end-1),1] + pmean[2:end,1])/2)
plot(year_x , sptrue[:,1] .- gravity_fac,"r:", label="True soil")

plot(year_x,pmean[:,1],"b", label="All leaf")
#plot(year_x[1:(end-1)],(pmean[1:(end-1),1] + pmean[2:end,1])/2)
plot(year_x , spmean[:,1] .- gravity_fac,":" , color="orange",label="All soil")
=#

figure()

plot(year_x, sptrue[:,1],color="k",linewidth=3, label="Truth")
for i in 1:3
	plot(year_x, spmean[:,i],color=colors_list[i],linewidth=3, label=label_list[i])
end

for i in 1:3
	fill_between(year_x, sp25[:,i], sp75[:,i] , color=colors_list[i],alpha=0.33)
end

legend()
xlabel("Time (years)")
ylabel("Soil water potential (MPa)")
#xlim(2011, 2013.5)


swmean = CSV.read("smc_means.csv",DataFrame, header=0);
sw25 = CSV.read("smc_lower.csv",DataFrame,header=0) ;
sw75 = CSV.read("smc_upper.csv",DataFrame,header=0) ;
swtrue = CSV.read("smc_true.csv",DataFrame,header=0) ;

figure()

plot(year_x, swtrue[:,1],color="k",linewidth=3, label="Truth")
for i in 1:3
	plot(year_x, swmean[:,i],color=colors_list[i],linewidth=3, label=label_list[i])
end

for i in 1:3
	fill_between(year_x, sw25[:,i], sw75[:,i] , color=colors_list[i],alpha=0.33)
end

legend()
xlabel("Time (years)")
ylabel("Soil moisture")
#xlim(2011, 2013.5)

=#

#need to start at 0.3 soil moisture in year 3