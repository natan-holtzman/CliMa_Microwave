#include("../get_home_dir.jl")

using DataFrames
using CSV
using PyPlot
using Random
using StatsBase
using Statistics

#include(string(PROJECT_HOME,"/assim_code/time_averaging.jl"))

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 1;
rcParams["font.size"] = 24;

pmean = CSV.read("dlwp3/DL_means2.csv",DataFrame, header=0);
p25 = CSV.read("dlwp3/DL_lower2.csv",DataFrame,header=0);
p75 = CSV.read("dlwp3/DL_upper2.csv",DataFrame,header=0);

ptrue = CSV.read("dlwp3/DL_true2.csv",DataFrame,header=0);


colors_list = ["tab:blue","tab:green","tab:orange","tab:red"];
label_list = ["All obs.", "1 AM/PM obs.", "6 AM/PM obs."];

year_x = collect(1:24);

#=
figure()
plot(year_x, ptrue[:,1],color="k",linewidth=3, label="Truth")
for i in 1:3
	plot(year_x, pmean[:,i],color=colors_list[i],linewidth=3, label=label_list[i])
end
for i in 1:3
	fill_between(year_x, p25[:,i], p75[:,i] , color=colors_list[i],alpha=0.33)
end
=#

figure()
plot(year_x, (ptrue[:,1] .- mean(ptrue[:,1])) / std(ptrue[:,1]),color="k",linewidth=3, label="Truth")
for i in 1:3
	plot(year_x, (pmean[:,i] .- mean(pmean[:,i])) / std(pmean[:,i]),color=colors_list[i],linewidth=3, label=label_list[i])
end

legend()
xlabel("Time of day")
ylabel("Leaf water potential, diurnal mean (MPa)")


#need to start at 0.3 soil moisture in year 3