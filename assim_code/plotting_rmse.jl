#include("../get_home_dir.jl")

using DataFrames
using CSV
using PyPlot
using Random
using StatsBase

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 2;
rcParams["font.size"] = 24;
rcParams["mathtext.default"] = "regular";


#leafRMSE =  [0.04030680922286562, 0.1114426954163669, 0.19222915929366224]
#etRMSE = [5.11686931714865e-05, 3.992463510637098e-05, 5.691494713989705e-05]*18.02/1000*60*60*24;
#rzRMSE = [0.0030338018951135853, 0.003058330701815462, 0.003075149421206762];

leafRMSE = [0.04030680922286563, 0.11144269541636687, 0.19222915929366224];
etRMSE = [5.11686931714865e-05, 3.992463510637098e-05, 5.691494713989705e-05];
rzRMSE = [0.0012556391476988972, 0.003711247616340684, 0.0016033746735104248];


mylabels = ["All", "1 AM/PM", "6 AM/PM"];

#=
figure()
subplot(131)
bar([1,2,3],leafRMSE)
xticks([1,2,3], mylabels)
title("Leaf water potential (MPa)")

subplot(132)
bar([1,2,3],rzRMSE)
xticks([1,2,3], mylabels)
title("Column-average soil moisture")

subplot(133)
bar([1,2,3],etRMSE)
xticks([1,2,3], mylabels)
title(L"Daily ET $(mm/m^2/s)$")
ticklabel_format(style="plain",axis="y")

#tight_layout()
=#

gsRMSE = [0.0009707089849704303, 0.000809590695874667, 0.0012329717941626317]
leafscaleRMSE = [0.02781028047973449, 0.01783460829728346, 0.03495551737372624]
betaRMSE = [0.01874038306283515, 0.013533028689382927, 0.023987804499910826]

figure()
subplot(231)
bar([1,2,3],leafRMSE)
xticks([1,2,3], mylabels)
title("Leaf water potential (MPa)")

subplot(232)
bar([1,2,3],leafscaleRMSE)
xticks([1,2,3], mylabels)
title("LWP divided by P63")
ticklabel_format(style="plain",axis="y")

subplot(233)
bar([1,2,3],betaRMSE)
xticks([1,2,3], mylabels)
title("Beta factor")
ticklabel_format(style="plain",axis="y")

subplot(234)
bar([1,2,3],gsRMSE)
xticks([1,2,3], mylabels)
title(L"Stomatal conductance $(mol/m^2/s)$")

subplot(235)
bar([1,2,3],etRMSE)
xticks([1,2,3], mylabels)
title(L"Daily ET $(mm/m^2/s)$")
ticklabel_format(style="plain",axis="y")

leafmean = convert(Array, CSV.read("leaf_mean_tab.csv",DataFrame,header=false));
leafscale = convert(Array, CSV.read("leaf_scale_tab.csv",DataFrame,header=false));

figure()
for i in 1:3
	plot(leafmean[i,:], 1 ./ leafscale[i,:],"o",label=mylabels[i])
end
plot(leafmean[1,end], 1 ./ leafscale[1,end],"k+",label="True model",markersize=25,markeredgewidth=3)
xlabel("Mean LWP (MPa)")
ylabel("Stomatal P63 (MPa)")
legend()
