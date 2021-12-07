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