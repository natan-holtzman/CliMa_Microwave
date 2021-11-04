
using DataFrames
using CSV
using PyPlot
using Random
using StatsBase

rawfile = readlines("log_2x3.log")[101:end];

function procline(L)
    return [parse(Float64, x) for x in split(L[2:(length(L)-1)], ", ")];
end
	
g1 = hcat([procline(x) for x in rawfile[27:27:end]]...);

prior_min = [10, 0.01, 0.1,  1e-6, 500, 0.75, 0.75];
prior_max = [120,0.75,  100, 2e-5, 3000, 10, 8];

true_val = [22, 0.33, 15, 1e-5, 800, 1.5, 1];
par_names = ["Vcmax", "Stomatal margin", "Kmax_plant", "Kmax_soil", "Soil depth", "Kplant location", "Kplant shape"];


rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 1;
rcParams["font.size"] = 20;

#=
j = 1
fig, ax_all = subplots(2,4,figsize=(12,8))
for ax in vec(ax_all)[1:7]
    ax.plot(exp.(g1[j,:]))
    ax.set_ylim(prior_min[j],prior_max[j])
    ax.set_title(par_names[j])
    ax.set_xticks([])
    ax.hlines(true_val[j], 0, length(g1[j,:]), color="red")
    global j += 1
end
tight_layout()
ax_all[2,4].axis("off")
=#



include("mironov.jl");
include("tau_omega_funs.jl");
include("../simulation_code/rebuild_sim_Scrit_drain_plugin_ET_setsoil.jl");
include("time_averaging.jl")

N = 48*365
istart = 48*365*2 - 52*48 #+ 230*48
soil0 = 0.39;

function run_sim_2(vcmax_par::FT, k_frac::FT, k_plant::FT, k_soil::FT, z_soil::FT, weibB::FT, weibC::FT)
	return convert_sim(run_sim(vcmax_par, k_frac, k_plant, k_soil, z_soil, istart, N, soil0, FT(0), 0, weibB, weibC));
end

sim_res1 = run_sim_2(FT(22),FT(0.33), FT(15.0),FT(1e-5), FT(800),FT(1.5),FT(1));

noise_std = 1.3;
noise_var = noise_std^2

trueTB = get_TB(sim_res1, mean(sim_res1[3],dims=2), 0.25, 1, 1/50);
obsH = trueTB[1] + noise_std*randn(Float64, length(trueTB[1]));
obsV = trueTB[2] + noise_std*randn(Float64, length(trueTB[1]));

partest = convert(Array{FT},exp.(g1[8:14,1]));
sim_res2 = run_sim_2(partest...);

