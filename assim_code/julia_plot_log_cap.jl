include("../get_home_dir.jl")

using DataFrames
using CSV
using PyPlot
using Random
using StatsBase

file1 = readlines("julia_test_1.log")[101:end];
file2 = readlines("julia_test_2.log")[101:end];
file3 = readlines("julia_test_3.log")[101:end];


function procline(L)
    return [parse(Float64, x) for x in split(L[2:(length(L)-1)], ", ")];
end
	
g1 = hcat([procline(x) for x in file1[27:27:end]]...);
g2 = hcat([procline(x) for x in file2[27:27:end]]...);
g3 = hcat([procline(x) for x in file3[27:27:end]]...);



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
    ax.plot(exp.(g1[j,:]),color="blue")
	ax.plot(exp.(g2[j,:]),color="orange")
	ax.plot(exp.(g3[j,:]),color="purple")

	
    ax.set_ylim(prior_min[j],prior_max[j])
    ax.set_title(par_names[j])
    ax.set_xticks([])
    ax.hlines(true_val[j], 0, length(g1[j,:]), color="black",alpha=0.75)
    global j += 1
end
tight_layout()
ax_all[2,4].axis("off")
=#

df_raw = CSV.read("../data/moflux_land_data_newnames_7.csv", DataFrame);


include("mironov.jl");
include("tau_omega_funs.jl");
#include("../simulation_code/sim_cap_nov4.jl");


#include("../simulation_code/new_capacitance/sim_vary_new_stomata.jl");
include("time_averaging.jl")

N = 48*365
istart = 48*365*2 - 52*48 #+ 230*48
soil0 = 0.39;

include("../simulation_code/rebuild_sim_Scrit_drain_plugin_ET_setsoil.jl");
function run_sim_old(vcmax_par::FT, k_frac::FT, k_plant::FT, k_soil::FT, z_soil::FT, weibB::FT, weibC::FT)
	return convert_sim(run_sim(vcmax_par, k_frac, k_plant, k_soil, z_soil, istart, N, soil0, FT(0), 0, weibB, weibC, df_raw));
end

mypars = convert(Array{FT}, [22, 0.33, 1, 5e-6, 800, 1.5, 1])

#sim_res0 = run_sim_old(mypars...);

include("../simulation_code/new_capacitance/sim_vary_new_stomata.jl");
function run_sim_new(vcmax_par::FT, k_frac::FT, k_plant::FT, k_soil::FT, z_soil::FT, weibB::FT, weibC::FT, scheme::Int, storage_mult::FT)
	return convert_sim(run_sim_vary(vcmax_par, k_frac, weibB, weibC, k_plant, k_soil, z_soil, istart, N, soil0, storage_mult, 1e-5, scheme, df_raw));
end

#sim_res0_x = run_sim_new(mypars...,1, FT(1));

sim_res0_c0 = run_sim_new(mypars...,3, FT(0.1));

sim_res0_c1 = run_sim_new(mypars...,3, FT(1));

sim_res0_c10 = run_sim_new(mypars...,3, FT(10));

#=
figure()
plot(get_diurnal(mean(sim_res0_c1[3],dims=2),24)); plot(get_diurnal(mean(sim_res0_c1[4],dims=2), 24)); plot(get_diurnal(sim_res0_c1[5], 24))

figure()
plot(get_diurnal(mean(sim_res0_x[3],dims=2),24)); plot(get_diurnal(mean(sim_res0_x[4],dims=2), 24)); plot(get_diurnal(sim_res0_x[5], 24))
=#

figure()
plot(sim_res0_c0[1].leafpot); plot(sim_res0_c1[1].leafpot); plot(sim_res0_c10[1].leafpot)

#=
plot(mean(sim_res0_c1[3],dims=2)[2:24:end]); plot(mean(sim_res0_c1[4],dims=2)[2:24:end]); plot(sim_res0_c1[5][2:24:end])
plot(mean(sim_res0_c1[3],dims=2)[13:24:end]); plot(mean(sim_res0_c1[4],dims=2)[13:24:end]); plot(sim_res0_c1[5][13:24:end])
plot(mean(sim_res0_c1[3],dims=2)[13:24:end] - sim_res0_c1[5][13:24:end])
=#

#=
noise_std = 1.3;
noise_var = noise_std^2

#trueTB = get_TB(sim_res0, mean(sim_res0[3],dims=2), 0.25, 1, 1/50);
#obsH = trueTB[1] + noise_std*randn(Float64, length(trueTB[1]));
#obsV = trueTB[2] + noise_std*randn(Float64, length(trueTB[1]));


partest = convert(Array{FT},exp.(g1[1:7,end]));
sim_res1 = run_sim_2(partest...);
#retTB = get_TB(sim_res1, mean(sim_res1[3],dims=2), g1[17:end,end]...);

partest = convert(Array{FT},exp.(g2[1:7,end]));
sim_res2 = run_sim_2(partest...);

partest = convert(Array{FT},exp.(g3[1:7,end]));
sim_res3 = run_sim_2(partest...);


#leaf, stem, trunk
#rzsm, ET
function format_out(sim_res)
	leaf = mean(sim_res[3],dims=2);
	branch = mean(sim_res[4],dims=2);
	trunk = sim_res[5];
	rzsm = get_daily(sim_res[2][:,1], 24);
	et = get_daily(sim_res[1].ETmod, 24);
	
	return leaf[2:24:end],leaf[13:24:end],branch[2:24:end],branch[13:24:end],trunk[2:24:end],trunk[13:24:end],rzsm,et
end

outvar_true = format_out(sim_res0);
outvar_1 = format_out(sim_res1);
outvar_2 = format_out(sim_res2);
outvar_3 = format_out(sim_res3);

out_names = ["Leaf","Leaf", "Branch","Branch","Trunk", "Trunk","RZSM","ET"];

j = 1
fig, ax_all = subplots(2,4,figsize=(12,8))
for ax in vec(ax_all)[1:8]
	ax.plot(outvar_true[j],color="black")
    ax.plot(outvar_1[j],color="blue",alpha=0.75)
	ax.plot(outvar_2[j],color="orange",alpha=0.75)
    ax.plot(outvar_3[j],color="purple",alpha=0.75)

    ax.set_title(out_names[j])
    #ax.set_xticks([])
    global j += 1
end
tight_layout()
=#



