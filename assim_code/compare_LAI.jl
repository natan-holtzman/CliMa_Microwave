include("../get_home_dir.jl")

using DataFrames
using CSV
using PyPlot
using Random
using StatsBase

df_raw = CSV.read("../data/moflux_land_data_newnames_7.csv", DataFrame);

include("time_averaging.jl")

N = 48*365*2
istart = 48*365*2 - 52*48 #+ 230*48
#istart = 48*365*3 - 52*48 #+ 230*48

soil0 = 0.39;

#=
include("../simulation_code/sim_vary_new_stomata.jl");
function run_sim_new(vcmax_par::FT, k_frac::FT, k_plant::FT, k_soil::FT, z_soil::FT, weibB::FT, weibC::FT, storage_mult::FT)
	return convert_sim(run_sim_vary(vcmax_par, k_frac, weibB, weibC, k_plant, k_soil, z_soil, istart, N, soil0, storage_mult, 1e-5, 3, df_raw));
end
=#

include("../simulation_code/sim_vary_new_stomata_med.jl");
function run_sim_new(vcmax_par::FT, k_frac::FT, k_plant::FT, k_soil::FT, z_soil::FT, weibB::FT, weibC::FT, storage_mult::FT)
	return convert_sim(run_sim_vary(vcmax_par, k_frac, weibB, weibC, k_plant, k_soil, z_soil, istart, N, soil0, storage_mult, 1e-5, 3, df_raw, FT(600)));
end


mypars1 = convert(Array{FT}, [22, 0.2, 1, 1e-5, 500, 4.0, 2, 1]);
mypars_k0 = 1*mypars1; mypars_k0[3] *= 0.5;
mypars_k2 = 1*mypars1; mypars_k2[3] *= 2;
mypars_v0 = 1*mypars1; mypars_v0[8] *= 0.2;# mypars_v0[3] /= 0.2;
mypars_v2 = 1*mypars1; mypars_v2[8] *= 5; #mypars_v2[3] /= 5;
#mypars_c0 = 1*mypars1; mypars_c0[7] *= 0.5;
#mypars_c2 = 1*mypars1; mypars_c2[7] *= 2;
#mypars_b0 = 1*mypars1; mypars_b0[6] *= 0.5;
#mypars_b2 = 1*mypars1; mypars_b2[6] *= 2;


#par_list = [mypars1, mypars_k0, mypars_k2, mypars_v0,mypars_v2, mypars_c0, mypars_c2, mypars_b0, mypars_b2];
par_list = [mypars1, mypars_k0, mypars_k2, mypars_v0,mypars_v2];

sim_res_all = [run_sim_new(x...) for x in par_list];

surf_list = hcat([get_daily(x[2][:,1],24) for x in sim_res_all]...);
canopy_list = hcat([mean(x[3],dims=2) for x in sim_res_all]...);

figure()
plot(surf_list)

figure()
plot(canopy_list[2:24:end,:])

figure()
plot(canopy_list[13:24:end,:])

function getR2(a,b)
	return 1 - sum((a-b).^2) / sum((b .- mean(b)).^2)
end

diurnal_y1 = hcat([get_diurnal(x[1].leafpot[1:(365*24)],24) for x in sim_res_all]...);
diurnal_y2 = hcat([get_diurnal(x[1].leafpot[(365*24+1):end],24) for x in sim_res_all]...);

