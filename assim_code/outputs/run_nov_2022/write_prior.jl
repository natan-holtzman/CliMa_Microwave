include("../../../get_home_dir.jl")
using Pkg
Pkg.activate(string(PROJECT_HOME,"/feb_j171"));
#Pkg.activate(string(PROJECT_HOME,"/land_v011_env"));


using DataFrames
using CSV
#using Plots
using Random
using StatsBase
using Distributions

prior_min = [10, 0.75, 0.1, 500, 0.01,  0.01,1,0.5e-7,0.3,0.01,1.1,0.125];
prior_max = [300,15,  50, 3000, 10, 10,120,    20e-6,0.5,5,1.9,8];

lp_min = log.(prior_min);
lp_max = log.(prior_max);
lp_range = lp_max - lp_min;

prior_dist = Product(Uniform.(lp_min, lp_max));

prior_params = [rand(prior_dist) for j in 1:120];
println(size(prior_params))

CSV.write("prior_sample_inlog.csv", DataFrame(prior_params,:auto));

