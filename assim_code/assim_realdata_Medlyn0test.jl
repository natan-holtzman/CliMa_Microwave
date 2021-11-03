
using DataFrames
using CSV
#using Plots
using Random
using StatsBase

include("mironov.jl");

#include("../simulation_code/rebuild_sim_medlyn_0g0.jl");
include("../simulation_code/rebuild_sim_medlyn_0g0_plugin_ET.jl");


include("time_averaging.jl")

#N = 48*3
#istart = 48*365*2 - 52*48 + 230*48


N = 48*365*2
istart = 48*365*2 - 52*48 #+ 230*48
soil0 = 0.39;

sim_res1 = convert_sim(run_sim(FT(22),FT(1.5), FT(0.15), FT(1e-4), FT(2.5), FT(2.0), FT(700),FT(0.45),istart,N, soil0, FT(0.05), FT(300), FT(2)));

obsET = sum(reshape(sim_res1[1].LE/44200, (24,:)),dims=1)[1,:]/24;
obsSMC = sum(reshape(sim_res1[1].SMC, (24,:)),dims=1)[1,:]/24;

using Distributions
using LabelledArrays
using LinearAlgebra

####################


function log_p_nolabel(a,  ET_var, SMC_var, oldET, oldSMC)

	v = LVector(logVcmax=a[1], logScrit=a[2], logKmaxPlant=a[3], logKmaxSoil=a[4], logBsoil=a[5],
			 logP20=a[6], logZsoil=a[7], logNsoil=a[8], logSlope=a[9], logG1=a[10], logWeibC=a[11]);

	p = 0.0;
	p += logpdf(Uniform(log(10), log(120)), v.logVcmax);
	p += logpdf(Uniform(log(0.1), log(100)), v.logKmaxPlant);
	p += logpdf(Uniform(log(0.075), log(0.3)), v.logScrit);
	#p += logpdf(Uniform(log(0.75), log(5)), v.logPcrit);
        p += logpdf(Uniform(log(1e-5), log(2e-3)), v.logKmaxSoil);
        p += logpdf(Uniform(log(1), log(8)), v.logBsoil);
        p += logpdf(Uniform(log(0.25), log(5)), v.logP20);
        p += logpdf(Uniform(log(0.4), log(0.6)), v.logNsoil);
        p += logpdf(Uniform(log(250), log(3000)), v.logZsoil);
        p += logpdf(Uniform(log(0.005), log(0.5)), v.logSlope);
 p += logpdf(Uniform(log(30), log(650)), v.logG1);
 p += logpdf(Uniform(log(0.75), log(8)), v.logWeibC);
		
	p0 = p+0 #assumes both previous and proposed parameters are in prior range 
        p0 += logpdf(MvNormal(obsET, sqrt(ET_var)), oldET);
        p0 += logpdf(MvNormal(obsSMC,sqrt(SMC_var)), oldSMC);


	try		
		sim_res = convert_sim(run_sim(FT(exp(v.logVcmax)),FT(exp(v.logScrit)), 
									  FT(exp(v.logKmaxPlant)), FT(exp(v.logKmaxSoil)), 
									  FT(exp(v.logBsoil)), FT(exp(v.logP20)), 
									FT(exp(v.logZsoil)), FT(exp(v.logNsoil)),
						 istart, N, soil0, FT(exp(v.logSlope)),
							FT(exp(v.logG1)), FT(exp(v.logWeibC))));

	if isnan(mean(sim_res[1].leafpot))
		return -Inf, 0,0,p0,0

	else
	
		simET = sum(reshape(sim_res[1].ETmod, (24,:)),dims=1)[1,:]/24;
		simSMC = sum(reshape(sim_res[2][:,1], (24,:)),dims=1)[1,:]/24;

		p += logpdf(MvNormal(obsET, sqrt(ET_var)), simET);
		p += logpdf(MvNormal(obsSMC,sqrt(SMC_var)), simSMC);

		return p, simET, simSMC, p0, sim_res
	end

	catch err_occurred
                return -Inf,0,0, p0,0
        end

end


function log_p_err(error_var, new_error_var, pred_ts, obs_ts)
	p0 = 0;
	p1 = 0;
	
	p0 += logpdf(truncated(Cauchy(), 0, Inf), error_var);
	p1 += logpdf(truncated(Cauchy(), 0, Inf), new_error_var);

	p0 -= 0.5*sum( (pred_ts-obs_ts) .^2 ./ error_var .+ log(error_var)) 
        p1 -= 0.5*sum( (pred_ts-obs_ts) .^2 ./ new_error_var .+ log(new_error_var))
	

return p0, p1
end


using PyPlot

init_name = "init_calib_Medlyn0.csv"
ipar_LL = Array(CSV.read(init_name, DataFrame));

NPAR = 11;

#=
best_i = argmax(ipar_LL[:,end]);
a01 = ipar_LL[1,1:NPAR];
ll0, et0, smc0, p_obs, leaf0 =  log_p_nolabel(a01,  0.001^2, 0.05^2, obsET, obsSMC);

a01 = ipar_LL[best_i,1:NPAR];
ll1, et1, smc1, p_obs, leaf1 =  log_p_nolabel(a01,  0.001^2, 0.05^2, obsET, obsSMC);

worst_i = argmin(ipar_LL[:,end]);
a01 = ipar_LL[worst_i,1:NPAR];
ll2, et2, smc2, p_obs, leaf2 =  log_p_nolabel(a01,  0.001^2, 0.05^2, obsET, obsSMC);

figure()
plot(obsET,"k",label="Obs"); plot(et0,color="green",label="OK"); plot(et1,color="blue",label="Best"); plot(et2,color="red", label="Worst")
legend()
xlabel("Day of year")
ylabel("ET (mol/m2/s)")

figure()
plot(obsSMC,"k",label="Obs"); plot(smc0,color="green",label="OK"); plot(smc1,color="blue",label="Best"); plot(smc2,color="red", label="Worst")
legend()
xlabel("Day of year")
ylabel("Surface SMC")
=#

a01 =log.([30,0.2,20,5e-4,2.5,2.0,1000,0.45,0.005,600,2])
#potential of -0.17 eq to 
a01 =log.([30,0.2,20,5e-4,2.5,2.0,1000,0.45,0.005,600,1.5]);


ll1, et1, smc1, p_obs, simres1 =  log_p_nolabel(a01,  0.001^2, 0.05^2, obsET, obsSMC);
sim_rz_day = simres1[2][:,4][1:24:end];
et_ratio = obsET ./ et1;
lai_daily=simres1[1].LAI_modis[1:24:end];
et_ratio[lai_daily .< 5] .= NaN;
vpd_noon=simres1[1].VPD[12:24:end];

using GLM
df1 = DataFrame(logratio = log.(et_ratio)[.! isnan.(et_ratio)], vpd=vpd_noon[.! isnan.(et_ratio)], rz=sim_rz_day[.! isnan.(et_ratio)], lai=lai_daily[.! isnan.(et_ratio)]);
mod1 = lm(@formula(logratio~vpd),df1[66:end,:]);
df1.ResVPD = df1.logratio - predict(mod1,df1);
rs13 = (0.13-0.05)/(0.45-0.05);
psi_sat = FT(2 / rs13^(-2.5));


df1.rzpot = ((df1.rz .- 0.05)/(0.45-0.05)) .^ (-2.5) * -psi_sat;
#=

pot_arr = collect(-0.35:0.01:-0.01);
weib = exp.(-(pot_arr/-0.17).^1.5);
plot(df1.rzpot, df1.ResVPD .- 0.34,"o")
plot(pot_arr,log.(weib));

rz_arr = (pot_arr/ -psi_sat) .^ (1/-2.5) * 0.4 .+ 0.05;

#%%
plot(df1.rzpot, exp.(df1.ResVPD .- 0.34),"o")
plot(pot_arr,weib);

plot(df1.rz, exp.(df1.ResVPD .- 0.34),"o")
plot(rz_arr,weib);
=#





