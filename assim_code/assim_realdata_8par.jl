
using DataFrames
using CSV
#using Plots
using Random
using StatsBase

include("mironov.jl");

include("../simulation_code/rebuild_sim_N.jl");




function avg2(x)
if (typeof(x[1])==Float64) | (typeof(x[1])==Float32)
        return (x[1:2:end] + x[2:2:end])/2;
else
        return x[1:2:end]
end

end


function avg2_2d(x)

return (x[1:2:end,:] + x[2:2:end,:])/2;

end


function avg2_df(x)

ncol = Integer(size(x)[2]);
nrow = Integer(size(x)[1]/2);
ans = deepcopy(x[1:nrow,:]);

for j in 1:ncol
#println(x[1:10,j])
#println(avg2(x[1:10,j]))
        ans[:,j] = avg2(x[:,j]);
end

return ans;

end

function convert_sim(x)
a = avg2_df(x[1]);
b = avg2_2d(x[2]);
c = avg2_2d(x[3]);
d = avg2_2d(x[4]);

return a,b,c,d
end




#N = 48*3
#istart = 48*365*2 - 52*48 + 230*48


N = 48*365
istart = 48*365*2 - 52*48 #+ 230*48
soil0 = 0.39;

sim_res1 = convert_sim(run_sim(FT(22),FT(1.5), FT(4.0), FT(1e-4), FT(2.5), FT(2.0), FT(700),FT(0.45),istart,N, soil0));

obsET = sum(reshape(sim_res1[1].LE/44200, (24,:)),dims=1)[1,:]/24;
obsSMC = sum(reshape(sim_res1[1].SMC, (24,:)),dims=1)[1,:]/24;

using Distributions
using LabelledArrays
using LinearAlgebra

####################


function log_p_nolabel(a,  ET_var, SMC_var, oldET, oldSMC)

	v = LVector(logVcmax=a[1], logPcrit=a[2], logKmaxPlant=a[3], logKmaxSoil=a[4], logBsoil=a[5], logP20=a[6], logZsoil=a[7], logNsoil=a[8]);

	p = 0.0;
	p += logpdf(Uniform(log(10), log(120)), v.logVcmax);
	p += logpdf(Uniform(log(0.1), log(100)), v.logKmaxPlant);
	p += logpdf(Uniform(log(0.75), log(5)), v.logPcrit);
        p += logpdf(Uniform(log(1e-5), log(2e-3)), v.logKmaxSoil);
        p += logpdf(Uniform(log(1), log(8)), v.logBsoil);
        p += logpdf(Uniform(log(0.25), log(5)), v.logP20);
        p += logpdf(Uniform(log(0.4), log(0.6)), v.logNsoil);
        p += logpdf(Uniform(log(500), log(3000)), v.logZsoil);

	p0 = p+0 #assumes both previous and proposed parameters are in prior range 
        p0 += logpdf(MvNormal(obsET, sqrt(ET_var)), oldET);
        p0 += logpdf(MvNormal(obsSMC,sqrt(SMC_var)), oldSMC);




	try		
		sim_res = convert_sim(run_sim(FT(exp(v.logVcmax)),FT(exp(v.logPcrit)), FT(exp(v.logKmaxPlant)), FT(exp(v.logKmaxSoil)), FT(v.logBsoil), FT(v.logP20), FT(exp(v.logZsoil)), 
			  FT(exp(v.logNsoil)), istart, N, soil0));

	if isnan(mean(sim_res[1].leafpot))
		return -Inf, 0,0,p0

	else
	
		simET = sum(reshape(sim_res[1].ETmod, (24,:)),dims=1)[1,:]/24;
		simSMC = sum(reshape(sim_res[2][:,1], (24,:)),dims=1)[1,:]/24;

		p += logpdf(MvNormal(obsET, sqrt(ET_var)), simET);
		p += logpdf(MvNormal(obsSMC,sqrt(SMC_var)), simSMC);

		return p, simET, simSMC, p0
	end

	catch err_occurred
                return -Inf,0,0, p0
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



function runAMH(x_init, niter, burnlen)

k0 = exp(rand(Uniform(log(0.01),log(0.75))));
alpha0 = exp(rand(Uniform(log(0.25),log(2))));
beta0 = exp(rand(Uniform(log(1/100),log(1/25))));

#pars = log.([FT(22),FT(1.5), FT(4.0), FT(1e-4), FT(2.5), FT(2.0), FT(700),FT(0.45)]);
#err0 = 1.7;


pars = Array(x_init);

smc_err0 = 0.05^2;
ET_err0 = 0.001^2;

nPar= length(pars);

ll0,et0,smc0,ll0prev = log_p_nolabel(pars,ET_err0,smc_err0,obsET,obsSMC);
mycov = I*0.01/length(pars);
#npar, ll0, llp, npar, 4
chain = zeros(niter, nPar*2+4);

for j in 1:niter
	if j % 2 == 0
		println(j)
	end

        if j % 25 == 0
                println(chain[j-1,:])
        end

	chain[j,1:nPar] = pars;
	chain[j,nPar*2+1] = ll0

	if j > burnlen
		mycov = (1-10^-6)*cov(chain[(j-burnlen):j,1:nPar])  + 10^-6*I;
	end


    parsP = rand(MvNormal(pars, mycov),1);
    llP,etP,smcP,ll0prev = log_p_nolabel(parsP,ET_err0 ,smc_err0,et0,smc0 );
   chain[j,(nPar+1):(nPar*2)] = parsP;
	chain[j,nPar*2+2] = llP
    mhratio = exp(llP - ll0prev);
    randI = rand();
    if randI < mhratio
		ll0 = llP*1;
		pars = parsP[:,1]*1;
		et0 = etP*1;
		smc0 = smcP*1;
	end
	
	chain[j,nPar*2+3] = smc_err0
	 chain[j,nPar*2+4] = ET_err0 

	smc_err_prop = rand(LogNormal(log(smc_err0), 0.01));
	p0_new, p1_new = log_p_err(smc_err0, smc_err_prop, smc0,obsSMC);
	mhratio_err = exp(p1_new - p0_new) * smc_err_prop/smc_err0;
	randIe = rand();
	if randIe < mhratio_err
		smc_err0 = smc_err_prop;
	end

        ET_err_prop = rand(LogNormal(log(ET_err0), 0.01));
        p0_new, p1_new = log_p_err(ET_err0, ET_err_prop,et0, obsET);
        mhratio_err = exp(p1_new - p0_new) * ET_err_prop/ET_err0;
        randIe = rand();
        if randIe < mhratio_err
                ET_err0 = ET_err_prop;
        end

end

return chain
end

function init_works()

k0 = exp(rand(Uniform(log(0.01),log(0.75))));
alpha0 = exp(rand(Uniform(log(0.25),log(2))));
beta0 = exp(rand(Uniform(log(1/100),log(1/25))));

smc_err0 = 0.05^2;
ET_err0 = 0.001^2;

x0 = LVector(logVcmax=rand(Uniform(log(10), log(120))),
			 logPcrit=rand(Uniform(log(0.75), log(5))),
			 logKmaxPlant=rand(Uniform(log(0.1), log(100))),
			 logKmaxSoil=rand(Uniform(log(1e-5), log(2e-3))), 
 logBsoil=rand(Uniform(log(1), log(8))),
 logP20=rand(Uniform(log(0.25), log(5))),
 logZsoil=rand(Uniform(log(500), log(3000))),
 logNsoil=rand(Uniform(log(0.4), log(0.6)))
			 );

ll_init,a2,a3,a4 = log_p_nolabel(Array(x0),ET_err0,smc_err0,obsET,obsSMC);

if isnan(ll_init)

while isnan(ll_init)



x0 = LVector(logVcmax=rand(Uniform(log(10), log(120))),
                         logPcrit=rand(Uniform(log(0.75), log(5))),
                         logKmaxPlant=rand(Uniform(log(0.1), log(100))),
                         logKmaxSoil=rand(Uniform(log(1e-5), log(2e-3))),
 logBsoil=rand(Uniform(log(1), log(8))),
 logP20=rand(Uniform(log(0.25), log(5))),
 logZsoil=rand(Uniform(log(500), log(3000))),
 logNsoil=rand(Uniform(log(0.4), log(0.6)))
                         );


ll_init,a2,a3,a4 = log_p_nolabel(Array(x0),ET_err0,smc_err0,obsET,obsSMC);

end
end

return x0, ll_init

end


init_name = "init_par_LL_calib_8par.csv"

if isfile(init_name)
	ipar_LL = Array(CSV.read(init_name, DataFrame))
else


ipar_LL = zeros(100,9);

for j in 1:100
println(j)
xi, LLi = init_works();
ipar_LL[j,1:8] = Array(xi);
ipar_LL[j,9] = LLi;
end

CSV.write(init_name, DataFrame(ipar_LL));
end

a01 = ipar_LL[argmax(ipar_LL[:,9]),1:8];

c1 = runAMH(a01, 10000, 500);

dfc = DataFrame(c1);

CSV.write("post_ET_SMC_newvar_8par.csv",dfc);


