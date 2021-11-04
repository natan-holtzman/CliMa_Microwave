using StatsBase
#using Plots
using PyPlot


#navg = 2

function avg2(x,navg)
	if (typeof(x[1])==Float64) | (typeof(x[1])==Float32)
			return mean(reshape(x,(navg,:)),dims=1)[1,:];
	else
			return x[1:navg:end]
	end
end



function get_diurnal(x,navg)
	if (typeof(x[1])==Float64) | (typeof(x[1])==Float32)
			return mean(reshape(x,(navg,:)),dims=2)[:,1];
	else
			return x[1:navg:end]
	end
end



function avg2_2d(x,navg)
return return mean(reshape(x,(navg,:,size(x)[2])),dims=1)[1,:,:];
end

function avg2_df(x,navg)

ncol = Integer(size(x)[2]);
nrow = Integer(size(x)[1]/navg);
ans = deepcopy(x[1:nrow,:]);

for j in 1:ncol
#println(x[1:10,j])
#println(avg2(x[1:10,j]))
        ans[:,j] = avg2(x[:,j]);
end

return ans;

end

function convert_sim(x)
a = avg2_df(x[1],2);
b = avg2_2d(x[2],2);
c = avg2_2d(x[3],2);
d = avg2_2d(x[4],2);

return a,b,c,d
end


#=
N = 48*3
istart = 48*365*2 - 52*48 + 230*48
soil_init1 = 0.2;
=#

N = 48*365
istart = 48*365*2 - 52*48 #+ 230*48
soil_init1 = 0.38;




#include("rebuild_sim_ij.jl")
include("sim_vary_new_stomata.jl")
println("No capacitance")
sim_res1 = run_sim(FT(22),FT(1.5), FT(20), FT(1e-4), FT(2.5), FT(2.0), FT(700),FT(0.45),istart,N, soil_init1,FT(1),FT(1e-5),1);

#sim_res1_time = @timed run_sim(FT(22),FT(1.5), FT(20), FT(1e-4), FT(2.5), FT(2.0), FT(700),FT(0.45),istart,N, soil_init1);


#include("rebuild_sim_ij_nss.jl")
println("Old capacitance")

sim_res2 = run_sim(FT(22),FT(1.5), FT(20), FT(1e-4), FT(2.5), FT(2.0), FT(700),FT(0.45),istart,N, soil_init1,FT(1),FT(1e-5),2);
#sim_res2_time = @timed run_sim(FT(22),FT(1.5), FT(20), FT(1e-4), FT(2.5), FT(2.0), FT(700),FT(0.45),istart,N, soil_init1,1, 1e-5);


#include("rebuild_sim_ij_par_vol.jl")
println("New capacitance")

sim_res3 = run_sim(FT(22),FT(1.5), FT(20), FT(1e-4), FT(2.5), FT(2.0), FT(700),FT(0.45),istart,N, soil_init1,FT(0.5),FT(1e-5),3);

sim_res4_time = @timed run_sim(FT(22),FT(1.5), FT(20), FT(1e-4), FT(2.5), FT(2.0), FT(700),FT(0.45),istart,N, soil_init1,FT(2),FT(1e-5),3);

sim_res4 = sim_res4_time.value;

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 1;
rcParams["font.size"] = 24;


xarr = collect(1:N)/48

figure(figsize=(12,6))
plot(xarr, sim_res1[1].leafpot, label="No storage")
plot(xarr, sim_res2[1].leafpot, label="Old capacitance")
plot(xarr, sim_res3[1].leafpot, label="New capacitance, small storage")
plot(xarr, sim_res4[1].leafpot, label="New capacitance, large storage")
xlabel("Time (days)")
ylabel("Leaf water potential (MPa)")
legend()
show()

#=
figure(figsize=(12,6))
plot(xarr, sim_res1[1].ETmod,"--" , label="No storage")
plot(xarr, sim_res2[1].ETmod,"--", label="Old capacitance")
plot(xarr, sim_res3[1].ETmod,"--", label="New capacitance, small storage")
plot(xarr, sim_res4[1].ETmod,"--", label="New capacitance, large storage")
xlabel("Time (hours)")
ylabel("ET (mol/m2/s)")
legend()
show()
=#

#=
dayrange = collect(1:Int(N/48))

figure(figsize=(12,6))
plot(dayrange, avg2(sim_res1[1].ETmod,48), label="No storage")
plot(dayrange, avg2(sim_res2[1].ETmod,48), label="Old capacitance")
plot(dayrange, avg2(sim_res3[1].ETmod,48), label="New capacitance, small storage")
plot(dayrange, avg2(sim_res4[1].ETmod,48), label="New capacitance, large storage")
xlabel("Time (days)")
ylabel("ET (mol/m2/s)")
legend()
show()

diurn_range = collect(1:48)/2

figure(figsize=(12,6))
plot(diurn_range, get_diurnal(sim_res1[1].ETmod,48), label="No storage")
plot(diurn_range, get_diurnal(sim_res2[1].ETmod,48), ":" ,label="Old capacitance")
plot(diurn_range, get_diurnal(sim_res3[1].ETmod,48), label="New capacitance, small storage")
plot(diurn_range, get_diurnal(sim_res4[1].ETmod,48), ":", label="New capacitance, large storage")
xlabel("Time (hours)")
ylabel("ET (mol/m2/s)")
legend()
show()



figure(figsize=(12,6))
plot(diurn_range, get_diurnal(sim_res1[1].leafpot,48), label="No storage")
plot(diurn_range, get_diurnal(sim_res2[1].leafpot,48),label="Old capacitance")
plot(diurn_range, get_diurnal(sim_res3[1].leafpot,48), label="New capacitance, small storage")
plot(diurn_range, get_diurnal(sim_res4[1].leafpot,48), label="New capacitance, large storage")
xlabel("Time (hours)")
ylabel("Leaf water potential (mol/m2/s)")
legend()
show()
=#