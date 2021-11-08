
using DataFrames
using CSV
using PyPlot
using Random
using StatsBase

prior_min = [10, 0.01, 0.1,  1e-6, 500, 0.75, 0.75];
prior_max = [120,0.75,  100, 2e-5, 3000, 10, 8];

#true_val = [22, 0.33, 15, 1e-5, 800, 1.5, 1];
par_names = ["Vcmax", "Stomatal margin", "Kmax_plant", "Kmax_soil", "Soil depth", "Kplant location", "Kplant shape"];


rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 1;
rcParams["font.size"] = 20;

out_folder = "outputs/realET/";

g1 = Array(CSV.read(string(out_folder,"post_par.csv")));

j = 1
fig, ax_all = subplots(2,4,figsize=(12,8))
for ax in vec(ax_all)[1:7]
    ax.plot(exp.(g1[j,:]),color="blue")
#ax.plot(exp.(g2[j,:]),color="orange")
#ax.plot(exp.(g3[j,:]),color="purple")

	
    ax.set_ylim(prior_min[j],prior_max[j])
    ax.set_title(par_names[j])
    ax.set_xticks([])
#    ax.hlines(true_val[j], 0, length(g1[j,:]), color="black",alpha=0.75)
    global j += 1
end
tight_layout()
ax_all[2,4].axis("off")
savefig("chain.png")


#leaf, stem, trunk
#rzsm, ET
out_names = ["Leaf pre-dawn","Leaf mid-day", "Branch pre-dawn","Branch mid-day","Trunk pre-dawn", "Trunk mid-day","RZSM","ET"];

leaftab = Array(CSV.read(string(out_folder,"postLeaf.csv")));
branchtab = Array(CSV.read(string(out_folder,"postBranch.csv")));
trunktab = Array(CSV.read(string(out_folder,"postTrunk.csv")));
RZSMtab = Array(CSV.read(string(out_folder,"postRZ.csv")));
ETtab = Array(CSV.read(string(out_folder,"postET.csv")));

outvar = [leaftab[2:24:end,:],leaftab[13:24:end,:], branchtab[2:24:end,:], branchtab[13:24:end,:],
	  trunktab[2:24:end,:], trunktab[13:24:end,:], RZSMtab[:,:], ETtab[:,:]]


j = 1
fig, ax_all = subplots(2,4,figsize=(12,8))
for ax in vec(ax_all)[1:8]
    ax.plot(outvar[j][:,300:400],color="grey",alpha=0.1)
    ax.set_title(out_names[j])
    #ax.set_xticks([])
    global j += 1
end
tight_layout()
ax.plot(ETtab[:,end],color="red")

savefig("post_vars.png")

