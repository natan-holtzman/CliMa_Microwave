import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


prior_min = [10, 0.01, 0.1,  1e-6, 500, 0.75, 0.75,0.01];
prior_max = [120,0.75,  100, 2e-5, 3000, 10, 8,10];

#true_val = [22, 0.33, 15, 1e-5, 800, 1.5, 1];
par_names = ["Vcmax", "Stomatal margin", "Kmax_plant", "Kmax_soil", "Soil depth", "Kplant location", "Kplant shape","Volume factor"];


plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["font.size"] = 16;

out_folder = "outputs/realET_cap/";

g1 = np.array(pd.read_csv(out_folder+"post_par.csv"));

j = 0
fig, ax_all = plt.subplots(2,4,figsize=(12,8))
for ax in ax_all.ravel()[:]:
    ax.plot(np.exp(g1[:,j]),color="blue")
#ax.plot(exp.(g2[j,:]),color="orange")
#ax.plot(exp.(g3[j,:]),color="purple")

	
    ax.set_ylim((prior_min[j]),(prior_max[j]))
    ax.set_title(par_names[j])
    ax.set_xticks([])
#    ax.hlines(true_val[j], 0, length(g1[j,:]), color="black",alpha=0.75)
    j += 1
#end
plt.tight_layout()
#ax = ax_all[1,3]
#ax.plot(np.sqrt(g1[:,-1]))
#ax.set_title("ET RMSE")

plt.savefig("chainC.png")


#leaf, stem, trunk
#rzsm, ET
out_names = ["Leaf pre-dawn","Leaf mid-day", "Branch pre-dawn","Branch mid-day","Trunk pre-dawn", "Trunk mid-day","RZSM","ET"];

leaftab = np.array(pd.read_csv(out_folder+"postLeaf.csv"));
branchtab = np.array(pd.read_csv(out_folder+"postBranch.csv"));
trunktab =  np.array(pd.read_csv(out_folder+"postTrunk.csv"));
RZSMtab =  np.array(pd.read_csv(out_folder+"postRZ.csv"));
ETtab =  np.array(pd.read_csv(out_folder+"postET.csv"));

outvar = [leaftab[1::24,:],leaftab[12::24,:], branchtab[1::24,:], branchtab[12::24,:],
	  trunktab[1::24,:], trunktab[12::24,:], RZSMtab[:,:], ETtab[:,:]]


j = 0
fig, ax_all = plt.subplots(2,4,figsize=(12,8))
for ax in ax_all.ravel():
    ax.plot(outvar[j][:,-100:-1],color="grey",alpha=0.1)
    ax.set_title(out_names[j])
    #ax.set_xticks([])
    j += 1
#end
plt.tight_layout()
ax.plot(ETtab[:,-1],color="red")

plt.savefig("post_varsC.png")

