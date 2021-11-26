import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));


prior_min = [10, 0.01, 0.1,  1e-6, 500, 0.75, 0.75,0.01,100];
prior_max = [120,0.75,  100, 2e-5, 3000, 10, 8,10,1000];

true_val = [22, 0.33, 2, 1e-5, 600, 4, 3, 1, 600];
par_names = ["Vcmax", "Stomatal margin", "Kmax_plant", "Kmax_soil", "Soil depth", "Kplant location", "Kplant shape","Volume factor","Medlyn g1"];


plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["font.size"] = 16;

out_folder = "outputs/obs2x3_med/";

g1 = np.array(pd.read_csv(out_folder+"post_par.csv"));

j = 0
fig, ax_all = plt.subplots(3,3,figsize=(12,8))
for ax in ax_all.ravel()[:]:
    ax.plot(np.exp(g1[:,j]),color="blue")
#ax.plot(exp.(g2[j,:]),color="orange")
#ax.plot(exp.(g3[j,:]),color="purple")

	
    ax.set_ylim((prior_min[j]),(prior_max[j]))
    ax.set_title(par_names[j])
    ax.set_xticks([])
    ax.hlines(true_val[j], 0, len(g1[:,j]), color="black",alpha=0.75)
    j += 1
#end
plt.tight_layout()
#ax = ax_all[1,3]
#ax.plot(np.sqrt(g1[:,-1]))
#ax.set_title("ET RMSE")

plt.savefig("chain2x3_med.png")


#leaf, stem, trunk
#rzsm, ET
out_names = ["Leaf pre-dawn","Leaf diurnal","RZSM","ET"];

leaftab = np.array(pd.read_csv(out_folder+"postLeaf.csv"));
branchtab = np.array(pd.read_csv(out_folder+"postBranch.csv"));
trunktab =  np.array(pd.read_csv(out_folder+"postTrunk.csv"));
RZSMtab =  np.array(pd.read_csv(out_folder+"postRZ.csv"));
ETtab =  np.array(pd.read_csv(out_folder+"postET.csv"));

def get_diurnal(x,nstep):
    y = np.reshape(x, (-1, nstep, x.shape[-1]))
    return np.mean(y, axis=0)	

outvar = [leaftab[1::24,:],get_diurnal(leaftab,24),
	  RZSMtab[:,:], ETtab[:,:]]


j = 0
fig, ax_all = plt.subplots(2,2,figsize=(12,8))
for ax in ax_all.ravel():
    ax.plot(outvar[j][:,-100:-1],color="grey",alpha=0.1)
    ax.plot(outvar[j][:,-1],"r")
    ax.set_title(out_names[j])
    #ax.set_xticks([])
    j += 1
#end
plt.tight_layout()
#ax.plot(ETtab[:,-1],color="red")

plt.savefig("post_vars2x3d_med.png")
