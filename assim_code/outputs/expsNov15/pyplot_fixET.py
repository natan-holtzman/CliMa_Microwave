import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));


prior_min = [ 0.1,  1e-6, 500, 0.75, 0.75,0.01];
prior_max = [ 100, 2e-5, 3000, 10, 8,10];


#FT(2),FT(1e-5), FT(1000),FT(4.0),FT(2),FT(1));
true_val = [ 2, 1e-5, 1000, 4, 2, 1];
par_names = ["Kmax_plant", "Kmax_soil", "Soil depth", "Kplant location", "Kplant shape","Volume factor"];


plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["font.size"] = 16;

#for each run in FixET, get the parameters
out_folder = "fixET_air/";

subdir_list = ["obsall", "obs1", "obs6", "obs1only"]

par_post = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"/post_par.csv"));
    par_post.append(g1)

par_post = np.array(par_post) #4 x 5000 x npar

obs_names = ["All", "1 AM/PM", "6 AM/PM", "1 AM"]

j = 0
fig, ax_all = plt.subplots(2,3,figsize=(12,8))
for ax in ax_all.ravel():
    for iobs in range(4):
    	ax.plot(np.exp(par_post[iobs,:,j]),label=obs_names[iobs],alpha=0.5)
    if j==0:
        ax.legend()
	
    ax.set_ylim((prior_min[j]),(prior_max[j]))
    ax.set_title(par_names[j])
    ax.set_xticks([])
    ax.hlines(true_val[j], 0, len(g1[:,j]), color="black")
    j += 1
#end
#ax_all[1,3].axis("off")

plt.tight_layout()
#ax = ax_all[1,3]
#ax.plot(np.sqrt(g1[:,-1]))
#ax.set_title("ET RMSE")

plt.savefig("chain_fixET.png")

par_post2 = np.exp(par_post[:,2500::25,:])

j = 0
fig, ax_all = plt.subplots(2,3,figsize=(12,8))
for ax in ax_all.ravel():
    for iobs in range(4):
        ax.hist(par_post2[iobs,:,j],label=obs_names[iobs], alpha=0.5)
    if j==0:
        ax.legend()

    #ax.set_ylim((prior_min[j]),(prior_max[j]))
    ax.set_title(par_names[j])
    ax.set_yticks([])
    ax.vlines(true_val[j], 0, ax.get_ylim()[1], color="black")
    j += 1

plt.savefig("hist_fixET.png")


leafpost = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"/postLeaf.csv"));
    leafpost.append(g1)

#leaf_post = np.array(leaf_post)

#out_names = ["Leaf pre-dawn","Leaf diurnal","RZSM"];

#leaftab = np.array(pd.read_csv("postLeaf.csv"));
#branchtab = np.array(pd.read_csv(out_folder+"postBranch.csv"));
#trunktab =  np.array(pd.read_csv(out_folder+"postTrunk.csv"));
#RZSMtab =  np.array(pd.read_csv(out_folder+"postRZ.csv"));
#ETtab =  np.array(pd.read_csv(out_folder+"postET.csv"));





def get_diurnal_2d(x,nstep):
    y = np.reshape(x, (-1, nstep, x.shape[-1]))
    return np.mean(y, axis=0)	

def get_diurnal_1d(x,nstep):
    y = np.reshape(x, (-1, nstep))
    return np.mean(y, axis=0)

#outvar = [leaftab[1::24,:],get_diurnal(leaftab,24),
#	  RZSMtab[:,:] ]


fig, ax = plt.subplots(2,4,figsize=(12,5))
for i in range(4):
    ax[0,i].plot(leafpost[i][1::24,-100:-1],color="grey",alpha=0.1)
    ax[0,i].plot(leafpost[i][1::24,-1],"r")
    ax[1,i].plot(get_diurnal_2d(leafpost[i][:,-100:-1],24),color="grey",alpha=0.1)
    ax[1,i].plot(get_diurnal_1d(leafpost[i][:,-1],24),"r")
    
    ax[0,i].set_title(obs_names[i])
    ax[1,i].set_xlabel("Time of day")
    ax[0,i].set_xlabel("Day of year")
    #ax.set_xticks([])
ax[0,0].set_ylabel("1 AM LWP (MPa)")
ax[1,0].set_ylabel("Diurnal mean LWP (MPa)")
plt.tight_layout()
#ax.plot(ETtab[:,-1],color="red")

plt.savefig("lwp_post_fixET.png")

