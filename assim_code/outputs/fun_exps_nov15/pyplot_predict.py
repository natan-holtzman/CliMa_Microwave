import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));


#prior_min = [ 0.1,  1e-6, 500, 0.75, 0.75,0.01];
#prior_max = [ 100, 2e-5, 3000, 10, 8,10];


#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));

prior_min = [10, 0.01, 0.1,  1e-7, 500, 0.75, 0.75, 0.1,100];
prior_max = [120,0.9,  50, 2e-5, 3000, 10, 8, 10,1000];

true_val = [31,0.5, 2, 0.4e-5,800, 5, 2, 1, 506];
par_names = ["Vcmax", "Stomatal margin", "Kmax_plant", "Kmax_soil", "Soil depth", "Kplant location", "Kplant shape","Volume factor","Medlyn g1"];


plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["font.size"] = 16;

#for each run in FixET, get the parameters
out_folder = "./";

subdir_list = ["oAll/", "o1AMPM/", "o6AMPM/", "o1AM/"]


obs_names = ["All", "1 AM/PM", "6 AM/PM", "1 AM"]


leafpost = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"postLeaf.csv"));
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
    ax[0,i].plot(leafpost[i][1::24,:-1],color="grey",alpha=0.5)
    ax[0,i].plot(leafpost[i][1::24,-1],"r")
    ax[1,i].plot(get_diurnal_2d(leafpost[i][:,:-1],24),color="grey",alpha=0.5)
    ax[1,i].plot(get_diurnal_1d(leafpost[i][:,-1],24),"r")
    
    ax[0,i].set_title(obs_names[i])
    ax[1,i].set_xlabel("Time of day")
    ax[0,i].set_xlabel("Day of year")
    #ax.set_xticks([])
ax[0,0].set_ylabel("1 AM LWP (MPa)")
ax[1,0].set_ylabel("Diurnal mean\nLWP (MPa)")
plt.tight_layout()
#ax.plot(ETtab[:,-1],color="red")

plt.savefig("lwp_pred.png")


leafpost = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"postTrunk.csv"));
    leafpost.append(g1)


fig, ax = plt.subplots(2,4,figsize=(12,5))
for i in range(4):
    ax[0,i].plot(leafpost[i][1::24,:-1],color="grey",alpha=0.5)
    ax[0,i].plot(leafpost[i][1::24,-1],"r")
    ax[1,i].plot(get_diurnal_2d(leafpost[i][:,:-1],24),color="grey",alpha=0.5)
    ax[1,i].plot(get_diurnal_1d(leafpost[i][:,-1],24),"r")

    ax[0,i].set_title(obs_names[i])
    ax[1,i].set_xlabel("Time of day")
    ax[0,i].set_xlabel("Day of year")
    #ax.set_xticks([])
ax[0,0].set_ylabel("1 AM LWP (MPa)")
ax[1,0].set_ylabel("Diurnal mean\nLWP (MPa)")
plt.tight_layout()
#ax.plot(ETtab[:,-1],color="red")

plt.savefig("trunk_pred.png")





leafpost = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"postET.csv"));
    leafpost.append(g1)

rzpost = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"postRZ.csv"));
    rzpost.append(g1)

fig, ax = plt.subplots(2,4,figsize=(12,5))
for i in range(4):
    ax[0,i].plot(leafpost[i][:,:-1],color="grey",alpha=0.5)
    ax[0,i].plot(leafpost[i][:,-1],"r")
    #ax[1,i].plot(get_diurnal_2d(leafpost[i][:,-100:-1],24),color="grey",alpha=0.1)
    #ax[1,i].plot(get_diurnal_1d(leafpost[i][:,-1],24),"r")
    ax[1,i].plot(rzpost[i][:,:-1],color="grey",alpha=0.5)
    ax[1,i].plot(rzpost[i][:,-1],"r")
    ax[0,i].set_title(obs_names[i])
    #ax[1,i].set_xlabel("Time of day")
    ax[1,i].set_xlabel("Day of year")
    #ax.set_xticks([])
ax[0,0].set_ylabel("Daily ET\n(mol/s/m2)")
ax[1,0].set_ylabel("Root zone\nsoil moisture")
plt.tight_layout()
#ax.plot(ETtab[:,-1],color="red")

plt.savefig("etrz_pred.png")





