import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));


#prior_min = [ 0.1,  1e-6, 500, 0.75, 0.75,0.01];
#prior_max = [ 100, 2e-5, 3000, 10, 8,10];


#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));
#prior_min = [10, 0.01, 0.1,  1e-6, 500, 0.75, 0.75, 0.1,100];
#prior_max = [120,0.75,  100, 2e-5, 3000, 10, 8, 10,1000];


#true_val = [22,0.33, 2, 1e-5,600, 4, 3, 1, 600];
par_names = ["Vcmax", "Stomatal margin", "Kmax_plant", "Kmax_soil", "Soil depth", "Kplant location", "Kplant shape","Volume factor","Medlyn g1"];


#plt.rcParams["lines.linewidth"] = 1;
#plt.rcParams["font.size"] = 16;

#for each run in FixET, get the parameters
out_folder = "./";

subdir_list = ["oAll_c2/", "o1AMPM_c2/","o6AMPM_c2/"]
subdir_list2 = ["oAll_c3/", "o1AMPM_c3/", "o6AMPM_c3/"]
subdir_list3 = ["oAll_c1/", "o1AMPM_c1/", "o6AMPM_c1/"]


obs_names = ["All", "1 AM/PM", "6 AM/PM"]

d24 = np.arange(24*365*12) % (24*365)
summer_24 = (d24 >= 150)*(d24 < 275)
summer_24 = d24 >= 0

d1 = np.arange(365*12) % 365 
summer_1 = (d1 >= 150)*(d1 < 275)
summer_1 = d1 >= 0


leafpost = []
for i in range(3):
    g0 = np.array(pd.read_csv(out_folder+subdir_list3[i]+"postLeaf.csv"))[:,:-1];  
    g1 = np.array(pd.read_csv(out_folder+subdir_list[i]+"postLeaf.csv"))[:,:-1];
    g2 = np.array(pd.read_csv(out_folder+subdir_list2[i]+"postLeaf.csv")); 
    g3 = np.concatenate((g0,g1,g2),axis=1)
    
    good_run = np.mean(g3==0,0) != 1
    g3 = g3[:,good_run]
    leafpost.append(g3[:,:])

#leaf_post = np.array(leaf_post)

#out_names = ["Leaf pre-dawn","Leaf diurnal","RZSM"];

#leaftab = np.array(pd.read_csv("postLeaf.csv"));
#branchtab = np.array(pd.read_csv(out_folder+"postBranch.csv"));
#trunktab =  np.array(pd.read_csv(out_folder+"postTrunk.csv"));
#RZSMtab =  np.array(pd.read_csv(out_folder+"postRZ.csv"));
#ETtab =  np.array(pd.read_csv(out_folder+"postET.csv"));

def getcor(x,y):
    return np.corrcoef(x,y)[0,1]

def getRMSE(x,y):
    return np.sqrt(np.mean((x-y)**2))/np.mean(y)*100

def getbias(pred,obs):
    return np.mean(pred-obs)/np.mean(obs)*100

def getR2(pred,obs):
    return 1 - np.mean((pred-obs)**2)/np.var(obs)


def wrapfunc(func,tab):
    truthval = tab[:,-1]
    posterior = tab[:,:-1]
    allstats = [func(posterior[:,i],truthval) for i in range(posterior.shape[1])]
    return np.nanmedian(allstats)

def wrapR2(tab):
    truthval = tab[:,-1]
    posterior = tab[:,:-1]
    allstats = [getR2(posterior[:,i],truthval) for i in range(posterior.shape[1])]
    return np.array(allstats)

def get_diurnal_2d(x,nstep):
    y = np.reshape(x, (-1, nstep, x.shape[-1]))
    return np.mean(y, axis=0)	

def get_diurnal_1d(x,nstep):
    y = np.reshape(x, (-1, nstep))
    return np.mean(y, axis=0)
def get_daily_2d(x,nstep):
    y = np.reshape(x, (-1, nstep, x.shape[-1]))
    return np.mean(y, axis=1)
#outvar = [leaftab[1::24,:],get_diurnal(leaftab,24),
#	  RZSMtab[:,:] ]

def make_stats_tab(tab_list):
    ans = np.zeros((4,4))
    ans[0,:] = [wrapfunc(getcor,x) for x in tab_list]
    ans[1,:] = [wrapfunc(getbias,x) for x in tab_list]
    ans[2,:] = [wrapfunc(getRMSE,x) for x in tab_list]
    ans[3,:] = [wrapfunc(getR2,x) for x in tab_list]
    return ans


obs_names = ["All", "1 AM/PM", "6 AM/PM", "1 AM"]

j = 0

colors_list = ["tab:blue","tab:green","tab:orange","tab:red"]



print("Leaf")
#leaftab = [wrapR2(x) for x in leafpost]
#print(leaftab)

#fig,ax = plt.subplots(figsize=(12,8))
#plt.subplot(2,2,1

#plt.figure(figsize=(12,8))

#daily_ser = [get_daily_2d(x,24) for x in leafpost]
daily_ser = [x[4::24,:] for x in leafpost]

daily_mean = [np.median(x[:,:-1],axis=1) for x in daily_ser]
daily_25 = [np.quantile(x[:,:-1],0.25,axis=1) for x in daily_ser]
daily_75 = [np.quantile(x[:,:-1],0.75,axis=1) for x in daily_ser]

np.savetxt("lwp_means.csv", np.transpose(np.array(daily_mean)))
np.savetxt("lwp_lower.csv", np.transpose(np.array(daily_25)))
np.savetxt("lwp_upper.csv", np.transpose(np.array(daily_75)))

np.savetxt("lwp_true.csv",np.transpose(daily_ser[0][:,-1]))



#bigtab = pd.DataFrame(np.concatenate([leaftab,leaftab_predawn ,trunktab,trunktab_predawn, ETtab, leaftab4]))
#bigtab.to_csv("post_stats_nov28summer.csv")


