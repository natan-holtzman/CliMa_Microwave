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

def getcor(x,y):
    return np.corrcoef(x,y)[0,1]

def getRMSE(x,y):
    return np.sqrt(np.mean((x-y)**2))/np.mean(y)*100

def getbias(pred,obs):
    return np.mean(pred-obs)/np.mean(obs)*100

def wrapfunc(func,tab):
    truthval = tab[:,-1]
    posterior = tab[:,:-1]
    allstats = [func(posterior[:,i],truthval) for i in range(posterior.shape[1])]
    return np.nanmedian(allstats)


def get_diurnal_2d(x,nstep):
    y = np.reshape(x, (-1, nstep, x.shape[-1]))
    return np.mean(y, axis=0)	

def get_diurnal_1d(x,nstep):
    y = np.reshape(x, (-1, nstep))
    return np.mean(y, axis=0)

#outvar = [leaftab[1::24,:],get_diurnal(leaftab,24),
#	  RZSMtab[:,:] ]

def make_stats_tab(tab_list):
    ans = np.zeros((3,4))
    ans[0,:] = [wrapfunc(getcor,x) for x in tab_list]
    ans[1,:] = [wrapfunc(getbias,x) for x in tab_list]
    ans[2,:] = [wrapfunc(getRMSE,x) for x in tab_list]
    return ans


print("Leaf")
leaftab = make_stats_tab(leafpost)
print(leaftab)

def make_stats_tab_1AM(tab_list):
    ans = np.zeros((3,4))
    ans[0,:] = [wrapfunc(getcor,x[1::24,:]) for x in tab_list]
    ans[1,:] = [wrapfunc(getbias,x[1::24,:]) for x in tab_list]
    ans[2,:] = [wrapfunc(getRMSE,x[1::24,:]) for x in tab_list]
    return ans

leaftab_predawn =make_stats_tab_1AM(leafpost)
print(leaftab_predawn)

leafpost = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"postTrunk.csv"));
    leafpost.append(g1)

print("Trunk")
trunktab = make_stats_tab(leafpost)
print(trunktab)

trunktab_predawn = make_stats_tab_1AM(leafpost)


leafpost = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"postET.csv"));
    leafpost.append(g1)

print("ET")
ETtab = make_stats_tab(leafpost)
print(ETtab)



rzpost = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"postRZ.csv"));
    rzpost.append(g1)

print("Root zone soil")
leaftab4 = make_stats_tab(rzpost)
print(leaftab4)

bigtab = pd.DataFrame(np.concatenate([leaftab,leaftab_predawn ,trunktab,trunktab_predawn, ETtab, leaftab4]))
bigtab.to_csv("post_stats_rel2.csv")


