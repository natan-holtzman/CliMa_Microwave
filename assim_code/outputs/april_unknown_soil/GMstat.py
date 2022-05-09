import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));
#plt.rcParams["lines.linewidth"] = 1;
#plt.rcParams["font.size"] = 22;

print("ALL YEARS")
#prior_min = [ 0.1,  1e-6, 500, 0.75, 0.75,0.01];
#prior_max = [ 100, 2e-5, 3000, 10, 8,10];


#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));
#prior_min = [10, 0.01, 0.1,  1e-6, 500, 0.75, 0.75, 0.1,100];
#prior_max = [120,0.75,  100, 2e-5, 3000, 10, 8, 10,1000];
out_folder = "./";


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
#summer_24 = (d24 >= 150)*(d24 < 275)
summer_24 = d24 >= 0

d1 = np.arange(365*12) % 365 
#summer_1 = (d1 >= 150)*(d1 < 275)
summer_1 = d1 >= 0

y1 = np.floor(np.arange(365*12) / 365)
#summer_1 = ((y1 < 6)+(y1 > 9)) * summer_1

y24 = np.floor(np.arange(24*365*12)/(24*365))
#summer_24 = ((y24<6)+(y24>9)) * summer_24


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
    return np.sqrt(np.mean((x-y)**2))

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

def meanR2(tab):
    truthval = tab[:,-1]
    postmean = np.nanmean(tab[:,:-1],1)
    return getRMSE(postmean, truthval)

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


def wrapRMSE(tab):
    truthval = tab[:,-1]
    posterior = tab[:,:-1]
    allstats = [getRMSE(posterior[:,i],truthval) for i in range(posterior.shape[1])]
    return np.array(allstats)


obs_names = ["All", "1 AM/PM", "6 AM/PM", "1 AM"]

j = 0

colors_list = ["tab:blue","tab:green","tab:orange","tab:red"]

leafpars = []


for i in range(3):
    p0 = np.array(pd.read_csv(out_folder+subdir_list[i]+"post_par.csv"))[-4000:,:12]
    p1 = np.array(pd.read_csv(out_folder+subdir_list2[i]+"post_par.csv"))[-4000:,:12]
    p2 = np.array(pd.read_csv(out_folder+subdir_list3[i]+"post_par.csv"))[-4000:,:12]
    p3 = [p0,p1,p2]
    leafpars.append(p3)

#leafpars is a list of 3 items
#each item corresponds to an observation scenerio
#within each item, there is a list of 3 chains
#within each chain, there is a Lx7 array, with dims samples x parameters 

def getGM(x):
    z = np.array(x)
    J,L,nP = z.shape
    chain_mean = np.mean(z,axis=1) #dims J x nP
    grand_mean = np.mean(z,axis=(0,1))
    B = L*np.var(chain_mean,axis=0,ddof=1)
    s_squared = np.var(z,axis=1,ddof=1)
    W = np.mean(s_squared,axis=0)
    R = ((L-1)/L*W + B/L) / W
    return R
print("Gelman Rubin statistic")
for i in range(3):
    print(getGM(leafpars[i]))


