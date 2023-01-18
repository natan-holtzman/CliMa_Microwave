import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));
plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["font.size"] = 12;
plt.rcParams["mathtext.default"] = "regular"

df_raw = pd.read_csv("../../../data/moflux_fluxnet_data_nov2022_lef.csv");
df_raw = df_raw.iloc[:(24*365*13)]
#dry_year = [2005, 2012, 2013, 2014]
summer_24 = np.ones(len(df_raw)) == 1
summer_1 = summer_24[::24]

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

#subdir_list = ["oAll_c2/", "o1AMPM_c2/","o6AMPM_c2/", "o1and6_c2/", "o16offset_c2/"]
#subdir_list2 = ["oAll_c3/", "o1AMPM_c3/", "o6AMPM_c3/", "o1and6_c3/", "o16offset_c3/"]
#subdir_list3 = ["oAll_c1/", "o1AMPM_c1/", "o6AMPM_c1/", "o1and6_c1/", "o16offset_c1/"]

obs_types = ['oAll','o1AMPM','o6AMPM',"o1and6","o16offset","o1AMPM_all"];



obs_names = ["All", "1 AM/PM", "6 AM/PM","1+6 sync.","1+6 offset"]

#d24 = np.arange(24*365*12) % (24*365)
#summer_24 = (d24 >= 150)*(d24 < 275)
#summer_24 = d24 >= 0

#d1 = np.arange(365*12) % 365 
#summer_1 = (d1 >= 150)*(d1 < 275)
#summer_1 = d1 >= 0

#y1 = np.floor(np.arange(365*12) / 365) + 2005
#summer_1 = ((y1 < 6)+(y1 > 9)) * summer_1
#summer_1 = (y1 == 6)*summer_1

#y24 = np.floor(np.arange(24*365*12)/(24*365)) + 2005

#summer_1 = (((y1 == 2005) + (y1 == 2012) + (y1 == 6)) > 0) * summer_1
#summer_24 = (((y24 == 2005) + (y24 == 2012) + (y24 == 6)) > 0) * summer_24


#summer_1 = (summer_1==0)
#summer_24 = (summer_24==0)


#summer_24 = ((y24<6)+(y24>9)) * summer_24
#summer_24 = (y24 == 6)*summer_24

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

def meanR2_med(tab):
    truthval = tab[:,-1]
    postmean = np.nanmedian(tab[:,:-1],1)
    return getRMSE(postmean, truthval)

def meanR2_dist(tab):
    diffs = tab[:,:-1] - tab[:,-1].reshape(-1,1)
    return np.sqrt(np.nanmean(diffs**2,0))

def meanR2_dist_unbiased(tab):
    std_post = tab[:,:-1]*1
    std_post -= np.nanmean(std_post,0).reshape(1,-1)
    std_post /= np.nanstd(std_post,0).reshape(1,-1)
    std_post *= np.nanstd(tab[:,-1])
    std_post += np.nanmean(tab[:,-1]) 
    diffs = np.reshape(std_post[:,:-1] - tab[:,-1].reshape(-1,1), (-1,1))
    
    return np.sqrt(np.nanmean(diffs**2))


def meanR2(tab):
    return meanR2_dist(tab)


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


def get_anom(x):
    xclim = np.mean(x.reshape(-1,365),0)
    xstd = np.std(x.reshape(-1,365),0)
    zanom = (x - np.tile(xclim,13))/np.tile(xstd,13)
    return zanom

def do_analysis(tabx):
    lwp_means = [np.mean(tabx[:,j]) for j in range(tabx.shape[1])]
    lwp_relative = [tabx[:,j]/np.mean(tabx[:,j]) for j in range(tabx.shape[1])]
    relative_clim = [np.mean(x.reshape(-1,365),0) for x in lwp_relative]
    relative_anom = [get_anom(x) for x in lwp_relative]
    return [np.array(lwp_means),meanR2_dist(np.array(relative_clim).T),
            meanR2_dist(np.array(relative_anom).T)]


#def do_analysis(listx):
#    to_unwrap = [do_analysis_base(x) for x in listx]
#    return [np.stack([x[i] for x in to_unwrap],1) for i in range(3)]

j = 0

mpa2mm = 10**6 / 9.8
grav_pot = 13.5*1000/mpa2mm;

allpars = []

fname = "postLeaf.csv"

for i in range(5):
    #datalist = []
    #normlist = []
    for chainI in range(1,4):
        if i == 1 and chainI == 99:
            pass
        else:
            #g0 = np.array(pd.read_csv(out_folder+obs_types[i]+"_c"+str(chainI)+"/"+fname));
            #g0 = get_daily_2d(g0,24)[summer_1,:]
            p0 = np.array(pd.read_csv(out_folder+obs_types[i]+"_c"+str(chainI)+"/"+"post_par.csv"))[6000::100,:13]
 
     #   datalist.append(g0[:,:-1])
        allpars.append(p0)
   # datalist.append(g0[:,-1].reshape(-1,1))
   # normlist.append(g0[:,-1].reshape(-1,1))
   # data_all = np.concatenate(datalist,axis=1)
   # norm_all = np.concatenate(normlist,axis=1)   
    
  #  leafpost.append(data_all)
  #  normpost.append(norm_all)

#normpost.append(g0[:,-1].reshape(-1,1))
#normpost = np.concatenate(normpost,axis=1)
#print(normpost.shape)

    #leafpost.append(get_daily_2d(g3,24)[summer_1,:])
    #leafpost_midnight.append(g3[3::24,:])
    #leafpost_noon.append(g3[15::24,:])    


#print("Leaf daily mean")
#leaftab = [meanR2(x) for x in leafpost]
#print(leaftab)
#print(np.mean(leafpost[1][:,-1]))


#print([x.shape for x in normpost])
#LWPerrs = np.array([meanR2(x) for x in normpost])


#LWPmean = np.abs(np.mean(normpost[0][:,-1]))

#print(leaftab)
#print(np.mean(normpost[1][:,-1]))
np.save("pars_nov25.npy",np.concatenate(allpars,axis=0))
