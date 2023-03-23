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
summer_24 = np.array(df_raw.YEAR) ==2007
summer_1 = summer_24[::24]

print("ALL YEARS")

out_folder = "./";


#true_val = [22,0.33, 2, 1e-5,600, 4, 3, 1, 600];
par_names = ["Vcmax", "Stomatal margin", "Kmax_plant", "Kmax_soil", "Soil depth", "Kplant location", "Kplant shape","Volume factor","Medlyn g1"];


out_folder = "./";

obs_types = ['oAll','o1AMPM','o6AMPM',"o1and6","o16offset","o1AMPM_all"];



obs_names = ["All", "1 AM/PM", "6 AM/PM","1+6 sync.","1+6 offset"]

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

def getmean(tab):
    return np.nanmean(tab,0)

def get_diurnal_2d(x,nstep):
    y = np.reshape(x, (-1, nstep, x.shape[-1]))
    return np.mean(y, axis=0)	

def get_diurnal_1d(x,nstep):
    y = np.reshape(x, (-1, nstep))
    return np.mean(y, axis=0)
def get_daily_2d(x,nstep):
    y = np.reshape(x, (-1, nstep, x.shape[-1]))
    return y[:,6,:] - y[:,14]
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


j = 0

mpa2mm = 10**6 / 9.8
grav_pot = 13.5*1000/mpa2mm;

leafpost = []
normpost = []

fname = "postLeaf.csv"

for i in range(5):
    datalist = []
    normlist = []
    for chainI in range(1,4):
        if i == 1 and chainI == 99:
            pass
        else:
            g0 = np.array(pd.read_csv(out_folder+obs_types[i]+"_c"+str(chainI)+"/"+fname));
            g0 = get_daily_2d(g0,24)[summer_1,:]
            p0 = np.array(pd.read_csv(out_folder+obs_types[i]+"_c"+str(chainI)+"/"+"post_par.csv"))[6000::100,:13]
 
        datalist.append(g0[:,:-1])
        normlist.append((g0[:,:-1]+grav_pot) * np.exp(p0[:,11]).reshape(1,40) - grav_pot )
    datalist.append(g0[:,-1].reshape(-1,1))
    normlist.append(g0[:,-1].reshape(-1,1))
    data_all = np.concatenate(datalist,axis=1)
    norm_all = np.concatenate(normlist,axis=1)   
    
    leafpost.append(data_all)
    normpost.append(norm_all)

    #leafpost.append(get_daily_2d(g3,24)[summer_1,:])
    #leafpost_midnight.append(g3[3::24,:])
    #leafpost_noon.append(g3[15::24,:])    


print("Leaf daily mean")
leaftab = [meanR2(x) for x in leafpost]
#print(leaftab)
#print(np.mean(leafpost[1][:,-1]))


print("Leaf normalized")
LWPerrs = np.array([meanR2(x) for x in normpost])
LWPmean = np.abs(np.mean(normpost[0][:,-1]))

#print(leaftab)
#print(np.mean(normpost[1][:,-1]))

def do_compare(fname,daily):
	leafpost = []

	for i in range(5):
		datalist = []
		for chainI in range(1,4):
			if i == 1 and chainI == 99:
				pass
			else:
				g0 = np.array(pd.read_csv(out_folder+obs_types[i]+"_c"+str(chainI)+"/"+fname));
				if daily:
					g0 = get_daily_2d(g0,24)[summer_1,:]
				else:
					g0 = g0[summer_1,:]
			datalist.append(g0[:,:-1])
		datalist.append(g0[:,-1].reshape(-1,1))
		data_all = np.concatenate(datalist,axis=1)
		leafpost.append(data_all)
		#g3[:,np.mean(g3==0,0) > 0.1] = np.nan

	ETtab = np.array([meanR2(x) for x in leafpost])
	#np.savetxt(outname,ETtab)#print(ETtab)
	#print(np.mean(leafpost[1][:,-1]))
	return ETtab,np.mean(g0[:,-1])


fname = "postET.csv"
print(fname)
ETerrs,ETmean = do_compare(fname,1)
ETerrs *= 18.02/1000*60*60*24;
ETmean *= 18.02/1000*60*60*24


fname = "postGPP.csv"
print(fname)
GPPerrs,GPPmean = do_compare(fname,1)

fname = "postGSW.csv"
print(fname)
GSWerrs,GSWmean = do_compare(fname,1)

fname = "postSWS.csv"
print(fname)
SMCerrs,SMCmean = do_compare(fname,1)


all_errs = [LWPerrs,SMCerrs,ETerrs,GPPerrs];
all_means = [LWPmean, SMCmean, ETmean, GPPmean]
print(all_means)

np.save("d52_rmse_2007.npy",np.array(all_errs))


