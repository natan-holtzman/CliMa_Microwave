import numpy as np
import pandas as pd

data_folder = "./"

out_folder = "test_stats/"

obs_types = ['oAll','o1AMPM','o6AMPM',"o1and6","o16offset","o1AMPM_all"];

obs_names = ["All", "1 AM/PM", "6 AM/PM","1+6 sync.","1+6 offset"]

#inputs should be:
#a day aggregator (daily mean, amplitude, hourly, etc.)
#a period chooser (all, summer, dry summer, 2007)
#and an accuracy metric

def do_all_stats(day_agg, period_select, metric, outfname):

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
            
            g0 = np.array(pd.read_csv(data_folder+obs_types[i]+"_c"+str(chainI)+"/"+fname));
            g0 = day_agg(g0)[period_select,:]
            p0 = np.array(pd.read_csv(data_folder+obs_types[i]+"_c"+str(chainI)+"/"+"post_par.csv"))[6000::100,:13]

            datalist.append(g0[:,:-1])
            normlist.append((g0[:,:-1]+grav_pot) * np.exp(p0[:,11]).reshape(1,40) - grav_pot )
        datalist.append(g0[:,-1].reshape(-1,1))
        normlist.append(g0[:,-1].reshape(-1,1))
        data_all = np.concatenate(datalist,axis=1)
        norm_all = np.concatenate(normlist,axis=1)   
        
        leafpost.append(data_all)
        normpost.append(norm_all)

    #print("Leaf daily mean")
    #leaftab = [metric(x) for x in leafpost]

    #print("Leaf normalized")
    LWPerrs = np.array([metric(x) for x in normpost])

    def do_compare(fname,daily):
        leafpost = []

        for i in range(5):
            datalist = []
            for chainI in range(1,4):
            
                g0 = np.array(pd.read_csv(data_folder+obs_types[i]+"_c"+str(chainI)+"/"+fname));
                g0 = day_agg(g0)[period_select,:]
                datalist.append(g0[:,:-1])
            datalist.append(g0[:,-1].reshape(-1,1))
            data_all = np.concatenate(datalist,axis=1)
            leafpost.append(data_all)
        ETtab = np.array([metric(x) for x in leafpost])
        return ETtab

    fname = "postET.csv"
    #print(fname)
    ETerrs = do_compare(fname,1)
    ETerrs *= 18.02/1000*60*60*24;

    fname = "postGPP.csv"
    #print(fname)
    GPPerrs = do_compare(fname,1)

    #fname = "postGSW.csv"
    #print(fname)
    #GSWerrs = do_compare(fname,1)

    fname = "postSWS.csv"
    #print(fname)
    SMCerrs = do_compare(fname,1)

    all_errs = [LWPerrs,SMCerrs,ETerrs,GPPerrs];

    np.save(out_folder+outfname,np.array(all_errs))
    print("Saving "+outfname)

df_raw = pd.read_csv("../../../data/moflux_fluxnet_data_nov2022_lef.csv");
df_raw = df_raw.iloc[:(24*365*13)]

year_is_2007 = np.array(df_raw.YEAR) == 2007

dry_year = [2005, 2012, 2013, 2014]
year_is_dry = np.array(df_raw.YEAR.isin(dry_year))
is_summer = np.array((df_raw.Day >= 152)&(df_raw.Day <= 273))
dry_summers = (is_summer & year_is_dry) == 1

allyears = np.array(df_raw.YEAR >= 0)

#Daily RMSE over all years
#Daily RMSE over dry summers
#Hourly RMSE over all years
#Hourly RMSE over all summers

def get_daily_2d(x,nstep):
    y = np.reshape(x, (-1, nstep, x.shape[-1]))
    return np.mean(y, axis=1)

def get_daily1(x):
    return get_daily_2d(x,24)
    
def get_hourly1(x):
    return x

def get_rmse_dist(tab):
    diffs = tab[:,:-1] - tab[:,-1].reshape(-1,1)
    return np.sqrt(np.nanmean(diffs**2,0))

do_all_stats(get_daily1, allyears[::24], get_rmse_dist, "all_daily_RMSE.npy")
do_all_stats(get_daily1, dry_summers[::24], get_rmse_dist, "dry_daily_RMSE.npy")
do_all_stats(get_hourly1, allyears, get_rmse_dist, "all_hourly_RMSE.npy")
do_all_stats(get_hourly1, is_summer, get_rmse_dist, "summer_hourly_RMSE.npy")

#5 AM RMSE in 2007
#Diurnal RMSE in 2007
#Daily cor in 2007
#Daily mean in 2007

def get_mean_tab(tabx):
    lwp_means = [np.mean(tabx[:,j]) for j in range(tabx.shape[1])]
    return np.array(lwp_means)

def getcor(tab):
    return np.array([np.corrcoef(tab[:,j],tab[:,-1])[0,1] for j in range(tab.shape[1]-1)])

def get_5am(x,nstep):
    y = np.reshape(x, (-1, nstep, x.shape[-1]))
    return y[:,6,:]


def get_diurnal_amp(x,nstep):
    y = np.reshape(x, (-1, nstep, x.shape[-1]))
    return y[:,6,:] - y[:,14]

do_all_stats(get_daily1, year_is_2007[::24], get_mean_tab, "y2007_daily_mean.npy")
do_all_stats(get_daily1, year_is_2007[::24], getcor, "y2007_daily_cor.npy")
do_all_stats(get_5am, year_is_2007[::24], get_rmse_dist, "y2007_5am_RMSE.npy")
do_all_stats(get_diurnal_amp, year_is_2007[::24], get_mean_tab, "y2007_diurnal_RMSE.npy")
