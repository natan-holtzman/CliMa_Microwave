# -*- coding: utf-8 -*-
"""
Created on Sun May  2 14:17:24 2021

@author: natan
"""



import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import statsmodels.api as sm
#%%
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 14
fig_size[1] = 7

plt.rcParams['font.size']=28
#%%
fluxnet_folder = r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\AMF_US-MOz_FLUXNET_SUBSET_2004-2019_3-5"
rawdata = pd.read_csv(fluxnet_folder+"/AMF_US-MOz_FLUXNET_SUBSET_HH_2004-2019_3-5.csv",comment='#')
#rawdata = pd.read_csv('../../../AMF_US-MOz_BASE-BADM_8-5/AMF_US-MOz_BASE_HH_8-5.csv',comment='#')



#rawdata = alldata[20016:].copy()

timestamps = [datetime.datetime.strptime(str(x), '%Y%m%d%H%M') for x in rawdata['TIMESTAMP_START']]
hours = [(x - datetime.datetime(2000,1,1,0,0)).total_seconds()/(60*60.0) for x in timestamps]
hours_ind = np.array(hours) - hours[0]
rawdata['DateTime'] = timestamps

#%%
years = np.array([x.year for x in rawdata['DateTime']])
#%%
days = np.array([x.day_of_year for x in rawdata['DateTime']])
hours = np.array([x.hour + x.minute/60 for x in rawdata['DateTime']])
#hours = np.array([x.hour for x in rawdata['DateTime']])
hour_of_year = days*24+hours
#%%
hours_only = np.array([x.hour for x in rawdata['DateTime']])
minutes = np.array([x.minute for x in rawdata['DateTime']])
#%%
# le_arr = pd.DataFrame({'LE':rawdata['LE_1_1_1'],'hoy':hour_of_year})
# le_arr['LE'].loc[le_arr['LE'] == -9999] = np.nan 
# le_hoy = le_arr.groupby('hoy').mean().reset_index()

# le_hoy2 = pd.merge(le_arr,le_hoy,how='left',on='hoy')

#%%

#%%
def fillbad_mean(x):
    good_avg = np.mean(x[x > -9999])
    x2 = 1*x
    x2[x == -9999] = good_avg
    return x2

def fillbad(x):
    x_nanmask = 1*x
    x_nanmask[x == -9999] = np.nan
    le_arr = pd.DataFrame({'LE':x_nanmask,'hoy':hour_of_year})
    le_hoy = le_arr.groupby('hoy').mean()
    le_hoy2 = pd.merge(le_arr,le_hoy,how='left',on='hoy')
    
    gap_filler = np.array(le_hoy2['LE_y'])
    x2 = 1*x
    x2[x == -9999] = gap_filler[x == -9999]
    return x2

def fillbad_interp(x):
    xrange = np.arange(len(x))
    x2 = np.interp(xrange, xrange[x > -9999], x[x > -9999])
    return x2
#%%
def fillna(x):
    xrange = np.arange(len(x))
    x2 = np.interp(xrange, xrange[np.isfinite(x)], x[np.isfinite(x)])
    return x2
#%%

#%%
#lai2 = pd.read_csv("lai_spline_lag.csv",parse_dates = [6])
lai2 = pd.read_csv("lai_seasonal_lag.csv",parse_dates = [6])

rawdata = pd.merge(rawdata, lai2[['DateTime','LAI_fitted','LAI_mean_site']], on='DateTime',how='left')
#%%
rawdata["YEAR"] = years
#%%

rawdata_isnan = rawdata == -9999
rawdata_isnan["YEAR"] = years
count_year = rawdata_isnan.groupby("YEAR").mean()
#%%
lai = fillna(rawdata['LAI_fitted']) 

tair = fillbad(np.array(rawdata['TA_F']))
tair[tair < -20] = -20
tair[tair > 40]  = 40

rain = fillbad(np.array(rawdata['P_F']))
solar_umol_in = fillbad(np.array(rawdata['PPFD_IN'])) 

solar_umol_out =  fillbad(np.array(rawdata['PPFD_OUT']))

solar_umol_net = solar_umol_in - solar_umol_out

#plt.plot(solar_umol_in[24:360*48:48])
#plt.plot(solar_umol_net[24:360*48:48])


#par_wm2 = solar_umol * 1e-6 * 6.022e23 * 4e-19

#netrad = fillbad(np.array(rawdata['NETRAD_1_1_1']))

lw_out = fillbad(np.array(rawdata['LW_OUT']))


relhum = fillbad(np.array(rawdata['RH']))
wind = fillbad(np.array(rawdata['WS_F']))
le = fillbad(np.array(rawdata['LE_CORR']))

smc = fillbad(np.array(rawdata['SWC_F_MDS_1']))

ts = fillbad(np.array(rawdata['TS_F_MDS_1']))


p_atm = fillbad(np.array(rawdata['PA_F']))
#%%
#sw_in = fillbad(np.array(rawdata['SW_IN_1_1_1']))
#sw_out = fillbad(np.array(rawdata['SW_OUT_1_1_1']))

#ustar = fillbad(np.array(rawdata['USTAR_1_1_1']))


#h2o = fillbad(np.array(rawdata['H2O_1_1_1']))
#co2 = fillbad(np.array(rawdata['CO2_1_1_1']))


# lw = fillbad(np.array(rawdata['LW_IN_1_1_1']))

# sw[par_wm2==0] = 0
#%%
nee = fillbad(np.array(rawdata['NEE_VUT_REF']))
gpp_n = fillbad(np.array(rawdata['GPP_NT_VUT_REF']))
gpp_d = fillbad(np.array(rawdata['GPP_DT_VUT_REF']))

#%%
SatVP = 6.1094*np.exp(17.625*tair/ (tair+ 243.04))/10  #kpa
VPD = SatVP * (1 - relhum/100)

#%%
obs = pd.read_csv('../../../predawn_potential_2020.csv',skiprows=[1])
obs_time = [datetime.datetime.strptime(
    str(obs['Year'][i])+'_'+str(obs[' DOY'][i]), '%Y_%j')
    for i in range(len(obs))]
obs['DateTime'] = [x+datetime.timedelta(hours=6) for x in obs_time]
common_names = pd.unique(obs[' Species_Common_Name'])
obsWO = obs.loc[obs[' Species_Common_Name']=='white oak']
obsWOm = obsWO.groupby('DateTime').mean()
y2data = pd.merge(rawdata,obsWOm,on='DateTime',how='left')
#daypot = np.array(y2data[' PLWP'])[12::48]
#%%
oakpot = np.array(y2data[' PLWP'])
#%%
mid_day = pd.read_csv("../../../MDLWP_MASTER.csv",skiprows=[1])
obs_time_MD = [datetime.datetime.strptime(
    str(mid_day['Year'][i])+'_'+str(mid_day['DOY'][i]), '%Y_%j')
    for i in range(len(mid_day))]


mid_day['DateTime'] = [x+datetime.timedelta(hours=6) for x in obs_time_MD]
obsWO_MD = mid_day.loc[mid_day['Species_Symbol']=='QUAL']
obsWO_MDmean = obsWO_MD.groupby('DateTime').mean()
#%%
y2data = pd.merge(y2data,obsWO_MDmean,on='DateTime',how='left')
#%%
oakpot_MD = np.array(y2data['MLWP'])
#%%
#plt.plot(oakpot,"o")
#plt.plot(oakpot_MD,"o")
# plt.plot(par_wm2[:500])
# plt.plot(netrad[:500])
#%%
# mydf = pd.DataFrame({'DateTime':np.array(rawdata['DateTime']),
#                      'Day':days, 'Hour':hours_only, 'Minu':minutes,
#                       'PPFD_in':solar_umol_in, 'PPFD_out':solar_umol_out,
#                       'SW_in': sw_in, 'SW_out':sw_out,'Rnet':netrad,
#                      'T_SOIL':ts,'T_AIR':tair, 'P_ATM':p_atm,
#                      'VPD':VPD*10, 'WIND':wind, 'LW_OUT':lw_out,'RAIN':rain,
#                      'RelHum':relhum/100,'USTAR': ustar,
#                      'H2O':h2o, 'CO2':co2,'P_ATM':patm,
#                      'LE':le, 'SMC':smc/100, 'NEE':nee,'LAI_site':lai,
#                      'LAI_modis':laiM, 'Oak_Psi':oakpot
#                      }).iloc[100000:]
#%%
mydf_all = pd.DataFrame({'DateTime':np.array(rawdata['DateTime']),
                      'Day':days, 'Hour':hours_only, 'Minu':minutes,"YEAR":years,
                      'PPFD_in':solar_umol_in, 
                      'T_SOIL':ts,'T_AIR':tair, 'P_ATM':p_atm,
                      'WIND':wind, 'LW_OUT':lw_out,'RAIN':rain,
                      'RelHum':relhum/100,
                      'LE':le, 'SMC':smc/100, 'NEE':nee,
                      'GPP_night':gpp_n, "GPP_day":gpp_d,
                      'LAI_modis':lai, 'LWP_predawn':oakpot,
                      'LWP_midday':oakpot_MD,
                      })
#%%
df_skipyear = mydf_all.loc[(mydf_all["YEAR"] != 2004) & (mydf_all["YEAR"] != 2011)].reset_index()
#%%
def df_avg(df,N):
    newlen = int(len(df)/N);
    newdf = df.iloc[:newlen,:].copy()
    coltypes = df.dtypes
    for j in range(len(df.columns)):
        if coltypes[j] == 'float64':
            newcol = np.nanmean(np.reshape(np.array(df.iloc[:,j]), (-1,N)), axis=1)
        else:
            newcol = np.array(df.iloc[:,j])[::N]
        newdf.iloc[:,j] = newcol
    return newdf
#%%
df_24 = df_avg(df_skipyear,2)
#%%
df_24.to_csv('moflux_fluxnet_data.csv')
#%%
