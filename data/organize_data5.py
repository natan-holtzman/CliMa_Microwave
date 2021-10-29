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
alldata = pd.read_csv('AMF_US-MOz_BASE-BADM_8-5/AMF_US-MOz_BASE_HH_8-5.csv',comment='#')
rawdata = alldata[20016:].copy()
#%%
timestamps = [datetime.datetime.strptime(str(x), '%Y%m%d%H%M') for x in rawdata['TIMESTAMP_START']]
hours = [(x - datetime.datetime(2000,1,1,0,0)).total_seconds()/(60*60.0) for x in timestamps]
hours_ind = np.array(hours) - hours[0]
#%%
rawdata['DateTime'] = timestamps
#%%

def meanfilt(amvod,N=15):
    vodmed = np.zeros(len(amvod))
    for i in np.arange(N, len(vodmed)-N):
        vodmed[i] = np.nanmean(amvod[i-N:i+N])
    
    vodmed[:N] = np.nanmean(amvod[:N])
    vodmed[-N:] = np.nanmean(amvod[-N:])
    return vodmed
#%%

laidata = pd.read_csv('modis_lai_point_mo.csv')
laidata['DateTime'] = [datetime.datetime.strptime(x[:-2], '%Y_%m_%d') for x in laidata['system:index']]
#%%
#plt.plot(laidata['DateTime'].loc[laidata['FparLai_QC'] < 50],laidata['Lai'].loc[laidata['FparLai_QC'] < 50],'o')
#plt.plot(laidata['DateTime'].loc[laidata['FparLai_QC'] >= 50],laidata['Lai'].loc[laidata['FparLai_QC'] >= 50],'o')
#%%
laidata['Lai'].loc[laidata['FparLai_QC'] >= 50] = np.nan


lai_avhrr = pd.read_csv('moflux_lai_avhrr.csv')

lai_avhrr['DateTime'] = [datetime.datetime.strptime(x.split('_')[0], '%Y%m%d') for x in lai_avhrr['system:index']]

lai2 = pd.merge(lai_avhrr, laidata,how='left',on='DateTime')

lai2['LAI_filt_A'] = meanfilt(np.array(lai2['LAI']),15)/1000
lai2['LAI_filt_M'] = meanfilt(np.array(lai2['Lai']),15)/10

#%%
lai_site = pd.read_csv('moflux_lai_site.csv')
#%%
lai_site_time = [datetime.datetime.strptime(
    str(lai_site['Year'].iloc[i])+'_'+str(lai_site['DOY'].iloc[i]), '%Y_%j')
    for i in range(len(lai_site))]
lai_site['DateTime'] = lai_site_time
#%%
lai_mat = np.array(lai_site.loc[:,'LAI_SE':'LAI_NW'])
lai_mean = np.nanmean(lai_mat,1)
lai_site['LAI_mean_site'] = lai_mean
#%%
lai2 = pd.merge(lai2, lai_site[['DateTime','LAI_mean_site']],how='left',on='DateTime')
#%%
plt.plot(lai2['DateTime'],lai2['LAI_filt_A'])
plt.plot(lai2['DateTime'],lai2['LAI_filt_M'])
plt.plot(lai2['DateTime'],lai2['LAI_mean_site'],'o')

plt.xlim(datetime.datetime(2007,1,1),datetime.datetime(2010,1,1))
#%%
def window_interp(x,y):
    z = x*0
    pts_to_interp = np.where(np.isfinite(y))[0]
    for i in np.arange(0,len(pts_to_interp)-5):
        x0 = x[pts_to_interp[i]:pts_to_interp[i+5]]
        y0 = y[pts_to_interp[i]:pts_to_interp[i+5]]
        lmi = sm.OLS(y0,sm.add_constant(x0),missing='drop').fit()
        lmpred = lmi.predict(sm.add_constant(x0))
        z[pts_to_interp[i]:pts_to_interp[i+5]] += lmpred/5
    x0 = x[0:pts_to_interp[5]]
    y0 = y[0:pts_to_interp[5]]
    lmi = sm.OLS(y0,sm.add_constant(x0),missing='drop').fit()
    lmpred = lmi.predict(sm.add_constant(x0))
    z[0:pts_to_interp[5]] = lmpred
    return z
#%%
def season_interp(x,y):
    z = x*1 + np.nan
    hyears = np.arange(0,len(z),182)
    for i in range(len(hyears)-1):
        x0 = x[hyears[i]:hyears[i+1]]
        y0 = y[hyears[i]:hyears[i+1]]
        if np.sum(np.isfinite(y0)) >= 3:
            lmi = sm.OLS(y0,sm.add_constant(x0),missing='drop').fit()
            lmpred = lmi.predict(sm.add_constant(x0))
            z[hyears[i]:hyears[i+1]] = lmpred
    return z
#%%
mod2site = window_interp(np.array(lai2['LAI_filt_M']),np.array(lai2['LAI_mean_site']) )
am2site = window_interp(np.array(lai2['LAI_filt_A']),np.array(lai2['LAI_mean_site']) )
#%%
plt.plot(lai2['DateTime'],mod2site,label='Fitted MODIS')
plt.plot(lai2['DateTime'],am2site,label='Fitted AVHRR')
plt.plot(lai2['DateTime'],lai2['LAI_mean_site'],'o',label='Site')

plt.xlim(datetime.datetime(2007,1,1),datetime.datetime(2010,1,1))
#%%
days = np.array([x.day_of_year for x in rawdata['DateTime']])
hours = np.array([x.hour + x.minute/60 for x in rawdata['DateTime']])
#hours = np.array([x.hour for x in rawdata['DateTime']])
hour_of_year = days*24+hours
#%%
hours_only = np.array([x.hour for x in rawdata['DateTime']])
minutes = np.array([x.minute for x in rawdata['DateTime']])
#%%
lai2['LAI_fitted'] = mod2site
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
rawdata = pd.merge(rawdata, lai2[['DateTime','LAI_fitted','LAI_filt_M']], on='DateTime',how='left')
#%%

#%%
#%%
lai = fillna(rawdata['LAI_fitted']) 
laiM = fillna(rawdata['LAI_filt_M'])

tair = fillbad(np.array(rawdata['TA_1_1_1']))
tair[tair < -20] = -20
tair[tair > 40]  = 40

rain = fillbad(np.array(rawdata['P_1_1_1']))
solar_umol_in = fillbad(np.array(rawdata['PPFD_IN_1_1_1'])) 

solar_umol_out =  fillbad(np.array(rawdata['PPFD_OUT_1_1_1']))

solar_umol_net = solar_umol_in - solar_umol_out

#plt.plot(solar_umol_in[24:360*48:48])
#plt.plot(solar_umol_net[24:360*48:48])


#par_wm2 = solar_umol * 1e-6 * 6.022e23 * 4e-19

netrad = fillbad(np.array(rawdata['NETRAD_1_1_1']))

lw_out = fillbad(np.array(rawdata['LW_OUT_1_1_1']))


relhum = fillbad(np.array(rawdata['RH_1_1_1']))
wind = fillbad(np.array(rawdata['WS_1_1_1']))
le = fillbad(np.array(rawdata['LE_1_1_1']))

smc = fillbad(np.array(rawdata['SWC_1_1_1']))

ts = fillbad(np.array(rawdata['TS_1_1_1']))


p_atm = fillbad(np.array(rawdata['PA_1_1_1']))
#%%
sw_in = fillbad(np.array(rawdata['SW_IN_1_1_1']))
sw_out = fillbad(np.array(rawdata['SW_OUT_1_1_1']))

ustar = fillbad(np.array(rawdata['USTAR_1_1_1']))


h2o = fillbad(np.array(rawdata['H2O_1_1_1']))
co2 = fillbad(np.array(rawdata['CO2_1_1_1']))

patm = fillbad(np.array(rawdata['PA_1_1_1']))

# lw = fillbad(np.array(rawdata['LW_IN_1_1_1']))

# sw[par_wm2==0] = 0
#%%
nee = fillbad(np.array(rawdata['NEE_PI_F']))
#%%
SatVP = 6.1094*np.exp(17.625*tair/ (tair+ 243.04))/10  #kpa
VPD = SatVP * (1 - relhum/100)

#%%
obs = pd.read_csv('predawn_potential_2018.csv',skiprows=[1])
obs_time = [datetime.datetime.strptime(
    str(obs['Year'][i])+'_'+str(obs[' DOY'][i]), '%Y_%j')
    for i in range(len(obs))]
obs['DateTime'] = [x+datetime.timedelta(hours=6) for x in obs_time]
common_names = pd.unique(obs[' Species_Common_Name'])
obsWO = obs.loc[obs[' Species_Common_Name']==common_names[5]]
obsWOm = obsWO.groupby('DateTime').mean()
y2data = pd.merge(rawdata,obsWOm,on='DateTime',how='left')
#daypot = np.array(y2data[' PLWP'])[12::48]
#%%
oakpot = np.array(y2data[' PLWP'])


# plt.plot(par_wm2[:500])
# plt.plot(netrad[:500])
#%%
mydf = pd.DataFrame({'DateTime':np.array(rawdata['DateTime']),
                     'Day':days, 'Hour':hours_only, 'Minu':minutes,
                      'PPFD_in':solar_umol_in, 'PPFD_out':solar_umol_out,
                      'SW_in': sw_in, 'SW_out':sw_out,'Rnet':netrad,
                     'T_SOIL':ts,'T_AIR':tair, 'P_ATM':p_atm,
                     'VPD':VPD*10, 'WIND':wind, 'LW_OUT':lw_out,'RAIN':rain,
                     'RelHum':relhum/100,'USTAR': ustar,
                     'H2O':h2o, 'CO2':co2,'P_ATM':patm,
                     'LE':le, 'SMC':smc/100, 'NEE':nee,'LAI_site':lai,
                     'LAI_modis':laiM, 'Oak_Psi':oakpot
                     }).iloc[100000:]
#%%
# fapar_noon = 1- (solar_umol_out / solar_umol_in)[24::48]
# lai_noon = lai[24::48]
# plt.plot(lai_noon[lai_noon > 1], fapar_noon[lai_noon > 1],'.')
# #%%
# n_per_year = 365*24*2
# df2 = mydf.iloc[n_per_year*3:n_per_year*4]
#%%
mydf.to_csv('moflux_land_data_newnames_7_pt2.csv')
#%%
# rain_only = rain*(tair>0)
# snow = rain*(tair <= 0)
# #%%
# slist = []
# s = 0.30
# for i in range(48*100):
#     slist.append(s)
#     s = s - 0.002*(s/0.5)**4
# #%%
# b = 5
# c = 3
# prange = np.arange(-10,-0.01,0.01)
# kval = np.exp(-(-prange/b)**c)
# #%%
# plt.plot(prange,kval)
# #%%
# h_planck = 6.63e-34
# c_light = 3e8
# n_avo = 6.022e23

# photo_data = pd.read_csv('par_basis.csv')
# #%%
# single_photon_energy_W = h_planck * c_light / (photo_data['WL_nm'] * 1e-9)
# band_energy_W = photo_data['Energy_mW_m2_nm']/1000 * photo_data['dWL_nm']
# band_n_photon = band_energy_W / single_photon_energy_W
# band_mol_photon = band_n_photon/n_avo
# total_umol_photon = np.sum(band_mol_photon)*1e6
#%%
x = np.arange(10)
y = np.arange(10)
plt.plot(x,y,'o')
plt.xlabel('$\pi_{tlp}$')
