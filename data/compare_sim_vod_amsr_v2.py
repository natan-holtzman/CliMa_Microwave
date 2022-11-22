# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 14:08:14 2021

@author: natan
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import scipy.optimize
#%%
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 14
fig_size[1] = 7

plt.rcParams['font.size']=20

import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%b %Y')
#%%
#simdata = pd.read_csv("sim_leafpot_nov30_2.csv",parse_dates=[0])
simdata = pd.read_csv("sim_leafpot_mar2022v2.csv",parse_dates=[1])

def setna(x):
    y = x*1
    y[y<0] = np.nan
    return y
#%%

amsrDtwain = pd.read_csv("marktwain_6year.csv",parse_dates=[0], skiprows=8)
#amsrAtwain = pd.read_csv("amsrAtwain.csv",parse_dates=[0], skiprows=8)

amsrDtwain["VODt"] = setna(amsrDtwain[' mean_LPRM_AMSRE_D_SOILM3_002_opt_depth_x'])
#amsrAtwain["VODt"] = setna(amsrAtwain[' mean_LPRM_AMSRE_A_SOILM3_002_opt_depth_x'])
#%%

amsrDtwain = pd.read_csv("amsr_mof_6year.csv",parse_dates=[0], skiprows=8)
#amsrAtwain = pd.read_csv("amsrAmoflux.csv",parse_dates=[0], skiprows=8)

amsrDtwain["VODt"] = setna(amsrDtwain[' mean_LPRM_AMSRE_D_SOILM3_002_opt_depth_x'])
#amsrAtwain["VODt"] = setna(amsrAtwain[' mean_LPRM_AMSRE_A_SOILM3_002_opt_depth_x'])
#%%
#ascending is 1:30 PM, descending is 1:30 AM

amsrDtwain["DateTime"] = amsrDtwain["time"] + datetime.timedelta(hours=1)
#amsrAtwain["DateTime"] = amsrAtwain["time"] + datetime.timedelta(hours=13.5)

alldata = pd.merge(simdata, amsrDtwain, on="DateTime",how="left")

alldata["Year"] = [x.year for x in alldata["DateTime"]]

morningdata = alldata[(alldata["Year"] >= 2005) * (alldata["Year"] <= 2010)].reset_index()
#%%
morningdata = morningdata.iloc[1::24,:]




#%%
#plt.plot(morningdata["DateTime"], morningdata["simVOD"])

#vodguess = (1 + 0.08*morningdata["Psi"]) * (0.55+0.01*morningdata["LAI"])
plt.figure(figsize=(10,10))
plt.subplot(3,1,1)
plt.plot(morningdata["DateTime"], morningdata["leafpot"],"b")
plt.subplot(3,1,2)
plt.plot(morningdata["DateTime"], morningdata["LAI_modis"],"g")
plt.subplot(3,1,3)
plt.plot(morningdata["DateTime"], morningdata["VODt"],"ro")
#%%
plt.plot(morningdata["LAI_modis"], morningdata["VODt"],"o")
#%%
lai_all = np.array(morningdata["LAI_modis"])
psi_all = np.array(morningdata["leafpot"])
vod_all = np.array(morningdata["VODt"])


temp_all = np.array(morningdata["T_AIR"])
#%%
doy = np.arange(len(vod_all)) % 365
#vod_all[vod_all < 0.75] = np.nan
#vod_all[temp_all <= 5] = np.nan
#vod_all[doy<100] = np.nan
#vod_all[doy>310] = np.nan
#%%
vod_smooth = 0*vod_all
for i in range(7, len(vod_all)-7):
    vod_smooth[i] = np.nanmean(vod_all[(i-7):(i+8)])
vod_smooth[:7] = vod_smooth[7]
vod_smooth[-7:] = vod_smooth[-7]

# myint = scipy.interpolate.UnivariateSpline(x=np.arange(len(vod_all))[np.isfinite(vod_all)],
#                                            y=vod_all[np.isfinite(vod_all)],
#                                            s = 3)

# vod_smooth = myint(np.arange(len(vod_all)))
#%%
lai_smooth = 0*vod_all
for i in range(7, len(vod_all)-7):
    lai_smooth[i] = np.nanmean(lai_all[(i-7):(i+8)])
lai_smooth[:7] = lai_smooth[7]
lai_smooth[-7:] = lai_smooth[-7]
#%%
plt.plot(vod_all[:],"o")
plt.plot(vod_smooth[:])



#%%
plt.plot(psi_all[lai_all > 3], vod_smooth[lai_all > 3],"o")
plt.plot([-2.5,-0.5],[0.9,0.98])
plt.plot([-2.5,0],[0.9,1.0])

#%%
lai_high = np.mean(lai_all[lai_all > 3])

#so when LAI is around lai_high, VOD has this relationship
psi_factor = 1 + 0.04*psi_all
#%%
myrem = vod_smooth / psi_factor
#%%
plt.plot(lai_all, myrem,"o")
#%%
psi_amp = -0.15/np.min(psi_all)
lai_amp = 0.15/(np.max(lai_all) - np.min(lai_all))
#%%
#lai_fac = 0.15*(lai-min(lai))/(max(lai)-min(lai)) + 0.85
# lai_int = 0.85-lai_amp*np.min(lai_all)
# #%%
# lai_int = 0.5-lai_amp*np.min(lai_all)
# #%%
# laifac = lai_int + lai_amp*lai_all
# psifac = 1 + psi_amp*psi_all
#%%
laifac = 0.82 + 0.051*lai_smooth
psifac = 1 + 0.067*psi_all
#%%
# laifac = 0.82 + 0.085*(lai_all-0.85)
# psifac = 1 + 0.085*psi_all
#%%
#plt.plot(vod_all,".")
plt.plot(morningdata.DateTime,vod_smooth,label="AMSR-E")
plt.plot(morningdata.DateTime,laifac*psifac,label="Model")
plt.ylabel("VOD")
plt.xlabel("Time")
plt.legend()
#%%
plt.plot(laifac*psifac, vod_smooth,"o")
plt.plot([0.8,1],[0.8,1])
#%%
plt.plot(vod_all,"ko")
plt.plot(laifac*psifac,"r",linewidth=3)
#%%
#psi_all_log = -np.log(-psi_all) - 2
#%%
def myfun(x):
    predvod = (1+x[0]*psi_all)*(x[1]+x[2]*lai_all)
    return np.nansum(np.abs(predvod-vod_all))
myopt = scipy.optimize.minimize(myfun,x0=np.array([0.067,0.82,0.051]),
                                bounds = ((0,1),(0,2),(0,1)))
#%%
laifac = myopt.x[1] + myopt.x[2]*lai_all
psifac = 1 + myopt.x[0]*psi_all

plt.plot(vod_all,".",label="Obs")
plt.plot(laifac*psifac,label="Model")
plt.xlabel("Time (days since 1 Jan 2005)")
plt.ylabel("VOD")
plt.legend()
#%%

