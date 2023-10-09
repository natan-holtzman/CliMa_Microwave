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
#ascending is 1:30 PM, descending is 1:30 AM

amsrDtwain["DateTime"] = amsrDtwain["time"] + datetime.timedelta(hours=1)
#amsrAtwain["DateTime"] = amsrAtwain["time"] + datetime.timedelta(hours=13.5)

alldata = pd.merge(simdata, amsrDtwain, on="DateTime",how="left")

alldata["Year"] = [x.year for x in alldata["DateTime"]]

morningdata = alldata[(alldata["Year"] >= 2005) * (alldata["Year"] <= 2010)].reset_index()
#%%
morningdata = morningdata.iloc[1::24,:]

#%%
lai_all = np.array(morningdata["LAI_modis"])
psi_all = np.array(morningdata["leafpot"])
vod_all = np.array(morningdata["VODt"])


temp_all = np.array(morningdata["T_AIR"])
#%%
doy = np.arange(len(vod_all)) % 365

#%%
vod_smooth = 0*vod_all
for i in range(7, len(vod_all)-7):
    vod_smooth[i] = np.nanmean(vod_all[(i-7):(i+8)])
vod_smooth[:7] = vod_smooth[7]
vod_smooth[-7:] = vod_smooth[-7]

#%%
laifac = 0.82 + 0.051*lai_all
psifac = 1 + 0.067*psi_all
#%%
#plt.plot(vod_all,".")
plt.figure()
plt.plot(morningdata.DateTime,vod_smooth,label="LPRM")
plt.plot(morningdata.DateTime,laifac*psifac,label="CliMA")
plt.ylabel("VOD")
plt.xlabel("Time")
plt.legend()
#%%
