# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 14:52:29 2021

@author: natan
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import scipy.interpolate
import scipy.optimize
#%%
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 14
fig_size[1] = 7

plt.rcParams['font.size']=20

import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%b %Y')
#%%

#alldata = pd.read_csv('moflux_land_data_skipyear_hourly2.csv').iloc[(24*365*9):(24*365*12),:]
alldata = pd.read_csv('moflux_land_data_skipyear_hourly2.csv')
#%%
allET_mols = np.clip(np.array(alldata['LE'])/44200, 0, 1000)
#W/m2 to mol/m2/s to kg/m2/s = mm/s to mm/s to mm/hr, integrate over hours = mm to m
daily_cum_ET = np.cumsum(allET_mols)[::24]*18.02/1000 * (60*60)  / 1000

cor_p = np.array(alldata['RAIN'])

#per half hour -> per hour -> integrate over hours
daily_cum_P = np.cumsum(cor_p)[::24] / 1000 * 2
#daily_cum_P[:3000] *= 4.7/3.7
#%%
plt.plot(daily_cum_ET)
plt.plot(daily_cum_P)
#%%
plt.plot(daily_cum_P - daily_cum_ET)
#%%
plt.plot(np.diff(daily_cum_ET[::365]),"o-")
plt.plot(np.diff(daily_cum_P[::365]),"o-")
#%%
plwp = np.array(alldata["LWP_predawn"])[6::24]

daily_ET = np.diff(daily_cum_ET)
daily_P = np.diff(daily_cum_P)

def soil_water_balance(limit,pmult):
    ans = np.zeros(len(daily_cum_P)-1)
    x = min(0,limit)
    for i in range(len(daily_cum_P)-1):
        ans[i] = x
        x = min(limit, x + daily_P[i]*pmult - daily_ET[i])
    return ans

#%%
z = soil_water_balance(0,1)
plt.plot(z)
plt.twinx()
plt.plot(-np.log(-plwp),"ko")
#%%
plt.plot(z[:3000], -np.log(-plwp[1:])[:3000], "o")
#%%
z2 = 0.45 + z/2

sarr = np.arange(0.19,0.45,0.005)
a = 1
n = 1.8
#a = 1.5
#n = 1.6
m = 1-1/n
parr = (((sarr - 0.067) / (0.45 - 0.067)) ** (-1/m) -1) ** (1/n) / a
#%%
plt.plot(z2, plwp[1:], "o")
plt.plot(sarr,-parr)
#%%
alpha = 0.5
nsoil = 1.4
msoil = 1-1/nsoil
s_pred = ((-plwp * alpha) ** nsoil + 1) ** (-msoil)
#%%
plt.plot(s_pred,"o")
plt.plot((z2-0.067)/(0.45-0.067))
#%%
z2rel = (z2-0.067)/(0.45-0.067)
#%%
plt.plot(z2rel[:400],plwp[:400],"o")
p_arr = np.arange(-2,0,0.01)
alpha = 0.5
nsoil = 1.6
msoil = 1-1/nsoil
s_arr = ((-p_arr * alpha) ** nsoil + 1) ** (-msoil)
plt.plot(s_arr,p_arr)

alpha = 1.5
nsoil = 1.3
msoil = 1-1/nsoil
s_arr = ((-p_arr * alpha) ** nsoil + 1) ** (-msoil)
plt.plot(s_arr,p_arr)
#%%
x = z2rel[:400]
y = plwp[:400]
x = x[np.isfinite(y)]
y = y[np.isfinite(y)]
#%%
def myerr(pars):
    n,a = pars
    m = (1-1/n)
    ypred_raw = (x ** (-1/m) -1) ** (1/n)
    #ypred_scale = -1.0/(0.81 ** (-1/m) -1) ** (1/n)
    return np.mean((y- ypred_raw/a)**2)
myopt = scipy.optimize.minimize(myerr, x0=np.array([1.4,15]), bounds=((0,10),(10,200)))
#%%
plt.plot(x,y,"o")
n,mscale = myopt.x
m = (1-1/n)*mscale
s_arr = np.arange(0.75,1,0.01)
ypred_raw = (s_arr ** (-1/m) -1) ** (1/n)
ypred_scale = -1.0/(0.81 ** (-1/m) -1) ** (1/n)
plt.plot(s_arr,ypred_scale*ypred_raw)
#%%
plt.plot(z2rel[2000:3000],plwp[2000:3000],"o")
p_arr = np.arange(-2,0,0.01)
for nsoil in np.arange(1.2,3.2,0.2):
    #nsoil = 1.6
    msoil = 1-1/nsoil
    s_arr = np.arange(0.63,1,0.01)
    ypred_raw = (s_arr ** (-1/msoil) -1) ** (1/nsoil)
    ypred_scale = -1.5/(0.8 ** (-1/msoil) -1) ** (1/nsoil)
    p_arr = ypred_scale*ypred_raw
    plt.plot(s_arr,p_arr,"k",alpha=0.5)
#%%
plt.plot(np.log(z2rel),-np.log(-plwp[1:]),"o")
#%%
import statsmodels.api as sm
#%%
myOLS = sm.OLS(-np.log(-plwp[1:]), sm.add_constant(np.log(z2rel[:])),missing="drop").fit()
#%%
plt.plot(z2rel[:], plwp[1:],"o")

sarr2 = np.arange(0.4,1,0.01)
plt.plot(sarr2, -np.exp(-1.3) * sarr2**(-2.7))
#%%
# rev_P = daily_P[-1::-1]
# rev_ET = daily_ET[-1::-1]

# def soil_water_balance_rev(limit,pmult):
#     ans = np.zeros(len(daily_cum_P)-1)
#     x = min(0,limit)
#     for i in range(len(daily_cum_P)-1):
#         ans[i] = x
#         x = min(limit, x - rev_P[i]*pmult + rev_ET[i])
#     return ans[-1::-1]
# #%%
# z = soil_water_balance(0,1.9)
# zR = soil_water_balance_rev(0,1.9)
# plt.plot(z)
# plt.plot(zR)
# plt.twinx()
# plt.plot(-np.log(-plwp),"ko")