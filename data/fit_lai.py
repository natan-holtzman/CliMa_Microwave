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
import scipy.interpolate
#%%
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 14
fig_size[1] = 7

plt.rcParams['font.size']=28
#%%
#rawdata = pd.read_csv('../../../AMF_US-MOz_BASE-BADM_8-5/AMF_US-MOz_BASE_HH_8-5.csv',comment='#')
#rawdata = alldata[20016:].copy()
#%%
#timestamps = [datetime.datetime.strptime(str(x), '%Y%m%d%H%M') for x in rawdata['TIMESTAMP_START']]
#hours = [(x - datetime.datetime(2000,1,1,0,0)).total_seconds()/(60*60.0) for x in timestamps]
#hours_ind = np.array(hours) - hours[0]
#%%
#rawdata['DateTime'] = timestamps
#%%

#%%

laidata = pd.read_csv('../../../modis_lai_point_mo.csv')
laidata['DateTime'] = [datetime.datetime.strptime(x[:-2], '%Y_%m_%d') for x in laidata['system:index']]
#%%
#plt.plot(laidata['DateTime'].loc[laidata['FparLai_QC'] < 50],laidata['Lai'].loc[laidata['FparLai_QC'] < 50],'o')
#plt.plot(laidata['DateTime'].loc[laidata['FparLai_QC'] >= 50],laidata['Lai'].loc[laidata['FparLai_QC'] >= 50],'o')
#%%
laidata['Lai'].loc[laidata['FparLai_QC'] >= 50] = np.nan


lai_avhrr = pd.read_csv('../../../moflux_lai_avhrr.csv')

lai_avhrr['DateTime'] = [datetime.datetime.strptime(x.split('_')[0], '%Y%m%d') for x in lai_avhrr['system:index']]

lai2 = pd.merge(lai_avhrr, laidata,how='left',on='DateTime')

#%%
lai_site = pd.read_csv('../../../moflux_lai_site.csv')
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
day_since_2000 = np.array([(x - datetime.datetime(2000,1,1,0,0)).total_seconds()/(60*60*24) for x in lai2["DateTime"]])
#interpolate MODIS to days when have site observations
modis_arr = np.array(lai2["Lai"])
site_arr = np.array(lai2["LAI_mean_site"])
modis_int = np.interp(day_since_2000,
                      day_since_2000[np.isfinite(modis_arr)], modis_arr[np.isfinite(modis_arr)])
#%%
df_int = pd.DataFrame({"MODIS":modis_int/10,"day":day_since_2000})
df_site = pd.DataFrame({"site":site_arr,"day":day_since_2000})
df_both = pd.merge(df_int, df_site, on="day",how="left")
#%%
#%%
df_both["DOY"] = df_both["day"] % 365
#%%
plt.figure()
plt.scatter(df_both["MODIS"],df_both["site"],c=df_both["DOY"])
#%%
mylag = 5
#%%
x1 = np.array(df_both["MODIS"])[:-mylag]
y1 = np.array(df_both["site"])[mylag:]
x1 = x1[np.isfinite(y1)]
y1 = y1[np.isfinite(y1)]
x1 = np.array(list(x1) + [0])
y1 = np.array(list(y1) + [0])
#%%
sdf = pd.DataFrame({"x":x1,"y":y1}).sort_values("x")
x2 = np.array(sdf["x"])
y2 = np.array(sdf["y"])
#x1[x1 <= 0.25] = np.nan
#%%
#plt.figure()
#plt.plot(np.log(x1)[:-mylag], np.log(np.array(df_both["site"])[mylag:]),"o")
myint = scipy.interpolate.UnivariateSpline(x=x2, y=y2)
#%%
xnew  = myint(np.array(df_both["MODIS"])[:-mylag])
xrange = np.arange(0,7,0.05)
xsamp  = myint(xrange)

#%%

plt.figure()
plt.scatter(np.array(df_both["MODIS"])[:-mylag],
            np.array(df_both["site"])[mylag:],c=np.array(df_both["DOY"])[mylag:])
plt.plot(xrange,xsamp,"k")
#%%
#%%
# plt.figure()
# plt.plot(np.arange(len(xnew))+mylag, xnew)
# plt.plot(df_both["site"],"o")
#%%
plt.figure()
plt.plot(xnew,np.array(df_both["site"])[mylag:],"o")
plt.plot([0,4],[0,4])

#%%
print(np.nanmean(np.abs(xnew - np.array(df_both["site"])[mylag:])))
#%%
lagged_fit = np.zeros(len(lai2))
lagged_fit[mylag:] = xnew
lagged_fit[:mylag] = lagged_fit[mylag]
lai2["LAI_fitted"] = lagged_fit
#%%
plt.figure()
plt.plot(lai2["DateTime"],lai2["LAI_fitted"])
plt.plot(lai2["DateTime"],lai2["LAI_mean_site"],"o")
#%%
lai2.to_csv("lai_spline_lag.csv")