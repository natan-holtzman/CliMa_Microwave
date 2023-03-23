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

laidata = pd.read_csv('modis_lai_point_mo.csv')
laidata['DateTime'] = [datetime.datetime.strptime(x[:-2], '%Y_%m_%d') for x in laidata['system:index']]
#%%
#plt.plot(laidata['DateTime'].loc[laidata['FparLai_QC'] < 50],laidata['Lai'].loc[laidata['FparLai_QC'] < 50],'o')
#plt.plot(laidata['DateTime'].loc[laidata['FparLai_QC'] >= 50],laidata['Lai'].loc[laidata['FparLai_QC'] >= 50],'o')
#%%
#laidata['Lai'].loc[laidata['FparLai_QC'] >= 50] = np.nan


lai_avhrr = pd.read_csv('moflux_lai_avhrr.csv')

lai_avhrr['DateTime'] = [datetime.datetime.strptime(x.split('_')[0], '%Y%m%d') for x in lai_avhrr['system:index']]

lai2 = pd.merge(lai_avhrr, laidata,how='left',on='DateTime')

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
plt.plot(lai_site['DateTime'],lai_site["LAI_SE"],"o-")
plt.plot(lai_site['DateTime'],lai_site["LAI_S"],"o-")
plt.plot(lai_site['DateTime'],lai_site["LAI_SW"],"o-")
plt.plot(lai_site['DateTime'],lai_site["LAI_W"],"o-")
plt.plot(lai_site['DateTime'],lai_site["LAI_NW"],"o-")
#%%
plt.plot(lai_site['DateTime'],lai_site["LAI_mean_site"],"o-")

#%%
lai2 = pd.merge(lai2, lai_site[['DateTime','LAI_mean_site']],how='left',on='DateTime')
#%%
day_since_2000 = np.array([(x - datetime.datetime(2000,1,1,0,0)).total_seconds()/(60*60*24) for x in lai2["DateTime"]])
#interpolate MODIS to days when have site observations
modis_arr = np.array(lai2["Lai"])
avhrr_arr = np.array(lai2["LAI"]/1000)

site_arr = np.array(lai2["LAI_mean_site"])
modis_int = np.interp(day_since_2000,
                      day_since_2000[np.isfinite(modis_arr)], modis_arr[np.isfinite(modis_arr)])
avhrr_int = np.interp(day_since_2000,
                      day_since_2000[np.isfinite(avhrr_arr)], avhrr_arr[np.isfinite(avhrr_arr)])


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
plt.plot([0,5],[0,5],"k")
#%%
plt.plot(df_both["day"],3.5/5.7*(df_both["MODIS"] - 0.51)+0.92)
#plt.plot(df_both["day"],(df_both["MODIS"]-1)*4/5.5+1)
plt.plot(df_both["day"],df_both["site"])
#plt.xlim()
#%%
#plt.plot(4/5.5*(df_both["MODIS"] - 1)+1, df_both["site"],"o")
plt.plot(3.5/5.7*(df_both["MODIS"] - 0.51)+0.92, df_both["site"],"o")
plt.plot([0,5],[0,5])
#%%
x1 = np.array(df_both["site"])[20:-20]
y0 = 0.55*(np.array(df_both["MODIS"])-1)+1
#%%
laglist = np.arange(-10,11)
rmse_list = []
for mylag in laglist:
    y1 = y0[(20+mylag):(-20+mylag)]
    rmse_list.append(np.sqrt(np.nanmean((x1-y1)**2)))
#%%
spring = np.array((df_both["DOY"] >= 100) * (df_both["DOY"] < 175))
summer = np.array((df_both["DOY"] >= 175) * (df_both["DOY"] < 250))
fall = np.array((df_both["DOY"] >= 250))
#%%
plt.plot(df_both["DOY"], y0,".")
plt.plot(df_both["DOY"], df_both["site"],"o")
#%%
y1a = y0[20:-20]
y1b = y0[10:-30]
doy2 = np.array(df_both["DOY"])[20:-20]# > 250
#%%
myX = np.stack((y1a,y1b,y1a*(doy2 >= 180),y1b*(doy2 >= 180),doy2 >= 180),1)
mymod = sm.OLS(x1, sm.add_constant(myX),missing='drop').fit()
#%%
mypred = mymod.predict(sm.add_constant(myX))
plt.plot(mypred,x1,"o")
plt.plot([0,4],[0,4])
#%%
ypred = 1*y0
ypred[20:-20] = mypred
#%%
plt.plot(df_both["day"], ypred)
plt.plot(df_both["day"], df_both["site"],"o")
#%%
plt.plot(y1a[doy2<170],x1[doy2<170],"o")
#plt.plot(y1a[doy2 >= 250],x1[doy2 >= 250],"o")
plt.plot(y1b[doy2 >= 170],x1[doy2 >= 170],"o")
plt.plot([0,4],[0,4])
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
lai2["LAI_fitted"] = ypred

#%%
plt.figure()
plt.plot(lai2["DateTime"],lai2["LAI_fitted"])
plt.plot(lai2["DateTime"],lai2["LAI_mean_site"],"o")
#%%
lai2.to_csv("lai_seasonal_lag.csv")