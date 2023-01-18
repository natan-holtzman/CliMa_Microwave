# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 10:40:52 2022

@author: natan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
import scipy.optimize
#%%
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 14
fig_size[1] = 7

plt.rcParams['font.size']=28
#%%
alldata = pd.read_csv('moflux_land_data_newLAI.csv')
psi = np.array(alldata["LWP_predawn"])
smc = np.array(alldata["SMC"])
no_data = np.isnan(psi)
smc = smc[np.logical_not(no_data)]
psi = psi[np.logical_not(no_data)]
years = np.array(alldata["YEAR"])[np.logical_not(no_data)]
#%%
yearlist = np.unique(years)
#%%
ti = yearlist[7]
x = smc[years==ti]
y = psi[years==ti]
xR = (x-0.067)/(0.55-0.067)
#%%
def pot_err(params):
    s0, s1, n, alpha = params
    relsoil = np.clip((x-s0)/(s1-s0),0.01,0.99)
    m = 1 - 1/n
    predpot = -1/alpha * (relsoil ** (-1/m) - 1) ** (1/n)
    return np.mean(np.abs(predpot - y))
#%%
myopt = scipy.optimize.minimize(pot_err,np.array([0.05,0.55,1.2,100]),
                                bounds = ((0,0.1),(0.4,0.6),(1.05,2),(1,1000)))
#%%
s0, s1, n, alpha = myopt.x
#relsoil = np.clip((x-s0)/(s1above),0.01,0.99)
m = 1 - 1/n
#predpot = -1/alpha * (relsoil ** (-1/m) - 1) ** (1/n)
#%%
srange = np.arange(np.min(x)-0.01,np.max(x),0.005)
relsoil = np.clip((srange-s0)/(s1-s0),0.01,0.99)
m = 1 - 1/n
predpot = -1/alpha * (relsoil ** (-1/m) - 1) ** (1/n)

#%%
plt.plot(x,y,"o")
plt.plot(srange,predpot)
#%%
def pot_err(params):
    n, alpha = params
    m = 1 - 1/n
    predpot = -1/alpha * (xR ** (-1/m) - 1) ** (1/n)
    return np.mean(np.abs(predpot - y))
#%%
myopt = scipy.optimize.minimize(pot_err,np.array([1.3,25]),
                                bounds = ((1.05,2),(1,1000)))
n, alpha = myopt.x
m = 1 - 1/n
srange = np.arange(np.min(x)-0.01,np.max(x)+0.01,0.005)
relsoil = np.clip((srange-0.067)/(0.55-0.067),0.01,0.99)
predpot = -1/alpha * (relsoil ** (-1/m) - 1) ** (1/n)
#%%
nlist = np.arange(1.1,2,0.01)

bestpar = []

for ti in range(len(yearlist)):
    x = smc[years==yearlist[ti]]
    y = psi[years==yearlist[ti]]
    xR = np.clip((x-0.067)/(0.55-0.067),0.01,0.99)
    
    corlist = []
    for n in nlist:
        m = 1 - 1/n
        predpot_base = -(xR ** (-1/m) - 1) ** (1/n)
        corlist.append(np.corrcoef(predpot_base,y)[0,1])
    bestN = nlist[np.argmax(corlist)]
    n = bestN
    m = 1 - 1/n
    predpot_base = -(xR ** (-1/m) - 1) ** (1/n)
    alphaI = np.std(predpot_base)/np.std(y)
    bestpar.append([bestN,alphaI])
bestpar = np.array(bestpar)
#%%
ti = 6

x = smc[years==yearlist[ti]]
y = psi[years==yearlist[ti]]
srange = np.arange(np.min(x)-0.01,np.max(x),0.005)
relsoil = np.clip((srange-s0)/(s1-s0),0.01,0.99)
xR = np.clip((x-0.067)/(0.55-0.067),0.01,0.99)
n,alphaI = bestpar[ti]
n = 1.3
m = 1 - 1/n
predpot_base = -(xR ** (-1/m) - 1) ** (1/n)
alphaI = np.std(predpot_base)/np.std(y)
predpot_range = -(relsoil ** (-1/m) - 1) ** (1/n)
plt.figure()
plt.plot(x,y,"o")
plt.plot(srange,predpot_range/alphaI)
    #plt.plot([0,-2],[0,-2])