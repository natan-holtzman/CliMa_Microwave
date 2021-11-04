# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 17:39:53 2021

@author: natan
"""



import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import statsmodels.api as sm
#%%
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 15
fig_size[1] = 10

plt.rcParams['font.size']=20
#%%
prior_min = [10, 0.075, 0.1, 1e-5, 1, 0.25, 500, 0.4, 0.005, 30, 0.75];
prior_max = [120, 0.3, 100, 2e-3, 8, 5, 3000, 0.6, 0.5, 1200, 8];

#true_val = [22, 1.5, 4, 1e-4, 2.5, 2.0, 700,0.45];
par_names = ["Vcmax", "Scrit", "Kmax_plant",
             "Kmax_soil", "B_soil", "Psi_sat", "Z_soil", "n_soil",
             "Slope","g1","WeibC"]
#%%
g1 = np.array(pd.read_csv("post_files_nov2/post_calib_medlyn0.csv"))[:,:12]

#g1 = np.array(pd.read_csv("post_allobs_errflex_8par.csv"))[:,:8]
#%%

#g1 = np.array(pd.read_csv("init_par_LL.csv"))[:,:8]

#%%
j = 0
fig, ax_all = plt.subplots(3,4)
for ax in ax_all.ravel()[:11]:
    ax.plot(np.exp(g1[:,j]))
    ax.set_ylim(prior_min[j],prior_max[j])
    ax.set_title(par_names[j])
    ax.set_xticks([])
    #ax.hlines(true_val[j], 0, 10000, color='red')
    j += 1
plt.tight_layout()
ax_all[2,3].axis('off')



#%%
j = 0
fig, ax_all = plt.subplots(3,3)
for ax in ax_all.ravel()[:8]:
    ax.hist(np.exp(g1[7500::25,j]))
    ax.set_xlim(prior_min[j],prior_max[j])
    #ax.vlines(true_val[j], 0, ax.get_ylim()[1], color='red')
    ax.set_title(par_names[j])
    ax.set_yticks([])
    j += 1
plt.tight_layout()
ax_all[2,2].axis('off')
#%%
trueB = 2.5
trueP20 = 2.0

postB = np.exp(g1[7500::25,4])
postP20 = np.exp(g1[7500::25,5])

B_dist = (postB - trueB)**2 + (postP20 - trueP20)**2
best_soil = np.argmin(B_dist)

#%%
mycov = np.cov(g1[7500::25,:].T)
detPrior = 0.0007051221222183998
detPost = np.linalg.det(mycov)

detPrior_noK = 0.000200395992236997

pmask = np.ones(8); pmask[2] = 0; pmask = pmask==1

mycov = np.cov(g1[7500::25,pmask].T)
detPost_noK = np.linalg.det(mycov)
#%%
et_daily = np.array(pd.read_csv("post_files_nov2/postET_scrit.csv"))[:,200:]
smc_file = np.array(pd.read_csv("post_files_nov2/postSMC_scrit.csv"))[:,200:]

#%%
fig_size[0] = 12
fig_size[1] = 8

fig, ax_all = plt.subplots(2,1)
ax = ax_all[0]
ax.plot(smc_file[:,:-1],color='grey',alpha=0.1)
ax.plot(smc_file[:,-1],color='red')
ax.set_title("Root zone soil moisture")

ax = ax_all[1]
ax.plot(et_daily[:,:-1],color='grey',alpha=0.1)
ax.plot(et_daily[:,-1],color='red')
ax.set_title("Daily ET (mol/s/m2)")

plt.tight_layout()
#%%
fig_size[0] = 12
fig_size[1] = 8

fig, ax_all = plt.subplots(2,2)
ax = ax_all[0,0]
ax.plot(pot_predawn[:,best_soil],color='blue',alpha=0.5)
ax.plot(pot_predawn[:,-1],color='red')
ax.set_title("Pre-dawn LWP (MPa)")

ax = ax_all[0,1]
ax.plot(pot_midday[:,best_soil],color='blue',alpha=0.5)
ax.plot(pot_midday[:,-1],color='red')
ax.set_title("Mid-day LWP (MPa)")

ax = ax_all[1,0]
ax.plot(rzsm_file[:,best_soil],color='blue',alpha=0.5)
ax.plot(rzsm_file[:,-1],color='red')
ax.set_title("Root zone soil moisture")

ax = ax_all[1,1]
ax.plot(et_daily[:,best_soil],color='blue',alpha=0.5)
ax.plot(et_daily[:,-1],color='red')
ax.set_title("Daily ET (mol/s/m2)")

plt.tight_layout()
#%%
def getstats(xtab):
    diffs = xtab[:,:-1] - xtab[:,-1].reshape((-1,1))
    rmse = np.sqrt(np.mean(diffs**2,axis=0)) / np.abs(np.mean(xtab[:,-1]))
    
    bias = np.mean(diffs, axis=0) / np.abs(np.mean(xtab[:,-1]))
    
    cor = np.corrcoef(xtab.T)[-1,:-1]
    
    #rsq = 1-np.mean(diffs**2, axis=0) / np.mean(xtab[:,-1]**2)
    
    return rmse, bias, cor #, rsq
#%%
pd_stats = getstats(pot_predawn)
md_stats = getstats(pot_midday)

rzsm_stats = getstats(rzsm_file)  
et_stats = getstats(et_daily)  
#%%
stat_tab = pd.DataFrame()
stat_tab["Variable"] = ["PreDawnLWP", "MidDayLWP", "RZSM", "ET"]
stat_tab["RMSE"] = [np.median(x[0]) for x in [pd_stats, md_stats, rzsm_stats, et_stats]]
stat_tab["Bias"] = [np.median(x[1]) for x in [pd_stats, md_stats, rzsm_stats, et_stats]]
stat_tab["Correlation"] = [np.median(x[2]) for x in [pd_stats, md_stats, rzsm_stats, et_stats]]


#%%
#pot_file = np.array(pd.read_csv("post_prior/psi_0post.csv"))

#pot_file = np.array(pd.read_csv("post_prior/psi_prior.csv"))

#pot_file = np.array(pd.read_csv("posts_oct2/psi_8par_3day.csv"))

pot_file = np.array(pd.read_csv("post_fix/psi_fix_3day.csv"))

#pot_file = (pot_file[::2,:] + pot_file[1::2,:])/2
#%%
pot_cycle = np.mean(np.array(pot_file).reshape((24,-1,101), order="F"), axis=1)
#%%
plt.plot(pot_cycle[:,:-1],color='grey',alpha=0.5)
plt.plot(pot_cycle[:,-1],color='red',linewidth=4)
#plt.ylim([-2,0])
plt.xlabel("Hour of day")
plt.ylabel("LWP (MPa)")
#%%
pot_cycle_norm = (pot_cycle-np.mean(pot_cycle,0))/np.std(pot_cycle,0)

#%%
plt.plot(pot_cycle_norm[:,:-1],color='grey',alpha=0.5)
plt.plot(pot_cycle_norm[:,-1],color='red',linewidth=4)
plt.xlabel("Hour of day")
plt.ylabel("Normalized LWP (MPa)")
#%%
