# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 21:45:22 2022

@author: natan
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats
#%%
#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));
plt.rcParams["lines.linewidth"] = 2;
plt.rcParams["mathtext.default"] = "regular"

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 15
fig_size[1] = 7

plt.rcParams['font.size']=15
#%%
obs_names = ["All", "1 AM/PM", "6 AM/PM","1+6 sync.","1+6 offset"]

#%%
means = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\means_nov25.npy")
truemeans = means[:,-1]
retmeans = means[:,:-1].reshape(4,5,120)
#%%
clim_cor = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\clim_cor_nov25.npy").reshape(4,5,120)
anom_cor = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\anom_cor_nov25.npy").reshape(4,5,120)
#%%
plt.figure();
plt.violinplot(retmeans[1,:,:].T,showmedians=True);
plt.axhline(y=truemeans[1],color='k');
#%%
plt.figure()
plt.violinplot(clim_cor[3,:,:].T,showmedians=True)
#%%
plt.figure()
plt.violinplot(anom_cor[3,:,:].T,showmedians=True)
#%%
plt.figure()
plt.boxplot(anom_cor[3,:,:].T,showfliers=False)
#%%
titles = ["Leaf water potential","Soil moisture","ET","GPP"]
units = ["MPa","$m^3/m^3$","mm/day","$\mu mol/m^2/s$"]
#%%
varj = 0
plt.figure(figsize=(9,9));
plt.subplot(3,1,1)
plt.boxplot(retmeans[varj,:,:].T,showfliers=False);
plt.axhline(y=truemeans[varj],color='blue');
plt.xticks([],[])
plt.ylabel(units[varj])
plt.title("Overall mean")

plt.subplot(3,1,2)
plt.boxplot(clim_cor[varj,:,:].T,showfliers=False)
plt.xticks([],[])
plt.title("Climatology correlation")
plt.ylim(0.8,1)

plt.subplot(3,1,3)
plt.boxplot(anom_cor[varj,:,:].T,showfliers=False)
plt.xticks(range(1,6),obs_names)
plt.title("Anomaly correlation")
plt.ylim(0.8,1)


plt.suptitle(titles[varj],fontsize=22)

