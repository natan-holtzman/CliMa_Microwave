# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 15:54:52 2022

@author: natan
"""

import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
import scipy.stats
#%%
#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));
plt.rcParams["lines.linewidth"] = 2;
plt.rcParams["mathtext.default"] = "regular"

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 5

plt.rcParams['font.size']=15
#%%
colors_list = ["#E69F00","#56B4E9","#009E73","#CC79A7","#0072B2"]

#colors_list = ["tab:blue","tab:green","tab:orange","tab:red","tab:purple","black"]
obs_names = ["HOURLY", "1 AM/PM", "6 AM/PM","1+6","1+6 offset","True model"]


#%%
mean_2007 = np.load("test_stats/y2007_meanvalue.npy")
bias_2007 = (mean_2007[:,:,:-1] - mean_2007[:,:,np.array([-1])]).reshape(4,5,120) 
#%%
am5_2007 = np.load("test_stats/y2007_5am_rmse.npy")
amp_2007 = np.load("test_stats/y2007_diurnal_rmse.npy")
#%%
cor_2007 = np.load("test_stats/y2007_daily_cor.npy")

summer_mase = np.load("test_stats/summer_hourly_mase.npy")
#%%

def add_sig(ax,pvalues,centers):

    min0,max0 = ax.get_ylim()
    ax.set_ylim(min0,max0+0.2*(max0-min0))

    newtix = ["0"]
    for scenario in range(1,4):
        if pvalues[scenario] < 0.05:
            if centers[scenario] < centers[0]:
                sig_lab = "(-)"
            else:
                sig_lab = "(+)"
        else:
            sig_lab = "(ns)"
        newtix.append(sig_lab)
    for ti in range(1,4):
        ax.text(ti+1,max0,newtix[ti],horizontalalignment="center",fontsize=18)
        

def colored_box(ax,data):
    vpi = ax.boxplot(np.transpose(data[:-1,:]),showfliers=False,patch_artist=True,notch=True)
    for patch, color in zip(vpi["boxes"], colors_list):
        patch.set_color(color)
        patch.set_alpha(0.75)
    for med in vpi["medians"]:
        med.set_color("black")
    utest = [scipy.stats.mannwhitneyu(data[0,:], data[j,:]).pvalue for j in range(4)]
    smeds = [np.nanmedian(data[j,:]) for j in range(4)]
    add_sig(ax,utest,smeds)

    #ttest = [scipy.stats.ttest_ind(data[0,:], data[j,:], equal_var=False).pvalue for j in range(4)]
    #smeans = [np.nanmean(data[j,:]) for j in range(4)]
    #add_sig(ax,ttest,smeans)


#%%
fig, ax_all = plt.subplots(2,2,figsize=(10,8))
ax = ax_all[0,0]
#plt.boxplot(am5_2007[0,:,:].T)
colored_box(ax,am5_2007[0,:,:])
mymax = ax.get_ylim()[1]
#ax.set_ylim(0,mymax)
#ax.set_title("(a)",loc="left",fontsize=28)
ax.text(-0.1,1.1,"(a)",fontsize=24,transform=ax.transAxes)
#ax.set_title("RMSE of pre-dawn $LWP^o$ (MPa)")
ax.set_title("RMSE of pre-dawn $\mathit{\psi_l^o}$ (MPa)")

ax.set_xticks(range(1,5),obs_names[:4])

ax = ax_all[0,1]
#plt.boxplot(am5_2007[0,:,:].T)
colored_box(ax,amp_2007[0,:,:])
mymax = ax.get_ylim()[1]
#ax.set_ylim(0,mymax)
#ax.set_title("(c)",loc="left",fontsize=28)
ax.text(-0.1,1.1,"(b)",fontsize=24,transform=ax.transAxes)
ax.set_title("RMSE of diurnal $\Delta \mathit{\psi_l^o}$ (MPa)")
ax.set_xticks(range(1,5),obs_names[:4])

#ax.set_xticks(range(1,6),obs_names[:-1])
ax = ax_all[1,0]
#plt.boxplot(am5_2007[0,:,:].T)
colored_box(ax,bias_2007[2,:,:])
#colored_box(ax,rmse_all[1,:,:])
ax.axhline(0,color='grey')
ax.set_title("Bias of daily ET (mm/day)")
#ax.set_title("(b)",loc="left",fontsize=28)
ax.text(-0.1,1.1,"(c)",fontsize=24,transform=ax.transAxes)
ax.set_xticks(range(1,5),obs_names[:4])

ax = ax_all[1,1]
#plt.boxplot(am5_2007[0,:,:].T)
colored_box(ax,cor_2007[2,:,:])
mymax = ax.get_ylim()[1]
#ax.set_ylim(0,mymax)
ax.set_title("Correlation of daily ET")
#ax.set_title("(d)",loc="left",fontsize=28)
ax.text(-0.1,1.1,"(d)",fontsize=24,transform=ax.transAxes)
ax.set_xticks(range(1,5),obs_names[:4])

fig.suptitle("Errors in assimilation year",fontsize=24)

# for coli in range(5):
#     ax.plot([],[],color=colors_list[coli], label=obs_names[coli],linewidth=4,alpha=0.75)
# fig.legend(loc="center",bbox_to_anchor=(0.5, -0.025),ncols=5,title="Observation scenario")

fig.tight_layout()
#%%
mean_all = np.load("test_stats/all_meanvalue.npy")
bias_all = (mean_all[:,:,:-1] - mean_all[:,:,np.array([-1])]).reshape(4,5,120) 

rmse3hr_all = np.load("test_stats/all_3hr_rmse.npy")

add_ubrmse = np.sqrt(rmse3hr_all**2 - bias_all**2)
#%%