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
ubrmse_all = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\all_ubrmse_nov30.npy")
rmse_all = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\all_rmse_nov25.npy")

mean_full = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\means_nov25.npy")
bias_all = (mean_full[:,:-1] - mean_full[:,np.array([-1])]).reshape(4,5,120) 

mean_2007 = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\mean2007.npy")
bias_2007 = (mean_2007[:,:-1] - mean_2007[:,np.array([-1])]).reshape(4,5,120) 


#std_2007 = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\std_rmse_2007.npy")
am5_2007 = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\am5_rmse_2007.npy")
#am5_2007 = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\am5_rmse_all.npy")

#cv_2007 = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\cv_rmse_2007.npy")
amp_2007 = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\d52_rmse_2007.npy")

#amp_2007 = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\d52_rmse_all.npy")

#%%
anom_cor = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\anom_cor_nov25.npy").reshape(4,5,120)
all_cor = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\all_cor_nov28.npy").reshape(4,5,120)
cor_2007 = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\cor_2007.npy").reshape(4,5,120)

all_mae = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\mase_dec31_hourly.npy")#.reshape(4,5,120)
summer_mae = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\mase_dec31_hourly_summer.npy")#.reshape(4,5,120)

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
ax.set_title("RMSE of pre-dawn $LWP^o$ (MPa)")
ax.set_xticks(range(1,5),obs_names[:4])

ax = ax_all[0,1]
#plt.boxplot(am5_2007[0,:,:].T)
colored_box(ax,amp_2007[0,:,:])
mymax = ax.get_ylim()[1]
#ax.set_ylim(0,mymax)
#ax.set_title("(c)",loc="left",fontsize=28)
ax.text(-0.1,1.1,"(b)",fontsize=24,transform=ax.transAxes)
ax.set_title("RMSE of diurnal $\Delta LWP^o$ (MPa)")
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

rmse3hr_all = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\rmse_3hr_jan18.npy")
cor3hr_all = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\cor_3hr_jan18.npy")


add_ubrmse = np.sqrt(rmse3hr_all**2 - bias_all**2)
#%%
# plt.rcParams['font.size']=20

# fig, ax = plt.subplots(1,1,figsize=(10,8))
# colored_box(ax,add_ubrmse[1,:,:])
# #plt.ylim(0,0.075)
# plt.xticks(range(1,6),obs_names[:5])
# plt.xlabel("Observation scenario")
# plt.title("Additive ubRMSE of soil moisture $(m^3/m^3)$")
# #%%

# fig, ax = plt.subplots(1,1,figsize=(10,8))
# colored_box(ax,all_cor[3,:,:])
# plt.ylim(0.5,1)
# plt.xticks(range(1,6),obs_names[:5])
# plt.xlabel("Observation scenario")
# plt.title("Correlation of soil moisture")
# #%%
# stds = [0.7823546800694551, 0.03949884314393721, 4.014208876396588, 9.024972530415111]
# stds_summer = [0.9456777104965098, 0.03855216081718514, 4.9815970376196645, 10.770501028458469]

# fig, ax = plt.subplots(1,1,figsize=(10,8))
# colored_box(ax,summer_mae[2,:,:]/stds_summer[2])
# #plt.ylim(0.5,1)
# plt.xticks(range(1,6),obs_names[:5])
# plt.xlabel("Observation scenario")
# plt.title("MASE of ET")
# #%%
# plt.rcParams['font.size']=20

# fig, ax = plt.subplots(1,1,figsize=(10,8))
# colored_box(ax,all_cor[2,:,:])
# #plt.ylim(0,0.075)
# plt.xticks(range(1,6),obs_names[:5])
# plt.xlabel("Observation scenario")
# plt.title("Additive ubRMSE of ET (mm/day)$")
#%%
# plt.plot(std_2007[0,:,:].reshape(600,-1),am5_2007[0,:,:].reshape(600,-1),'o')
# #%%
# plt.plot(am5_2007[0,:,:].reshape(600,-1),rmse_all[2,:,:].reshape(600,-1),'o')
# #%%
# plt.plot(std_2007[0,:,:].reshape(600,-1),ubrmse_all[2,:,:].reshape(600,-1),'o')