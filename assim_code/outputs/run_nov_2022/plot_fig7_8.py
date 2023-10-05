import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats
#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));
plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["mathtext.default"] = "regular"

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 18
fig_size[1] = 9

plt.rcParams['font.size']=15


out_folder = "./";


#obs_types = ['oAll','o1AMPM','o6AMPM',"o1and6","o16offset"];

obs_names = ["HOURLY", "1 AM/PM", "6 AM/PM","1+6"]

all_means = [0.8737915976892229, 0.3780452778386911, 2.395900386656546, 6.2669428043052715]
dry_means = [1.6770836191526692, 0.3329708971276553, 3.7649526834622287, 8.490219724548774]

def exclude_outliers(x):
    q25 = np.quantile(x,0.25)
    q75 = np.quantile(x,0.75)
    iqr = q75 - q25
    high_bd = q75 + 1.5*iqr
    low_bd = q25 - 1.5*iqr
    return x[(x > low_bd)*(x < high_bd)]

#titles = [r"$LWP^o$","Soil moisture","ET","GPP"]
titles = [r"$\mathit{\psi_l^o}$","Column-averaged\nsoil moisture","ET","GPP"]

units = ["MPa","$m^3/m^3$","mm/day","$\mu mol/m^2/s$"]
#colors_list = ["tab:blue","tab:green","tab:orange","tab:red","tab:purple"]

colors_list = ["#E69F00","#56B4E9","#009E73","#CC79A7","#0072B2"]

subplot_letters = "(a),(b),(c),(d)".split(",")

def makeplot(err_file,means_list,big_title):
    
    all_errs = np.load(err_file)
    
    #xcoords = range(1,6)
    
    j = 0
    fig, ax_all = plt.subplots(2,2,figsize=(10,8))
    for ax in ax_all.ravel():#[:7]:
        meds = [np.median(all_errs[j][scenario]) for scenario in range(5)]
        ures = [1]
        for scenario in range(1,4):
            utest = scipy.stats.mannwhitneyu(all_errs[j][0], all_errs[j][scenario])
            ures.append(utest.pvalue)
        # ures = [1]
        # for scenario in range(1,5):
        #     utest = scipy.stats.mannwhitneyu(all_errs[j][0], all_errs[j][scenario])
        #     ures.append(utest.pvalue)
        #print(ures)
        
        vpi = ax.boxplot(np.transpose(all_errs[j,:-1,:]),showfliers=False,patch_artist=True,notch=True)
        for patch, color in zip(vpi["boxes"], colors_list):
            patch.set_color(color)
            patch.set_alpha(0.75)
        for med in vpi["medians"]:
            med.set_color("black")
            #else:
            #    vpi[element].set_color(colors_list)
            
        max0 = ax.get_ylim()[1]
        ax.set_ylim(0,max0*1.2)

        newtix = ["0"]
        for scenario in range(1,4):
            if ures[scenario] < 0.05:
                if meds[scenario] < meds[0]:
                    sig_lab = "(-)"
                else:
                    sig_lab = "(+)"
            else:
                sig_lab = "(ns)"
            newtix.append(sig_lab)
        for ti in range(1,4):
            ax.text(ti+1,max0,newtix[ti],horizontalalignment="center",fontsize=18)
            
            
        ax.set_title(titles[j]+" (" + units[j] + ")")
        #ax.set_ylabel(units[j])
        mymax = ax.get_ylim()[1]
        ax.set_ylim(0,mymax)
        ylim_units = ax.get_ylim()
        #ax.plot([-1,6],[np.median(all_errs[j,0,:])]*2,"--", color="tab:blue")
        #ax.set_xlim(0.25,4.75)
        
        #ax.set_xticks(np.arange(1,6),obs_names,fontsize=10)
        #ax.set_xticks([],[])
        ax.set_xticks(range(1,5),obs_names[:4])
        ax.text(-0.1,1.1,subplot_letters[j],fontsize=24,transform=ax.transAxes)

            #ax.text(scenario+1,0,ttoplot,horizontalalignment="center")
        
        #ax2 = ax.twinx()
        #ax2.set_ylim(100 * np.array(ylim_units) / means_list[j])
        #ax2.set_ylabel("Percent")
        
    
    #ticklabel_format(style="plain",axis="y")
        j += 1
        
    
    # for coli in range(4):
    #     ax2.plot([],[],color=colors_list[coli], label=obs_names[coli],linewidth=4,alpha=0.75)
    # fig.legend(loc="center",bbox_to_anchor=(0.5, -0.025),ncols=5,title="Observation scenario")
    
    fig.suptitle(big_title,fontsize=24)
    
    plt.tight_layout()

#plt.savefig("err_dry_sept19f.png")
filepref = "test_stats/"

makeplot(filepref+"all_daily_rmse.npy",all_means,"RMSE over 13 years")
makeplot(filepref+"dry_daily_rmse.npy",dry_means,"RMSE over four dry summers")
makeplot(filepref+"all_hourly_rmse.npy",dry_means,"Hourly RMSE")

#%%
#prior_errs = np.load(r"C:\Users\natan\OneDrive - Stanford\Documents\moflux_docs\transfer_files_nov22\rmse_prior.npy")#.reshape(4,5,120)
all_errs_daily = np.load(filepref+"all_daily_rmse.npy")#.reshape(4,5,120)
dry_errs_daily = np.load(filepref+"dry_daily_rmse.npy")#.reshape(4,5,120)
#%%
et_daily_med = np.median(all_errs_daily[2,:,:],axis=1)
et_dry_med = np.median(dry_errs_daily[2,:,:],axis=1)
#%%
print(et_daily_med/et_daily_med[3])
print(et_dry_med/et_dry_med[3])
#%%
gpp_daily_med = np.median(all_errs_daily[3,:,:],axis=1)
gpp_dry_med = np.median(dry_errs_daily[3,:,:],axis=1)
print(gpp_daily_med/gpp_daily_med[3])
print(gpp_dry_med/gpp_dry_med[3])



