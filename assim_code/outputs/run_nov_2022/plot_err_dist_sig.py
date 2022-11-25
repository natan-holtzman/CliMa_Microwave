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

obs_names = ["All", "1 AM/PM", "6 AM/PM","1+6 sync.","1+6 offset"]

all_means = [0.8737915976892229, 0.3780452778386911, 2.395900386656546, 6.2669428043052715]
dry_means = [1.6770836191526692, 0.3329708971276553, 3.7649526834622287, 8.490219724548774]



def makeplot(err_file,means_list):
    
    all_errs = np.load(err_file)
    
    
    titles = ["Leaf water potential","Soil moisture","ET","GPP"]
    units = ["MPa","$m^3/m^3$","mm/day","$\mu mol/m^2/s$"]
    colors_list = ["tab:blue","tab:green","tab:orange","tab:red","tab:purple"]
    #xcoords = range(1,6)
    
    j = 0
    fig, ax_all = plt.subplots(2,2,figsize=(8,5))
    for ax in ax_all.ravel():#[:7]:
        meds = [np.median(all_errs[j][scenario]) for scenario in range(5)]
        ures = [1]
        for scenario in range(1,5):
            utest = scipy.stats.mannwhitneyu(all_errs[j][0], all_errs[j][scenario])
            ures.append(utest.pvalue)
        #print(ures)
        
        vpi = ax.violinplot(np.transpose(all_errs[j,:,:]),showmedians=True)
        for element in vpi.keys():
            if element=="bodies":
                for patch, color in zip(vpi[element], colors_list):
                    patch.set_color(color)
            else:
                vpi[element].set_color(colors_list)
        ax.set_title(titles[j])
        ax.set_ylabel(units[j])
        mymax = ax.get_ylim()[1]
        if j==0:
            mymax = np.max(all_errs[j,0,:])*1.33 #because of extreme outliers in LWP
        ax.set_ylim(0,mymax)
        ylim_units = ax.get_ylim()
        ax.plot([-1,6],[np.median(all_errs[j,0,:])]*2,"--", color="tab:blue")
        ax.set_xlim(0.25,5.75)
        
        newtix = []
        for scenario in range(1,5):
            if ures[scenario] < 0.05:
                if meds[scenario] < meds[0]:
                    sig_lab = "-*"
                else:
                    sig_lab = "+*"
            else:
                sig_lab = "ns"
            newtix.append(sig_lab)
        ax.set_xticks(np.arange(2,6),newtix)
    
            #ax.text(scenario+1,0,ttoplot,horizontalalignment="center")
        
        ax2 = ax.twinx()
        ax2.set_ylim(100 * np.array(ylim_units) / means_list[j])
        ax2.set_ylabel("Percent")
        
    
    #ticklabel_format(style="plain",axis="y")
        j += 1
    
    plt.tight_layout()

#plt.savefig("err_dry_sept19f.png")
makeplot("all_rmse_nov25.npy",all_means)
makeplot("dry_rmse_nov25.npy",dry_means)
makeplot("dry_rmse_nov25_hourly.npy",dry_means)
makeplot("dry_ubrmse_nov25.npy",dry_means)

