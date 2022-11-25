import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));
plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["font.size"] = 12;
plt.rcParams["mathtext.default"] = "regular"


print("ALL YEARS")
#prior_min = [ 0.1,  1e-6, 500, 0.75, 0.75,0.01];
#prior_max = [ 100, 2e-5, 3000, 10, 8,10];


#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));
#prior_min = [10, 0.01, 0.1,  1e-6, 500, 0.75, 0.75, 0.1,100];
#prior_max = [120,0.75,  100, 2e-5, 3000, 10, 8, 10,1000];
out_folder = "./";


#true_val = [22,0.33, 2, 1e-5,600, 4, 3, 1, 600];
par_names = ["Vcmax", "Stomatal margin", "Kmax_plant", "Kmax_soil", "Soil depth", "Kplant location", "Kplant shape","Volume factor","Medlyn g1"];


#plt.rcParams["lines.linewidth"] = 1;
#plt.rcParams["font.size"] = 16;

#for each run in FixET, get the parameters
out_folder = "./";

#subdir_list = ["oAll_c2/", "o1AMPM_c2/","o6AMPM_c2/", "o1and6_c2/", "o16offset_c2/"]
#subdir_list2 = ["oAll_c3/", "o1AMPM_c3/", "o6AMPM_c3/", "o1and6_c3/", "o16offset_c3/"]
#subdir_list3 = ["oAll_c1/", "o1AMPM_c1/", "o6AMPM_c1/", "o1and6_c1/", "o16offset_c1/"]

obs_types = ['oAll','o1AMPM','o6AMPM',"o1and6","o16offset"];



obs_names = ["All", "1 AM/PM", "6 AM/PM","1+6 sync.","1+6 offset"]


#all_errs = np.load("all_rmse.npy")
all_errs = np.load("dry_rmse.npy")

all_means = [0.8737915976892229, 0.3780452778386911, 2.395900386656546, 6.2669428043052715]
dry_means = [1.6770836191526692, 0.3329708971276553, 3.7649526834622287, 8.490219724548774]


titles = ["Leaf water potential","Soil moisture","ET","GPP"]
units = ["MPa","$m^3/m^3$","mm/day","$\mu mol/m^2/s$"]
colors_list = ["tab:blue","tab:green","tab:orange","tab:red","tab:purple"]
xcoords = range(1,6)
j = 0
fig, ax_all = plt.subplots(2,2,figsize=(8,5))
for ax in ax_all.ravel():#[:7]:
    vpi = ax.violinplot(np.transpose(all_errs[j,:,:]),showmedians=True)
    for element in vpi.keys():
        if element=="bodies":
            for patch, color in zip(vpi[element], colors_list):
                patch.set_color(color)
        else:
            vpi[element].set_color(colors_list)
    ax.set_xticks([],[])
    ax.set_title(titles[j])
    ax.set_ylabel(units[j])
    mymax = ax.get_ylim()[1]
    #if j==0:
    #    mymax = 0.6  #because of extreme outliers in LWP
    ax.set_ylim(0,mymax)
    ylim_units = ax.get_ylim()
    ax.plot([-1,6],[np.median(all_errs[j,0,:])]*2,"--", color="tab:blue")
    ax.set_xlim(0.25,5.75)
    ax2 = ax.twinx()
    ax2.set_ylim(100 * np.array(ylim_units) / all_means[j])
    ax2.set_ylabel("Percent")

#ticklabel_format(style="plain",axis="y")
    j += 1

plt.tight_layout()

plt.savefig("err_dry_sept19f.png")


