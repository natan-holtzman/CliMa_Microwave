import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats
#import pandas as pd
#%%

full_par_names = ["Vcmax $(\mu mol/m^2/s)$", "*Stomatal P63 (-MPa)", "*Max. xylem\nconductance $(mol/m^2/s/MPa)$",  "Soil depth (mm)", "*Xylem P63 (-MPa)",
 "*Max. water storage\nvolume $(kg/m^2)$","Ball-Berry g1","*Soil conductivity $(\mu m/s)$","Soil boundary condition","Root profile ($m^{-1}$)",
"Van Genuchten n","Water potential\nof dry soil (MPa)"];

# par_names_no_star = ["$V_{cmax}$", r"$P63_{\beta}$",
#              "$k_{plant}$","Z","$S$",
#              "V","$g1$","$k_{soil}$",
#              "$s_{lower}$",r"$\alpha_{root}$",
#              "n",r"$\psi_{ref}$"]


par_names = ["$V_{cmax}$", r"$P63_{\beta}^o$",
              "$k_{plant}^o$","Z",r"$P63_x^o$",
              "$V^o$","$g_1$","$k_{soil}^o$",
              "$s_{lower}$",r"$\alpha_{root}$",
              "n",r"$\psi_{ref}$"]


# par_names = ["$V_{cmax}$", r"$\overline{P63_{\beta}}$",
#              r"$\overline{k_{plant}}$","Z",r"$\overline{P63_x}$",
#              r"$\overline{V}$","$g_1$",r"$\overline{k_{soil}}$",
#              "$s_{lower}$",r"$\alpha_{root}$",
#              "n",r"$\psi_{ref}$"]

vod_names = [r"$a^o$","b","c","$\omega$"]

#order_index = [2,5,4,10,3,0,6,1,7,8,9,11,12]

order_index = [[0,6,1],[2,5,4],[3,8,9],[10,11,7]]


plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["font.size"] = 12;
plt.rcParams["mathtext.default"] = "regular"
#plt.rcParams['text.usetex'] = True;
#for each run in FixET, get the parameters

allpars = np.load("all_param_dec3.npy")[:-1,:,:]

#allpars[:,:,4] = np.log(np.exp(allpars[:,:,4])-np.exp(allpars[:,:,1]))
allpars[:,:,5] += np.log(12) - np.log(11)
#allpars[:,:,-4] /= np.exp(allpars[:,:,11])


#prior_min = [10, 0.75, 0.1, 500, 0.01+0.75,  0.01*12,1,0.1,0.3,0.01,1,1.1];
#prior_max = [120,15,  50, 3000, 10+15, 10*12,120, 20,0.5,5,8,1.9];

prior_min = [10, 0.75, 0.1, 500, 0.01+0.75,  0.01*12,1,0.5e-1,0.3,0.01,1.1,0.125];
prior_max = [300,15,  50, 3000, 10+15, 10*12,120,    20,0.5,5,1.9,8];


true_val = [90, 3, 10, 2000, 4, 1.0*12, None,0.4,
         0.4,2,1.5,1]


obs_names = ["HOURLY", "1 AM/PM", "6 AM/PM", "1+6"]#[:3]
#obs_names = np.array(obs_names)[np.array([0,2,4])]

j = 0

#colors_list = ["tab:blue","tab:green","tab:orange","tab:red","tab:purple"]#[:3]
colors_list = ["#E69F00","#56B4E9","#009E73","#CC79A7","#0072B2"]
#colors_list = np.array(colors_list)[np.array([0,2,4])]

def single_panel_plain(ax,data,realval,title):
    #vpi = ax.violinplot(np.transpose(data),showmedians=True)
    vpi = ax.boxplot(np.transpose(data),showfliers=False,patch_artist=True,notch=True)

    for patch, color in zip(vpi["boxes"], colors_list):
        patch.set_color(color)
        patch.set_alpha(0.75)
    for med in vpi["medians"]:
        med.set_color("black")
            
    if realval:
        ax.hlines(realval, 0, 6, color="black")

    ax.set_xticks([])
    ax.set_xlim(0.25,4.75) 
    ax.set_title(title,fontsize=24)
    
    min0,max0 = ax.get_ylim()
    ax.set_ylim(min0-0.1*(max0-min0),max0+0.2*(max0-min0))
    for ti in range(1,4):
        ures = scipy.stats.mannwhitneyu(data[0,:],data[ti,:])
        med_diff = np.median(data[ti,:]) - np.median(data[0,:])
        if ures.pvalue < 0.05 and med_diff > 0:
            towrite = "(+)"
        elif ures.pvalue < 0.05 and med_diff <= 0:
            towrite = "(-)"
        else:
            towrite = "(ns)"
        ax.text(ti+1,max0,towrite,horizontalalignment="center",fontsize=18)




vod_index = [[-1],[-4,-3,-2]]
vod_min = [0,0,0,0]
vod_max = [0.75,2,0.2,0.5]
true_vod = [0.067, 0.82, 0.051, 0.05]

#prior_min += vod_min
#prior_max += vod_max
#true_val += true_vod
#order_index += vod_index 
#par_names += vod_names
order_index = vod_index + order_index 

true_val += true_vod
par_names += vod_names

#%%
j0 = 0
fig, ax_all = plt.subplots(6,3,figsize=(12,18),dpi=150)
for row in range(len(order_index)):
    for column in range(len(order_index[row])):
        print(row,column)
        ax = ax_all[row,column]
        j = order_index[row][column]
        # if row == 1 and column >= 1:
        #     pass
        if row <= 1:
            single_panel_plain(ax,data=allpars[:,:,j],
                         realval=true_val[j],
                         title=par_names[j])
        else:
            single_panel_plain(ax,data=np.exp(allpars[:,:,j]),
                         realval=true_val[j],
                         title=par_names[j])
        if row == 5:
            ax.set_xticks(range(1,5),obs_names,fontsize=12,rotation=320)
        #j0 += 1

ax_all[0,1].axis("off")
ax_all[0,2].axis("off")

last_ax = ax_all[-1,-1]
#last_ax.axis("off")
for coli in range(4):
    last_ax.plot([],[],color=colors_list[coli], label=obs_names[coli],linewidth=4,alpha=0.75)
lg = fig.legend(loc="center",bbox_to_anchor=(0.67, 0.9),ncols=2,fontsize=24)
lg.set_title(title="Observation scenario", prop = {'size':24})

plt.tight_layout(w_pad = 2, h_pad=2)
#%%

