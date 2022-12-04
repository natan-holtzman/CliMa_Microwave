import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd


full_par_names = ["Vcmax $(\mu mol/m^2/s)$", "*Stomatal P63 (-MPa)", "*Max. xylem\nconductance $(mol/m^2/s/MPa)$",  "Soil depth (mm)", "*Xylem P63 (-MPa)",
 "*Max. water storage\nvolume $(kg/m^2)$","Ball-Berry g1","*Soil conductivity $(\mu m/s)$","Soil boundary condition","Root profile ($m^{-1}$)",
"Van Genuchten n","Water potential\nof dry soil (MPa)"];

par_names = ["Vcmax", r"$P^{\beta}_{63}$",
             "$K^{xyl}_{max}$","$Z_s$","$P^{xyl}_{63}$",
             "$V_{max}$","$g_1$","$k_s$",
             "$s_{lower}$","$d_{root}$",
             "n",r"$\psi_{dry}$"]

vod_names = ["a","b","c","$\omega$"]

#order_index = [2,5,4,10,3,0,6,1,7,8,9,11,12]

order_index = [[2,5,4],[0,6,1],[3,8,9],[10,11,7]]


plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["font.size"] = 12;
plt.rcParams["mathtext.default"] = "regular"
#plt.rcParams['text.usetex'] = True;
#for each run in FixET, get the parameters

allpars = np.load("all_param_dec3.npy")

#prior_min = [10, 0.75, 0.1, 500, 0.01+0.75,  0.01*12,1,0.1,0.3,0.01,1,1.1];
#prior_max = [120,15,  50, 3000, 10+15, 10*12,120, 20,0.5,5,8,1.9];

prior_min = [10, 0.75, 0.1, 500, 0.01+0.75,  0.01*12,1,0.5e-1,0.3,0.01,1.1,0.125];
prior_max = [300,15,  50, 3000, 10+15, 10*12,120,    20,0.5,5,1.9,8];


true_val = [90, 3, 10, 2000, 4, 1.0*12, None,0.4,
         0.4,2,1.5,1]


obs_names = ["Full diurnal", "1 AM/PM", "6 AM/PM", "1+6 sync.", "1+6 offset"]#[:3]
#obs_names = np.array(obs_names)[np.array([0,2,4])]

j = 0

colors_list = ["tab:blue","tab:green","tab:orange","tab:red","tab:purple"]#[:3]
#colors_list = np.array(colors_list)[np.array([0,2,4])]

def single_panel_log(ax,data,yrange,realval,title):
    vpi = ax.violinplot(np.transpose(data),showmedians=True)
    for element in vpi.keys():
        if element=="bodies":
            for patch, color in zip(vpi[element], colors_list):
                patch.set_color(color)
        else:
            vpi[element].set_color(colors_list)
            
    ax.set_ylim(np.log(yrange[0]), np.log(yrange[1]))
    ax.set_xticks([])
    if realval:
        ax.hlines(np.log(realval), 0, 6, color="black")
        ticklocs = np.array([yrange[0],realval,yrange[1]])
    else:
        ticklocs = np.array([yrange[0],yrange[1]])
    ax.set_xlim(0.25,5.75) 
    ax.set_yticks(np.log(ticklocs))
    ax.set_yticklabels( [str(x) for x in np.round(ticklocs,3)])
    ax.set_title(title,fontsize=12)

def single_panel_plain(ax,data,realval,title):
    vpi = ax.violinplot(np.transpose(data),showmedians=True)
    for element in vpi.keys():
        if element=="bodies":
            for patch, color in zip(vpi[element], colors_list):
                patch.set_color(color)
        else:
            vpi[element].set_color(colors_list)
    ax.hlines(realval, 0, 6, color="black")

    ax.set_xticks([])
    ax.set_xlim(0.25,5.75) 
    ax.set_title(title,fontsize=12)



vod_index = [[-4,-3,-2],[-1]]
vod_min = [0,0,0,0]
vod_max = [0.75,2,0.2,0.5]
true_vod = [0.067, 0.82, 0.051, 0.05]

#prior_min += vod_min
#prior_max += vod_max
true_val += true_vod
order_index += vod_index
par_names += vod_names

j0 = 0
fig, ax_all = plt.subplots(6,3,figsize=(10,14),dpi=200)
for row in range(len(order_index)):
    for column in range(len(order_index[row])):    
        j = order_index[row][column]
        ax = ax_all[row,column]
        if row <= 3:
            single_panel_log(ax,data=allpars[:,:,j],
                         yrange=[prior_min[j],prior_max[j]],
                         realval=true_val[j],
                         title=par_names[j])
        else:
            single_panel_plain(ax,data=allpars[:,:,j],
                         realval=true_val[j],
                         title=par_names[j])

ax_all[-1,1].axis("off")

last_ax = ax_all[-1,-1]
last_ax.axis("off")
for coli in range(5):
    last_ax.plot([],[],color=colors_list[coli], label=obs_names[coli],linewidth=4,alpha=0.75)
last_ax.legend(loc="center",borderpad=2,frameon=True,
               title="Observing scenario")

plt.tight_layout()

