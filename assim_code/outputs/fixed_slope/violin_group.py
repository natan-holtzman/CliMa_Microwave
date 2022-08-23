import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));


#prior_min = [ 0.1,  1e-6, 500, 0.75, 0.75,0.01];
#prior_max = [ 100, 2e-5, 3000, 10, 8,10];


#FT(31),FT(0.5), FT(2),FT(0.4e-5), FT(800),FT(5),FT(2),FT(1),FT((16-0.0*9.3)*sqrt(1000)));

#prior_min = [10, 0.01, 0.1, 500, 0.75,  0.01,100];
#prior_max = [120,0.9,  50, 3000, 10, 10,1000];
#sim_res1 = run_sim_2(FT(31),FT(0.5), FT(16), FT(2000),FT(5),FT(0.33),FT(506));

#true_val = [31,0.5, 16, 2000, 5, 0.33, 506];

#par_names = ["Vcmax $(\mu mol/m^2/s)$", "Stomatal P63 (-MPa)", "Max. xylem cond-\nuctance $(mol/m^2/s/MPa)$",  "Soil depth (mm)", "Xylem P63 (-MPa)",
# "Max. water storage\nvolume $(kg/m^2)$","Medlyn g1 $(Pa^{1/2})$","Soil conductivity $(\mu m/s)$","Soil boundary\ncondition","Root distrib. ($m^{-1}$)",
#"PLC steepness","Van Genuchen n"];

par_names = ["Vcmax $(\mu mol/m^2/s)$", "*Stomatal P63 (-MPa)", "*Max. xylem\nconductance $(mol/m^2/s/MPa)$",  "Soil depth (mm)", "*Xylem P63 (-MPa)",
 "*Max. water storage\nvolume $(kg/m^2)$","Ball-Berry g1","*Soil conductivity $(\mu m/s)$","Soil boundary condition","Root profile ($m^{-1}$)",
"PLC steepness","Van Genuchten n","Water potential\nof dry soil (MPa)"];


#order_index = [2,5,4,10,3,0,6,1,7,8,9,11,12]

order_index = [[2,5,4,10],[0,6,1],[3,8,9],[11,12,7]]


plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["font.size"] = 12;
plt.rcParams["mathtext.default"] = "regular"
#plt.rcParams['text.usetex'] = True;
#for each run in FixET, get the parameters
out_folder = "./";

subdir_list = ["oAll_c3/", "o1AMPM_c3/", "o6AMPM_c3/", "o1and6_c3/", "o16offset_c3/"]

par_post = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"post_par.csv"))[-4000::100,:13];
    #if x == "o1AMPM_c3":
    #    g1 *= np.nan
    par_post.append(g1)

par_post = np.array(par_post) #4 x 5000 x npar

subdir_list = ["oAll_c2/", "o1AMPM_c2/", "o6AMPM_c2/", "o1and6_c2/", "o16offset_c2/"]

par_post2 = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"post_par.csv"))[-4000::100,:13];
    par_post2.append(g1)

par_post2 = np.array(par_post2)

subdir_list = ["oAll_c1/", "o1AMPM_c1/", "o6AMPM_c1/", "o1and6_c1/", "o16offset_c1/"]

par_post3 = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"post_par.csv"))[-4000::100,:13];
    par_post3.append(g1)

par_post3 = np.array(par_post3)
print(par_post.shape)
print(par_post2.shape)
print(par_post3.shape)
#allpars = np.concatenate((par_post[:,6000::100,:12], par_post2[:,6000::100,:12], par_post3[:,6000::100,:12]), axis=1)
allpars = np.concatenate((par_post,par_post2,par_post3),axis=1)

#allpars = allpars[np.array([0,2,4]),:,:]


print(allpars.shape)
par_list  = np.mean(allpars,axis=1)
#np.savetxt("meanpars_raw.csv",par_list,delimiter=",")


#allpars[:,:,1] =  np.log((1-np.exp(allpars[:,:,1])) *np.exp(allpars[:,:,4]))
allpars[:,:,5] += np.log(12)
allpars[:,:,7] += np.log(1e6) 
allpars[:,:,4] = np.log(np.exp(allpars[:,:,1]) + np.exp(allpars[:,:,4]))
allpars[:,:,12] *= -1 #it's in terms of log mpa^-1, need to convert to log mpa

#prior_min = [10, 0.75, 0.1, 500, 0.01+0.75,  0.01*12,1,0.1,0.3,0.01,1,1.1];
#prior_max = [120,15,  50, 3000, 10+15, 10*12,120, 20,0.5,5,8,1.9];

prior_min = [10, 0.75, 0.1, 500, 0.01+0.75,  0.01*12,1,0.5e-1,0.3,0.01,1,1.1,0.125];
prior_max = [300,15,  50, 3000, 10+15, 10*12,120,    20,0.5,5,12,1.9,8];


true_val = [90, 3, 10, 2000, 4, 1.0*12, 300,0.4,
         0.4,2,4,1.5,1]


obs_names = ["Full diurnal", "1 AM/PM", "6 AM/PM", "1+6 sync.", "1+6 offset"]
#obs_names = np.array(obs_names)[np.array([0,2,4])]

j = 0

colors_list = ["tab:blue","tab:green","tab:orange","tab:red","tab:purple"]
#colors_list = np.array(colors_list)[np.array([0,2,4])]

j0 = 0
fig, ax_all = plt.subplots(4,4,figsize=(12,8))
for i in range(4):
    for i2 in range(len(order_index[i])):    
        j = order_index[i][i2]
        ax = ax_all[i,i2]
    #toplot = np.transpose(par_post[:,5500::90,j])
    #toplot2 = np.transpose(par_post2[:,5500::90,j])
    #toplot3 = np.transpose(par_post3[:,5500::90,j])
#toplot_both = np.conVcatenate((toplot,toplot2,toplot3),axis=0)
        vpi = ax.violinplot(np.transpose(allpars[:,:,j]),showmedians=True)
        for element in vpi.keys():
            if element=="bodies":
                for patch, color in zip(vpi[element], colors_list):
                    patch.set_color(color)
            else:
                vpi[element].set_color(colors_list)
        ax.set_ylim(np.log(prior_min[j]), np.log(prior_max[j]))
        ax.set_xticks([])
        if j != 6:
            ax.hlines(np.log(true_val[j]), 0, 6, color="black")
            ticklocs = np.array([prior_min[j],true_val[j],prior_max[j]])
        else:
            ticklocs = np.array([prior_min[j],prior_max[j]])
        ax.set_xlim(0.25,5.75) 
        ax.set_yticks(np.log(ticklocs))
        ax.set_yticklabels( [str(x) for x in np.round(ticklocs,3)])
        ax.set_title(par_names[j],fontsize=12)
    #ax.hlines(np.log(true_val[j]), 0, 4, color="black")


for i in range(1,4):
    last_ax = ax_all[i,-1]
    last_ax.set_xticks([])
    last_ax.set_yticks([])
    last_ax.axis("off")
last_ax = ax_all[1,-1]
for coli in range(5):
    last_ax.plot([],[],color=colors_list[coli], label=obs_names[coli],linewidth=4,alpha=0.75)
last_ax.legend(loc="center",borderpad=2,frameon=False)

#ax_all[1,3].axis("off")

plt.tight_layout()
plt.savefig("violin_june27c.png")

#plt.figure()
#errs = np.sqrt(np.concatenate((par_post[:,5500::90,-4], par_post2[:,5500::90,-4], par_post3[:,5500::90,-4]), axis=1))
#plt.boxplot(np.transpose(errs))
#plt.savefig("errs_est.png")

