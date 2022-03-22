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

par_names = ["Vcmax $(\mu mol/m^2/s)$", "Stomatal P63 (-MPa)", "Max. xylem cond-\nuctance $(mol/m^2/s/MPa)$",  "Rooting depth (mm)", "Xylem P63 (-MPa)",
 "Max. water storage\nvolume $(kg/m^2)$","Medlyn g1 $(Pa^{1/2})$","Soil conductivity $(\mu m/s)$","Soil boundary\ncondition","Root distrib. ($m^{-1}$)",
"PLC steepness"];

order_index = [2,5,4,10,3,0,6,1,7,8,9]

plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["font.size"] = 12;
plt.rcParams["mathtext.default"] = "regular"
#plt.rcParams['text.usetex'] = True;
#for each run in FixET, get the parameters
out_folder = "./";

subdir_list = ["oAll_c3/", "o1AMPM_c3/", "o6AMPM_c3/"]

par_post = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"post_par.csv"));
    par_post.append(g1)

par_post = np.array(par_post) #4 x 5000 x npar

subdir_list = ["oAll_c2/", "o1AMPM_c2/", "o6AMPM_c2/"]

par_post2 = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"post_par.csv"));
    par_post2.append(g1)

par_post2 = np.array(par_post2)

subdir_list = ["oAll_c1/", "o1AMPM_c1/", "o6AMPM_c1/"]

par_post3 = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"post_par.csv"));
    par_post3.append(g1)

par_post3 = np.array(par_post3)

allpars = np.concatenate((par_post[:,6000::100,:11], par_post2[:,6000::100,:11], par_post3[:,6000::100,:11]), axis=1)
print(allpars.shape)
par_list  = np.mean(allpars,axis=1)
#np.savetxt("meanpars_raw.csv",par_list,delimiter=",")


#allpars[:,:,1] =  np.log((1-np.exp(allpars[:,:,1])) *np.exp(allpars[:,:,4]))
allpars[:,:,5] += np.log(12)
allpars[:,:,7] += np.log(1e6) 
allpars[:,:,4] = np.log(np.exp(allpars[:,:,1]) + np.exp(allpars[:,:,4]))

#par_post[:,:,1] = np.log((1-np.exp(par_post[:,:,1])) *np.exp( par_post[:,:,4]))
#par_post2[:,:,1] = np.log((1-np.exp(par_post2[:,:,1])) *np.exp( par_post2[:,:,4]))
#par_post3[:,:,1] = np.log((1-np.exp(par_post3[:,:,1])) *np.exp( par_post3[:,:,4]))

#par_post[:,:,5] += np.log(12)
#par_post2[:,:,5] += np.log(12)
#par_post3[:,:,5] += np.log(12)

prior_min = [10, 0.75, 0.1, 500, 0.01+0.75,  0.01*12,100,0.1,0.3,0.01,1];
prior_max = [120,15,  50, 3000, 10+15, 10*12,1000,    4,0.5,5,8];


#true_val = [60, 0.5*5, 12, 2000, 5, 1.0*12, 250,5e-5, 5.0,2]
#pars0 = convert(Array{FT}, [90, 0.25, 10, 2000, 4, 1.0, 300, 0.4e-6,
#         0.4,2]);
true_val = [90, 3, 10, 2000, 4, 1.0*12, 300,0.4,
         0.4,2,4]

#print(par_post.shape)
#print(par_post2.shape)

obs_names = ["All", "1 AM/PM", "6 AM/PM", "1 AM"]

j = 0

colors_list = ["tab:blue","tab:green","tab:orange","tab:red"]


j0 = 0
fig, ax_all = plt.subplots(3,4,figsize=(12,6))
for ax in ax_all.ravel()[:11]:#[:7]:
    j = order_index[j0]
    #toplot = np.transpose(par_post[:,5500::90,j])
    #toplot2 = np.transpose(par_post2[:,5500::90,j])
    #toplot3 = np.transpose(par_post3[:,5500::90,j])
    #toplot_both = np.concatenate((toplot,toplot2,toplot3),axis=0)
    vpi = ax.violinplot(np.transpose(allpars[:,:,j]))
    for patch, color in zip(vpi['bodies'], colors_list[:3]):         
        patch.set_color(color)
    ax.set_ylim(np.log(prior_min[j]), np.log(prior_max[j]))
    #ax.set_yscale("log")
    ax.set_xticks([])
    ticklocs = np.array([prior_min[j],true_val[j],prior_max[j]])
    ax.set_yticks(np.log(ticklocs))
    ax.set_yticklabels( [str(x) for x in np.round(ticklocs,3)])
    #ax.set_ylim((prior_min[j]),(prior_max[j]))
    ax.set_title(par_names[j])
    ax.hlines(np.log(true_val[j]), 0, 4, color="black")
    j0 += 1

last_ax = ax_all.ravel()[-1]
for coli in range(3):
    last_ax.plot([],[],color=colors_list[coli], label=obs_names[coli])
last_ax.legend()
last_ax.set_xticks([])
last_ax.set_yticks([])


#ax_all[1,3].axis("off")

plt.tight_layout()
plt.savefig("violin_11parC2.png")

#plt.figure()
#errs = np.sqrt(np.concatenate((par_post[:,5500::90,-4], par_post2[:,5500::90,-4], par_post3[:,5500::90,-4]), axis=1))
#plt.boxplot(np.transpose(errs))
#plt.savefig("errs_est.png")

