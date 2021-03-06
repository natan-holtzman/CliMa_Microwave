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
 "Max. water storage\nvolume $(kg/m^2)$","Medlyn g1 $(Pa^{1/2})$"];

order_index = [2,5,4,3,0,6,1]

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

allpars = np.concatenate((par_post[:,5500::90,:7], par_post2[:,5500::90,:7], par_post3[:,5500::90,:7]), axis=1)
print(allpars.shape)
par_list  = np.mean(allpars,axis=1)
np.savetxt("meanpars_raw.csv",par_list,delimiter=",")

allparsV = np.concatenate((par_post[:,5500::90,-3:], par_post2[:,5500::90,-3:], par_post3[:,5500::90,-3:]), axis=1)
par_listV  = np.mean(allparsV,axis=1)
np.savetxt("Vmeanpars_raw.csv",par_listV,delimiter=",")


allpars[:,:,1] =  np.log((1-np.exp(allpars[:,:,1])) *np.exp(allpars[:,:,4]))
allpars[:,:,5] += np.log(12)
 

#par_post[:,:,1] = np.log((1-np.exp(par_post[:,:,1])) *np.exp( par_post[:,:,4]))
#par_post2[:,:,1] = np.log((1-np.exp(par_post2[:,:,1])) *np.exp( par_post2[:,:,4]))
#par_post3[:,:,1] = np.log((1-np.exp(par_post3[:,:,1])) *np.exp( par_post3[:,:,4]))

#par_post[:,:,5] += np.log(12)
#par_post2[:,:,5] += np.log(12)
#par_post3[:,:,5] += np.log(12)

prior_min = [10, 0.75*0.01, 0.1, 500, 0.75,  0.01*12,100];
prior_max = [120,10*0.9,  50, 3000, 10, 10*12,1000];
true_val = [31,0.5*5, 16, 2000, 5, 0.33*12, 506];

#print(par_post.shape)
#print(par_post2.shape)

obs_names = ["All", "1 AM/PM", "6 AM/PM", "1 AM"]

j = 0

colors_list = ["tab:blue","tab:green","tab:orange","tab:red"]


j0 = 0
fig, ax_all = plt.subplots(2,4,figsize=(12,6))
for ax in ax_all.ravel()[:7]:
    j = order_index[j0]
    #toplot = np.transpose(par_post[:,5500::90,j])
    #toplot2 = np.transpose(par_post2[:,5500::90,j])
    #toplot3 = np.transpose(par_post3[:,5500::90,j])
    #toplot_both = np.concatenate((toplot,toplot2,toplot3),axis=0)
    ax.boxplot(np.transpose(allpars[:,:,j]))
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
#ax_all[1,3].axis("off")

plt.tight_layout()
plt.savefig("box_pars_logscale_fullchainC.png")

plt.figure()
errs = np.sqrt(np.concatenate((par_post[:,5500::90,-4], par_post2[:,5500::90,-4], par_post3[:,5500::90,-4]), axis=1))
plt.boxplot(np.transpose(errs))
plt.savefig("errs_est.png")

