import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));


#prior_min = [ 0.1,  1e-6, 500, 0.75, 0.75,0.01];
#prior_max = [ 100, 2e-5, 3000, 10, 8,10];


#FT(31),FT(0.5), FT(2),FT(0.4e-5), FT(800),FT(5),FT(2),FT(1),FT((16-0.0*9.3)*sqrt(1000)));
prior_min = [10, 0.01, 0.1,  1e-7, 500, 0.75, 0.75, 0.1,100];
prior_max = [120,0.9,  50, 2e-5, 3000, 10, 8, 10,1000];

true_val = [31,0.5, 2, 0.4e-5,800, 5, 2, 1, 506];
par_names = ["Vcmax", "Stomatal margin", "Kmax_plant", "Kmax_soil", "Soil depth", "Kplant location", "Kplant shape","Volume factor","Medlyn g1"];

plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["font.size"] = 16;

#for each run in FixET, get the parameters
out_folder = "./";

subdir_list = ["oAll/", "o1AMPM/", "o6AMPM/", "o1AM/"]

par_post = []
for x in subdir_list:
    g1 = np.array(pd.read_csv(out_folder+x+"post_par.csv"));
    par_post.append(g1)

par_post = np.array(par_post) #4 x 5000 x npar

subdir_list2 = ["oAll_c2/", "o1AMPM_c2/", "o6AMPM_c2/", "o1AM_c2/"]

par_post2 = []
for x in subdir_list2:
    g1 = np.array(pd.read_csv(out_folder+x+"post_par.csv"));
    par_post2.append(g1)

par_post2 = np.array(par_post2)



obs_names = ["All", "1 AM/PM", "6 AM/PM", "1 AM"]

j = 0

colors_list = ["tab:blue","tab:green","tab:orange","tab:red"]

fig, ax_all = plt.subplots(3,3,figsize=(12,8))
for ax in ax_all.ravel():
    for iobs in range(4):
        ax.plot(np.exp(par_post[iobs,:,j]),label=obs_names[iobs],alpha=0.5,color=colors_list[iobs])
        ax.plot(np.exp(par_post2[iobs,:,j]),alpha=0.5,color=colors_list[iobs])
    if j==0:
        ax.legend()
	
    ax.set_ylim((prior_min[j]),(prior_max[j]))
    ax.set_title(par_names[j])
    ax.set_xticks([])
    ax.hlines(true_val[j], 0, len(g1[:,j]), color="black")
    j += 1
#end
#ax_all[1,3].axis("off")

plt.tight_layout()
#ax = ax_all[1,3]
#ax.plot(np.sqrt(g1[:,-1]))
#ax.set_title("ET RMSE")

plt.savefig("chain_coup.png")

par_post1b = np.exp(par_post[:,-1000::10,:])
par_post2b = np.exp(par_post2[:,-1000::10,:])


j = 0
fig, ax_all = plt.subplots(3,3,figsize=(12,8))
for ax in ax_all.ravel():
    for iobs in range(4):
        toplot = np.array(list(par_post1b[iobs,:,j]) + list(par_post2b[iobs,:,j]))
        ax.hist(toplot,label=obs_names[iobs], alpha=0.5,color=colors_list[iobs])
    if j==0:
        ax.legend()

    #ax.set_ylim((prior_min[j]),(prior_max[j]))
    ax.set_title(par_names[j])
    ax.set_yticks([])
    ax.vlines(true_val[j], 0, ax.get_ylim()[1], color="black")
    j += 1
plt.tight_layout()
plt.savefig("hist_coup.png")

