import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));
plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["font.size"] = 22;

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

subdir_list = ["oAll_c2/", "o1AMPM_c2/","o6AMPM_c2/"]
subdir_list2 = ["oAll_c3/", "o1AMPM_c3/", "o6AMPM_c3/"]
subdir_list3 = ["oAll_c1/", "o1AMPM_c1/", "o6AMPM_c1/"]





obs_names = ["All", "1 AM/PM", "6 AM/PM"]

d24 = np.arange(24*365*12) % (24*365)
#summer_24 = (d24 >= 150)*(d24 < 275)
summer_24 = d24 >= 0

d1 = np.arange(365*12) % 365 
#summer_1 = (d1 >= 150)*(d1 < 275)
summer_1 = d1 >= 0

y1 = np.floor(np.arange(365*12) / 365)
#summer_1 = ((y1 < 6)+(y1 > 9)) * summer_1

y24 = np.floor(np.arange(24*365*12)/(24*365))
#summer_24 = ((y24<6)+(y24>9)) * summer_24


#leaf_post = np.array(leaf_post)

#out_names = ["Leaf pre-dawn","Leaf diurnal","RZSM"];

#leaftab = np.array(pd.read_csv("postLeaf.csv"));
#branchtab = np.array(pd.read_csv(out_folder+"postBranch.csv"));
#trunktab =  np.array(pd.read_csv(out_folder+"postTrunk.csv"));
#RZSMtab =  np.array(pd.read_csv(out_folder+"postRZ.csv"));
#ETtab =  np.array(pd.read_csv(out_folder+"postET.csv"));

def getcor(x,y):
    return np.corrcoef(x,y)[0,1]

def getRMSE(x,y):
    return np.sqrt(np.mean((x-y)**2))

def getbias(pred,obs):
    return np.mean(pred-obs)/np.mean(obs)*100

def getR2(pred,obs):
    return 1 - np.mean((pred-obs)**2)/np.var(obs)


def wrapfunc(func,tab):
    truthval = tab[:,-1]
    posterior = tab[:,:-1]
    allstats = [func(posterior[:,i],truthval) for i in range(posterior.shape[1])]
    return np.nanmedian(allstats)

def wrapR2(tab):
    truthval = tab[:,-1]
    posterior = tab[:,:-1]
    allstats = [getR2(posterior[:,i],truthval) for i in range(posterior.shape[1])]
    return np.array(allstats)

def meanR2(tab):
    truthval = tab[:,-1]
    postmean = np.nanmean(tab[:,:-1],1)
    return getRMSE(postmean, truthval)

def get_diurnal_2d(x,nstep):
    y = np.reshape(x, (-1, nstep, x.shape[-1]))
    return np.mean(y, axis=0)	

def get_diurnal_1d(x,nstep):
    y = np.reshape(x, (-1, nstep))
    return np.mean(y, axis=0)
def get_daily_2d(x,nstep):
    y = np.reshape(x, (-1, nstep, x.shape[-1]))
    return np.mean(y, axis=1)
#outvar = [leaftab[1::24,:],get_diurnal(leaftab,24),
#	  RZSMtab[:,:] ]

def make_stats_tab(tab_list):
    ans = np.zeros((4,4))
    ans[0,:] = [wrapfunc(getcor,x) for x in tab_list]
    ans[1,:] = [wrapfunc(getbias,x) for x in tab_list]
    ans[2,:] = [wrapfunc(getRMSE,x) for x in tab_list]
    ans[3,:] = [wrapfunc(getR2,x) for x in tab_list]
    return ans


def wrapRMSE(tab):
    truthval = tab[:,-1]
    posterior = tab[:,:-1]
    allstats = [getRMSE(posterior[:,i],truthval) for i in range(posterior.shape[1])]
    return np.array(allstats)


obs_names = ["All", "1 AM/PM", "6 AM/PM", "1 AM"]

j = 0

colors_list = ["tab:blue","tab:green","tab:orange","tab:red"]

leafpars = []
leafpost = []
for i in range(3):
    g0 = np.array(pd.read_csv(out_folder+subdir_list[i]+"postLeaf.csv"))[:,:-1];
    print(g0.shape)
    g1 = np.array(pd.read_csv(out_folder+subdir_list2[i]+"postLeaf.csv"))[:,:-1];
    print(g1.shape)
    g2 = np.array(pd.read_csv(out_folder+subdir_list3[i]+"postLeaf.csv"));
    print(g2.shape)
    g3 = np.concatenate((g0,g1,g2),axis=1)
    print(g3.shape)
    g3[:,np.mean(g3==0,0)==1] = np.nan
    
    #leafpost.append(g3[5::24,:])
    leafpost.append(get_daily_2d(g3,24)[summer_1,:])


    p0 = np.array(pd.read_csv(out_folder+subdir_list[i]+"post_par.csv"))[5500::90,:7]
    p1 = np.array(pd.read_csv(out_folder+subdir_list2[i]+"post_par.csv"))[5500::90,:7]
    p2 = np.array(pd.read_csv(out_folder+subdir_list3[i]+"post_par.csv"))[5500::90,:7]
    p3 = np.concatenate((p0,p1,p2),axis=0)
    leafpars.append(p3[:,:])

print("Leaf")
leaftab = [meanR2(x) for x in leafpost]
print(leaftab)

alpha = 50
nsoil = 1.2
msoil = 1-1/nsoil

def eqsmc(psi):
    return ((-psi * alpha) ** nsoil + 1) ** (-msoil) * (0.45 - 0.067) + 0.067;


leafpost = []
for i in range(3):
    g0 = np.array(pd.read_csv(out_folder+subdir_list[i]+"postLeaf.csv"))[:,:-1];
    g1 = np.array(pd.read_csv(out_folder+subdir_list2[i]+"postLeaf.csv"))[:,:-1];
    g2 = np.array(pd.read_csv(out_folder+subdir_list3[i]+"postLeaf.csv"));
    g3 = np.concatenate((g0,g1,g2),axis=1)
    g3[:,np.mean(g3==0,0)==1] = np.nan
    g3[g3==0] = np.nan
    leafpost.append(eqsmc(g3[5::24,:]))
    #leafpost.append(get_daily_2d(g3,24)[summer_1,:])


#print("Backtrack RZ")
#leaftab = [meanR2(x) for x in leafpost]
#print(leaftab)














#print(leaftab)

#fig,ax = plt.subplots(figsize=(12,8))
#plt.subplot(2,2,1

#plt.figure(figsize=(12,8))
#plt.figure()
#plt.boxplot(leaftab)
#plt.xticks([1,2,3],obs_names)
#plt.title("RMSE of leaf water potential (MPa)")
#plt.savefig("leaf_RMSE.png")

leafpost = []
for i in range(3):
    g0 = np.array(pd.read_csv(out_folder+subdir_list[i]+"postET.csv"))[:,:-1];
    g1 = np.array(pd.read_csv(out_folder+subdir_list2[i]+"postET.csv"))[:,:-1];
    g2 = np.array(pd.read_csv(out_folder+subdir_list3[i]+"postET.csv"));
    g3 = np.concatenate((g0,g1,g2),axis=1)
    print(g3.shape)
    g3[:,np.mean(g3==0,0)==1] = np.nan

    leafpost.append(get_daily_2d(g3,24)[summer_1,:])

print("ET")
ETtab = [meanR2(x) for x in leafpost]
print(ETtab)

#print(trunktab)
#plt.figure()
#plt.boxplot(ETtab)
#plt.xticks([1,2,3],obs_names)
#plt.savefig("ET_errs.png")
#plt.ylim(-0.1,1.1)
#plt.title("ET")
#plt.savefig("ET_RMSE.png")

leafpost = []
for i in range(3):
    g0 = np.array(pd.read_csv(out_folder+subdir_list[i]+"postRZ.csv"))[:,:-1];
    g1 = np.array(pd.read_csv(out_folder+subdir_list2[i]+"postRZ.csv"))[:,:-1];
    g2 = np.array(pd.read_csv(out_folder+subdir_list3[i]+"postRZ.csv"));
    g3 = np.concatenate((g0,g1,g2),axis=1)
    g3[:,np.mean(g3==0,0)==1] = np.nan

    leafpost.append(g3[summer_1,:])


print("Root zone soil")
roottab = [meanR2(x) for x in leafpost]
print(roottab)



leafpost = []
for i in range(3):
    g0 = np.array(pd.read_csv(out_folder+subdir_list[i]+"postSWS.csv"))[:,:-1];
    g1 = np.array(pd.read_csv(out_folder+subdir_list2[i]+"postSWS.csv"))[:,:-1];
    g2 = np.array(pd.read_csv(out_folder+subdir_list3[i]+"postSWS.csv"))[:,:-1];
    g3 = np.concatenate((g0,g1,g2),axis=1)
    g3[:,np.mean(g3==0,0)==1] = np.nan
    g3 *= 2
    g3 /= (np.exp(leafpars[i][:,3])/1000)
    g4  =  np.array(pd.read_csv(out_folder+subdir_list3[i]+"postSWS.csv"))[:,-1].reshape((-1,1))
    g5 =  np.concatenate((g3, g4),axis=1)
    leafpost.append(g5[summer_1,:])


print("Column soil")
roottab = [meanR2(x) for x in leafpost]
print(roottab)

leafpost = []
for i in range(3):
    g0 = np.array(pd.read_csv(out_folder+subdir_list[i]+"postSWP.csv"))[:,:-1];
    g1 = np.array(pd.read_csv(out_folder+subdir_list2[i]+"postSWP.csv"))[:,:-1];
    g2 = np.array(pd.read_csv(out_folder+subdir_list3[i]+"postSWP.csv"))[:,:-1];
    g3 = np.concatenate((g0,g1,g2),axis=1)
    g3[:,np.mean(g3==0,0)==1] = np.nan
    g3 *= 2
    g3 /= (np.exp(leafpars[i][:,3])/1000)
    g4  =  np.array(pd.read_csv(out_folder+subdir_list3[i]+"postSWP.csv"))[:,-1].reshape((-1,1))
    g5 =  np.concatenate((g3, g4),axis=1)
    leafpost.append(g5[summer_1,:])



print("Column soil potential")
roottab = [meanR2(x) for x in leafpost]
print(roottab)



#plt.figure()
#plt.boxplot(roottab)
#plt.xticks([1,2,3],obs_names)
#plt.ylim(-0.1,1.1)
#plt.title("Root zone SMC")
#plt.savefig("RZ_RMSE.png")

#bigtab = pd.DataFrame(np.concatenate([leaftab,leaftab_predawn ,trunktab,trunktab_predawn, ETtab, leaftab4]))
#bigtab.to_csv("post_stats_nov28summer.csv")


