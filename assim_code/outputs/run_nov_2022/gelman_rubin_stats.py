import numpy as np
import pandas as pd


print("ALL YEARS")

out_folder = "./";


par_names = ["Vcmax", "Stomatal margin", "Kmax_plant", "Kmax_soil", "Soil depth", "Kplant location", "Kplant shape","Volume factor","Medlyn g1"];

#for each run in FixET, get the parameters
out_folder = "./";

subdir_list = ["oAll_c2/", "o1AMPM_c2/","o6AMPM_c2/", "o1and6_c2/", "o16offset_c2/"]
subdir_list2 = ["oAll_c3/", "o1AMPM_c3/", "o6AMPM_c3/", "o1and6_c3/", "o16offset_c3/"]
subdir_list3 = ["oAll_c1/", "o1AMPM_c1/", "o6AMPM_c1/", "o1and6_c1/","o16offset_c1/"]


leafpars = []

NPAR = 12
for i in range(5):
    p0 = np.array(pd.read_csv(out_folder+subdir_list[i]+"post_par.csv"))[-4000:,:NPAR]
    p1 = np.array(pd.read_csv(out_folder+subdir_list2[i]+"post_par.csv"))[-4000:,:NPAR]
    p2 = np.array(pd.read_csv(out_folder+subdir_list3[i]+"post_par.csv"))[-4000:,:NPAR]
    p3 = [p0,p1,p2]
    leafpars.append(p3)

#leafpars is a list of 3 items
#each item corresponds to an observation scenerio
#within each item, there is a list of 3 chains
#within each chain, there is a Lx7 array, with dims samples x parameters 

def getGM(x):
    z = np.array(x)
    J,L,nP = z.shape
    chain_mean = np.mean(z,axis=1) #dims J x nP
    grand_mean = np.mean(z,axis=(0,1))
    B = L*np.var(chain_mean,axis=0,ddof=1)
    s_squared = np.var(z,axis=1,ddof=1)
    W = np.mean(s_squared,axis=0)
    R = ((L-1)/L*W + B/L) / W
    return R
print("Gelman Rubin statistic")
for i in range(len(leafpars)):
    print(getGM(leafpars[i]))


