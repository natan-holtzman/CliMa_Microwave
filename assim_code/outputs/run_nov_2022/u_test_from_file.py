import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats

#sim_res1 = run_sim_2(FT(22),FT(0.33), FT(2),FT(1e-5), FT(600),FT(4.0),FT(3),FT(1),FT(600));
plt.rcParams["lines.linewidth"] = 1;
plt.rcParams["font.size"] = 12;
plt.rcParams["mathtext.default"] = "regular"


out_folder = "./";

all_errs = np.load("dry_rmse_nov2_hourly.npy")
print(all_errs.shape)
for i in range(4):
	print(np.nanmedian(all_errs[i][0]))
	for j in range(1,5):
		utest = scipy.stats.mannwhitneyu(all_errs[i][0], all_errs[i][j])
		print(np.nanmedian(all_errs[i][j]),", ", utest.pvalue)
	#alldiffs.append(d1)






