#include("../get_home_dir.jl")

using DataFrames
using CSV
using PyPlot
using Random
using StatsBase

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["lines.linewidth"] = 2;
rcParams["font.size"] = 12;
rcParams["mathtext.default"] = "regular";

pygui(true)

#=
#all
leafRMSE =  [0.022440029317634377, 0.2233237350733526, 0.16241515258946543]
etRMSE = [3.499848703698141e-05, 1.8282748027309184e-05, 0.00012361236818981321]*18.02/1000*60*60*24;
rzRMSE = [0.002419874171990943, 0.004836427076959402, 0.003798081303897795];

#old run
leafRMSE = [0.04030680922286563, 0.11144269541636687, 0.19222915929366224];
etRMSE = [5.11686931714865e-05, 3.992463510637098e-05, 5.691494713989705e-05]*18.02/1000*60*60*24;;
rzRMSE = [0.0012556391476988972, 0.003711247616340684, 0.0016033746735104248];



#summer only
leafRMSE = [0.027500784540943245, 0.31042667792524514, 0.2247603748223069];
etRMSE = [5.685735238025922e-05, 2.9178249797272345e-05, 0.00020069684599956946]*18.02/1000*60*60*24;
rzRMSE = [0.008622295202592903, 0.011245728719432036, 0.008830522894998668];

#smc is below 25th percentile
leafRMSE = [0.039628926693386926, 0.36747422274898844, 0.3023460875451147];
etRMSE = [5.117035264574388e-05, 2.6811949415476925e-05, 0.00010842262949138796]*18.02/1000*60*60*24;
ezRMSE = [0.006981421698789021, 0.009174284328462524, 0.008990465281692152];



#for run with 1 yr of obs 
leafRMSE = [0.3829476523393593, 0.46002236065688096, 0.4673672121505479];
etRMSE = [0.00013979589607411575, 0.00010461094524003548, 0.00019008095002526475]*18.02/1000*60*60*24;
rzRMSE = [0.015044365492075336, 0.014902333663947649, 0.01330902392320819];
gppRMSE = [0,0,0];



#of mean model
leafRMSE = [0.08945551524372355, 0.14251331000460887, 0.13241443754414872];
etRMSE = [1.5222648547911854e-05, 0.00013099436536534266, 9.21809614358626e-05]*18.02/1000*60*60*24;
rzRMSE = [0.003997532112235683, 0.008655108081845725, 0.005809986996467336];
gppRMSE = [0.21847112634602872, 0.4596499061391961, 0.6818102285952408];
=#
#of all posterior models
leafRMSE = [0.26133356838960564, 0.3074713890611421, 0.2918077098629171];
etRMSE = [9.842854882620204e-05, 0.00016921201488709768, 0.0001639950079942043]*18.02/1000*60*60*24;
rzRMSE = [0.011871646096340541, 0.014053328027630875, 0.012019206607290491];
gppRMSE = [0.9254825193504083, 0.9754004090810706, 1.1506586915475916];

#=
leafClim = 0.49159888520831974
etClim = 0.0005536780748759148*18.02/1000*60*60*24;
GPPClim = 1.108767203419914
SMCClim = 0.033217017084657637
=#

leafClim = 0.29
etClim = 0.0005485*18.02/1000*60*60*24;
GPPClim = 1.02
SMCClim = 0.02039


#hourly
#=
leafRMSE = [0.2700912852977157, 0.32065105501956126, 0.31745946727017765];
etRMSE = [0.00021427381041728977, 0.0003089649531727586, 0.0003283132486278844]*18.02/1000*60*60*24;
rzRMSE = [0.011871646096340541, 0.014053328027630875, 0.012019206607290491];
gppRMSE = [1.6918560659021655, 1.72694208229437, 2.360289397605617];
=#


#=
#of all posterior models in summer
leafRMSE = [0.35566616803265233, 0.41697165529764213, 0.40700052940551557];
etRMSE = [0.00015233797942265678, 0.0002756494984164114, 0.00026221773927576563]*18.02/1000*60*60*24;
rzRMSE = [0.012168110075194902, 0.01428246340606223, 0.013860139883244405];
gppRMSE = [1.1143159302378243, 1.228936877836732, 1.15036177618971];

#all models in summer of drought year
leafRMSE = [0.6231152544689271, 0.7787730785339397, 0.76263162013701];
etRMSE = [0.0001712816248330619, 0.0005417779507410804, 0.0003061501035671869]*18.02/1000*60*60*24;
rzRMSE = [0.009354288730134051, 0.01573702108409179, 0.011127186878137709];
gppRMSE = [0.7958373144571512, 1.1167087093931025, 0.8835488547589718];
=#

#with new soil drainage, only drought
leafRMSE = [0.3708963480381229, 0.7747645919304228, 1.3513359411628316];
etRMSE = [0.0006387851877550574, 0.000489038870705169, 0.000609516078627799]*18.02/1000*60*60*24;
rzRMSE =[0.052932234669433524, 0.03712856698811909, 0.06859739132784316];
gppRMSE = [1.7483021927417743, 2.2840484995809023, 1.9240131724189131];

#all years
leafRMSE = [0.11883196748196301, 0.2617287555126143, 0.5319219121716251];
etRMSE = [0.0003087018978885237, 0.00038174022588142445, 0.00039408562820244565]*18.02/1000*60*60*24;
rzRMSE =[0.022971679828772417, 0.020572573892492677, 0.03924394442779017];
gppRMSE = [1.3384870617159015, 2.1920721225574464, 1.8773748923191838];

mylabels = ["Full\ndiurnal", "1AM/PM", "6AM/PM"];


figure()
subplot(141)
bar([1,2,3],leafRMSE)
xticks([1,2,3], mylabels)
title("Leaf water\npotential (MPa)")
#plot([0.5,3.5],[1,1]*leafClim,"k--")

subplot(142)
bar([1,2,3],rzRMSE)
xticks([1,2,3], mylabels)
title("Column\nsoil moisture")
#plot([0.5,3.5],[1,1]*SMCClim,"k--")


subplot(143)
bar([1,2,3],etRMSE)
xticks([1,2,3], mylabels)
title(string("Evapotranpiration\n", L"$(mm/day)$"))
ticklabel_format(style="plain",axis="y")
#plot([0.5,3.5],[1,1]*etClim,"k--")


subplot(144)
bar([1,2,3],gppRMSE)
xticks([1,2,3], mylabels)
title(string("Gross primary\nproduction ", L"$(\mu mol/m^2/s)$"))
ticklabel_format(style="plain",axis="y")
#plot([0.5,3.5],[1,1]*GPPClim,"k--")


#tight_layout()

#=
gsRMSE = [0.0007360393180670877, 0.0006048789163496154, 0.001909214673064903]
leafscaleRMSE = [0.010286256839412344, 0.01665548366779598, 0.09018817219779075]
betaRMSE = [0.010342191353158734, 0.012702320594122641, 0.060695502174849196]


figure()
subplot(231)
bar([1,2,3],leafRMSE)
xticks([1,2,3], mylabels)
title("Leaf water potential (MPa)")

subplot(232)
bar([1,2,3],leafscaleRMSE)
xticks([1,2,3], mylabels)
title("LWP divided by P63")
ticklabel_format(style="plain",axis="y")

subplot(233)
bar([1,2,3],betaRMSE)
xticks([1,2,3], mylabels)
title("Beta factor")
ticklabel_format(style="plain",axis="y")

subplot(234)
bar([1,2,3],gsRMSE)
xticks([1,2,3], mylabels)
title(L"Stomatal conductance $(mol/m^2/s)$")

subplot(235)
bar([1,2,3],etRMSE)
xticks([1,2,3], mylabels)
title(L"Daily ET $(mm/m^2/s)$")
ticklabel_format(style="plain",axis="y")


leafmean = convert(Array, CSV.read("leaf_mean_tab2.csv",DataFrame,header=false));
leafscale = convert(Array, CSV.read("leaf_scale_tab2.csv",DataFrame,header=false));

figure()
for i in 1:3
	plot(leafmean[i,:], 1 ./ leafscale[i,:],"o",label=mylabels[i])
end
plot(leafmean[1,end], 1 ./ leafscale[1,end],"k+",label="True model",markersize=25,markeredgewidth=3)
xlabel("Mean LWP (MPa)")
ylabel("Stomatal P63 (MPa)")
legend()
=#

colors_list = ["tab:blue","tab:green","tab:orange","tab:red"];
label_list = ["All obs.", "1 AM/PM obs.", "6 AM/PM obs."];


#=
leafmeans = convert(Array, CSV.read("leaf_means.csv",DataFrame,header=0));
ETmeans = convert(Array, CSV.read("et_means.csv",DataFrame,header=0))*18.02/1000*60*60*24;
ETmeans[ETmeans .< 0] .= 0;

leafmean_clim = mean(reshape(leafmeans, (4,365,12)), dims=3)[:,:,1];
leafmean_anom = leafmeans - repeat(leafmean_clim, outer=(1,12));

ETmean_clim = mean(reshape(ETmeans, (4,365,12)), dims=3)[:,:,1];
ETmean_anom = ETmeans - repeat(ETmean_clim, outer=(1,12));
=#

#=
timeX = (1:(365*12))/365;

figure()
for i in 1:3
	plot(timeX, leafmeans[i,:],color=colors_list[i],linewidth=1, label=label_list[i])
end
plot(timeX, leafmeans[4,:],color="k",linewidth=1, label="True")
legend()

figure()
for i in 1:3
	plot(timeX, ETmeans[i,:]*18.02/1000*60*60*24,color=colors_list[i],linewidth=1, label=label_list[i])
end
plot(timeX, ETmeans[4,:]*18.02/1000*60*60*24,color="k",linewidth=1, label="True")
legend()
=#

#=
figure()
lmf = convert(Array,CSV.read("leaf_means_full.csv",DataFrame, header=0));
for i in 1:3
	plot(get_diurnal(lmf[i,:],24), color=colors_list[i], label=label_list[i])
end
plot(get_diurnal(lmf[4,:],24), color="k", label="True model")
legend()	
=#

