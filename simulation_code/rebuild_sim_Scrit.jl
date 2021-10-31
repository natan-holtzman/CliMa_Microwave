#sim_res = run_sim(FT(22),FT(1.5),FT(5),FT(10),FT(0.025*0.4),FT(32));

using DataFrames
using CSV
using Parameters
using Land
using Land.Photosynthesis
using Land.CanopyLayers
using Land.PlantHydraulics
using Land.SoilPlantAirContinuum
using Land.StomataModels
using Thermodynamics

using CLIMAParameters
using CLIMAParameters:AbstractEarthParameterSet
using CLIMAParameters.Planet: grav



const FT = Float32

struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

const FT = Float32

#=
#import Land.PlantHydraulics.xylem_p_crit
function xylem_p_crit2(
            vc::WeibullSingle{FT},
            f_st::FT = FT(1)
)
    @unpack b,c = vc;
    return -b * log( FT(100) ) ^ (1 / c) * f_st
end
=#





KG_H_2_MOL_S = FT(55.55 / 3600);
mpa2mm = FT(10^6/9.8);

deltaT = FT(60*30)


include("create_spac_no_angles.jl")
include("land_utils4.jl")

#using Plots

K_STEFAN = FT(Stefan());

#N      = 48*7 #only run for 4 days for quick testing
#istart = 48*365*1 - 52*48 + 165*48 #+ 180*48 - 3*48
#N = 48*365
#istart = 48*365*2 - 52*48 #+ 165*48
df_raw = CSV.read("../data/moflux_land_data_newnames_7.csv", DataFrame);

function run_sim(vcmax_par, s_crit, k_plant, k_soil, b_soil, p20, z_soil, n_soil, istart, N, smc0)

df = deepcopy(df_raw[istart+1:istart+N,:]);

df[!,"leafpot"] = zeros(N)
df[!,"leafvol"] = zeros(N)
df[!,"stempot"] = zeros(N)
df[!,"soilpot"] = zeros(N)
df[!,"ETmod"] = zeros(N)
df[!,"trunkflow"] = zeros(N)

df[!,"APARsun"] = zeros(N)
df[!,"Zenith"] = zeros(N)

df[!,"APARshade"] = zeros(N)

df[!,"APAR_ps"] = zeros(N)

df[!,"Ecrit"] = zeros(N)

df[!,"leafStore"] = zeros(N)

df[!,"fvis"] = zeros(N)
df[!,"glw"] = zeros(N)
df[!,"GPP"] = zeros(N)
df[!,"NPP"] = zeros(N)
df[!,"Pcrit"] = zeros(N)
df[!,"Runoff"] = zeros(N)


#beta_curve = WeibullSingle(FT(weibB_par), FT(weibC_par));
#sm1 = Land.StomataModels.ESMBallBerry(g0_par,g1_par);
#node = create_spac(sm1,vcmax_par,k_plant,FT(40));

#psi_sat = FT(p20 / 0.2^(-1*b_soil));

rs13 = (0.13-0.05)/(n_soil-0.05);
psi_sat = FT(p20 / rs13^(-1*b_soil));
p_crit = FT(psi_sat*((s_crit - 0.05)/(n_soil-0.05))^(-1*b_soil));

node = create_spac(OSMWang{FT}(),vcmax_par,k_plant,psi_sat, b_soil, z_soil, n_soil, FT(40));

smc_record = zeros(FT,N,length(node.swc));
psi_record_leaf = zeros(FT,N,length(node.plant_hs.leaves));
psi_record_branch = zeros(FT,N,length(node.plant_hs.leaves));


rbase =  Q10TD{FT}(0, 298.15, 1.7)


@unpack angles, can_opt, can_rad, canopy_rt, envirs, f_SL, ga, in_rad,
		latitude, leaves_rt, n_canopy, photo_set, plant_hs, plant_ps,
		rt_con, soil_opt, stomata_model, wl_set = node;
@unpack lidf, nAzi, nIncl = canopy_rt;
@unpack dWL, iPAR = wl_set;
in_rad_bak = deepcopy(in_rad);
nSL = nAzi * nIncl;
in_Erad = in_rad_bak.E_direct .+ in_rad_bak.E_diffuse;
#in_PPFD = sum( e2phot(dWL, in_Erad)[iPAR] ) * FT(1e6);
in_PPFD = sum( Land.CanopyLayers.e2phot(wl_set.WL, in_Erad/1000)[iPAR] .* dWL[iPAR] ) * FT(1e6);

# iterate through the weather data
for i in eachindex(df.Day)
	
	if i == 1
	
		w_soil = FT(smc0);		
		node.swc = ones(FT, length(node.swc))*w_soil;
		# update soil water matrices
		# need to adjust SWC to avoid problem in residual SWC at Niwot Ridge
		ψ_soil = soil_p_25_swc(plant_hs.roots[1].sh, w_soil);
		for root in plant_hs.roots
			root.p_ups = ψ_soil;
		end;
	end
	
	t_soil = FT(df.T_SOIL[1] + 273.15);


#=
	if mod(i, 48) == 0
		println(Int(i/48))
		println(node.plant_hs.leaves[length(node.plant_hs.leaves)].p_leaf)
	end	
=#


	update_LAI!(node, FT(df.LAI_modis[i]));
	

	# update PAR related information
	zenith = zenith_angle(latitude, FT(df.Day[i]), FT(df.Hour[i]),
						  FT(df.Minu[i]));
	zenith = min(88, zenith);
	angles.sza = zenith;
	in_rad.E_direct  .= in_rad_bak.E_direct  .* df.PPFD_in[i] ./ in_PPFD;
	in_rad.E_diffuse .= in_rad_bak.E_diffuse .* df.PPFD_in[i] ./ in_PPFD;
	canopy_geometry!(canopy_rt, angles, can_opt, rt_con);
	canopy_matrices!(leaves_rt, can_opt);
	short_wave!(canopy_rt, can_opt, can_rad, in_rad, soil_opt, rt_con);
	canopy_fluxes!(canopy_rt, can_opt, can_rad, in_rad, soil_opt,
				   leaves_rt, wl_set, rt_con);

				   
	#update_Weibull!(node, FT(weib_b), FT(3.0));
	#update_Kmax!(node, FT(k_plant));
	
	# update fluxes
	f_H₂O = 0;
	f_CO₂ = 0;
	for i_can in 1:n_canopy
		iEN = envirs[i_can];
		iHS = plant_hs.leaves[i_can];
		iPS = plant_ps[i_can];
		iRT = n_canopy + 1 - i_can;

		# update environmental conditions
		iEN.t_air = df.T_AIR[i] + 273.15;
		iEN.p_atm = df.P_ATM[i] * 1000;
		iEN.p_a   = iEN.p_atm * 4e-4;
		iEN.p_O₂  = iEN.p_atm * 0.209;
		iEN.p_sat = saturation_vapor_pressure( param_set, iEN.t_air , Liquid());
		iEN.vpd   = df.VPD[i] * 100;
		iEN.p_H₂O = iEN.p_sat - iEN.vpd;
		iEN.RH    = iEN.p_H₂O / iEN.p_sat;
		iEN.wind  = df.WIND[i];

		# prescribe leaf temperature
		_tl = (df.LW_OUT[i] / 0.97 / K_STEFAN ) ^ 0.25;
		iPS.T = _tl;
		update_leaf_TP!(photo_set, iPS, iHS, iEN);
		temperature_effects!(iHS, FT(_tl));

		# calculate the fraction of sunlit and shaded leaves
		f_view = (can_opt.Ps[iRT] + can_opt.Ps[iRT+1]) / 2;
		for iLF in 1:nSL
			iPS.APAR[iLF] = can_rad.absPAR_sun[(iRT-1)*nSL+iLF] *
							FT(1e6);
			iPS.LAIx[iLF] = f_view * f_SL[iLF];
		end;
		iPS.APAR[end] = can_rad.absPAR_shade[iRT] * FT(1e6);
		iPS.LAIx[end] = 1 - f_view;
		
	end
	
	
	glw_mean = 0	
	
	subIter = 15
	
	for i_can in 1:n_canopy
		
		iEN = envirs[i_can];
		iHS = plant_hs.leaves[i_can];
		iPS = plant_ps[i_can];
		iRT = n_canopy + 1 - i_can;
		
		
		iHS.p_crt = FT(-p_crit);
		iPS.ec    = critical_flow(iHS, iPS.ec);
        iPS.ec    = max(FT(0), iPS.ec);


		
		# iterate for 15 times to find steady state solution
		for subI in 1:subIter
		# calculate the photosynthetic rates
			gas_exchange_new!(photo_set, iPS, iEN, GswDrive());
			#update_gsw!(iPS, stomata_model, photo_set, iEN, FT(deltaT/subIter),FT(1));
			update_gsw!(iPS, stomata_model, photo_set, iEN, FT(deltaT/subIter));
			gsw_control!(photo_set, iPS, iEN);
		end
		
		#1/1000 = exp(-(p_crt_new / -vc.b) ^ vc.c)
		# log(1000) = (p_crt_new / -vc.b) ^ vc.c
		# log(1000)^(1/vc.c) * -vc.b = p_crt_new
		iHS.p_crt = FT(-iHS.vc.b * log(1000)^(1/iHS.vc.c));
		iPS.ec    = critical_flow(iHS, iPS.ec);
        iPS.ec    = max(FT(0), iPS.ec);
		
		# update the flow rates
		for iLF in 1:(nSL+1)
			f_CO₂ += iPS.An[iLF] * iPS.LAIx[iLF] * iPS.LA;
			f_H₂O += iPS.g_lw[iLF] * (iPS.p_sat - iEN.p_H₂O) / iEN.p_atm *
					 iPS.LAIx[iLF] * iPS.LA;
		end;
		
		glw_mean += sum(iPS.g_lw)/(nSL+1)/n_canopy;
	end;
		
		#update_Weibull!(node, FT(5.0), FT(5.0));
		
		# update flow profile and pressure history along the tree
	branch_total_flow = 0;
	for i_can in 1:n_canopy
		iEN = envirs[i_can];
		iLF = plant_hs.leaves[i_can];
		iPS = plant_ps[i_can];
		iST = plant_hs.branch[i_can];
		iLF.flow = sum(iPS.g_lw .* iPS.LAIx) *
				   (iPS.p_sat - iEN.p_H₂O) / iEN.p_atm;
		#iLF.flow = FT(max(0,df.LE[i])/44200)*node.ga/node.la
		iST.flow = iLF.flow * iPS.LA;
		branch_total_flow += iST.flow;
	end;
	plant_hs.trunk.flow = branch_total_flow; #sum([iST.flow for iST in plant_hs.branch]);
	for iRT in plant_hs.roots
		iRT.flow = plant_hs.trunk.flow / length(plant_hs.roots);
	end;
	
	#update_Kmax!(node, FT(3));
	#update_Weibull!(node, FT(5.7), FT(0.95));
	
	subIter2 = 4
	for subI2 in 1:subIter2
		
		do_soil!(node, FT(df.RAIN[i])/subIter2, deltaT/subIter2, k_soil)		
		
		for i_root in eachindex(node.plant_hs.roots)
			rootI = node.plant_hs.root_index_in_soil[i_root]
			node.plant_hs.roots[i_root].p_ups = soil_p_25_swc(node.plant_hs.roots[1].sh, node.swc[rootI])
		end
		
		pressure_profile!(node.plant_hs, SteadyStateMode(); update=false);

		# update canopy layer p_ups, which will be passed to each leaf
		for _i_can in eachindex(node.plant_hs.leaves)
			_iHS = node.plant_hs.leaves[_i_can];
			_iPS = node.plant_ps[_i_can];
			_iPS.p_ups = _iHS.p_ups;
		end
	end

	# calculate respiration from other components
	_r = photo_TD_from_set(rbase, t_soil);

	# update the data frame
	df.glw[i] = glw_mean
	df.NPP[i] = f_CO₂ / ga - _r;
	df.ETmod[i] = f_H₂O / ga;
	df.trunkflow[i] = plant_hs.trunk.flow
	df.stempot[i] = mean(node.plant_hs.trunk.p_element);
	df.leafpot[i] = node.plant_hs.leaves[length(node.plant_hs.leaves)].p_leaf;
	smc_record[i,:] = node.swc
	df.Ecrit[i] = node.plant_ps[length(node.plant_hs.leaves)].ec;
	df.Pcrit[i] = node.plant_hs.leaves[length(node.plant_hs.leaves)].p_crt;
	
	for ican in eachindex(node.plant_hs.leaves)
		psi_record_leaf[i,ican] = node.plant_hs.leaves[ican].p_leaf;
		psi_record_branch[i,ican] = mean(node.plant_hs.branch[ican].p_element);
	end

end;

return df, smc_record, psi_record_leaf, psi_record_branch

end

#=

etObs = sum(reshape(sim_res[1].LE/44200, (48,:)),dims=1)[1,:]/48;
etMod = sum(reshape(sim_res[1].ETmod, (48,:)),dims=1)[1,:]/48;
potMod = sum(reshape(sim_res[1].leafpot, (48,:)),dims=1)[1,:]/48;

plot(etObs); plot!(etMod)

psi_ex = -2.5:0.05:0; weib_ex = exp.(-1*(psi_ex / -1) .^ 5);
scatter(potMod[100:300], (etObs ./ etMod)[100:300]); plot!(psi_ex,weib_ex)

morningPot = sim_res[1].Oak_Psi[13:48:17520];
hasPot = .!(ismissing.(morningPot));
potObs = convert(Array{Float64}, morningPot[hasPot]);
plot(sim_res[1].leafpot[12:48:17520]); scatter!((1:365)[hasPot], potObs)
=#

