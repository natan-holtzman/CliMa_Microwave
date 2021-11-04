#sim_res = run_sim(FT(22),FT(1.5),FT(5),FT(10),FT(0.025*0.4),FT(32));

#sim_res1 = run_sim(FT(1.8e1), FT(8.2e-1), FT(1.5e1), FT(9.5e-4), FT(2.5), FT(2.0), FT(1.1e3), FT(4.2e-1), 48*365*1 - 52*48, 48*365, 0.35);


using StatsBase
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
include("cap_funs.jl")
include("create_node_newvol.jl")


#using Plots

K_STEFAN = FT(Stefan());

#N      = 48*7 #only run for 4 days for quick testing
#istart = 48*365*1 - 52*48 + 165*48 #+ 180*48 - 3*48
#N = 48*365
#istart = 48*365*2 - 52*48 #+ 165*48
df_raw = CSV.read("../../data/moflux_land_data_newnames_7.csv", DataFrame);

function run_sim(vcmax_par, p_crit, k_plant, k_soil, b_soil, p20, z_soil, n_soil, istart, N, smc0, storage_mult)

df = deepcopy(df_raw[istart+1:istart+N,:]);

df[!,"leafpot"] = zeros(N)
df[!,"leafpotStore"] = zeros(N)

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

df[!,"pl1"] = zeros(N)
df[!,"pl2"] = zeros(N)

node = create_moflux_node(vcmax_par, p_crit, k_plant, k_soil, b_soil, p20, z_soil, n_soil, smc0, storage_mult);


smc_record = zeros(FT,N,length(node.swc));
psi_record_leaf = zeros(FT,N,length(node.plant_hs.leaves));
psi_record_branch = zeros(FT,N,length(node.plant_hs.leaves));


@unpack angles, can_opt, can_rad, canopy_rt, envirs, f_SL, ga, in_rad,
		latitude, leaves_rt, n_canopy, photo_set, plant_hs, plant_ps,
		rt_con, soil_opt, stomata_model, wl_set, n_root = node;
@unpack lidf, nAzi, nIncl = canopy_rt;
@unpack dWL, iPAR = wl_set;

rbase =  Q10TD{FT}(0, 298.15, 1.7)

global hs0 = deepcopy(plant_hs);

in_rad_bak = deepcopy(in_rad);
nSL = nAzi * nIncl;
in_Erad = in_rad_bak.E_direct .+ in_rad_bak.E_diffuse;
#in_PPFD = sum( e2phot(dWL, in_Erad)[iPAR] ) * FT(1e6);
in_PPFD = sum( Land.CanopyLayers.e2phot(wl_set.WL, in_Erad/1000)[iPAR] .* dWL[iPAR] ) * FT(1e6);

# iterate through the weather data
for i in eachindex(df.Day)

	
	if i==2
		global hs1 = deepcopy(plant_hs);
	end
	
	if i==3
		global hs2 = deepcopy(plant_hs);
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
	
	df.pl1[i] = node.plant_hs.leaves[1].p_leaf*1

	
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
	
	df.pl2[i] = node.plant_hs.leaves[1].p_leaf*1
	
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
		
		
		iHS.p_crt = -iHS.vc.b * log(FT(1000)) ^ (1/iHS.vc.c);
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
	
	if i < 1
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
	end
		
	if i >= 1
	for i_can in 1:n_canopy
		iEN = envirs[i_can];
		iLF = plant_hs.leaves[i_can];
		iPS = plant_ps[i_can];
		iST = plant_hs.branch[i_can];
		iLF.q_out = sum(iPS.g_lw .* iPS.LAIx) *
				   (iPS.p_sat - iEN.p_H₂O) / iEN.p_atm #/ iPS.LAI ;
		#iLF.q_out = FT(max(0,df.LE[i])/44200)*node.ga/node.la;

	end;
	end
	
	subIter3 = 1
	subIter2 = 4
	
	dd1, dd2 = update_cap_mat!(plant_hs,deltaT);
	
	for subI2 in 1:subIter2
		
		
		if i >= 1
			#=
			for si3 in 1:subIter3
				try
					dd1, dd2 = update_cap_mat!(node,deltaT/subIter2/subIter3);
				catch err_update
					return df, smc_record, psi_record_leaf, psi_record_branch
				end
			end
			=#
			do_soil_nss!(node, FT(df.RAIN[i])/subIter2, deltaT/subIter2, k_soil)		
		end
		if i < 1
			pressure_profile!(node.plant_hs, SteadyStateMode(); update=false);
			do_soil!(node, FT(df.RAIN[i])/subIter2, deltaT/subIter2, k_soil)
		end
		
		for i_root in eachindex(node.plant_hs.roots)
			rootI = node.plant_hs.root_index_in_soil[i_root]
			node.plant_hs.roots[i_root].p_ups = soil_p_25_swc(node.plant_hs.roots[1].sh, node.swc[rootI])
		end
				
		
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
	df.trunkflow[i] = plant_hs.trunk.q_in
	df.stempot[i] = mean(node.plant_hs.trunk.p_element);
	df.leafpot[i] = node.plant_hs.leaves[length(node.plant_hs.leaves)].p_leaf;
	smc_record[i,:] = node.swc
	df.Ecrit[i] = node.plant_ps[length(node.plant_hs.leaves)].ec;
	df.Pcrit[i] = node.plant_hs.leaves[length(node.plant_hs.leaves)].p_crt;
	df.leafStore[i] = node.plant_hs.leaves[length(node.plant_hs.leaves)].v_storage;
	df.leafpotStore[i] = node.plant_hs.leaves[length(node.plant_hs.leaves)].p_storage;
	
	for ican in eachindex(node.plant_hs.leaves)
		psi_record_leaf[i,ican] = node.plant_hs.leaves[ican].p_leaf;
		psi_record_branch[i,ican] = mean(node.plant_hs.branch[ican].p_element);
	end

end;

return df, smc_record, psi_record_leaf, psi_record_branch, (hs0, hs1, hs2)

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

