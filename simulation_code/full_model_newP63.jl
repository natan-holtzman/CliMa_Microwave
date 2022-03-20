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

KG_H_2_MOL_S = FT(55.55 / 3600);
mpa2mm = FT(10^6/9.8);
K_STEFAN = FT(Stefan());

include(string(PROJECT_HOME,"/simulation_code/create_spac.jl"))

include(string(PROJECT_HOME,"/simulation_code/land_utils.jl"))
include(string(PROJECT_HOME,"/simulation_code/capacitance_funs.jl"))
include(string(PROJECT_HOME,"/simulation_code/initialize_node.jl"))

function run_sim_varyB(vcmax_par, beta_weibB, k_weibB, k_weibC, k_plant, k_soil, z_soil, istart, N,
 smc0, storage_mult, buffrate,scheme_number,df_raw,g1,deltaT,alpha,n,full_can, full_angles,
 drainage_smc, root_dist_par,canopy_pvslope,trunk_pvslope)

df = deepcopy(df_raw[istart:(istart+N-1),:]);

df[!,"leafpot"] = zeros(N)
df[!,"leafpotStore"] = zeros(N)

df[!,"VPD"] = zeros(N)

df[!,"stempot"] = zeros(N)
df[!,"soilpot"] = zeros(N)
df[!,"ETmod"] = zeros(N)
df[!,"trunk_in"] = zeros(N)
df[!,"trunk_out"] = zeros(N)
df[!,"leafStore"] = zeros(N)

df[!,"ColumnSMC"] = zeros(N)
df[!,"ColumnSWP"] = zeros(N)
df[!,"glw"] = zeros(N)
df[!,"GPP"] = zeros(N)
df[!,"NPP"] = zeros(N)
df[!,"Runoff"] = zeros(N)

node = create_moflux_node(vcmax_par, k_plant, z_soil, smc0, storage_mult,1,deltaT,alpha,n,
full_can,full_angles,root_dist_par,
canopy_pvslope, trunk_pvslope);

beta_curve = WeibullSingle(FT(beta_weibB), FT(k_weibC));
sm1 = Land.StomataModels.ESMMedlyn(FT(0),FT(g1));
#sm1 = Land.StomataModels.ESMBallBerry(FT(0),FT(g1));

node.stomata_model = sm1;

update_Weibull!(node, FT(k_weibB), FT(k_weibC));

@unpack angles, can_opt, can_rad, canopy_rt, envirs, f_SL, ga, in_rad,
		latitude, leaves_rt, n_canopy, photo_set, plant_hs, plant_ps,
		rt_con, soil_opt, stomata_model, wl_set, n_root = node;
@unpack lidf, nAzi, nIncl = canopy_rt;
@unpack dWL, iPAR = wl_set;
in_rad_bak = deepcopy(in_rad);
nSL = nAzi * nIncl;
in_Erad = in_rad_bak.E_direct .+ in_rad_bak.E_diffuse;
in_PPFD = sum( Land.CanopyLayers.e2phot(wl_set.WL, in_Erad/1000)[iPAR] .* dWL[iPAR] ) * FT(1e6);

maxvol = get_v_max(plant_hs);

smc_record = zeros(FT,N,length(node.swc));

psi_record_leaf = zeros(FT,N,length(node.plant_hs.leaves));
psi_record_branch = zeros(FT,N,length(node.plant_hs.leaves));
leaf_qout_record = zeros(FT,N,length(node.plant_hs.leaves));
apar_record = zeros(FT,N,length(node.plant_hs.leaves));
anet_record = zeros(FT,N,length(node.plant_hs.leaves));
gs_record = zeros(FT,N,length(node.plant_hs.leaves));

root_qin_record = zeros(FT,N,length(node.plant_hs.roots));

v_profile = zeros(FT,N,n_canopy*2+1+n_root);
p_profile = zeros(FT,N,n_canopy*2+1+n_root);
soil_p_profile = zeros(FT,N,n_root);
psi_record_trunk = zeros(FT,N);

rbase =  Photosynthesis.Q10TD{FT}(4, 298.15, 1.7)
drainage_pot = soil_p_25_swc(node.plant_hs.roots[1].sh, drainage_smc);

# iterate through the weather data
for i in eachindex(df.Day)
	t_soil = FT(df.T_SOIL[1] + 273.15);
	la_previous = FT(node.la);
	update_LAI!(node, FT(df.LAI_modis[i]));
	#update_LAI!(node, FT(2));

	if i > 1
		update_pv_leaf!(plant_hs, la_previous, node.la);
	end

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

	# update fluxes
	f_H₂O = 0;
	f_CO₂ = 0;
	gppi = 0;

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
		
		iEN.RH    = FT(df.RelHum[i]);
		iEN.p_H₂O = iEN.RH * iEN.p_sat;
		iEN.vpd   = (FT(1)-iEN.RH) * iEN.p_sat;
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
		
		# iterate for 15 times to find steady state solution
		for subI in 1:subIter
		# calculate the photosynthetic rates
			gas_exchange_new!(photo_set, iPS, iEN, GswDrive());
			update_gsw!(iPS, stomata_model, photo_set, iEN, FT(deltaT/subIter),FT(1));
			gsw_control!(photo_set, iPS, iEN);
		end
		
		# update the flow rates
		for iLF in 1:(nSL+1)
			f_CO₂ += iPS.An[iLF] * iPS.LAIx[iLF] * iPS.LA;
			gppi += iPS.Ag[iLF] * iPS.LAIx[iLF] * iPS.LA;
			f_H₂O += iPS.g_lw[iLF] * max(0,iPS.p_sat - iEN.p_H₂O) / iEN.p_atm *
					 iPS.LAIx[iLF] * iPS.LA;
		end;
		
		glw_mean += sum(iPS.g_lw)/(nSL+1)/n_canopy;
	end;
		
	# update flow profile and pressure history along the tree		
	for i_can in 1:n_canopy
		iEN = envirs[i_can];
		iLF = plant_hs.leaves[i_can];
		iPS = plant_ps[i_can];
		iST = plant_hs.branch[i_can];
		iLF.q_out = sum(iPS.g_lw .* iPS.LAIx) * #iPS.LA *
					max(0,iPS.p_sat - iEN.p_H₂O) / iEN.p_atm;
		#iLF.q_out = FT(max(0,df.LE[i])/44200)*node.ga/node.la;
		#uncomment the above line to instead force ET to be prescribed from the flux tower data
	end;
	
	try  #try to run the plant hydraulics. if an error is returned, exit gracefully
		dd1, dd2 = find_new_cap(plant_hs,deltaT);
		push_cap!(plant_hs,dd1,dd2)
		if mean(dd1 .> maxvol) > 0
			return df, smc_record, root_qin_record, leaf_qout_record, v_profile, p_profile,soil_p_profile, apar_record, anet_record, gs_record, node
		end
	catch err
		return df, smc_record, root_qin_record, leaf_qout_record, v_profile, p_profile,soil_p_profile, apar_record, anet_record, gs_record, node
	end

	next_with_rain = node.swc[1] + df.RAIN[i]/abs(1000*node.soil_bounds[2]);
	subIter2 = 4
	if next_with_rain > 0.5
		subIter2 = 20
	end

	total_runoff = 0
	for subI2 in 1:subIter2
		runoff_i = do_soil_nss_drain3!(node, FT(df.RAIN[i])/subIter2, deltaT/subIter2, k_soil, drainage_pot)		
		total_runoff += runoff_i/subI2
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
	
	# calculate respiration from other components
	_r = Photosynthesis.photo_TD_from_set(rbase, t_soil);

	for i_can in 1:n_canopy
		ratioI = xylem_k_ratio(beta_curve, plant_hs.leaves[i_can].p_leaf);
		update_VJR_leaf!(node.plant_ps[i_can], ratioI)
	end
	
	# update the data frame
	df.glw[i] = glw_mean
	df.GPP[i] = gppi / ga;
	df.NPP[i] = f_CO₂ / ga - _r;
	df.ETmod[i] = f_H₂O / ga;
	df.trunk_in[i] = plant_hs.trunk.q_in
	df.trunk_out[i] = plant_hs.trunk.q_out

	df.stempot[i] = mean(node.plant_hs.trunk.p_element);
	df.leafpot[i] = node.plant_hs.leaves[length(node.plant_hs.leaves)].p_leaf;
	smc_record[i,:] = node.swc

    df.ColumnSMC[i] = sum(node.swc .* -diff(node.soil_bounds)) / abs(node.soil_bounds[end])
    swpI = [soil_p_25_swc(node.plant_hs.roots[1].sh, x) for x in node.swc];
    df.ColumnSWP[i] = sum(swpI .* -diff(node.soil_bounds)) / abs(node.soil_bounds[end])
	df.Runoff[i] = total_runoff;
	df.leafStore[i] = node.plant_hs.leaves[length(node.plant_hs.leaves)].v_storage;
	df.leafpotStore[i] = node.plant_hs.leaves[length(node.plant_hs.leaves)].p_storage;
	df.VPD[i] = node.envirs[1].vpd
        	
	for ican in eachindex(node.plant_hs.leaves)
		psi_record_leaf[i,ican] = node.plant_hs.leaves[ican].p_leaf;
		psi_record_branch[i,ican] = mean(node.plant_hs.branch[ican].p_element);
	end
	psi_record_trunk[i] = mean(node.plant_hs.trunk.p_element);
	root_qin_record[i,:] = [x.q_in for x in node.plant_hs.roots];
	leaf_qout_record[i,:] = [x.q_out*x.area for x in node.plant_hs.leaves];
	apar_record[i,:] = [sum(x.APAR .* x.LAIx) for x in node.plant_ps];

	anet_record[i,:] = [sum(x.An .* x.LAIx) for x in node.plant_ps];
	gs_record[i,:] = [sum(x.g_lw .* x.LAIx) for x in node.plant_ps];

	p_profile[i,:] = get_p_prof2(node.plant_hs);	
	soil_p_profile[i,:] = swpI;	

	v_profile[i,:] = get_v_prof2(node.plant_hs);

end;

return df, smc_record, root_qin_record, leaf_qout_record,
 v_profile, p_profile,soil_p_profile, apar_record, anet_record,
  gs_record, node

end

