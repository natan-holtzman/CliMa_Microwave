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

KG_H_2_MOL_S = FT(55.55 / 3600);
mpa2mm = FT(10^6/9.8);
K_STEFAN = FT(Stefan());

#deltaT = FT(60*60)

include(string(PROJECT_HOME,"/simulation_code/create_spac_setsoil_10.jl"))
#include(string(PROJECT_HOME,"/simulation_code/create_spac_multi.jl"))

include(string(PROJECT_HOME,"/simulation_code/land_utils4.jl"))
include(string(PROJECT_HOME,"/simulation_code/cap_funs_inv.jl"))
include(string(PROJECT_HOME,"/simulation_code/create_node_newvol_varyPV_old.jl"))

function run_sim_varyB(vcmax_par, k_rel_crit, k_weibB, k_weibC, k_plant, k_soil, z_soil, istart, N,
 smc0, storage_mult, buffrate,scheme_number,df_raw,g1,deltaT,alpha,n,full_can, full_angles,
 drainage_smc, root_dist_par,canopy_pvslope,trunk_pvslope)

df = deepcopy(df_raw[istart:(istart+N-1),:]);

df[!,"leafpot"] = zeros(N)
df[!,"leafpotStore"] = zeros(N)

df[!,"VPD"] = zeros(N)

#df[!,"leafvol"] = zeros(N)
df[!,"stempot"] = zeros(N)
df[!,"soilpot"] = zeros(N)
df[!,"ETmod"] = zeros(N)
df[!,"trunk_in"] = zeros(N)
df[!,"trunk_out"] = zeros(N)

#df[!,"APARsun"] = zeros(N)
#df[!,"Zenith"] = zeros(N)

#df[!,"APARshade"] = zeros(N)

#df[!,"APAR_ps"] = zeros(N)

#df[!,"Ecrit"] = zeros(N)

df[!,"leafStore"] = zeros(N)

df[!,"ColumnSMC"] = zeros(N)
df[!,"ColumnSWP"] = zeros(N)
#df[!,"fvis"] = zeros(N)
df[!,"glw"] = zeros(N)
df[!,"GPP"] = zeros(N)
df[!,"NPP"] = zeros(N)
#df[!,"Pcrit"] = zeros(N)
df[!,"Runoff"] = zeros(N)

#df[!,"pl1"] = zeros(N)
#df[!,"pl2"] = zeros(N)

#if scheme_number == 3
node = create_moflux_node(vcmax_par, k_plant, z_soil, smc0, storage_mult,1,deltaT,alpha,n,
full_can,full_angles,root_dist_par,
canopy_pvslope, trunk_pvslope);
#else
#	node = create_spac(OSMWang{FT}(),vcmax_par,k_plant,z_soil, FT(40),alpha,n);
#end

beta_curve = WeibullSingle(FT(k_weibB*(1-k_rel_crit)), FT(k_weibC));
sm1 = Land.StomataModels.ESMMedlyn(FT(0),FT(g1));
#sm1 = Land.StomataModels.ESMBallBerry(FT(0),FT(g1));

node.stomata_model = sm1;

#k_weibB = 3;
#k_weibC = 2;
#k_rel_crit = FT(0.33);
update_Weibull!(node, FT(k_weibB), FT(k_weibC));

@unpack angles, can_opt, can_rad, canopy_rt, envirs, f_SL, ga, in_rad,
		latitude, leaves_rt, n_canopy, photo_set, plant_hs, plant_ps,
		rt_con, soil_opt, stomata_model, wl_set, n_root = node;
@unpack lidf, nAzi, nIncl = canopy_rt;
@unpack dWL, iPAR = wl_set;
in_rad_bak = deepcopy(in_rad);
nSL = nAzi * nIncl;
in_Erad = in_rad_bak.E_direct .+ in_rad_bak.E_diffuse;
#in_PPFD = sum( e2phot(dWL, in_Erad)[iPAR] ) * FT(1e6);
in_PPFD = sum( Land.CanopyLayers.e2phot(wl_set.WL, in_Erad/1000)[iPAR] .* dWL[iPAR] ) * FT(1e6);

maxvol = get_v_max(plant_hs);


if scheme_number == 2

	#node.plant_hs.branch[1].v_maximum = node.plant_hs.branch[2].v_maximum
	#node.plant_hs.branch[1].v_storage = node.plant_hs.branch[2].v_storage

	mypv = PVCurveLinear{Float32}(0.2f0, FT(buffrate));
	for i_can in 1:n_canopy
		plant_hs.leaves[i_can].pv = mypv
		plant_hs.branch[i_can].pv = mypv
	end

	plant_hs.trunk.pv = mypv;
		
	for i_root in 1:n_root
		plant_hs.roots[i_root].pv = mypv
	end


	kstore = storage_mult
	for i in 1:n_canopy
		plant_hs.leaves[i].v_maximum = FT(kstore) * plant_hs.leaves[i].v_maximum
		plant_hs.branch[i].v_maximum = FT(kstore) * plant_hs.branch[i].v_maximum
	end

	plant_hs.trunk.v_maximum = FT(kstore) * plant_hs.trunk.v_maximum

	for i in 1:n_root
		plant_hs.roots[i].v_maximum = FT(kstore) * plant_hs.roots[i].v_maximum
	end



	kinit = 0.95
	for i in 1:n_canopy
		plant_hs.leaves[i].v_storage = FT(kinit) * plant_hs.leaves[i].v_maximum
		plant_hs.branch[i].v_storage = FT(kinit) * plant_hs.branch[i].v_maximum
	end

	plant_hs.trunk.v_storage = FT(kinit) * plant_hs.trunk.v_maximum

	for i in 1:n_root
		plant_hs.roots[i].v_storage = FT(kinit) * plant_hs.roots[i].v_maximum
	end
	
	
	w_soil = FT(smc0);		
	node.swc = ones(FT, length(node.swc))*w_soil;
	# update soil water matrices
	# need to adjust SWC to avoid problem in residual SWC at Niwot Ridge
	ψ_soil = soil_p_25_swc(plant_hs.roots[1].sh, w_soil);
	for root in plant_hs.roots
		root.p_ups = ψ_soil;
	end;

	pot_init = 1+0.2*ψ_soil
	
	for i in 1:n_canopy
		plant_hs.leaves[i].v_storage = FT(pot_init) * plant_hs.leaves[i].v_maximum
		plant_hs.branch[i].v_storage = FT(pot_init) * plant_hs.branch[i].v_maximum
	end

	plant_hs.trunk.v_storage = FT(pot_init) * plant_hs.trunk.v_maximum

	for i in 1:n_root
		plant_hs.roots[i].v_storage = FT(pot_init) * plant_hs.roots[i].v_maximum
	end
	
	
end

if scheme_number == 1
	w_soil = FT(smc0);		
	node.swc = ones(FT, length(node.swc))*w_soil;
	# update soil water matrices
	# need to adjust SWC to avoid problem in residual SWC at Niwot Ridge
	ψ_soil = soil_p_25_swc(plant_hs.roots[1].sh, w_soil);
	for root in plant_hs.roots
		root.p_ups = ψ_soil;
	end;
end


smc_record = zeros(FT,N,length(node.swc));
smc_record2 = zeros(FT,N,length(node.swc));

psi_record_leaf = zeros(FT,N,length(node.plant_hs.leaves));
psi_record_branch = zeros(FT,N,length(node.plant_hs.leaves));
leaf_qin_record = zeros(FT,N,length(node.plant_hs.leaves));
leaf_qout_record = zeros(FT,N,length(node.plant_hs.leaves));
apar_record = zeros(FT,N,length(node.plant_hs.leaves));
anet_record = zeros(FT,N,length(node.plant_hs.leaves));
gs_record = zeros(FT,N,length(node.plant_hs.leaves));

root_qin_record = zeros(FT,N,length(node.plant_hs.roots));
root_qout_record = zeros(FT,N,length(node.plant_hs.roots));

v_profile = zeros(FT,N,n_canopy*2+1+n_root);
p_profile = zeros(FT,N,n_canopy*2+1+n_root);
soil_p_profile = zeros(FT,N,n_root);


psi_record_trunk = zeros(FT,N);


rbase =  Photosynthesis.Q10TD{FT}(4, 298.15, 1.7)
drainage_pot = soil_p_25_swc(node.plant_hs.roots[1].sh, drainage_smc);

# iterate through the weather data
for i in eachindex(df.Day)

#=
	if mod(i, (24*365)) == 1
		println(Int((i-1)/(24*365)))
		#println(node.plant_hs.leaves[length(node.plant_hs.leaves)].p_leaf)
	end	
=#

	t_soil = FT(df.T_SOIL[1] + 273.15);
	la_previous = FT(node.la);
	update_LAI!(node, FT(df.LAI_modis[i]));
	if i > 1
		update_pv_leaf!(plant_hs, la_previous, node.la);
	end
	#update_LAI!(node, FT(2.5));


	#update_LAI!(node, FT(df.LAI_modis[i]*0.5 + 1));
	#update_LAI!(node, FT(mean(df.LAI_modis*0.5 .+ 1)));
	

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
	
	#df.pl1[i] = node.plant_hs.leaves[1].p_leaf*1

	
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
		
		iEN.RH    = FT(df.RelHum[i]);
		iEN.p_H₂O = iEN.RH * iEN.p_sat;
		iEN.vpd   = (FT(1)-iEN.RH) * iEN.p_sat;
	
		#iEN.vpd   = df.VPD[i] * 100;
		#iEN.p_H₂O = iEN.p_sat - iEN.vpd;
		#iEN.RH    = iEN.p_H₂O / iEN.p_sat;
		
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
	
	#df.pl2[i] = node.plant_hs.leaves[1].p_leaf*1
	
	glw_mean = 0	
	
	subIter = 15
	
	for i_can in 1:n_canopy
		
		iEN = envirs[i_can];
		iHS = plant_hs.leaves[i_can];
		iPS = plant_ps[i_can];
		iRT = n_canopy + 1 - i_can;
		
		iHS.p_crt = FT(-iHS.vc.b * log(FT(1/k_rel_crit))^(1/iHS.vc.c));
		#iHS.p_crt = FT(-8);
		iPS.ec    = critical_flow(iHS, iPS.ec);
        iPS.ec    = max(FT(0), iPS.ec);

		# iterate for 15 times to find steady state solution
		for subI in 1:subIter
		# calculate the photosynthetic rates
			gas_exchange_new!(photo_set, iPS, iEN, GswDrive());
			#prognostic_gsw!(iPS, iEN, stomata_model, FT(1), deltaT/subIter);
			update_gsw!(iPS, stomata_model, photo_set, iEN, FT(deltaT/subIter),FT(1));
			#update_gsw!(iPS, stomata_model, photo_set, iEN, FT(deltaT/subIter));
			gsw_control!(photo_set, iPS, iEN);
		end
		
		iHS.p_crt = -iHS.vc.b * log(FT(1000)) ^ (1/iHS.vc.c);
		iPS.ec    = critical_flow(iHS, iPS.ec);
        iPS.ec    = max(FT(0), iPS.ec);
		
		# update the flow rates
		for iLF in 1:(nSL+1)
			f_CO₂ += iPS.An[iLF] * iPS.LAIx[iLF] * iPS.LA;
			f_H₂O += iPS.g_lw[iLF] * max(0,iPS.p_sat - iEN.p_H₂O) / iEN.p_atm *
					 iPS.LAIx[iLF] * iPS.LA;
		end;
		
		glw_mean += sum(iPS.g_lw)/(nSL+1)/n_canopy;
	end;
		
	# update flow profile and pressure history along the tree
	if scheme_number==1
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
		
	if (scheme_number==2) | (scheme_number==3)
		for i_can in 1:n_canopy
			iEN = envirs[i_can];
			iLF = plant_hs.leaves[i_can];
			iPS = plant_ps[i_can];
			iST = plant_hs.branch[i_can];
			iLF.q_out = sum(iPS.g_lw .* iPS.LAIx) *
					   max(0,iPS.p_sat - iEN.p_H₂O) / iEN.p_atm #/ iPS.LAI ;
			#iLF.q_out = FT(max(0,df.LE[i])/44200)*node.ga/node.la;

		end;
	end
	
	next_with_rain = node.swc[1] + df.RAIN[i]/abs(1000*node.soil_bounds[2]);
	subIter2 = 4
	if next_with_rain > 0.4
		subIter2 = 20
	end
	
	if scheme_number==3
		#dd1, dd2 = update_cap_mat!(plant_hs,deltaT);
		#=
		try
			dd1, dd2 = update_cap_mat!(plant_hs,deltaT);
		catch err
			return df, smc_record, psi_record_leaf, psi_record_branch, psi_record_trunk, node
		end
		=#
		
		try
			dd1, dd2 = find_new_cap(plant_hs,deltaT);
			push_cap!(plant_hs,dd1,dd2)
			if mean(dd1 .> maxvol) > 0
				return df, smc_record, root_qin_record, leaf_qout_record, v_profile, p_profile,soil_p_profile, apar_record, anet_record, gs_record, node
			end
#=
			if mean(dd1 .<= maxvol) == 1
				push_cap!(plant_hs,dd1,dd2)
			else
				cap_iter = 4
				for icap in 1:cap_iter
					dd1, dd2 = find_new_cap(plant_hs,deltaT/cap_iter);
					push_cap!(plant_hs,dd1,dd2)
				end
			end
=#
		catch err
			#return df, smc_record, psi_record_leaf, psi_record_branch, psi_record_trunk, node
			return df, smc_record, root_qin_record, leaf_qout_record, v_profile, p_profile,soil_p_profile, apar_record, anet_record, gs_record, smc_record2, node

			#=
			cap_iter = 4
			for icap in 1:cap_iter
				dd1, dd2 = find_new_cap(plant_hs,deltaT/cap_iter);
				push_cap!(plant_hs,dd1,dd2)
			end
			=#
		end
	end
	total_runoff = 0
	smcsol1 = 0*node.swc;
	for subI2 in 1:subIter2
		
		if scheme_number==1
			pressure_profile!(node.plant_hs, SteadyStateMode(); update=false);
			do_soil_drain!(node, FT(df.RAIN[i])/subIter2, deltaT/subIter2, k_soil,FT(0))
		end
		if scheme_number==2
			do_soil_nss_drain!(node, FT(df.RAIN[i])/subIter2, deltaT/subIter2, k_soil,FT(0))		
			update_PVF!(node.plant_hs,deltaT/subIter2);
		end
		if scheme_number==3
			runoff_i = do_soil_nss_drain3!(node, FT(df.RAIN[i])/subIter2, deltaT/subIter2, k_soil, drainage_pot)		
			#smcsol1 = soil_ode!(node, FT(df.RAIN[i])/subIter2, deltaT/subIter2, k_soil, slope_runoff)		
			total_runoff += runoff_i
		end
	end
	#total_runoff /= subIter2
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
	df.GPP[i] = f_CO₂ / ga;
	df.NPP[i] = f_CO₂ / ga - _r;
	df.ETmod[i] = f_H₂O / ga;
	df.trunk_in[i] = plant_hs.trunk.q_in
	df.trunk_out[i] = plant_hs.trunk.q_out

	df.stempot[i] = mean(node.plant_hs.trunk.p_element);
	df.leafpot[i] = node.plant_hs.leaves[length(node.plant_hs.leaves)].p_leaf;
	smc_record[i,:] = node.swc
	smc_record2[i,:] = smcsol1

    df.ColumnSMC[i] = sum(node.swc .* -diff(node.soil_bounds)) / abs(node.soil_bounds[end])
    swpI = [soil_p_25_swc(node.plant_hs.roots[1].sh, x) for x in node.swc];
    df.ColumnSWP[i] = sum(swpI .* -diff(node.soil_bounds)) / abs(node.soil_bounds[end])
	df.Runoff[i] = total_runoff;
	#df.Ecrit[i] = node.plant_ps[length(node.plant_hs.leaves)].ec;
	#df.Pcrit[i] = node.plant_hs.leaves[length(node.plant_hs.leaves)].p_crt;
	df.leafStore[i] = node.plant_hs.leaves[length(node.plant_hs.leaves)].v_storage;
	df.leafpotStore[i] = node.plant_hs.leaves[length(node.plant_hs.leaves)].p_storage;
	df.VPD[i] = node.envirs[1].vpd
        	
	for ican in eachindex(node.plant_hs.leaves)
		psi_record_leaf[i,ican] = node.plant_hs.leaves[ican].p_leaf;
		psi_record_branch[i,ican] = mean(node.plant_hs.branch[ican].p_element);
	end
	psi_record_trunk[i] = mean(node.plant_hs.trunk.p_element);
	root_qin_record[i,:] = [x.q_in for x in node.plant_hs.roots];
	#root_qout_record[i,:] = [x.q_out for x in node.plant_hs.roots];

	
	#leaf_qin_record[i,:] = [x.q_in*x.area for x in node.plant_hs.leaves];
	leaf_qout_record[i,:] = [x.q_out*x.area for x in node.plant_hs.leaves];
	apar_record[i,:] = [sum(x.APAR .* x.LAIx) for x in node.plant_ps];

	anet_record[i,:] = [sum(x.An .* x.LAIx) for x in node.plant_ps];
	gs_record[i,:] = [sum(x.g_lw .* x.LAIx) for x in node.plant_ps];

	p_profile[i,:] = get_p_prof2(node.plant_hs);	
	soil_p_profile[i,:] = swpI;	

	v_profile[i,:] = get_v_prof2(node.plant_hs);

end;

return df, smc_record, root_qin_record, leaf_qout_record, v_profile, p_profile,soil_p_profile, apar_record, anet_record, gs_record, smc_record2, node

end

