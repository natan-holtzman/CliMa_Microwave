
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

N = 48*365
istart = 48*365*2 - 52*48 #+ 230*48

const FT = Float32

struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

df_raw = CSV.read("../data/moflux_land_data_newnames_7.csv", DataFrame);

include("create_spac_setsoil.jl")
include("land_utils4.jl")


df = deepcopy(df_raw[istart+1:istart+N,:]);
i = argmax(df.T_AIR)# - 14*48;

deltaT = FT(30*60);
K_STEFAN = FT(Stefan());


node = create_spac(OSMWang{FT}(),FT(22),FT(5), FT(1000), FT(40));
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

_tl = (df.LW_OUT[i] / 0.97 / K_STEFAN ) ^ 0.25;			   
update_Weibull!(node, FT(1.5), FT(2.0));

function get_ET(hum, smois)


	w_soil = FT(smois);		
	node.swc = ones(FT, length(node.swc))*w_soil;
	ψ_soil = soil_p_25_swc(plant_hs.roots[1].sh, w_soil);
	for root in plant_hs.roots
		root.p_ups = ψ_soil;
	end;
	pressure_profile!(node.plant_hs, SteadyStateMode(); update=false);


	i_can = 1;
	iEN = envirs[i_can];
	iHS = plant_hs.leaves[i_can];
	iPS = plant_ps[i_can];
	iRT = n_canopy + 1 - i_can;

	
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
	


	# update environmental conditions
	iEN.t_air = df.T_AIR[i] + 273.15;
	iEN.p_atm = df.P_ATM[i] * 1000;
	iEN.p_a   = iEN.p_atm * 4e-4;
	iEN.p_O₂  = iEN.p_atm * 0.209;
	iEN.p_sat = saturation_vapor_pressure( param_set, iEN.t_air , Liquid());
	iEN.RH    = hum;
	iEN.p_H₂O = iEN.RH * iEN.p_sat;
	iEN.vpd   = (FT(1)-iEN.RH) * iEN.p_sat;
	iEN.wind  = df.WIND[i];

	iHS.p_crt = FT(-4);
	iPS.ec    = critical_flow(iHS, iPS.ec);
	iPS.ec    = max(FT(0), iPS.ec);

	subIter = 100;

	for subI in 1:subIter
	# calculate the photosynthetic rates
		gas_exchange_new!(photo_set, iPS, iEN, GswDrive());
		update_gsw!(iPS, stomata_model, photo_set, iEN, FT(deltaT/subIter));
		gsw_control!(photo_set, iPS, iEN);
	end

	fco2 = 0;
	fh2o = 0;

	for iLF in 1:(nSL+1)
		fco2 += iPS.An[iLF] * iPS.LAIx[iLF] * iPS.LA;
		fh2o += iPS.g_lw[iLF] * (iPS.p_sat - iEN.p_H₂O) / iEN.p_atm *
				 iPS.LAIx[iLF] * iPS.LA;
	end;
	
	return [fco2, fh2o, iEN.vpd]
end

hum_list = collect(0.23:0.025:0.995);
soil_list = collect(0.14:0.01:0.45);

a1 = get_ET(hum_list[1],soil_list[1])


outputs = zeros(length(hum_list),length(soil_list), 3);
for j in eachindex(hum_list)
	for k in eachindex(soil_list)
		outputs[j,k,:] = get_ET(hum_list[j],soil_list[k])
	end
end

vpd_list = outputs[:,1,3];

glw = outputs[:,:,2] ./ outputs[:,:,3];

#=
figure()
plot(vpd_list, glw[:,1],label="Low soil moisture")
plot(vpd_list, glw[:,end],label="High soil moisture")
xlabel("VPD")
ylabel("G_lw")
legend()

figure()
plot(soil_list, glw[1,:],label="High VPD")
plot(soil_list, glw[10,:],label="Low VPD")
xlabel("SMC")
ylabel("G_lw")
legend()
=#

istart2 = 48*365*1 - 52*48 #+ 230*48
df2 = deepcopy(df_raw[istart2+1:end,:]);
vpd_noon = df2.VPD[24:48:end];
vpd_noon[vpd_noon .< 0.5] .= 0.5;
lai_noon = df2.LAI_modis[24:48:end];
et_noon = df2.LE[24:48:end]/44200;
rad_in = df2.PPFD_in[24:48:end];
smc_noon = df2.SMC[24:48:end];
glw_noon = et_noon ./ vpd_noon;

figure()
plot(vpd_noon[(lai_noon .> 4) .&(smc_noon .< 0.2)], glw_noon[(lai_noon .> 4) .&(smc_noon .< 0.2)], "o")
plot(vpd_list/100, glw[:,7])
plot(vpd_list/100, glw[:,1])
xlabel("VPD")
ylabel("G_lw")

figure()
plot(vpd_noon[(lai_noon .> 4) .&(smc_noon .> 0.3)], glw_noon[(lai_noon .> 4) .&(smc_noon .> 0.3)], "o")
plot(vpd_list/100, glw[:,17])
plot(vpd_list/100, glw[:,end])
xlabel("VPD")
ylabel("G_lw")

#plot(vpd_noon[lai_noon .> 4], glw_noon[lai_noon .> 4], "o") 
#plot(vpd_noon[(lai_noon .> 4) .&(smc_noon .< 0.2)], glw_noon[(lai_noon .> 4) .&(smc_noon .< 0.2)], "o") 
#plot(vpd_noon[(lai_noon .> 4) .&(smc_noon .> 0.35)], glw_noon[(lai_noon .> 4) .&(smc_noon .> 0.35)], "o") 

#when it's dry, it behaves like has Kmax of 1
#when it's wet like has kmax of 5

