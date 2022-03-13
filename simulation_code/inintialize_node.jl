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

#this function creates a SPAC node and modifies it by editing the 
#storage volume and PV curves

function create_moflux_node(vcmax_par, k_plant, z_soil, smc0, storage_mult,
	nslice,deltaT,alpha,n,full_can,full_angles,root_dist_par,
	canopy_pvslope, trunk_pvslope)

node = create_spac(OSMWang{FT}(),vcmax_par,k_plant, z_soil, FT(40),
 nslice,alpha,n,full_can,full_angles,root_dist_par);

 
@unpack angles, can_opt, can_rad, canopy_rt, envirs, f_SL, ga, in_rad,
		latitude, leaves_rt, n_canopy, photo_set, plant_hs, plant_ps,
		rt_con, soil_opt, stomata_model, wl_set, n_root, soil_bounds = node;

canopyPV = PVCurveLinear{Float32}(canopy_pvslope, FT(1e-4));
trunkPV = PVCurveLinear{Float32}(trunk_pvslope, FT(1e-4));

for i_can in 1:n_canopy
	plant_hs.leaves[i_can].pv = canopyPV
	plant_hs.branch[i_can].pv = canopyPV
end

plant_hs.trunk.pv = trunkPV;
	
for i_root in 1:n_root
	plant_hs.roots[i_root].pv = trunkPV
end

total_water_vol = FT(storage_mult * 12 * node.ga * 1000/18.02); #in mol over the tree's area

#actually the branch volumes should be normalized by each branch's share of total branch length
for i in 1:n_canopy
	plant_hs.leaves[i].v_maximum = FT(0.1)*total_water_vol / n_canopy /  plant_hs.leaves[i].area
	plant_hs.branch[i].v_maximum[1] = FT(0.1)*total_water_vol / n_canopy
end

plant_hs.trunk.v_maximum[1] = FT(0.8)*total_water_vol

for i in 1:n_root
	plant_hs.roots[i].v_maximum[1] = FT(0.05)*total_water_vol / n_root #* plant_hs.roots[i].area
end

w_soil = FT(smc0);		
node.swc = ones(FT, length(node.swc))*w_soil;
ψ_soil = soil_p_25_swc(plant_hs.roots[1].sh, w_soil);
for root in plant_hs.roots
	root.p_ups = ψ_soil;
end;

pot_init = ψ_soil - 0.025

for i in 1:n_canopy
	plant_hs.leaves[i].v_storage = (FT(pot_init)*plant_hs.leaves[i].pv.slope + FT(1)) * plant_hs.leaves[i].v_maximum
	plant_hs.branch[i].v_storage = (FT(pot_init)*plant_hs.branch[i].pv.slope + FT(1)) * plant_hs.branch[i].v_maximum
end

plant_hs.trunk.v_storage = (FT(pot_init)*plant_hs.trunk.pv.slope + FT(1)) * plant_hs.trunk.v_maximum

for i in 1:n_root
	plant_hs.roots[i].v_storage = (FT(pot_init)*plant_hs.roots[i].pv.slope + FT(1)) * plant_hs.roots[i].v_maximum
end

#then calculate the steady state plant water storage
#and set plant state to that
update_pk_tree!(plant_hs)

volumes0 = get_v_prof2(plant_hs);
dmat0, dvec0 = create_deriv_mat(plant_hs)
set_vol!(plant_hs, -inv(dmat0)*dvec0)

update_pk_tree!(plant_hs)

return node

end
