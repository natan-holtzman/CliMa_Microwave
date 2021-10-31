

function update_gsw!(clayer::CanopyLayer{FT},
            sm::EmpiricalStomatalModel{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            envir::AirLayer{FT},
            Δt::FT,
			betafac::FT;
            τ::FT = FT(600)
) where {FT<:AbstractFloat} 
(
    # unpack values
    @unpack g_sw, n_leaf, ps = clayer;

    # calculate steady state values
    # assume τ = 10 minutes
    for _iLF in 1:n_leaf
        _gsw_ss = max(sm.g0,
                      stomatal_conductance(sm, clayer, envir, betafac, _iLF));
        g_sw[_iLF] += (_gsw_ss - g_sw[_iLF]) / τ * Δt;
    end;

    return nothing
)
end



function update_gsw!(clayer::CanopyLayer{FT},
            sm::OSMWang{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            envir::AirLayer{FT},
            Δt::FT;
            τ::FT = FT(1e-5)
) where {FT<:AbstractFloat} 
(
    # unpack values
    @unpack An, ec, g_bc, g_bw, g_lw, g_m, g_max, g_min, g_sw, n_leaf, p_sat,
            ps = clayer;
    @unpack p_atm, p_H₂O = envir;

    # update g_sw
    for iLF in 1:n_leaf
        _gsw = g_sw[iLF] .+ FT(0.001);
        _glw = 1 / ( 1/_gsw + 1/g_bw[iLF] );
        _gsc = _gsw / FT(1.6);
        _glc = 1 / ( 1/_gsc + 1/g_m[iLF] + 1/g_bc[iLF] );
        leaf_photosynthesis!(photo_set, ps, envir, GCO₂Mode(), _glc);
        _∂A   = ps.An - An[iLF];
        _e0   = g_lw[iLF] * (p_sat - p_H₂O) / p_atm;
        _e1   = _glw * (p_sat - p_H₂O) / p_atm;
        _∂E   = _e1 - _e0;
        _∂A∂E = _∂A / _∂E;
        _∂Θ∂E = FT(max(0.1, An[iLF]) / max(ec - _e0, 1e-7));

        # ensure that dgsw does not change too rapidly
        _Δgsw = (_∂A∂E - _∂Θ∂E) * τ * Δt;
        if _Δgsw > 0
            _Δgsw = min(_Δgsw, (g_max-g_sw[iLF]) / 4);
        else
            _Δgsw = max(_Δgsw, (g_min-g_sw[iLF]) / 4);
        end;
        g_sw[iLF] += _Δgsw;
    end;

    return nothing
)
end


function gas_exchange_nofluor!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            drive::GswDrive,
            ind::Int
) where {FT<:AbstractFloat}
    @unpack g_max, g_min = canopyi;

    # update the conductances
    if canopyi.g_sw[ind] > g_max
        canopyi.g_sw[ind] = g_max;
    end

    if canopyi.g_sw[ind] < g_min
        canopyi.g_sw[ind] = g_min;
    end

    canopyi.g_lw[ind] = 1 / ( 1 / canopyi.g_sw[ind] +
                              1 / canopyi.g_bw[ind] );
    canopyi.g_sc[ind] = canopyi.g_sw[ind] / FT(1.6);
    canopyi.g_lc[ind] = 1 / ( FT(1.6) / canopyi.g_sw[ind] +
                              1 / canopyi.g_m[ind] +
                              1 / canopyi.g_bc[ind] );

    # update the photosynthetic rates
    leaf_photosynthesis!(photo_set, canopyi.ps, envir, GCO₂Mode(),
                         canopyi.g_lc[ind]);
   
    canopyi.Ac[ind] = canopyi.ps.Ac;
    canopyi.Aj[ind] = canopyi.ps.Aj;
    canopyi.Ap[ind] = canopyi.ps.Ap;
    canopyi.Ag[ind] = canopyi.ps.Ag;
    canopyi.An[ind] = canopyi.ps.An;
    canopyi.ϕs[ind] = canopyi.ps.ϕs;

    # update the pressures
    canopyi.p_i[ind] = canopyi.ps.p_i;
    canopyi.p_s[ind] = canopyi.ps.p_s;

    return nothing
end




function gas_exchange_new!(
            photo_set::AbstractPhotoModelParaSet{FT},
            canopyi::CanopyLayer{FT},
            envir::AirLayer{FT},
            drive::GswDrive
) where {FT<:AbstractFloat}
    # update the conductances for each "leaf"
    for i in eachindex(canopyi.g_sw)
        canopyi.ps.APAR = canopyi.APAR[i];
        leaf_ETR!(photo_set, canopyi.ps);
        gas_exchange_nofluor!(photo_set, canopyi, envir, drive, i);
    end

    return nothing
end


function update_LAI!(node::SPACMono{FT}, lai::FT) where {FT<:AbstractFloat}
    # 1. Update the general information of the SPAC
    node.la = node.ga * lai;

    # 2. Update the LAI for canopy radiation model
    node.canopy_rt.LAI = lai;

    # 3. Update the LA and LAI for stomatal and photosynthesis models
    for _iPS in node.plant_ps
        _iPS.LA   = node.la / node.n_canopy;
        _iPS.LAI  = lai / node.n_canopy;
        _iPS.tLAI = lai;
    end

    # 4. Update the LA for hydraulic model
    #    TODO abstract this to avoid memory allocation
    for _iHS in node.plant_hs.leaves
        _iHS.area = node.la / node.n_canopy;
    end

    return nothing
end


function update_VJR!(node::SPACMono{FT}, ratio::FT) where {FT<:AbstractFloat}

    # TODO change the ratio accordingly to photo_set
    # TODO add another ratio V2J in photo_set
    # Update Vcmax25, Jmax25 (1.67 Vcmax), and Rd25 (0.015 Vcmax)
    for _iPS in node.plant_ps
        _iPS.ps.Vcmax   = _iPS.ps.Vcmax25WW * ratio;
        _iPS.ps.Vcmax25 = _iPS.ps.Vcmax25WW * ratio;
        _iPS.ps.Jmax    = _iPS.ps.Jmax25WW * ratio;
        _iPS.ps.Jmax25  = _iPS.ps.Jmax25WW * ratio;
        _iPS.ps.Rd      = _iPS.ps.Rd25WW * ratio;
        _iPS.ps.Rd25    = _iPS.ps.Rd25WW * ratio;
    end

    return nothing
	
end

function do_soil!(node::SPACMono{FT}, rain_in::FT, deltaT::FT, k_soil_max::FT)
	layer_thick = -1*diff(node.soil_bounds)*1000

	runoff = FT(0);
		
	for i_root in eachindex(node.plant_hs.roots)
		soil_layerI = node.plant_hs.root_index_in_soil[i_root]
		
		out_layer = node.plant_hs.roots[i_root].flow /KG_H_2_MOL_S * FT(1/node.ga*deltaT/(60*60));
		
		node.swc[soil_layerI] -= out_layer / layer_thick[soil_layerI]
		node.swc[soil_layerI] = max(node.plant_hs.roots[soil_layerI].sh.Θr,min(node.swc[soil_layerI],node.plant_hs.roots[soil_layerI].sh.Θs));
	end	

	smc_rain = node.swc[1] + rain_in /layer_thick[1]
	smc_rain = max(min(smc_rain, node.plant_hs.roots[1].sh.Θs),node.plant_hs.roots[1].sh.Θr);
	node.swc[1] = smc_rain
	
	k_soil = zeros(length(node.swc))
	p_soil_tot = zeros(length(node.swc))
	midpoints = zeros(length(node.swc))

	for i_soil in eachindex(node.swc)
		k_soil[i_soil] = soil_k_ratio_swc(node.plant_hs.roots[1].sh, node.swc[i_soil])
		midpoint = (node.soil_bounds[i_soil] + node.soil_bounds[i_soil+1]) /2 * 1000
		midpoints[i_soil] = midpoint;
		p_soil_tot[i_soil] = soil_p_25_swc(node.plant_hs.roots[1].sh, node.swc[i_soil])*mpa2mm + midpoint
	end

	k_soil *= k_soil_max

	p_diff = -1*diff(p_soil_tot) #in mm
	layer_len = -1*diff(midpoints)
	q_down = k_soil[1:length(k_soil)-1] .* p_diff ./ layer_len * deltaT #mm/s * mm / mm * s

	#drainage = k_soil*FT(0.15*6)*deltaT
	drainage = k_soil*FT(0.2)*deltaT
	#drainage = 0
	
	runoff += drainage[length(k_soil)]
	
	node.swc[1:(length(k_soil)-1)] -= q_down ./ layer_thick[1:(length(k_soil)-1)]
	node.swc[2:length(k_soil)] += q_down ./ layer_thick[2:length(k_soil)]
	
	node.swc[length(k_soil)] -= drainage[length(k_soil)] ./ layer_thick[length(k_soil)]
	#node.swc -= drainage ./ layer_thick
	
	for i_soil in eachindex(node.swc)
		node.swc[i_soil] = max(min(node.swc[i_soil], node.plant_hs.roots[1].sh.Θs),node.plant_hs.roots[1].sh.Θr);
	end
	
	#return runoff
end	


function do_soil_drain!(node::SPACMono{FT}, rain_in::FT, deltaT::FT, k_soil_max::FT, slope_index::FT)
	layer_thick = -1*diff(node.soil_bounds)*1000

	runoff = FT(0);
		
	for i_root in eachindex(node.plant_hs.roots)
		soil_layerI = node.plant_hs.root_index_in_soil[i_root]
		
		out_layer = node.plant_hs.roots[i_root].flow /KG_H_2_MOL_S * FT(1/node.ga*deltaT/(60*60));
		
		node.swc[soil_layerI] -= out_layer / layer_thick[soil_layerI]
		node.swc[soil_layerI] = max(node.plant_hs.roots[soil_layerI].sh.Θr,min(node.swc[soil_layerI],node.plant_hs.roots[soil_layerI].sh.Θs));
	end	

	smc_rain = node.swc[1] + rain_in /layer_thick[1]
	smc_rain = max(min(smc_rain, node.plant_hs.roots[1].sh.Θs),node.plant_hs.roots[1].sh.Θr);
	node.swc[1] = smc_rain
	
	k_soil = zeros(length(node.swc))
	p_soil_tot = zeros(length(node.swc))
	midpoints = zeros(length(node.swc))

	for i_soil in eachindex(node.swc)
		k_soil[i_soil] = soil_k_ratio_swc(node.plant_hs.roots[1].sh, node.swc[i_soil])
		midpoint = (node.soil_bounds[i_soil] + node.soil_bounds[i_soil+1]) /2 * 1000
		midpoints[i_soil] = midpoint;
		p_soil_tot[i_soil] = soil_p_25_swc(node.plant_hs.roots[1].sh, node.swc[i_soil])*mpa2mm + midpoint
	end

	k_soil *= k_soil_max

	p_diff = -1*diff(p_soil_tot) #in mm
	layer_len = -1*diff(midpoints)
	q_down = k_soil[1:length(k_soil)-1] .* p_diff ./ layer_len * deltaT #mm/s * mm / mm * s

	#drainage = k_soil*FT(0.15*6)*deltaT
	drainage = k_soil*slope_index*deltaT
	#drainage = 0
	
	runoff += drainage[length(k_soil)]
	
	node.swc[1:(length(k_soil)-1)] -= q_down ./ layer_thick[1:(length(k_soil)-1)]
	node.swc[2:length(k_soil)] += q_down ./ layer_thick[2:length(k_soil)]
	
	node.swc[length(k_soil)] -= drainage[length(k_soil)] ./ layer_thick[length(k_soil)]
	#node.swc -= drainage ./ layer_thick
	
	for i_soil in eachindex(node.swc)
		node.swc[i_soil] = max(min(node.swc[i_soil], node.plant_hs.roots[1].sh.Θs),node.plant_hs.roots[1].sh.Θr);
	end
	
	#return runoff
end	


function do_soil_nss!(node::SPACMono{FT}, rain_in::FT, deltaT::FT, k_soil_max::FT)
	layer_thick = -1*diff(node.soil_bounds)*1000

	runoff = FT(0);
		
	for i_root in eachindex(node.plant_hs.roots)
		soil_layerI = node.plant_hs.root_index_in_soil[i_root]
		
		out_layer = node.plant_hs.roots[i_root].q_in /KG_H_2_MOL_S * FT(1/node.ga*deltaT/(60*60));
		
		node.swc[soil_layerI] -= out_layer / layer_thick[soil_layerI]
		node.swc[soil_layerI] = max(node.plant_hs.roots[soil_layerI].sh.Θr,min(node.swc[soil_layerI],node.plant_hs.roots[soil_layerI].sh.Θs));
	end	

	smc_rain = node.swc[1] + rain_in /layer_thick[1]
	smc_rain = max(min(smc_rain, node.plant_hs.roots[1].sh.Θs),node.plant_hs.roots[1].sh.Θr);
	node.swc[1] = smc_rain
	
	k_soil = zeros(length(node.swc))
	p_soil_tot = zeros(length(node.swc))
	midpoints = zeros(length(node.swc))

	for i_soil in eachindex(node.swc)
		k_soil[i_soil] = soil_k_ratio_swc(node.plant_hs.roots[1].sh, node.swc[i_soil])
		midpoint = (node.soil_bounds[i_soil] + node.soil_bounds[i_soil+1]) /2 * 1000
		midpoints[i_soil] = midpoint;
		p_soil_tot[i_soil] = soil_p_25_swc(node.plant_hs.roots[1].sh, node.swc[i_soil])*mpa2mm + midpoint
	end

	k_soil *= k_soil_max

	p_diff = -1*diff(p_soil_tot) #in mm
	layer_len = -1*diff(midpoints)
	q_down = k_soil[1:length(k_soil)-1] .* p_diff ./ layer_len * deltaT #mm/s * mm / mm * s

	#drainage = k_soil*FT(0.15*6)*deltaT
	drainage = k_soil*FT(0.2)*deltaT
	#drainage = 0
	
	runoff += drainage[length(k_soil)]
	
	node.swc[1:(length(k_soil)-1)] -= q_down ./ layer_thick[1:(length(k_soil)-1)]
	node.swc[2:length(k_soil)] += q_down ./ layer_thick[2:length(k_soil)]
	
	node.swc[length(k_soil)] -= drainage[length(k_soil)] ./ layer_thick[length(k_soil)]
	#node.swc -= drainage ./ layer_thick
	
	for i_soil in eachindex(node.swc)
		node.swc[i_soil] = max(min(node.swc[i_soil], node.plant_hs.roots[1].sh.Θs),node.plant_hs.roots[1].sh.Θr);
	end
	
	#return runoff
end	

