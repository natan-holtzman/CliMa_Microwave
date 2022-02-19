function create_tree2(
#modified the default create_tree function to include the number of vertical slices as an argument
#which was by default set to 5
            z_root::FT,
            z_trunk::FT,
            z_canopy::FT,
            soil_bounds::Array{FT,1},
            air_bounds::Array{FT,1},
			N_subunit::Int
) where {FT<:AbstractFloat}
    # determine how many layers in roots and canopy
    #=
	_n_root  = 0;
	_r_index = Int[]
	for i in eachindex(soil_bounds)
        _z = soil_bounds[i];
        if _z > z_root
            _n_root += 1;
            push!(_r_index, i);
        else
            break
        end
    end
	=#
	#_n_root = 4;
	#_r_index = [4,5,6,7]
	
	#distribute roots in lower half of soil column
	#_r_index = [i for i in eachindex(soil_bounds[1:(end-1)]) if abs(soil_bounds[i]) >= abs(soil_bounds[end]/2)];
	_r_index = eachindex(soil_bounds[1:(end-1)]);

	_n_root = length(_r_index);

    _n_canopy = 0;
    _c_index  = Int[];
    for i in eachindex(air_bounds)
        _z0 = air_bounds[i];
        _z1 = air_bounds[i+1];
        if _z0 <= z_trunk < z_canopy <= _z1
            _n_canopy += 1;
            push!(_c_index, i);
            break
        elseif _z0 < z_trunk < _z1 < z_canopy
            _n_canopy += 1;
            push!(_c_index, i);
        elseif z_trunk <= _z0 < _z1 < z_canopy
            _n_canopy += 1;
            push!(_c_index, i);
        elseif z_trunk <= _z0 <= z_canopy < _z1
            _n_canopy += 1;
            push!(_c_index, i-1);
            break
        elseif z_canopy <= _z0 < _z1
            break
        end
    end

	
    # create evenly distributed root system for now
    _roots = RootHydraulics{FT}[];
    
	soil_thick = abs.(soil_bounds[2:end] - soil_bounds[1:(end-1)]);
	soil_mid = abs.(soil_bounds[2:end] + soil_bounds[1:(end-1)])/2;
	#root_dist_len = exp.(-1*abs(2/z_root)*soil_mid);
	root_dist_len = exp.(-2*soil_mid);
	root_dist_layer = root_dist_len .* soil_thick / sum(root_dist_len .* soil_thick);
	
	for i in _r_index
        _Δh = abs(soil_bounds[i+1] + soil_bounds[i]) / 2;
		dZ = abs(soil_bounds[i+1]) - abs(soil_bounds[i]);
		zfrac = dZ / abs(z_root);
		#=
		if i <= 4
			zfrac = 0;
		else
			zfrac = 1/4;
		end
		=#
        _rt = RootHydraulics{FT}(N=N_subunit,
                    #area=1*zfrac,
					area=FT(root_dist_layer[i]),
                    k_max=25,#*zfrac,
                    k_rhiz=5e14,#*zfrac,
                    Δh=_Δh);
        push!(_roots, _rt);
    end
	
	#=
	for i in _r_index
        _Δh = abs(soil_bounds[i+1] + soil_bounds[i]) / 2;
        _rt = RootHydraulics{FT}(N=N_subunit,
                    area=1/_n_root,
                    k_max=25/_n_root,
                    k_rhiz=5e14/_n_root,
                    Δh=_Δh);
        push!(_roots, _rt);
    end
	=#
		

    # create Trunk
    _k_trunk = z_canopy / z_trunk * 50;
    _trunk   = StemHydraulics{FT}(N=N_subunit,
                    k_max=_k_trunk,
                    Δh=z_trunk);

    # create Branches
    _k_branch = z_canopy / (z_canopy - z_trunk) * 50;
    _branch   = StemHydraulics{FT}[];
    for i in _c_index
        _Δh = (air_bounds[i] + max(air_bounds[i+1], z_trunk)) / 2 - z_trunk;
        _st = StemHydraulics{FT}(N=N_subunit,
                    area=1/_n_canopy,
                    k_max=_k_branch/_n_canopy,
                    Δh=_Δh);
        push!(_branch, _st);
    end

    # create leaves
    _leaves = [LeafHydraulics{FT}(N=N_subunit,area=1500/_n_canopy) for i in 1:_n_canopy];

    # return plant
    return TreeLikeOrganism{FT}(_n_root,
                                _n_canopy,
                                _roots,
                                _trunk,
                                _branch,
                                _leaves,
                                _r_index,
                                _c_index,
                                zeros(FT,_n_root),
                                zeros(FT,_n_root),
                                zeros(FT,_n_root))
end