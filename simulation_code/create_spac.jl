include(string(PROJECT_HOME,"/simulation_code/updatefuns.jl"))
include(string(PROJECT_HOME,"/simulation_code/create_hydraulics.jl"))

function make_layers(_totaldepth)
	linear_layers = collect(0:0.1:0.8);
	toadd = [0,0,1,2,3,4,5,6,7];
	quad_layers = cumsum(toadd)/sum(toadd);
	
	if _totaldepth < 0.8
		lower_layer_thick = (_totaldepth-0.1)/Int(floor((_totaldepth-0.1) / 0.1));
		layers = vcat(0, 0.1:lower_layer_thick:_totaldepth);
	end
	
	if _totaldepth >= 0.8
		layers = linear_layers + quad_layers*(_totaldepth-0.8);
	end
	return convert(Array{FT}, -1*layers)
end

function create_spac(
            sm::OptimizationStomatalModel = OSMWang{FT}(),
            vcmax::FT = FT(30),
            kmax::FT = FT(0.3),
		z_soil::FT = FT(1000),       
     chl::FT = FT(40),
	 N_slice::Int = 1,
	alpha::FT = FT(1.368),
	n::FT = FT(2.6257),
	full_can::Int = 1,
	full_angles::Int = 1,
	root_dist_par::FT = 2
) where {FT<:AbstractFloat} 
(
	_soil_hs = VanGenuchten{FT}(stype = "Ozark",
									α = alpha,
									n = n,
								   Θs = 0.55,
								   Θr = 0.067);

	_totaldepth = FT(z_soil/1000);
	_rootdepth = FT(z_soil/1000);
	_soil_bounds = make_layers(_totaldepth);
	
	_totaldepth = -1*_soil_bounds[end];
	_rootdepth = -1*_soil_bounds[end];

	if full_can == 1
		_air_bounds = convert(Array{FT},collect(0:1:20));
	else
		_air_bounds = convert(Array{FT},[0,5,9,12,16,18.5,20] );
	end;
	
	_tree_hs = create_tree2(FT(-1*_rootdepth), FT(9), FT(18.5), _soil_bounds,
						   _air_bounds, N_slice, root_dist_par);
	for _root in _tree_hs.roots
		_root.sh = deepcopy(_soil_hs);
	end;

	_node = SPACMono{FT}(soil_bounds = _soil_bounds,
						  air_bounds = _air_bounds,
							  z_root = FT(-1*_rootdepth),
							z_canopy = 18.5,
							plant_hs = _tree_hs,
								  ga = 413.223,
								  la = 1735.537,
							latitude = 38.74,
						   longitude = -92.2,
						   elevation = 219.4,
					   stomata_model = sm);
	_node.swc = collect(FT,0.4*ones(length(_soil_bounds)-1));
	
	n_canopy = length(_node.plant_hs.canopy_index_in_air);
	
	if full_angles == 1
		incl_grid =  collect(FT,5:10:85);
		azi_grid = collect(FT,5:10:355);
		incl_grid_bnd = [collect(0:10:80) collect(FT,10:10:90)];
	else
		incl_grid =  collect(FT,[45]);
		azi_grid = collect(FT,[180]);
		incl_grid_bnd = [collect(FT,[0]) collect(FT,[90])];
	end;
	
	nleaf = length(incl_grid)*length(azi_grid) + 1;
	
	_node.plant_ps = [CanopyLayer{FT}(n_leaf=nleaf) for
										  i in 1:n_canopy];
	
	_node.canopy_rt = Canopy4RT{FT}(nLayer=n_canopy, LAI=_node.la/_node.ga, litab  = incl_grid, litab_bnd = incl_grid_bnd, lazitab = azi_grid);	

	#"RT dimensions"
	_node.rt_dim = create_rt_dims(_node.canopy_rt, _node.wl_set);
	#"CanopyRads container"
	_node.can_rad = create_canopy_rads(FT, _node.rt_dim);
	#"CanopyOpticals container"
	_node.can_opt = create_canopy_opticals(FT, _node.rt_dim);
	#"Array of LeafBios container"
	_node.leaves_rt = [create_leaf_bios(FT, _node.rt_dim) for
										 i in 1:n_canopy];

	#"RT container"
	_node.rt_con = create_rt_cache(FT, _node.rt_dim);
	#"Container for sunlit leaf area fraction in each layer"
	_node.f_SL = repeat(_node.canopy_rt.lidf, outer=[ _node.canopy_rt.nAzi ]) /
						_node.canopy_rt.nAzi;
	
	for _iPS in _node.plant_ps
		_iPS.g_min   = 0.001;
		_iPS.g_min25 = 0.001;
		_iPS.g_max   = 0.5; 
		_iPS.g_max25 = 0.5;
	end;
	
	_node.canopy_rt.Ω       = 0.69;
	_node.canopy_rt.clump_a = 0.69;
	
	update_LAI!(_node, FT(4.2));
	update_Weibull!(_node, FT(5.703), FT(0.953));

    # update fitted Vcmax, Kmax, and Chlrophyll content
    update_VJRWW!(_node, vcmax);

	#update_Kmax_ratio!(_node, kmax, convert(Array{FT},[0.25,8,0.25,0.25]));
	update_Kmax_ratio!(_node, kmax, convert(Array{FT},[2,8,4,2]));
	update_Cab!(_node, chl);
    initialize_spac_canopy!(_node);
	
    return _node
)
end
