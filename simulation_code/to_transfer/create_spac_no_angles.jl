include("updatefuns.jl")
include("create_hydraulics_layered.jl")

function create_spac(
            sm::AbstractStomatalModel = OSMWang{FT}(),
            vcmax::FT = FT(30),
            kmax::FT = FT(0.3),
			psi_sat::FT = FT(0.014),
			b_soil::FT = FT(5.3),
		z_soil::FT = FT(0.7),  
	porosity::FT = FT(0.45),     
     chl::FT = FT(40),
	 N_slice::Int = 1
) where {FT<:AbstractFloat} 
(

	#_soil_hs = create_soil_VC(VanGenuchten{FT}(), "Silt Loam");
	
	
	_soil_hs = VanGenuchten{FT}(stype = "Ozark",
									α = 1.368,
									n = 2.6257,
								   Θs = 0.45,
								   Θr = 0.067);
						   
								   
	#_soil_hs = BrooksCorey{FT}(stype= "Ozark", Θs = FT(porosity), Θr = FT(0.05), ϕs = psi_sat, b = b_soil);

	_totaldepth = FT(z_soil/1000);
	_rootdepth = FT(z_soil/1000);
	_soil_bounds = collect(FT,[0,-0.1,-0.25,-0.45,-z_soil/1000]);
	#_soil_bounds = collect(FT,[0,-0.1,-0.35,-0.7]);
	
	#_totaldepth = FT(0.85);
	#_rootdepth = FT(0.85);
	#_soil_bounds = collect(FT,[0,-0.1,-0.3,-0.55,-0.85]);

	_tree_hs = create_tree2(FT(-1*_rootdepth), FT(9), FT(18.5), _soil_bounds,
						   collect(FT,0:5:20), N_slice);
	for _root in _tree_hs.roots
		_root.sh = deepcopy(_soil_hs);
	end;

	_node = SPACMono{FT}(soil_bounds = _soil_bounds,
						  air_bounds = collect(FT,0:5:20),
							  z_root = FT(-1*_rootdepth),
							z_canopy = 18.5,
							plant_hs = _tree_hs,
								  ga = 413.223,
								  la = 1735.537,
							latitude = 38.74,
						   longitude = -92.2,
						   elevation = 219.4,
					   stomata_model = sm);
	_node.swc = collect(FT,[0.4,0.4,0.4,0.4]);
	
	n_canopy = length(_node.plant_hs.canopy_index_in_air);
	
	#=
	incl_grid =  collect(FT,[22.5,45+22.5]);
	azi_grid = collect(FT,[90,270]);
	incl_grid_bnd = [collect(FT,[0,45]) collect(FT,[45,90])];
	=#
	
	incl_grid =  collect(FT,[45]);
	azi_grid = collect(FT,[180]);
	incl_grid_bnd = [collect(FT,[0]) collect(FT,[90])];
	
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
		_iPS.g_max   = 0.08; #0.50;
		_iPS.g_max25 = 0.08; #0.50;
	end;
	
	
	_node.canopy_rt.Ω       = 0.69;
	_node.canopy_rt.clump_a = 0.69;
	update_LAI!(_node, FT(4.2));
	update_Weibull!(_node, FT(5.703), FT(0.953));

	#update_Weibull!(_node, FT(4.0), FT(3.0));
	#for _leaf in _node.plant_hs.leaves
        #_leaf.vc.b = 1.0;
        #_leaf.vc.c = 3.0;
    #end;
	
	
	
	
	#for _leaf in _node.plant_hs.leaves
    #    _leaf.vc.b = 2.0;
        #_leaf.vc.c = c;
    #end;
	


    # update fitted Vcmax, Kmax, and Chlrophyll content
    update_VJRWW!(_node, vcmax);
    update_Kmax!(_node, kmax);
    update_Cab!(_node, chl);
    initialize_spac_canopy!(_node);
	
    return _node
)
end
