using DifferentialEquations

mpa2mm = FT(10^6/9.8); #mPa of water pressure to mm of water height

#functions to update the pressure and conductivity of a plant, after its storage volumes have already been updated
function update_pk_leaf!(tissue::LeafHydraulics{FT})
	tissue.p_storage = (tissue.v_storage ./ tissue.v_maximum .- FT(1)) / tissue.pv.slope
	tissue.p_element[1] = tissue.p_storage
	tissue.p_leaf = tissue.p_storage
	tissue.p_ups = tissue.p_storage[1]
	tissue.k_element = exp.(-(tissue.p_element ./ -tissue.vc.b) .^ tissue.vc.c) * tissue.k_sla
end
	
function update_pk_nonleaf!(tissue::Land.PlantHydraulics.AbstractHydraulicOrgan{FT})
	tissue.p_storage = (tissue.v_storage ./ tissue.v_maximum .- FT(1)) / tissue.pv.slope
	tissue.p_element = tissue.p_storage
	tissue.k_element = exp.(-(tissue.p_element ./ -tissue.vc.b) .^ tissue.vc.c) * tissue.k_max
end
	
function update_pk_tree!(plant_hs)

	for can_i in 1:plant_hs.n_canopy
		iLeaf = plant_hs.leaves[can_i]		
		update_pk_leaf!(iLeaf)
			
		iBranch = plant_hs.branch[can_i]
		update_pk_nonleaf!(iBranch)
	end

	update_pk_nonleaf!(plant_hs.trunk)
	
	for root_i in 1:plant_hs.n_root
		iRoot = plant_hs.roots[root_i]
		update_pk_nonleaf!(iRoot)
	end

end

#conserve the total leaf water volume when changing LAI
#alternatively, this water could be removed or added to branches
function update_pv_leaf!(plant_hs, lai_old, lai_new)
	for can_i in 1:plant_hs.n_canopy
		iLeaf = plant_hs.leaves[can_i]	
		iLeaf.v_storage *= (lai_old / lai_new)
		update_pk_leaf!(iLeaf)
	end
end

#The capacitance model solves a set of ODEs at each timestep
#described by dv/dt = A*v + b
#this function fills in the matrix A and vector b
#using the state of the plant
#note that it is designed to work with a system where each component only has one vertical slice

function create_deriv_mat(plant_hs)

	n_canopy= plant_hs.n_canopy
	n_tissue = plant_hs.n_canopy*2 + 1 + plant_hs.n_root

	ans_lhs = zeros(n_tissue, n_tissue);
	ans_rhs = zeros(n_tissue)

	i = 1
	for can_i in 1:plant_hs.n_canopy
		iLeaf = plant_hs.leaves[can_i]	
		iBranch = plant_hs.branch[can_i]
		
		leaf_k_total = iLeaf.k_element[1]*iLeaf.area
		leaf_vmax_total = iLeaf.v_maximum*iLeaf.area
		
		kbl = 2/(1/leaf_k_total + 1/iBranch.k_element[1])
		#kbl = leaf_k_total

		
		dzbl = iBranch.??h/2 * 1000/mpa2mm;
		#dzbl = 0;
		
		leaf_coef_leaf_in = -kbl/iLeaf.pv.slope/leaf_vmax_total
		ans_lhs[i,i] += leaf_coef_leaf_in
		ans_lhs[i+n_canopy,i] += -leaf_coef_leaf_in
		
		branch_coef_leaf_in = kbl/iBranch.pv.slope/iBranch.v_maximum[1]

		ans_lhs[i,i+n_canopy] += branch_coef_leaf_in
		ans_lhs[i+n_canopy,i+n_canopy] += -branch_coef_leaf_in
		
		#constant terms in pressure difference flow into leaf
		rhs_term1 = kbl*(-1/iBranch.pv.slope + 1/iLeaf.pv.slope)
		ans_rhs[i] += rhs_term1
		ans_rhs[i+n_canopy] += -rhs_term1
		
		ans_rhs[i] += -dzbl*kbl - iLeaf.q_out*iLeaf.area
		ans_rhs[i+n_canopy] += dzbl*kbl
		
		i = i+1
	end
	
	m1= ans_lhs*1;
	#println(ans_lhs)

	trunk = plant_hs.trunk
	itrunk = plant_hs.n_canopy*2+1;

	for can_i in 1:plant_hs.n_canopy
		iBranch = plant_hs.branch[can_i]
		ktb = 2/(1/iBranch.k_element[1] + 1/trunk.k_element[1])
		#ktb = iBranch.k_element[1] 		
		dztb = (iBranch.??h/2 + trunk.??h/2) * 1000/mpa2mm
		#dztb = iBranch.??h* 1000/mpa2mm
		
		branch_coef_branch_in = -ktb/iBranch.pv.slope/iBranch.v_maximum[1]
		ans_lhs[i,i] += branch_coef_branch_in
		ans_lhs[itrunk,i] += -branch_coef_branch_in
		
		trunk_coef_branch_in = ktb/trunk.pv.slope/trunk.v_maximum[1]
		ans_lhs[i,itrunk] += trunk_coef_branch_in
		ans_lhs[itrunk,itrunk] += -trunk_coef_branch_in
		
		rhs_term1 = ktb*(-1/trunk.pv.slope + 1/iBranch.pv.slope)
		ans_rhs[i] += rhs_term1
		ans_rhs[itrunk] += -rhs_term1
		
		ans_rhs[i] += -dztb*ktb
		ans_rhs[itrunk] += dztb*ktb

		i = i+1
		
	end

	m2 = ans_lhs*1;
	#println(ans_lhs)

	
	i = itrunk+1
	
	total_k = sum([x.area/x.??h for x in plant_hs.roots]);

	for root_i in 1:plant_hs.n_root

		iRoot = plant_hs.roots[root_i]
		area_frac = iRoot.area/iRoot.??h/total_k  #total root area is 1		
		krt = 2/(1/(iRoot.k_element[1]*area_frac) + 1/trunk.k_element[1])
		#krt = trunk.k_element[1]

		dzrt = (iRoot.??h/2 + trunk.??h/2)*1000/mpa2mm
		#dzrt = trunk.??h/2*1000/mpa2mm
		
		trunk_coef_trunk_in = -krt/trunk.pv.slope/trunk.v_maximum[1]
		ans_lhs[itrunk,itrunk] += trunk_coef_trunk_in
		ans_lhs[i,itrunk] += -trunk_coef_trunk_in
		
		root_coef_trunk_in = krt/iRoot.pv.slope/iRoot.v_maximum[1]
		ans_lhs[itrunk,i] += root_coef_trunk_in
		ans_lhs[i,i] += -root_coef_trunk_in
		
		rhs_term1 = krt*(-1/iRoot.pv.slope + 1/trunk.pv.slope)
		ans_rhs[itrunk] += rhs_term1
		ans_rhs[i] += -rhs_term1
		
		ans_rhs[itrunk] += -dzrt*krt
		ans_rhs[i] += dzrt*krt
		
		ksr = iRoot.k_element[1]*area_frac
		dzsr = iRoot.??h/2*1000/mpa2mm
		#dzsr = iRoot.??h*1000/mpa2mm
		
		root_coef_root_in = -ksr/iRoot.pv.slope/iRoot.v_maximum[1]
		ans_lhs[i,i] += root_coef_root_in
		ans_rhs[i] += -dzsr*ksr + iRoot.p_ups*ksr + ksr/iRoot.pv.slope
		i = i+1

	end

	#println(ans_lhs)
	
	return ans_lhs, ans_rhs#, m1, m2


end

using LinearAlgebra

#this function solves the system of linear ODEs
#dx/dt = deriv_mat*x + deriv_const
#given the initial condition of x
#it returns the value of x after the time interval deltaT
#it also returns the integrals of the components of x over that time interval
function solve_odes(deriv_mat, deriv_const, init_cond, deltaT)

	e_val = eigvals(deriv_mat);
	e_vec = eigvecs(deriv_mat);

	inv_evec = inv(e_vec);
	inv_deriv = inv(deriv_mat);
	cprime = -inv_deriv * deriv_const;
	mateI = e_vec * diagm(exp.(e_val*deltaT)) * inv_evec
	new_x = cprime + mateI * (init_cond-cprime)
	int_x = cprime*deltaT + (inv_deriv * mateI - inv_deriv) * (init_cond-cprime)

	return new_x, int_x
end

#=
#slower way of solving ODEs
function solve_odes_exp(deriv_mat, deriv_const, init_cond, deltaT)
	inv_deriv = inv(deriv_mat);
	cprime = -inv_deriv * deriv_const;
	mateI = exp(deltaT*deriv_mat)
	new_x = cprime + mateI * (init_cond-cprime)
	int_x = cprime*deltaT + (inv_deriv * mateI - inv_deriv) * (init_cond-cprime)
	return new_x, int_x
end
=#

#this function takes a plant hydraulic system and returns the pressure of each leaf,branch,trunk,and root
#note that it is designed to work with a system where each component only has one vertical slice
function get_p_prof2(plant_hs)
	root_p = [x.p_element[1] for x in plant_hs.roots]
	trunk_p = plant_hs.trunk.p_element[1];
	branch_p = [x.p_element[1] for x in plant_hs.branch]
	leaf_p = [x.p_element[1] for x in plant_hs.leaves]
	
	ncan = plant_hs.n_canopy
	nroot = plant_hs.n_root
	n_tissue = ncan*2 + 1 + nroot

	p_list = vcat(leaf_p,branch_p,trunk_p,root_p)
	
	#=
	p_list = zeros(FT, n_tissue)
	p_list[1:ncan] = leaf_p;
	p_list[(ncan+1):(ncan*2)] = branch_p
	p_list[ncan*2+1] = trunk_p
	p_list[(ncan*2+2):end] = root_p
	=#
	
	return p_list
end

#same as previous function but returning current storage volume instead of pressure
function get_v_prof2(plant_hs)
	root_p = [x.v_storage[1] for x in plant_hs.roots]
	trunk_p = plant_hs.trunk.v_storage[1];
	branch_p = [x.v_storage[1] for x in plant_hs.branch]
	leaf_p = [x.v_storage[1]*x.area for x in plant_hs.leaves]

	ncan = plant_hs.n_canopy
	nroot = plant_hs.n_root
	n_tissue = ncan*2 + 1 + nroot

	p_list = vcat(leaf_p,branch_p,trunk_p,root_p)
	
	#=
	p_list = zeros(FT, n_tissue)
	p_list[1:ncan] = leaf_p;
	p_list[(ncan+1):(ncan*2)] = branch_p
	p_list[ncan*2+1] = trunk_p
	p_list[(ncan*2+2):end] = root_p
	=#
	return p_list
end

#sets volumes of each component to the values listed in newvols
function set_vol!(plant_hs, newvals)
	for ican in 1:plant_hs.n_canopy
		plant_hs.leaves[ican].v_storage = newvals[ican]/plant_hs.leaves[ican].area;
		plant_hs.branch[ican].v_storage[1] = newvals[ican+plant_hs.n_canopy];
	end
	
	plant_hs.trunk.v_storage[1] = newvals[plant_hs.n_canopy*2+1];
	
	for iroot in 1:plant_hs.n_root
		plant_hs.roots[iroot].v_storage[1] = newvals[plant_hs.n_canopy*2+1+iroot];
	end
end

#combines the above functions to advance the capacitance model by a time step
function update_cap_mat!(plant_hs, deltaT)

	volumes = get_v_prof2(plant_hs);
	dmat, dvec = create_deriv_mat(plant_hs)

	newvals, newints = solve_odes(dmat, dvec, volumes, deltaT)
	#newvals, newints = solve_odes_exp(dmat, dvec, volumes, deltaT)
	
	set_vol!(plant_hs,newvals)

	#compute the average value of root water uptake over the time step
	#for each root
	for iroot in 1:plant_hs.n_root
		rootI = plant_hs.roots[iroot]
		m = 1/rootI.pv.slope;
		ksr = rootI.k_element[1]*rootI.area #total root area is 1
		ps = rootI.p_ups
		dzr = rootI.??h*1000/mpa2mm
		vroot = rootI.v_maximum[1]
		integral_rootvol = newints[plant_hs.n_canopy*2+1+iroot];
		rootI.q_in = ksr*(ps*deltaT - m/vroot*integral_rootvol + m*deltaT - dzr*deltaT)/deltaT
	end
		
	update_pk_tree!(plant_hs)

	return newvals,newints
	
end

function find_new_cap(plant_hs, deltaT)

	volumes = get_v_prof2(plant_hs);
	dmat, dvec = create_deriv_mat(plant_hs)

	newvals, newints = solve_odes(dmat, dvec, volumes, deltaT)

	return newvals,newints
end

function push_cap!(plant_hs, newvals, newints)
	set_vol!(plant_hs,newvals)
	
	total_k = sum([x.area/x.??h for x in plant_hs.roots]);

	
	for iroot in 1:plant_hs.n_root
		rootI = plant_hs.roots[iroot]
		m = 1/rootI.pv.slope;
		ksr = rootI.k_element[1]*rootI.area/rootI.??h/total_k #total root area is 1
		ps = rootI.p_ups
		dzr = rootI.??h/2*1000/mpa2mm
		vroot = rootI.v_maximum[1]
		integral_rootvol = newints[plant_hs.n_canopy*2+1+iroot];
		
		
		rootI.q_in = ksr*(ps*deltaT - m/vroot*integral_rootvol + m*deltaT - dzr*deltaT)/deltaT
	end
	
	update_pk_tree!(plant_hs)
end


function get_v_max(plant_hs)
	root_p = [x.v_maximum[1] for x in plant_hs.roots]
	trunk_p = plant_hs.trunk.v_maximum[1];
	branch_p = [x.v_maximum[1] for x in plant_hs.branch]
	leaf_p = [x.v_maximum[1]*x.area for x in plant_hs.leaves]
	
	p_list = vcat(leaf_p,branch_p,trunk_p,root_p)
	
	return p_list
end

function soil_ode!(node::SPACMono{FT}, rain_in::FT, deltaT::FT, k_soil_max::FT, slope_index::FT)
	layer_thick = -1*diff(node.soil_bounds)*1000
	
	#fluxes_inout = -[x.q_in for x in node.plant_hs.roots] /KG_H_2_MOL_S * FT(1/node.ga*deltaT/(60*60));
	fluxes_inout = -[x.q_in for x in node.plant_hs.roots] *FT(18.02/1000)* FT(1/node.ga) #mol/s over basal to kg/s over ground;
	#node.swc[soil_layerI] = max(node.plant_hs.roots[i_root].sh.??r,min(node.swc[soil_layerI],node.plant_hs.roots[i_root].sh.??s));
	#fluxes_inout = zeros(FT,8);
	max_rain = (node.plant_hs.roots[1].sh.??s - node.swc[1])*layer_thick[1]/(60*60)
	rain_rate = rain_in/(60*60);
	#fluxes_inout[1] += rain_in/(60*60)*(node.swc[1] < node.plant_hs.roots[1].sh.??s) #mm/hour to mm/s;
	fluxes_inout[1] += min(rain_rate, max_rain)
	#fluxes_inout *= 0;

	wrc = node.plant_hs.roots[1].sh;
	midpoints = (node.soil_bounds[2:end] + node.soil_bounds[1:(end-1)]) /2 * 1000;
	layer_len = -1*diff(midpoints);

	function soil_derivative(smc_vec,p,t)
		k_soil = k_soil_max/(60*60) * [soil_k_ratio_swc(wrc, xi) for xi in smc_vec];
		p_soil_tot = [soil_p_25_swc(wrc,xi) for xi in smc_vec]*mpa2mm .+ midpoints;
		p_diff = -1*diff(p_soil_tot) #in mm
		q_down = k_soil[1:(end-1)] .* p_diff ./ layer_len * deltaT #mm/s * mm / mm * s
		free_drainage = k_soil[end]*slope_index*deltaT
		ans = fluxes_inout ./ layer_thick;
		ans[1:(end-1)] -= q_down ./ layer_thick[1:(end-1)];
		ans[2:end] += q_down ./ layer_thick[2:end];
		ans[end] -= free_drainage / layer_thick[end];

		#shortfall_to_sat = node.plant_hs.roots[1].sh.??s .- smc_vec;
		#only_negative = ans .* (ans .< 0);
		#ans[smc_vec .>= node.plant_hs.roots[1].sh.??s] = only_negative[smc_vec .>= node.plant_hs.roots[1].sh.??s]
		return ans
	end
	u0 = node.swc;
	tspan = (FT(0.0),deltaT)
	problem1 = ODEProblem(soil_derivative,u0,tspan);
	sol1 = solve(problem1);
	finalsol = sol1.u[end];
	runoff = 0;
	
	for i_soil in eachindex(node.swc)
		node.swc[i_soil] = max(min(finalsol[i_soil], node.plant_hs.roots[1].sh.??s),node.plant_hs.roots[1].sh.??r);
	end
	
	return finalsol #in units of mm
end	

