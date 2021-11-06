mpa2mm = FT(10^6/9.8);

function update_pk_leaf!(tissue::LeafHydraulics{FT})
	tissue.p_storage = (tissue.v_storage ./ tissue.v_maximum .- FT(1)) * tissue.pv.slope
	tissue.p_element[1] = tissue.p_storage
	tissue.p_leaf = tissue.p_storage
	tissue.p_ups = tissue.p_storage[1]
	tissue.k_element = exp.(-(tissue.p_element ./ -tissue.vc.b) .^ tissue.vc.c) * tissue.k_sla
end
	
function update_pk_nonleaf!(tissue::AbstractHydraulicOrgan{FT})
	tissue.p_storage = (tissue.v_storage ./ tissue.v_maximum .- FT(1)) * tissue.pv.slope
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


function create_deriv_mat(plant_hs)

	n_canopy= plant_hs.n_canopy
	n_tissue = plant_hs.n_canopy*2 + 1 + plant_hs.n_root

	ans_lhs = zeros(n_tissue, n_tissue);
	ans_rhs = zeros(n_tissue)

	i = 1 #row index of plant organ in the output matrix
	for can_i in 1:plant_hs.n_canopy
		iLeaf = plant_hs.leaves[can_i]	
		iBranch = plant_hs.branch[can_i]
		
		leaf_k_total = iLeaf.k_element[1]*iLeaf.area
		leaf_vmax_total = iLeaf.v_maximum*iLeaf.area
		
		kbl = 2/(1/leaf_k_total + 1/iBranch.k_element[1])
		#kbl = leaf_k_total

		dzbl = iBranch.Δh/2 * 1000/mpa2mm;
		#dzbl = 0;
		
		leaf_coef_leaf_in = -kbl/iLeaf.pv.slope/leaf_vmax_total
		ans_lhs[i,i] += leaf_coef_leaf_in
		ans_lhs[i+n_canopy,i] += -leaf_coef_leaf_in
		
		branch_coef_leaf_in = kbl/iBranch.pv.slope/iBranch.v_maximum[1]

		ans_lhs[i,i+n_canopy] += branch_coef_leaf_in
		ans_lhs[i+n_canopy,i+n_canopy] += -branch_coef_leaf_in
		
		#constant terms in pressure difference flow into leaf
		rhs_term1 = kbl*(-1/iBranch.pv.slope + 1/iLeaf.pv.slope)
		#ans_rhs[i] += rhs_term1
		#ans_rhs[i+n_canopy] += -rhs_term1
		
		#flow due to gravity
		ans_rhs[i] += -dzbl*kbl
		ans_rhs[i+n_canopy] += dzbl*kbl
		
		#flow out of leaf due to ET
		ans_rhs[i] += -iLeaf.q_out*iLeaf.area
		
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
		dztb = (iBranch.Δh/2 + trunk.Δh/2) * 1000/mpa2mm
		#dztb = iBranch.Δh* 1000/mpa2mm
		
		branch_coef_branch_in = -ktb/iBranch.pv.slope/iBranch.v_maximum[1]
		ans_lhs[i,i] += branch_coef_branch_in
		ans_lhs[itrunk,i] += -branch_coef_branch_in
		
		trunk_coef_branch_in = ktb/trunk.pv.slope/trunk.v_maximum[1]
		ans_lhs[i,itrunk] += trunk_coef_branch_in
		ans_lhs[itrunk,itrunk] += -trunk_coef_branch_in
		
		rhs_term1 = ktb*(-1/trunk.pv.slope + 1/iBranch.pv.slope)
		#ans_rhs[itrunk] += rhs_term1
		#ans_rhs[i] += -rhs_term1
		
		ans_rhs[i] += -dztb*ktb
		ans_rhs[itrunk] += dztb*ktb

		i = i+1
		
	end

	m2 = ans_lhs*1;
	#println(ans_lhs)

	
	i = itrunk+1

	for root_i in 1:plant_hs.n_root

		iRoot = plant_hs.roots[root_i]
		krt = 2/(1/iRoot.k_element[1] + 1/trunk.k_element[1])
		#krt = trunk.k_element[1]

		dzrt = (iRoot.Δh/2 + trunk.Δh/2)*1000/mpa2mm
		#dzrt = trunk.Δh/2*1000/mpa2mm
		
		trunk_coef_trunk_in = -krt/trunk.pv.slope/trunk.v_maximum[1]
		ans_lhs[itrunk,itrunk] += trunk_coef_trunk_in
		ans_lhs[i,itrunk] += -trunk_coef_trunk_in
		
		root_coef_trunk_in = krt/iRoot.pv.slope/iRoot.v_maximum[1]
		ans_lhs[itrunk,i] += root_coef_trunk_in
		ans_lhs[i,i] += -root_coef_trunk_in
		
		ans_rhs[itrunk] += -dzrt*krt
		ans_rhs[i] += dzrt*krt
		
		rhs_term1 = krt*(-1/iRoot.pv.slope + 1/trunk.pv.slope)
		#ans_rhs[i] += rhs_term1
		#ans_rhs[itrunk] += -rhs_term1
		
		ksr = iRoot.k_element[1]
		dzsr = iRoot.Δh/2*1000/mpa2mm
		#dzsr = iRoot.Δh*1000/mpa2mm
		
		root_coef_root_in = -ksr/iRoot.pv.slope/iRoot.v_maximum[1]
		ans_lhs[i,i] += root_coef_root_in
		ans_rhs[i] += -dzsr*ksr + iRoot.p_ups*ksr + ksr/iRoot.pv.slope
		i = i+1

	end
	
	return ans_lhs, ans_rhs
end

using LinearAlgebra


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
function solve_odes_exp(deriv_mat, deriv_const, init_cond, deltaT)
	inv_deriv = inv(deriv_mat);
	cprime = -inv_deriv * deriv_const;
	mateI = exp(deltaT*deriv_mat)
	new_x = cprime + mateI * (init_cond-cprime)
	int_x = cprime*deltaT + (inv_deriv * mateI - inv_deriv) * (init_cond-cprime)
	return new_x, int_x
end
=#

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


function update_cap_mat!(plant_hs, deltaT)

	volumes = get_v_prof2(plant_hs);
	dmat, dvec = create_deriv_mat(plant_hs)

	newvals, newints = solve_odes(dmat, dvec, volumes, deltaT)
	#newvals, newints = solve_odes_exp(dmat, dvec, volumes, deltaT)
	
	set_vol!(plant_hs,newvals)

	for iroot in 1:plant_hs.n_root
		rootI = plant_hs.roots[iroot]
		m = rootI.pv.slope;
		ksr = rootI.k_element[1]
		ps = rootI.p_ups
		dzr = rootI.Δh*1000/mpa2mm
		vroot = rootI.v_maximum[1]
		integral_rootvol = newints[plant_hs.n_canopy*2+1+iroot];
		rootI.q_in = ksr*(ps*deltaT - m/vroot*integral_rootvol + m*deltaT - dzr*deltaT)/deltaT
	end
		
	update_pk_tree!(plant_hs)

	return newvals,newints
	
end



