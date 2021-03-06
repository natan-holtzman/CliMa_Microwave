#these functions are copied, and slightly modified, from the code included with the GMD paper

###############################################################################
#
# As the land model is modularized, there could be multiple places that conatin
#     exact the same information. Thus, these information needs to be synced.
#     This type of information includes leaf area, hydraulic traits, and
#     photosynthetic traits.
#
# The following function is meant for leaf area
#
###############################################################################
"""
    update_LAI!(node::SPACMono{FT}, lai::FT) where {FT<:AbstractFloat}
Update SPAC leaf area index
"""
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




"""
    update_Kmax!(node::SPACMono{FT}, kmax::FT) where {FT<:AbstractFloat}
Update the maximal hydraulic conductance for SPAC
TODO need to abstractize this to PlantHydraulics.jl
"""
function update_Kmax0!(node::SPACMono{FT}, kmax::FT) where {FT<:AbstractFloat}
    # Root:Stem:Leaf conductance ratio is set as 1:2:2.
    #    For trees, roots: 2X, trunk: 8X, branch: 8X, leaves: 4X
    #    For palms, roots: 2X, trunk: 4X, leaves: 4X
    #    For grasses, roots: 2X, leaves: 2X

    # 1. update the hydraulic conductance in Roots
    for root in node.plant_hs.roots
        root.k_max      = 2 * kmax / node.plant_hs.n_root;
        root.k_element .= root.k_max * root.N;
    end

    # 2. If is a tree
    if typeof(node.plant_hs) <: TreeLikeOrganism
        node.plant_hs.trunk.k_max      = 8 * kmax;
        node.plant_hs.trunk.k_element .= node.plant_hs.trunk.k_max *
                                         node.plant_hs.trunk.N;
        for stem in node.plant_hs.branch
            stem.k_max      = 8 * kmax / node.plant_hs.n_canopy;
            stem.k_element .= stem.k_max * stem.N;
        end
        for leaf in node.plant_hs.leaves
            leaf.k_sla      = 4 * kmax / node.la;
            leaf.k_element .= leaf.k_sla * leaf.N;
        end
    end

    # 3. If is a palm
    if typeof(node.plant_hs) <: PalmLikeOrganism
        node.plant_hs.trunk.k_max      = 4 * kmax;
        node.plant_hs.trunk.k_element .= node.plant_hs.trunk.k_max *
                                         node.plant_hs.trunk.N;
        for leaf in node.plant_hs.leaves
            leaf.k_sla      = 4 * kmax / node.la;
            leaf.k_element .= leaf.k_sla * leaf.N;
        end
    end

    # 4. If is a grass
    if typeof(node.plant_hs) <: GrassLikeOrganism
        for leaf in node.plant_hs.leaves
            leaf.k_sla      = 2 * kmax / node.la;
            leaf.k_element .= leaf.k_sla * leaf.N;
        end
    end

    return nothing
end

function update_Kmax_ratio!(node::SPACMono{FT}, kmax::FT, ratio::Array{FT}) where {FT<:AbstractFloat}
    # Root:Stem:Leaf conductance ratio is set as 1:2:2.
    #    For trees, roots: 2X, trunk: 8X, branch: 8X, leaves: 4X
    #    For palms, roots: 2X, trunk: 4X, leaves: 4X
    #    For grasses, roots: 2X, leaves: 2X

	resist_sum = sum(1 ./ ratio);
	ratio_norm = ratio * resist_sum;
	
    # 1. update the hydraulic conductance in Roots
    for root in node.plant_hs.roots
        root.k_max      = ratio_norm[1] * kmax / node.plant_hs.n_root;
        root.k_element .= root.k_max * root.N;
    end

    # 2. If is a tree
    if typeof(node.plant_hs) <: TreeLikeOrganism
        node.plant_hs.trunk.k_max      = ratio_norm[2] * kmax;
        node.plant_hs.trunk.k_element .= node.plant_hs.trunk.k_max *
                                         node.plant_hs.trunk.N;
        for stem in node.plant_hs.branch
            stem.k_max      = ratio_norm[3] * kmax / node.plant_hs.n_canopy;
            stem.k_element .= stem.k_max * stem.N;
        end
        for leaf in node.plant_hs.leaves
            leaf.k_sla      = ratio_norm[4] * kmax / node.plant_hs.n_canopy / node.la;
            leaf.k_element .= leaf.k_sla * leaf.N;
        end
    end

    return nothing
end



"""
    update_VJRWW!(node::SPACMono{FT}, vcmax::FT) where {FT<:AbstractFloat}
Update Vcmax25(WW), Jmax25(WW), and Rd25(WW) from a given Vcmax25
"""
function update_VJRWW!(node::SPACMono{FT}, vcmax::FT) where {FT<:AbstractFloat}
    # TODO change the ratio accordingly to photo_set
    # TODO add another ratio V2J in photo_set
    # Update Vcmax25, Jmax25 (1.67 Vcmax), and Rd25 (0.015 Vcmax)
    for _iPS in node.plant_ps
        _iPS.ps.Vcmax     = vcmax;
        _iPS.ps.Vcmax25   = vcmax;
        _iPS.ps.Vcmax25WW = vcmax;
        _iPS.ps.Jmax      = vcmax * 1.67;
        _iPS.ps.Jmax25    = vcmax * 1.67;
        _iPS.ps.Jmax25WW  = vcmax * 1.67;
        _iPS.ps.Rd        = vcmax * 0.015;
        _iPS.ps.Rd25      = vcmax * 0.015;
        _iPS.ps.Rd25WW    = vcmax * 0.015;
    end

    return nothing
end




"""
    update_VJR!(node::SPACMono{FT}, ratio::FT) where {FT<:AbstractFloat}
Update effective Vcmax25, Jmax25, and Rd25 from a given tune factor
"""
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


function update_VJR_leaf!(iPS::CanopyLayer{FT}, ratio::FT) where {FT<:AbstractFloat}
    # TODO change the ratio accordingly to photo_set
    # TODO add another ratio V2J in photo_set
    # Update Vcmax25, Jmax25 (1.67 Vcmax), and Rd25 (0.015 Vcmax)

	iPS.ps.Vcmax   = iPS.ps.Vcmax25WW * ratio;
	iPS.ps.Vcmax25 = iPS.ps.Vcmax25WW * ratio;
	iPS.ps.Jmax    = iPS.ps.Jmax25WW * ratio;
	iPS.ps.Jmax25  = iPS.ps.Jmax25WW * ratio;
	iPS.ps.Rd      = iPS.ps.Rd25WW * ratio;
	iPS.ps.Rd25    = iPS.ps.Rd25WW * ratio;

    return nothing
end





"""
    update_Cab!(node::SPACMono{FT}, cab::FT) where {FT<:AbstractFloat}
Update leaf chlorophyll content, and then rerun `fluspect!`
"""
function update_Cab!(node::SPACMono{FT}, cab::FT) where {FT<:AbstractFloat}
    for _leaf in node.leaves_rt
        _leaf.Cab = cab;
        fluspect!(_leaf, node.wl_set);
    end

    return nothing
end




"""
    update_Weibull!(node::SPACMono{FT}, b::FT, c::FT) where {FT<:AbstractFloat}
Update Weibull B and C for all components
"""
function update_Weibull!(node::SPACMono{FT}, b::FT) where {FT<:AbstractFloat}
    # 1. update B and C for roots
    for _root in node.plant_hs.roots
        _root.vc.b = b;
    end

    # 2. update B and C for trunk (if exists)
    if !(typeof(node.plant_hs) <: GrassLikeOrganism)
        node.plant_hs.trunk.vc.b = b;
    end

    # 3. update B and C for branches (if exist)
    if typeof(node.plant_hs) <: TreeLikeOrganism
        for _stem in node.plant_hs.branch
            _stem.vc.b = b;
        end
    end

    # 4. update B and C for leaves
    for _leaf in node.plant_hs.leaves
        _leaf.vc.b = b;
    end

    return nothing
end


function update_Weibull!(node::SPACMono{FT}, b::FT, c::FT) where {FT<:AbstractFloat}
    # 1. update B and C for roots
    for _root in node.plant_hs.roots
        _root.vc.b = b;
        _root.vc.c = c;
    end

    # 2. update B and C for trunk (if exists)
    if !(typeof(node.plant_hs) <: GrassLikeOrganism)
        node.plant_hs.trunk.vc.b = b;
        node.plant_hs.trunk.vc.c = c;
    end

    # 3. update B and C for branches (if exist)
    if typeof(node.plant_hs) <: TreeLikeOrganism
        for _stem in node.plant_hs.branch
            _stem.vc.b = b;
            _stem.vc.c = c;
        end
    end

    # 4. update B and C for leaves
    for _leaf in node.plant_hs.leaves
        _leaf.vc.b = b;
        _leaf.vc.c = c;
    end

    return nothing
end


function update_Kmax!(node::SPACMono{FT}, kmax::FT) where {FT<:AbstractFloat}
    # Root:Stem:Leaf conductance ratio is set as 1:2:2.
    #    For trees, roots: 2X, trunk: 8X, branch: 8X, leaves: 4X
    #    For palms, roots: 2X, trunk: 4X, leaves: 4X
    #    For grasses, roots: 2X, leaves: 2X

    # 1. update the hydraulic conductance in Roots
    for root in node.plant_hs.roots
        root.k_max      = 2 * kmax / node.plant_hs.n_root;
        root.k_element .= root.k_max * root.N;
    end

    # 2. If is a tree
    if typeof(node.plant_hs) <: TreeLikeOrganism
        node.plant_hs.trunk.k_max      = 8 * kmax;
        node.plant_hs.trunk.k_element .= node.plant_hs.trunk.k_max *
                                         node.plant_hs.trunk.N;
        for stem in node.plant_hs.branch
            stem.k_max      = 8 * kmax / node.plant_hs.n_canopy;
            stem.k_element .= stem.k_max * stem.N;
        end
        for leaf in node.plant_hs.leaves
            leaf.k_sla      = 4 * kmax / node.la;
            leaf.k_element .= leaf.k_sla * leaf.N;
        end
    end

    # 3. If is a palm
    if typeof(node.plant_hs) <: PalmLikeOrganism
        node.plant_hs.trunk.k_max      = 4 * kmax;
        node.plant_hs.trunk.k_element .= node.plant_hs.trunk.k_max *
                                         node.plant_hs.trunk.N;
        for leaf in node.plant_hs.leaves
            leaf.k_sla      = 4 * kmax / node.la;
            leaf.k_element .= leaf.k_sla * leaf.N;
        end
    end

    # 4. If is a grass
    if typeof(node.plant_hs) <: GrassLikeOrganism
        for leaf in node.plant_hs.leaves
            leaf.k_sla      = 2 * kmax / node.la;
            leaf.k_element .= leaf.k_sla * leaf.N;
        end
    end

    return nothing
end
