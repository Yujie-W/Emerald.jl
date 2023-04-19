
function canopy_geometry(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT}) where {FT}
    DIM_CANOPY = length(state.chl);
    DIM_WL = length(config.WLSET.Λ);

    # compute leaf optical properties and save it to cache variables
    _cache_leaf = leaf_optical_properties.((config,), (state,), eachindex(state.chl), eachindex(config.WLSET.Λ)');

    return CanopyCache{FT,DIM_CANOPY,DIM_WL,DIM_CANOPY*DIM_WL}(
                leaf_optics = _cache_leaf,
    )
end
