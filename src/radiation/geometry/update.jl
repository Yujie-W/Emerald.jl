
function canopy_radiation_cache(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT,DIM_CANOPY}, sza::FT, p_incl::SVector{DIM_INCL,FT}) where {FT,DIM_CANOPY,DIM_INCL}
    DIM_AZI = dim_azi(config);
    DIM_WL = dim_wl(config);

    # compute leaf optical properties and save it to cache variables
    _cache_leaf = leaf_optical_properties.((config,), (state,), eachindex(state.chl), eachindex(config.WLSET.Λ)');
    (_coefs, _abs_fs, _f_sunlit) = sun_geometry(config, state, sza, p_incl);

    return CanopyRadiationCache{FT,DIM_AZI,DIM_CANOPY,DIM_INCL,DIM_WL,DIM_AZI*DIM_INCL,DIM_CANOPY*DIM_WL}(
                leaf_optics      = _cache_leaf,
                abs_fs           = _abs_fs,
                extinction_coefs = _coefs,
                f_sunlit         = _f_sunlit,
    )
end
