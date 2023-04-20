
function canopy_radiation_cache(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT,DIM_CANOPY}, sza::FT, p_incl::SVector{DIM_INCL,FT}) where {FT,DIM_CANOPY,DIM_INCL}
    DIM_AZI = dim_azi(config);
    DIM_WL = dim_wl(config);

    # compute leaf optical properties and save it to cache variables
    _leaf_optics = leaf_optical_properties(config, state);
    (_coefs, _abs_fs, _f_sunlit) = sun_geometry(config, state, sza, p_incl);
    (_ρ_lw, _τ_lw) = longwave_coefs(config, state, _coefs);
    (_r_lw, _t_lw) = effective_longwave_coefs(config, _ρ_lw, _τ_lw);
    # _sun_scatter = sun_scatter_coefs.(state.δlai, _leaf_optics, (_coefs,));




    # compute soil sw albedo
    # effective_sun_scatter_coefs(view(_sun_scatter,:,1), FT(0.1));




    return CanopyRadiationCache{FT,DIM_AZI,DIM_CANOPY,DIM_CANOPY+1,DIM_INCL,DIM_WL,DIM_AZI*DIM_INCL}(
                leaf_optics      = _leaf_optics,
                abs_fs           = _abs_fs,
                extinction_coefs = _coefs,
                f_sunlit         = _f_sunlit,
                #scatter_coefs    = _sun_scatter,
                ρ_lw             = _ρ_lw,
                τ_lw             = _τ_lw,
                r_lw             = _r_lw,
                t_lw             = _t_lw,
    )
end
