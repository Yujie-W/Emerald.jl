
function canopy_radiation_cache(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT,DIM_CANOPY}, sza::FT, p_incl::SVector{DIM_INCL,FT}) where {FT,DIM_CANOPY,DIM_INCL}
    DIM_AZI = dim_azi(config);
    DIM_WL = dim_wl(config);

    # compute leaf optical properties and save it to cache variables
    _leaf_optics = leaf_optical_properties(config, state);
    _coefs, _abs_fs, _f_sunlit = sun_geometry(config, state, sza, p_incl);
    _ρ_lw, _τ_lw = longwave_coefs(config, state, _coefs);
    _eff_ρ_lw, _eff_τ_lw = effective_longwave_coefs(config, _ρ_lw, _τ_lw);
    _coef_sw = shortwave_coefs(_leaf_optics, state.δlai, _coefs);




    # compute soil sw albedo
    _eff_coef_sw = effective_shortwave_coefs(_coef_sw, SVector{DIM_WL,FT}(ones(FT,DIM_WL) .* FT(0.1)));




    return CanopyRadiationCache{FT,DIM_AZI,DIM_CANOPY,DIM_CANOPY+1,DIM_INCL,DIM_WL,DIM_AZI*DIM_INCL}(
                leaf_optics  = _leaf_optics,
                abs_fs       = _abs_fs,
                ext_coefs    = _coefs,
                f_sunlit     = _f_sunlit,
                sw_coefs     = _coef_sw,
                ρ_lw         = _ρ_lw,
                τ_lw         = _τ_lw,
                eff_ρ_lw     = _eff_ρ_lw,
                eff_τ_lw     = _eff_τ_lw,
                eff_sw_coefs = _eff_coef_sw,
    )
end
