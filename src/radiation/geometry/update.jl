
function canopy_radiation_cache(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT}, sza::FT, p_incl::SVector{DIM_INCL,FT}) where {FT,DIM_INCL}
    DIM_AZI = dim_azi(config);
    DIM_CANOPY = dim_canopy(state);
    DIM_WL = dim_wl(config);

    # compute leaf optical properties and save it to cache variables
    _leaf_optics = leaf_optical_properties(config, state);
    _coefs, _abs_fs, _f_sunlit = sun_geometry(config, state, sza, p_incl);
    _ρ_lw, _τ_lw = longwave_coefs(config, state, _coefs);
    _eff_ρ_lw, _eff_τ_lw = effective_longwave_coefs(config, _ρ_lw, _τ_lw);
    _coef_sw = shortwave_coefs(_leaf_optics, state.δlai, _coefs);




    # compute soil sw albedo
    _soil_albedo = USE_STATIC_ARRAY ? SVector{DIM_WL,FT}(ones(FT,DIM_WL) .* FT(0.1)) : ones(FT,DIM_WL) .* FT(0.1);
    _eff_coef_sw = effective_shortwave_coefs(_coef_sw, _soil_albedo);




    if USE_STATIC_ARRAY
        VT_CANOPY          = SVector{DIM_CANOPY,FT};
        VT_CANOPY_1        = SVector{DIM_CANOPY+1,FT};
        VT_WL_TVVVV        = SVector{DIM_WL,Tuple{VT_CANOPY_1,VT_CANOPY_1,VT_CANOPY,VT_CANOPY}};
        VT_WL_VT_CANOPY_T3 = SVector{DIM_WL,SVector{DIM_CANOPY,NTuple{3,FT}}};
        VT_WL_VT_CANOPY_T5 = SVector{DIM_WL,SVector{DIM_CANOPY,NTuple{5,FT}}};
        MT_AZI_INCL        = SMatrix{DIM_INCL,DIM_AZI,FT,DIM_AZI*DIM_INCL};
    else
        VT_CANOPY          = Vector{FT};
        VT_CANOPY_1        = Vector{FT};
        VT_WL_TVVVV        = Vector{Tuple{VT_CANOPY_1,VT_CANOPY_1,VT_CANOPY,VT_CANOPY}};
        VT_WL_VT_CANOPY_T3 = Vector{Vector{NTuple{3,FT}}};
        VT_WL_VT_CANOPY_T5 = Vector{Vector{NTuple{5,FT}}};
        MT_AZI_INCL        = Matrix{FT};
    end;

    return CanopyRadiationCache{FT,VT_CANOPY,VT_CANOPY_1,VT_WL_TVVVV,VT_WL_VT_CANOPY_T3,VT_WL_VT_CANOPY_T5,MT_AZI_INCL}(
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
