"""

$(TYPEDEF)

Structure for all configurations

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct CanopyRadiationCache{FT,DIM_AZI,DIM_CANOPY,DIM_INCL,DIM_WL,SIZE_AZI_INCL,SIZE_CANOPY_WL}
    # Leaf optics
    "Leaf reflectance, transmittance, and f_ppar"
    leaf_optics::SMatrix{DIM_CANOPY,DIM_WL,NTuple{3,FT},SIZE_CANOPY_WL}

    # Canopy geometry
    "Absolute value of fs"
    abs_fs::SMatrix{DIM_INCL,DIM_AZI,FT,SIZE_AZI_INCL}
    "Extinction coefficients based on leaf angle distribution (_ks,_bf,_ci)"
    extinction_coefs::NTuple{3,FT}
    "Fraction of sunlit leaves"
    f_sunlit::SVector{DIM_CANOPY,FT}
end
