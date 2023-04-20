"""

$(TYPEDEF)

Structure for all configurations

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct CanopyRadiationCache{FT,DIM_AZI,DIM_CANOPY,DIM_CANOPY_1,DIM_INCL,DIM_WL,SIZE_AZI_INCL}
    # Leaf optical properties
    "Leaf reflectance, transmittance, and f_ppar"
    leaf_optics::SVector{DIM_WL,SVector{DIM_CANOPY,NTuple{3,FT}}}

    # # Canopy optical properties
    "Absolute value of fs (conversion factor for angles from solar at different inclination and azimuth angles)"
    abs_fs::SMatrix{DIM_INCL,DIM_AZI,FT,SIZE_AZI_INCL}
    "Extinction coefficients based on leaf angle distribution (_ks,_bf,_ci)"
    extinction_coefs::NTuple{3,FT}
    "Fraction of sunlit leaves"
    f_sunlit::SVector{DIM_CANOPY,FT}
    # "Canopy radiation scattering coefficients"
    # scatter_coefs::SMatrix{DIM_CANOPY,DIM_WL,NTuple{5,FT},SIZE_CANOPY_WL}
    "Effective longwave radiation reflectance based on level LAI"
    ρ_lw::SVector{DIM_CANOPY,FT}
    "Effective longwave radiation transmittance based on level LAI"
    τ_lw::SVector{DIM_CANOPY,FT}

    # # Effective optical properties after accounting for scattering
    "Effective longwave radiation reflectance based on level LAI after accounting for scattering"
    r_lw::SVector{DIM_CANOPY_1,FT}
    "Effective longwave radiation transmittance based on level LAI after accounting for scattering"
    t_lw::SVector{DIM_CANOPY,FT}
end
