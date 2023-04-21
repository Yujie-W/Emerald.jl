"""

$(TYPEDEF)

Structure for all canopy radiation cache variables

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct CanopyRadiationCache{FT,VT_CANOPY,VT_CANOPY_1,VT_WL_TVVVV,VT_WL_VT_CANOPY_T3,VT_WL_VT_CANOPY_T5,MT_AZI_INCL}
    # Leaf optical properties
    "Leaf reflectance, transmittance, and f_ppar"
    leaf_optics::VT_WL_VT_CANOPY_T3

    # # Canopy optical properties
    "Absolute value of fs (conversion factor for angles from solar at different inclination and azimuth angles)"
    abs_fs::MT_AZI_INCL
    "Extinction coefficients based on leaf angle distribution (_ks,_bf,_ci)"
    ext_coefs::NTuple{3,FT}
    "Fraction of sunlit leaves"
    f_sunlit::VT_CANOPY
    "Canopy shortwave coefficients based on level LAI"
    sw_coefs::VT_WL_VT_CANOPY_T5
    "Effective longwave radiation reflectance based on level LAI"
    ρ_lw::VT_CANOPY
    "Effective longwave radiation transmittance based on level LAI"
    τ_lw::VT_CANOPY

    # # Effective optical properties after accounting for scattering
    "Effective longwave radiation reflectance based on level LAI after accounting for scattering"
    eff_ρ_lw::VT_CANOPY_1
    "Effective longwave radiation transmittance based on level LAI after accounting for scattering"
    eff_τ_lw::VT_CANOPY
    "Effective shortwave radiation coefficients based on level LAI after accounting for scattering"
    eff_sw_coefs::VT_WL_TVVVV
end
