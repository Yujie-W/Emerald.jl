"""

$(TYPEDEF)

Structure for all configurations

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct EmeraldConfiguration{FT,DIM_AXES,DIM_AZI,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}
    # On/off features
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool = true
    "Whether to use steady state plant hydraulic system"
    STEADY_STATE_HS::Bool = true
    "Whether to convert energy to photons when computing fluorescence"
    Φ_PHOTON::Bool = true

    # Embedded structures
    "Canopy structure"
    CAN::CanopyStructure{FT,DIM_AZI,DIM_INCL}
    "Universal constants"
    CONST::UniversalConstants{FT}
    "GSV soil albedo model matrix"
    GSV::GSVSoilAlbedo{FT,DIM_AXES,DIM_WL}
    "Leaf hyperspectral absorption coefficients"
    LHA::HyperspectralAbsorption{FT,DIM_WL}
    "Reference hyperspectral radiation profiles"
    RAD::HyperspectralRadiation{FT,DIM_WL}
    "Wavelength sets"
    WLSET::WaveLengthSet{FT,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}
end

"""

    EmeraldConfiguration{FT}(
                dset::String = LAND_2021;
                gsv_axes::Int = 4,
                n_azi::Int = 36,
                n_incl::Int = 9,
                wl_nir = (700, 2500),
                wl_par = (400, 750),
                wl_sif = (640, 850),
                wl_sife = (400, 750)
    ) where {FT}

Construct an EmeraldConfiguration struct, given
- `dset` Dataset to read the wavelength information
- `gsv_axes` Number of GSV model axes to use to derive hyperspectral soil albedo
- `n_azi` Number of azimuth angles
- `n_incl` Number of leaf inclination angles
- `wl_nir` NIR wavelength range (used for hyperspectral soil albedo)
- `wl_par` PAR wavelength range
- `wl_sif` SIF wavelength range
- `wl_sife` SIF excitation wavelength range

# Examples
```julia
using Emerald.EmeraldLand.NameSpace;
config = NameSpace.EmeraldConfiguration{Float64}();
```
"""
EmeraldConfiguration{FT}(
            dset::String = LAND_2021;
            gsv_axes::Int = 4,
            n_azi::Int = 36,
            n_incl::Int = 9,
            wl_nir = (700, 2500),
            wl_par = (400, 750),
            wl_sif = (640, 850),
            wl_sife = (400, 750)
) where {FT} = (
    _cst = UniversalConstants{FT}();
    _can = CanopyStructure{FT}(n_azi = n_azi, n_incl = n_incl);
    _gsv = GSVSoilAlbedo{FT}(dset; n = gsv_axes);
    _lha = HyperspectralAbsorption{FT}(dset);
    _rad = HyperspectralRadiation{FT}(dset);
    _wls = WaveLengthSet{FT}(dset; wl_nir = wl_nir, wl_par = wl_par, wl_sif = wl_sif, wl_sife = wl_sife);

    (DIM_AXES,DIM_WL) = dims(_gsv);
    (DIM_AZI,DIM_INCL) = dims(_can);
    (DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL) = dims(_wls);

    return EmeraldConfiguration{FT,DIM_AXES,DIM_AZI,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}(
                CAN   = _can,
                CONST = _cst,
                GSV   = _gsv,
                LHA   = _lha,
                RAD   = _rad,
                WLSET = _wls,
    )
);

dims(::EmeraldConfiguration{FT,DIM_AXES,DIM_AZI,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}) where {FT,DIM_AXES,DIM_AZI,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL} = (
    return DIM_AXES,DIM_AZI,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL
);
