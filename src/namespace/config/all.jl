"""

$(TYPEDEF)

Structure for all configurations

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct EmeraldConfiguration{FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL,SIZE_GSV_WL}
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
    GSV::GSVSoilAlbedo{FT,DIM_GSV,DIM_WL,SIZE_GSV_WL}
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
                DIM_AZI::Int = 36,
                DIM_GSV::Int = 4,
                DIM_INCL::Int = 9,
                wl_nir = (700, 2500),
                wl_par = (400, 750),
                wl_sif = (640, 850),
                wl_sife = (400, 750)
    ) where {FT}

Construct an EmeraldConfiguration struct, given
- `dset` Dataset to read the wavelength information
- `DIM_AZI` Number of azimuth angles
- `DIM_GSV` Number of GSV model axes to use to derive hyperspectral soil albedo
- `DIM_INCL` Number of leaf inclination angles
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
            DIM_AZI::Int = 36,
            DIM_GSV::Int = 4,
            DIM_INCL::Int = 9,
            wl_nir = (700, 2500),
            wl_par = (400, 750),
            wl_sif = (640, 850),
            wl_sife = (400, 750)
) where {FT} = (
    _cst = UniversalConstants{FT}();
    _can = CanopyStructure{FT}(DIM_AZI = DIM_AZI, DIM_INCL = DIM_INCL);
    _gsv = GSVSoilAlbedo{FT}(dset; DIM_GSV = DIM_GSV);
    _lha = HyperspectralAbsorption{FT}(dset);
    _rad = HyperspectralRadiation{FT}(dset);
    _wls = WaveLengthSet{FT}(dset; wl_nir = wl_nir, wl_par = wl_par, wl_sif = wl_sif, wl_sife = wl_sife);

    DIM_NIR = dim_nir(_wls);
    DIM_PAR = dim_par(_wls);
    DIM_SIF = dim_sif(_wls);
    DIM_SIFE = dim_sife(_wls);
    DIM_WL = dim_wl(_wls);

    return EmeraldConfiguration{FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL,DIM_GSV*DIM_WL}(
                CAN   = _can,
                CONST = _cst,
                GSV   = _gsv,
                LHA   = _lha,
                RAD   = _rad,
                WLSET = _wls,
    )
);

dim_azi(::EmeraldConfiguration{FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}) where {FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL} = DIM_AZI;
dim_gsv(::EmeraldConfiguration{FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}) where {FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL} = DIM_GSV;
dim_incl(::EmeraldConfiguration{FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}) where {FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL} = DIM_INCL;
dim_nir(::EmeraldConfiguration{FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}) where {FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL} = DIM_NIR;
dim_par(::EmeraldConfiguration{FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}) where {FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL} = DIM_PAR;
dim_sif(::EmeraldConfiguration{FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}) where {FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL} = DIM_SIF;
dim_sife(::EmeraldConfiguration{FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}) where {FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL} = DIM_SIFE;
dim_wl(::EmeraldConfiguration{FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}) where {FT,DIM_AZI,DIM_GSV,DIM_INCL,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL} = DIM_WL;
