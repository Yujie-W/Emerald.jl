"""

$(TYPEDEF)

Structure for all configurations

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct EmeraldConfiguration{FT,IT_NIR,IT_PAR,IT_SIF,IT_SIFE,VT_AZI,VT_INCL,VT_PAR,VT_SIF,VT_SIFE,VT_WL,MT_GSV_WL}
    # On/off features
    "Whether APAR absorbed by carotenoid is counted as PPAR"
    APAR_CAR::Bool = true
    "Whether to use hyperspectral soil albedo"
    HYPER_SOIL::Bool = false
    "Whether to use steady state plant hydraulic system"
    STEADY_STATE_HS::Bool = true
    "Whether to convert energy to photons when computing fluorescence"
    Φ_PHOTON::Bool = true

    # Embedded structures
    "Canopy structure"
    CAN::CanopyStructure{FT,VT_AZI,VT_INCL}
    "Universal constants"
    CONST::UniversalConstants{FT}
    "GSV soil albedo model matrix"
    GSV::GSVSoilAlbedo{FT,MT_GSV_WL}
    "Leaf hyperspectral absorption coefficients"
    LHA::HyperspectralAbsorption{FT,VT_WL}
    "Reference hyperspectral radiation profiles"
    RAD::HyperspectralRadiation{FT,VT_WL}
    "Wavelength sets"
    WLSET::WaveLengthSet{FT,IT_NIR,IT_PAR,IT_SIF,IT_SIFE,VT_PAR,VT_SIF,VT_SIFE,VT_WL}
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

    if USE_STATIC_ARRAY
        IT_NIR    = SVector{DIM_NIR,Int};
        IT_PAR    = SVector{DIM_PAR,Int};
        IT_SIF    = SVector{DIM_SIF,Int};
        IT_SIFE   = SVector{DIM_SIFE,Int};
        VT_AZI    = SVector{DIM_AZI,FT};
        VT_INCL   = SVector{DIM_INCL,FT};
        VT_PAR    = SVector{DIM_PAR,FT};
        VT_SIF    = SVector{DIM_SIF,FT};
        VT_SIFE   = SVector{DIM_SIFE,FT};
        VT_WL     = SVector{DIM_WL,FT};
        MT_GSV_WL = SMatrix{DIM_WL,DIM_GSV,FT,DIM_GSV*DIM_WL};
    else
        IT_NIR    = Vector{Int};
        IT_PAR    = Vector{Int};
        IT_SIF    = Vector{Int};
        IT_SIFE   = Vector{Int};
        VT_AZI    = Vector{FT};
        VT_INCL   = Vector{FT};
        VT_PAR    = Vector{FT};
        VT_SIF    = Vector{FT};
        VT_SIFE   = Vector{FT};
        VT_WL     = Vector{FT};
        MT_GSV_WL = Matrix{FT};
    end;

    return EmeraldConfiguration{FT,IT_NIR,IT_PAR,IT_SIF,IT_SIFE,VT_AZI,VT_INCL,VT_PAR,VT_SIF,VT_SIFE,VT_WL,MT_GSV_WL}(
                CAN   = _can,
                CONST = _cst,
                GSV   = _gsv,
                LHA   = _lha,
                RAD   = _rad,
                WLSET = _wls,
    )
);

dim_azi(var::EmeraldConfiguration) = dim_azi(var.CAN);
dim_gsv(var::EmeraldConfiguration) = dim_gsv(var.GSV);
dim_incl(var::EmeraldConfiguration) = dim_incl(var.CAN);
dim_nir(var::EmeraldConfiguration) = dim_nir(var.WLSET);
dim_par(var::EmeraldConfiguration) = dim_par(var.WLSET);
dim_sif(var::EmeraldConfiguration) = dim_sif(var.WLSET);
dim_sife(var::EmeraldConfiguration) = dim_sife(var.WLSET);
dim_wl(var::EmeraldConfiguration) = dim_wl(var.WLSET);
