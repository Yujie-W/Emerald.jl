"""

$(TYPEDEF)

Structure for all configurations

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct EmeraldConfiguration{FT,DIM_AXES,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}
    "Universal constants"
    CONST::UniversalConstants{FT}
    "GSV soil albedo model matrix"
    GSV::GSVSoilAlbedo{FT,DIM_AXES,DIM_WL}
    "Leaf hyperspectral absorption coefficients"
    LHA::HyperspectralAbsorption{FT,DIM_WL}
    "Reference hyperspectral radiation profiles"
    RAD_REF::HyperspectralRadiation{FT,DIM_WL}
    "Wavelength sets"
    WLSET::WaveLengthSet{FT,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}
end

EmeraldConfiguration{FT}(dset::String; gsv_axes::Int = 4, wl_nir = (700, 2500), wl_par = (400, 750), wl_sif = (640, 850), wl_sife = (400, 750)) where {FT} = (
    _cst = UniversalConstants{FT}();
    _gsv = GSVSoilAlbedo{FT}(dset, gsv_axes);
    _lha = HyperspectralAbsorption{FT}(dset);
    _rad = HyperspectralRadiation{FT}(dset);
    _wls = WaveLengthSet{FT}(dset; wl_nir = wl_nir, wl_par = wl_par, wl_sif = wl_sif, wl_sife = wl_sife);

    (DIM_AXES,DIM_WL) = dims(_gsv);
    (DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL) = dims(_wls);

    return EmeraldConfiguration{FT,DIM_AXES,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}(
                CONST   = _cst,
                GSV     = _gsv,
                LHA     = _lha,
                RAD_REF = _rad,
                WLSET   = _wls,
    )
);
