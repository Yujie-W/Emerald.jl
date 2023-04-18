"""

$(TYPEDEF)

Immutable structure that stores wave length information.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct WaveLengthSet{FT,DIM_NIR,DIM_PAR,DIM_SIF,DIM_SIFE,DIM_WL}
    # Constants
    "Wavelength (bins) `[nm]`"
    Λ::SVector{DIM_WL,FT}
    "Lower boundary wavelength `[nm]`"
    Λ_LOWER::SVector{DIM_WL,FT}
    "Upper boundary wavelength `[nm]`"
    Λ_UPPER::SVector{DIM_WL,FT}

    # Indices
    "Indicies of Λ_NIR in Λ"
    IΛ_NIR::SVector{DIM_NIR,Int}
    "Indicies of Λ_PAR in Λ"
    IΛ_PAR::SVector{DIM_PAR,Int}
    "Indicies of Λ_SIF in Λ"
    IΛ_SIF::SVector{DIM_SIF,Int}
    "Indicies of Λ_SIFE in Λ"
    IΛ_SIFE::SVector{DIM_SIFE,Int}

    # Constants based on the ones above
    "Differential wavelength `[nm]`"
    ΔΛ::SVector{DIM_WL,FT} = Λ_UPPER .- Λ_LOWER
    "Differential wavelength for PAR `[nm]`"
    ΔΛ_PAR::SVector{DIM_PAR,FT} = ΔΛ[IΛ_PAR]
    "Differential wavelength for SIF excitation `[nm]`"
    ΔΛ_SIFE::SVector{DIM_SIFE,FT} = ΔΛ[IΛ_SIFE]
    "Wavelength bins for PAR `[nm]`"
    Λ_PAR::SVector{DIM_PAR,FT} = Λ[IΛ_PAR]
    "Wavelength bins for SIF `[nm]`"
    Λ_SIF::SVector{DIM_SIF,FT} = Λ[IΛ_SIF]
    "Wavelength bins for SIF excitation `[nm]`"
    Λ_SIFE::SVector{DIM_SIFE,FT} = Λ[IΛ_SIFE]
end

WaveLengthSet{FT}(dset::String; wl_nir = (700, 2500), wl_par = (400, 750), wl_sif = (640, 850), wl_sife = (400, 750)) where {FT} =  (
    _λ = read_nc(dset, "WL");
    _dim_nir  = length(findall(wl_nir[1]  .<= _λ .<= wl_nir[2]));
    _dim_par  = length(findall(wl_par[1]  .<= _λ .<= wl_par[2]));
    _dim_sif  = length(findall(wl_sif[1]  .<= _λ .<= wl_sif[2]));
    _dim_sife = length(findall(wl_sife[1] .<= _λ .<= wl_sife[2]));
    _dim_wl   = length(_λ);

    return WaveLengthSet{FT,_dim_nir,_dim_par,_dim_sif,_dim_sife,_dim_wl}(
                Λ       = _λ,
                Λ_LOWER = read_nc(dset, "WL_LOWER"),
                Λ_UPPER = read_nc(dset, "WL_UPPER"),
                IΛ_NIR  = findall(wl_nir[1]  .<= _λ .<= wl_nir[2]),
                IΛ_PAR  = findall(wl_par[1]  .<= _λ .<= wl_par[2]),
                IΛ_SIF  = findall(wl_sif[1]  .<= _λ .<= wl_sif[2]),
                IΛ_SIFE = findall(wl_sife[1] .<= _λ .<= wl_sife[2]),
    )
);
