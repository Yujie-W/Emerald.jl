"""

$(TYPEDEF)

Immutable structure that stores wave length information.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct WaveLengthSet{FT,IT_NIR,IT_PAR,IT_SIF,IT_SIFE,VT_PAR,VT_SIF,VT_SIFE,VT_WL}
    # Constants
    "Wavelength (bins) `[nm]`"
    Λ::VT_WL
    "Lower boundary wavelength `[nm]`"
    Λ_LOWER::VT_WL
    "Upper boundary wavelength `[nm]`"
    Λ_UPPER::VT_WL

    # Indices
    "Indicies of Λ_NIR in Λ"
    IΛ_NIR::IT_NIR
    "Indicies of Λ_PAR in Λ"
    IΛ_PAR::IT_PAR
    "Indicies of Λ_SIF in Λ"
    IΛ_SIF::IT_SIF
    "Indicies of Λ_SIFE in Λ"
    IΛ_SIFE::IT_SIFE

    # Constants based on the ones above
    "Differential wavelength `[nm]`"
    ΔΛ::VT_WL = Λ_UPPER .- Λ_LOWER
    "Differential wavelength for PAR `[nm]`"
    ΔΛ_PAR::VT_PAR = ΔΛ[IΛ_PAR]
    "Differential wavelength for SIF excitation `[nm]`"
    ΔΛ_SIFE::VT_SIFE = ΔΛ[IΛ_SIFE]
    "Wavelength bins for PAR `[nm]`"
    Λ_PAR::VT_PAR = Λ[IΛ_PAR]
    "Wavelength bins for SIF `[nm]`"
    Λ_SIF::VT_SIF = Λ[IΛ_SIF]
    "Wavelength bins for SIF excitation `[nm]`"
    Λ_SIFE::VT_SIFE = Λ[IΛ_SIFE]
end

WaveLengthSet{FT}(dset::String = LAND_2021; wl_nir = (700, 2500), wl_par = (400, 750), wl_sif = (640, 850), wl_sife = (400, 750)) where {FT} =  (
    _λ = read_nc(dset, "WL");
    DIM_NIR  = length(findall(wl_nir[1]  .<= _λ .<= wl_nir[2]));
    DIM_PAR  = length(findall(wl_par[1]  .<= _λ .<= wl_par[2]));
    DIM_SIF  = length(findall(wl_sif[1]  .<= _λ .<= wl_sif[2]));
    DIM_SIFE = length(findall(wl_sife[1] .<= _λ .<= wl_sife[2]));
    DIM_WL   = length(_λ);

    if USE_STATIC_ARRAY
        IT_NIR  = SVector{DIM_NIR,Int};
        IT_PAR  = SVector{DIM_PAR,Int};
        IT_SIF  = SVector{DIM_SIF,Int};
        IT_SIFE = SVector{DIM_SIFE,Int};
        VT_PAR  = SVector{DIM_PAR,FT};
        VT_SIF  = SVector{DIM_SIF,FT};
        VT_SIFE = SVector{DIM_SIFE,FT};
        VT_WL   = SVector{DIM_WL,FT};
    else
        IT_NIR  = Vector{Int};
        IT_PAR  = Vector{Int};
        IT_SIF  = Vector{Int};
        IT_SIFE = Vector{Int};
        VT_PAR  = Vector{FT};
        VT_SIF  = Vector{FT};
        VT_SIFE = Vector{FT};
        VT_WL   = Vector{FT};
    end;

    return WaveLengthSet{FT,IT_NIR,IT_PAR,IT_SIF,IT_SIFE,VT_PAR,VT_SIF,VT_SIFE,VT_WL}(
                Λ       = _λ,
                Λ_LOWER = read_nc(dset, "WL_LOWER"),
                Λ_UPPER = read_nc(dset, "WL_UPPER"),
                IΛ_NIR  = findall(wl_nir[1]  .<= _λ .<= wl_nir[2]),
                IΛ_PAR  = findall(wl_par[1]  .<= _λ .<= wl_par[2]),
                IΛ_SIF  = findall(wl_sif[1]  .<= _λ .<= wl_sif[2]),
                IΛ_SIFE = findall(wl_sife[1] .<= _λ .<= wl_sife[2]),
    )
);

dim_nir(var::WaveLengthSet) = length(var.IΛ_NIR);
dim_par(var::WaveLengthSet) = length(var.IΛ_PAR);
dim_sif(var::WaveLengthSet) = length(var.IΛ_SIF);
dim_sife(var::WaveLengthSet) = length(var.IΛ_SIFE);
dim_wl(var::WaveLengthSet) = length(var.ΔΛ);
