"""

$(TYPEDEF)

Structure that stores hyperspectral radiation information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct HyperspectralRadiation{FT,VT_WL}
    "Diffuse radiation `[mW m⁻² nm⁻¹]`"
    E_DIF::VT_WL
    "Direct radiation `[mW m⁻² nm⁻¹]`"
    E_DIR::VT_WL
end

HyperspectralRadiation{FT}(dset::String = LAND_2021) where {FT} = (
    DIM_WL = size_nc(dset, "E_DIFF")[2][1];

    VT_WL = USE_STATIC_ARRAY ? SVector{DIM_WL,FT} : Vector{FT};

    HyperspectralRadiation{FT,VT_WL}(read_nc(dset, "E_DIFF"), read_nc(dset, "E_DIR"))
);

dim(var::HyperspectralRadiation) = length(var.E_DIF);
