"""

$(TYPEDEF)

Structure that stores hyperspectral radiation information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct HyperspectralRadiation{FT,DIM_WL}
    "Diffuse radiation `[mW m⁻² nm⁻¹]`"
    e_diffuse::SVector{DIM_WL,FT}
    "Direct radiation `[mW m⁻² nm⁻¹]`"
    e_direct::SVector{DIM_WL,FT}
end

HyperspectralRadiation{FT}(dset::String) where {FT} = (
    _dim_λ = size_nc(dset, "E_DIFF")[2][1];

    HyperspectralRadiation{FT,_dim_λ}(read_nc(dset, "E_DIFF"), read_nc(dset, "E_DIR"))
);
