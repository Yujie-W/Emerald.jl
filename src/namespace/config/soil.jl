"""

$(TYPEDEF)

Soil albedo matrix

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct GSVSoilAlbedo{FT,DIM_GSV,DIM_WL,SIZE_GSV_WL}
    "A matrix of characteristic curves"
    MAT_ρ::SMatrix{DIM_WL,DIM_GSV,FT,SIZE_GSV_WL}
end

GSVSoilAlbedo{FT}(dset::String = LAND_2021; DIM_GSV::Int = 4) where {FT} = (
    DIM_WL = size_nc(dset, "GSV_1")[2][1];
    _mat_ρ = FT[read_nc(dset, "GSV_1") read_nc(dset, "GSV_2") read_nc(dset, "GSV_3") read_nc(dset, "GSV_4")];

    return GSVSoilAlbedo{FT,DIM_GSV,DIM_WL,DIM_GSV*DIM_WL}(_mat_ρ)
);

dims(::GSVSoilAlbedo{FT,DIM_GSV,DIM_WL}) where {FT,DIM_GSV,DIM_WL} = (DIM_GSV,DIM_WL,);
