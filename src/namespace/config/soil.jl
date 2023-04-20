"""

$(TYPEDEF)

Soil albedo matrix

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct GSVSoilAlbedo{FT,MT_GSV_WL}
    "A matrix of characteristic curves"
    MAT_ρ::MT_GSV_WL
end

GSVSoilAlbedo{FT}(dset::String = LAND_2021; DIM_GSV::Int = 4) where {FT} = (
    DIM_WL = size_nc(dset, "GSV_1")[2][1];

    MT_GSV_WL = USE_STATIC_ARRAY ? SMatrix{DIM_WL,DIM_GSV,FT,DIM_GSV*DIM_WL} : Matrix{FT};

    _mat_ρ = FT[read_nc(dset, "GSV_1") read_nc(dset, "GSV_2") read_nc(dset, "GSV_3") read_nc(dset, "GSV_4")];

    return GSVSoilAlbedo{FT,MT_GSV_WL}(_mat_ρ)
);

dim_gsv(var::GSVSoilAlbedo) = size(var.MT_GSV_WL,2);
dim_wl(var::GSVSoilAlbedo) = size(var.MT_GSV_WL,1);
