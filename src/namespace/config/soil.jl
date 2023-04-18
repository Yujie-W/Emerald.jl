"""

$(TYPEDEF)

Soil albedo matrix

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct GSVSoilAlbedo{FT,DIM_AXES,DIM_WL}
    "A matrix of characteristic curves"
    MAT_ρ::NTuple{DIM_AXES,SVector{DIM_WL,FT}}
end

GSVSoilAlbedo{FT}(dset::String, n::Int = 4) where {FT} = (
    _dim_λ = size_nc(dset, "GSV_1")[2][1];
    _mat_ρ = FT[read_nc(dset, "GSV_1") read_nc(dset, "GSV_2") read_nc(dset, "GSV_3") read_nc(dset, "GSV_4")];
    _svecs = [SVector{_dim_λ,FT}(_mat_ρ[:,_i]) for _i in 1:n];

    return GSVSoilAlbedo{FT,n,_dim_λ}(NTuple{n,SVector{_dim_λ,FT}}(_svecs))
);
