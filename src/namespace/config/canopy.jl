"""

$(TYPEDEF)

Canopy structure.

# Fields

$(TYPEDFIELDS)

"""
struct CanopyStructure{FT,DIM_AZI,DIM_INCL}
    "Azimuth angles `[°]`"
    Θ_AZI::SVector{DIM_AZI,FT}
    "Inclination angles `[°]`"
    Θ_INCL::SVector{DIM_INCL,FT}
    "Bounds of inclination angles `[°]`"
    Θ_INCL_BNDS::NTuple{2,SVector{DIM_INCL,FT}}
end

CanopyStructure{FT}(; n_azi::Int = 360, n_incl::Int = 9) where {FT} = (
    _inc_min = collect(FT, range(0, 90; length=n_incl+1))[1:end-1];
    _inc_max = collect(FT, range(0, 90; length=n_incl+1))[2:end];
    _azis = collect(FT, range(0, 360; length=n_azi+1))[1:end-1] .+ 360 / n_azi / 2;
    _incs = [ (_inc_min[_i] + _inc_max[_i]) / 2 for _i in 1:n_incl ];

    return CanopyStructure{FT,n_azi,n_incl}(_azis, _incs, (_inc_min,_inc_max))
);

dims(::CanopyStructure{FT,DIM_AZI,DIM_INCL}) where {FT,DIM_AZI,DIM_INCL} = (DIM_AZI,DIM_INCL,);
