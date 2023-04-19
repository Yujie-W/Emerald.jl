"""

$(TYPEDEF)

Canopy structure.

# Fields

$(TYPEDFIELDS)

"""
struct CanopyStructure{FT,DIM_AZI,DIM_INCL}
    "Hot spot parameter"
    HOT_SPOT::FT
    "Azimuth angles `[°]`"
    Θ_AZI::SVector{DIM_AZI,FT}
    "Bounds of azimuth angles `[°]`"
    Θ_AZI_BNDS::NTuple{2,SVector{DIM_AZI,FT}}
    "Inclination angles `[°]`"
    Θ_INCL::SVector{DIM_INCL,FT}
    "Bounds of inclination angles `[°]`"
    Θ_INCL_BNDS::NTuple{2,SVector{DIM_INCL,FT}}
end

CanopyStructure{FT}(; hotspot::Number = 0.05, n_azi::Int = 360, n_incl::Int = 9) where {FT} = (
    _inc_min = collect(FT, range(0, 90; length=n_incl+1))[1:end-1];
    _inc_max = collect(FT, range(0, 90; length=n_incl+1))[2:end];
    _azi_min = collect(FT, range(0, 360; length=n_azi+1))[1:end-1];
    _azi_max = collect(FT, range(0, 360; length=n_azi+1))[2:end];
    _azis = [ (_azi_min[_i] + _azi_max[_i]) / 2 for _i in 1:n_azi ];
    _incs = [ (_inc_min[_i] + _inc_max[_i]) / 2 for _i in 1:n_incl ];

    return CanopyStructure{FT,n_azi,n_incl}(hotspot, _azis, (_azi_min,_azi_max), _incs, (_inc_min,_inc_max))
);

dims(::CanopyStructure{FT,DIM_AZI,DIM_INCL}) where {FT,DIM_AZI,DIM_INCL} = (DIM_AZI,DIM_INCL,);
