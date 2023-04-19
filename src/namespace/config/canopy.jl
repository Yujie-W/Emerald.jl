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

CanopyStructure{FT}(; hotspot::Number = 0.05, DIM_AZI::Int = 360, DIM_INCL::Int = 9) where {FT} = (
    _inc_min = collect(FT, range(0, 90; length=DIM_INCL+1))[1:end-1];
    _inc_max = collect(FT, range(0, 90; length=DIM_INCL+1))[2:end];
    _azi_min = collect(FT, range(0, 360; length=DIM_AZI+1))[1:end-1];
    _azi_max = collect(FT, range(0, 360; length=DIM_AZI+1))[2:end];
    _azis = [ (_azi_min[_i] + _azi_max[_i]) / 2 for _i in 1:DIM_AZI ];
    _incs = [ (_inc_min[_i] + _inc_max[_i]) / 2 for _i in 1:DIM_INCL ];

    return CanopyStructure{FT,DIM_AZI,DIM_INCL}(hotspot, _azis, (_azi_min,_azi_max), _incs, (_inc_min,_inc_max))
);

dim_azi(::CanopyStructure{FT,DIM_AZI,DIM_INCL}) where {FT,DIM_AZI,DIM_INCL} = DIM_AZI;
dim_incl(::CanopyStructure{FT,DIM_AZI,DIM_INCL}) where {FT,DIM_AZI,DIM_INCL} = DIM_INCL;
