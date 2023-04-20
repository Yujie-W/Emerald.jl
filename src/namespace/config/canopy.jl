"""

$(TYPEDEF)

Canopy structure.

# Fields

$(TYPEDFIELDS)

"""
struct CanopyStructure{FT,VT_AZI,VT_INCL}
    "Hot spot parameter"
    HOT_SPOT::FT
    "Azimuth angles `[°]`"
    Θ_AZI::VT_AZI
    "Bounds of azimuth angles `[°]`"
    Θ_AZI_BNDS::NTuple{2,VT_AZI}
    "Inclination angles `[°]`"
    Θ_INCL::VT_INCL
    "Bounds of inclination angles `[°]`"
    Θ_INCL_BNDS::NTuple{2,VT_INCL}
end

CanopyStructure{FT}(; DIM_AZI::Int = 360, DIM_INCL::Int = 9, hotspot::Number = 0.05) where {FT} = (
    _inc_min = collect(FT, range(0, 90; length=DIM_INCL+1))[1:end-1];
    _inc_max = collect(FT, range(0, 90; length=DIM_INCL+1))[2:end];
    _azi_min = collect(FT, range(0, 360; length=DIM_AZI+1))[1:end-1];
    _azi_max = collect(FT, range(0, 360; length=DIM_AZI+1))[2:end];
    _azis = [ (_azi_min[_i] + _azi_max[_i]) / 2 for _i in 1:DIM_AZI ];
    _incs = [ (_inc_min[_i] + _inc_max[_i]) / 2 for _i in 1:DIM_INCL ];

    if USE_STATIC_ARRAY
        VT_AZI  = SVector{DIM_AZI,FT};
        VT_INCL = SVector{DIM_INCL,FT};
    else
        VT_AZI  = Vector{FT};
        VT_INCL = Vector{FT};
    end;

    return CanopyStructure{FT,VT_AZI,VT_INCL}(hotspot, _azis, (_azi_min,_azi_max), _incs, (_inc_min,_inc_max))
);

dim_azi(var::CanopyStructure) = length(var.Θ_AZI);
dim_incl(var::CanopyStructure) = length(var.Θ_INCL);
