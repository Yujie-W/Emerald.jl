"""

$(TYPEDEF)

Structure for all configurations

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct CanopyCache{FT,DIM_CANOPY,DIM_WL,SIZE_CANOPY_WL}
    "Leaf reflectance, transmittance, and f_ppar"
    leaf_optics::SMatrix{DIM_CANOPY,DIM_WL,NTuple{3,FT},SIZE_CANOPY_WL}
end
