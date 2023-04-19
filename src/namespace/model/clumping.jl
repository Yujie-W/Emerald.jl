"""

Abstract type for clumping index models

"""
abstract type AbstractClumpingIndexModel{FT} end


"""

$(TYPEDEF)

Immutable structure that stores wave length information.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct ClumpingIndexPinty{FT} <:AbstractClumpingIndexModel{FT}
    "Clumping structure a"
    Ω_A::FT = 1
    "Clumping structure b"
    Ω_B::FT = 0
end
