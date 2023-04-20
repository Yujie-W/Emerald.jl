module CanopyRT

using QuadGK: quadgk
using SpecialFunctions: expint
using StaticArrays: SMatrix, SVector

using ..NameSpace: CanopyRadiationCache, ClumpingIndexPinty, EmeraldConfiguration, HyperspectralAbsorption, MultipleLayerSPACState, WaveLengthSet
using ..NameSpace: ρ_LEAF_LW, ρ_SOIL_LW, τ_LEAF_LW
using ..NameSpace: dim_azi, dim_canopy, dim_incl, dim_wl


include("geometry/clumping.jl")
include("geometry/extinction.jl")
include("geometry/geometry.jl")
include("geometry/optics.jl")
include("geometry/soil.jl")
include("geometry/update.jl")


end # module
