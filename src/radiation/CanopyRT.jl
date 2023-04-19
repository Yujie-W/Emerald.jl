module CanopyRT

using QuadGK: quadgk
using SpecialFunctions: expint
using StaticArrays: SVector

using ..NameSpace: CanopyCache, ClumpingIndexPinty, EmeraldConfiguration, HyperspectralAbsorption, MultipleLayerSPACState, WaveLengthSet
using ..NameSpace: dim_azi, dim_canopy


include("geometry/clumping.jl")
include("geometry/extinction.jl")
include("geometry/geometry.jl")
include("geometry/optics.jl")
include("geometry/update.jl")


end # module
