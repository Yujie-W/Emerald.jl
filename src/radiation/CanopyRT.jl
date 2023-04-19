module CanopyRT

using SpecialFunctions: expint

using ..NameSpace: CanopyCache, EmeraldConfiguration, HyperspectralAbsorption, MultipleLayerSPACState, WaveLengthSet


include("geometry/extinction.jl")
include("geometry/optics.jl")
include("geometry/shortwave.jl")
include("geometry/update.jl")


end # module
