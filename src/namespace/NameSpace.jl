module NameSpace

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using StaticArrays: SMatrix, SVector

using ..EmeraldIO.Netcdf: read_nc, size_nc


const LAND_2017 = artifact"land_model_spectrum_V3" * "/clima_land_spectra_2017.nc";
const LAND_2021 = artifact"land_model_spectrum_V3" * "/clima_land_spectra_2021.nc";


"""

Return the dimensions of the sturct

"""
function dims end


include("config/canopy.jl")
include("config/constant.jl")
include("config/leaf_bio.jl")
include("config/radiation.jl")
include("config/soil.jl")
include("config/wavelength.jl")

include("config/all.jl")

include("cache/canopy.jl")

include("state/spac.jl")


end # module
