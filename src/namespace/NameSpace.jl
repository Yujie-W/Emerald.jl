module NameSpace

using LazyArtifacts

using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using StaticArrays: SMatrix, SVector

using ..EmeraldIO.Netcdf: read_nc, size_nc


const LAND_2017 = artifact"land_model_spectrum_V3" * "/clima_land_spectra_2017.nc";
const LAND_2021 = artifact"land_model_spectrum_V3" * "/clima_land_spectra_2021.nc";


USE_STATIC_ARRAY = true;
use_static_arrays!(yesno::Bool = true) = (global USE_STATIC_ARRAY; USE_STATIC_ARRAY = yesno; return nothing);


include("config/canopy.jl")
include("config/constant.jl")
include("config/leaf_bio.jl")
include("config/radiation.jl")
include("config/soil.jl")
include("config/wavelength.jl")

include("config/all.jl")

include("model/clumping.jl")
include("model/soil_vc.jl")

include("state/spac.jl")

include("cache/canopy.jl")


end # module
