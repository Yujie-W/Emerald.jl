module Netcdf

using DataFrames: DataFrame
using NCDatasets: Dataset, defDim, defVar


include("netcdf/info.jl"  )
include("netcdf/read.jl")


end # module
