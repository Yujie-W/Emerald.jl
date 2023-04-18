"""

$(TYPEDEF)

Immutable struct that contains leaf biophysical traits used to run leaf reflection and transmittance.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct HyperspectralAbsorption{FT,DIM_WL}
    "Specific absorption coefficients of anthocynanin `[-]`"
    K_ANT::SVector{DIM_WL,FT}
    "Specific absorption coefficients of senescent material (brown pigments) `[-]`"
    K_BROWN::SVector{DIM_WL,FT}
    "Specific absorption coefficients of chlorophyll a and b `[-]`"
    K_CAB::SVector{DIM_WL,FT}
    "Specific absorption coefficients of violaxanthin carotenoid `[-]`"
    K_CAR_V::SVector{DIM_WL,FT}
    "Specific absorption coefficients of zeaxanthin carotenoid `[-]`"
    K_CAR_Z::SVector{DIM_WL,FT}
    "Specific absorption coefficients of carbon-based constituents `[-]`"
    K_CBC::SVector{DIM_WL,FT}
    "Specific absorption coefficients of water `[-]`"
    K_H₂O::SVector{DIM_WL,FT}
    "Specific absorption coefficients of dry matter `[-]`"
    K_LMA::SVector{DIM_WL,FT}
    "Specific absorption coefficients of protein `[-]`"
    K_PRO::SVector{DIM_WL,FT}
    "Specific absorption coefficients of PS I and II `[-]`"
    K_PS::SVector{DIM_WL,FT}
    "Refractive index `[-]`"
    NR::SVector{DIM_WL,FT}
end

HyperspectralAbsorption{FT}(dset::String) where {FT} = (
    _dim_λ = size_nc(dset, "K_ANT")[2][1];

    return HyperspectralAbsorption{FT,_dim_λ}(
                K_ANT   = read_nc(dset, "K_ANT") .* 1e-5,
                K_BROWN = read_nc(dset, "K_BROWN") .* 1e-5,
                K_CAB   = read_nc(dset, "K_CAB") .* 1e-5,
                K_CAR_V = read_nc(dset, "K_CAR_V") .* 1e-5,
                K_CAR_Z = read_nc(dset, "K_CAR_Z") .* 1e-5,
                K_CBC   = read_nc(dset, "K_CBC") .* 10,
                K_H₂O   = read_nc(dset, "K_H₂O") ./ 100 .* ρ_H₂O(),
                K_LMA   = read_nc(dset, "K_LMA") .* 10,
                K_PRO   = read_nc(dset, "K_PRO") .* 10,
                K_PS    = read_nc(dset, "K_PS"),
                NR      = read_nc(dset, "NR"),
    )
);

dims(::HyperspectralAbsorption{FT,DIM_WL}) where {FT,DIM_WL} = (DIM_WL,);
