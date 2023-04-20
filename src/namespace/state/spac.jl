"""

$(TYPEDEF)

Structure for all state variables in a multiple layer SPAC

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct MultipleLayerSPACState{FT,DIM_CANOPY}
    # Leaf biophysics used for canopy radiative transfer
    "Anthocyanin content `[kg m⁻²]`"
    ant::SVector{DIM_CANOPY,FT} = zeros(FT,DIM_CANOPY)
    "Senescent material (brown pigments) fraction `[kg m⁻²]`"
    brown::SVector{DIM_CANOPY,FT} = zeros(FT,DIM_CANOPY)
    "Leaf carotenoid contents `[kg m⁻²]`"
    car::SVector{DIM_CANOPY,FT} = ones(FT,DIM_CANOPY) .* 40 ./ 7
    "Carbon-based constituents in lma `[kg m⁻²]`"
    cbc::SVector{DIM_CANOPY,FT} = zeros(FT,DIM_CANOPY)
    "Leaf chlorophyll contents `[kg m⁻²]`"
    chl::SVector{DIM_CANOPY,FT} = ones(FT,DIM_CANOPY) .* 40
    "Zeaxanthin fraction in Carotenoid (1=all Zeaxanthin, 0=all Violaxanthin) `[-]`"
    f_zeax::SVector{DIM_CANOPY,FT} = zeros(FT,DIM_CANOPY)
    "Dry matter content (dry leaf mass per unit area) `[kg m⁻²]`"
    lma::SVector{DIM_CANOPY,FT} = ones(FT,DIM_CANOPY) .* FT(0.012)
    "Leaf mesophyll structural parameter that describes mesophyll reflectance and transmittance"
    mesophyll::SVector{DIM_CANOPY,FT} = ones(FT,DIM_CANOPY) .* FT(1.4)
    "Protein content in lma (pro = lma - cbc) `[kg m⁻²]`"
    pro::SVector{DIM_CANOPY,FT} = zeros(FT,DIM_CANOPY)
    "Leaf water content `[kg m⁻²]`"
    water::SVector{DIM_CANOPY,FT} = ones(FT,DIM_CANOPY) .* FT(0.06)

    # Canopy structure used for canopy radiative transfer
    "Clumping index model"
    ci::ClumpingIndexPinty{FT} = ClumpingIndexPinty{FT}()
    "Leaf area index"
    lai::FT = 3
    "Leaf area index distribution"
    δlai::SVector{DIM_CANOPY,FT} = ones(FT,DIM_CANOPY) .* lai ./ DIM_CANOPY
end

dim_canopy(::MultipleLayerSPACState{FT,DIM_CANOPY}) where {FT,DIM_CANOPY} = DIM_CANOPY;
