"""

$(TYPEDEF)

Structure for all state variables in a multiple layer SPAC

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct MultipleLayerSPACState{FT,SOIL_VC,VT_CANOPY,VT_SOIL}
    # Leaf biophysics used for canopy radiative transfer
    "Anthocyanin content `[kg m⁻²]`"
    ant::VT_CANOPY
    "Senescent material (brown pigments) fraction `[kg m⁻²]`"
    brown::VT_CANOPY
    "Leaf carotenoid contents `[kg m⁻²]`"
    car::VT_CANOPY
    "Carbon-based constituents in lma `[kg m⁻²]`"
    cbc::VT_CANOPY
    "Leaf chlorophyll contents `[kg m⁻²]`"
    chl::VT_CANOPY
    "Zeaxanthin fraction in Carotenoid (1=all Zeaxanthin, 0=all Violaxanthin) `[-]`"
    f_zeax::VT_CANOPY
    "Dry matter content (dry leaf mass per unit area) `[kg m⁻²]`"
    lma::VT_CANOPY
    "Leaf mesophyll structural parameter that describes mesophyll reflectance and transmittance"
    mesophyll::VT_CANOPY
    "Protein content in lma (pro = lma - cbc) `[kg m⁻²]`"
    pro::VT_CANOPY
    "Leaf water content `[kg m⁻²]`"
    water::VT_CANOPY

    # Canopy structure used for canopy radiative transfer
    "Clumping index model"
    ci::ClumpingIndexPinty{FT}
    "Leaf area index"
    lai::FT
    "Leaf area index distribution"
    δlai::VT_CANOPY

    # Soil parameters
    "Soil color"
    soil_color::Int
    "Soil vulnerability"
    soil_vc::SOIL_VC
    "Soil moisture"
    θ::VT_SOIL
end

MultipleLayerSPACState{FT}(; DIM_CANOPY::Int = 10, DIM_SOIL::Int = 4) where {FT} = (
    SOIL_VC_TYPE = VanGenuchten;

    if USE_STATIC_ARRAY
        VT_CANOPY = SVector{DIM_CANOPY,FT};
        VT_SOIL   = SVector{DIM_SOIL,FT};
        SOIL_VC   = SVector{DIM_SOIL,SOIL_VC_TYPE{FT}};
    else
        VT_CANOPY = Vector{FT};
        VT_SOIL   = Vector{FT};
        SOIL_VC   = Vector{SOIL_VC_TYPE{FT}};
    end;

    return MultipleLayerSPACState{FT,SOIL_VC,VT_CANOPY,VT_SOIL}(
                ant        = zeros(FT,DIM_CANOPY),
                brown      = zeros(FT,DIM_CANOPY),
                car        = ones(FT,DIM_CANOPY) .* 40 ./ 7 * 1e-5,
                cbc        = zeros(FT,DIM_CANOPY),
                chl        = ones(FT,DIM_CANOPY) .* 40 * 1e-5,
                f_zeax     = zeros(FT,DIM_CANOPY),
                lma        = ones(FT,DIM_CANOPY) .* FT(0.012) * 10,
                mesophyll  = ones(FT,DIM_CANOPY) .* FT(1.4),
                pro        = zeros(FT,DIM_CANOPY),
                water      = ones(FT,DIM_CANOPY) .* FT(0.06),
                ci         = ClumpingIndexPinty{FT}(),
                lai        = 3,
                δlai       = ones(FT,DIM_CANOPY) .* 3 ./ DIM_CANOPY,
                soil_color = 1,
                soil_vc    = [SOIL_VC_TYPE{FT}("Loam") for _ in 1:DIM_SOIL],
                θ          = ones(FT,DIM_SOIL) .* 0.2,
    )
);

dim_canopy(var::MultipleLayerSPACState) = length(var.ant);
dim_soil(var::MultipleLayerSPACState) = length(var.θ);
