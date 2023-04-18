Base.@kwdef struct MultipleLayerSPACState{FT,DIM_CANOPY}
    # Leaf biophysics used for canopy radiative transfer
    "Anthocyanin content `[μg cm⁻²]`"
    ant::NTuple{DIM_CANOPY,FT} = NTuple{DIM_CANOPY,FT}(zeros(FT,DIM_CANOPY))
    "Senescent material (brown pigments) fraction `[-]`"
    brown::NTuple{DIM_CANOPY,FT} = NTuple{DIM_CANOPY,FT}(zeros(FT,DIM_CANOPY))
    "Leaf carotenoid contents `[μg cm⁻²]`"
    car::NTuple{DIM_CANOPY,FT} = NTuple{DIM_CANOPY,FT}(ones(FT,DIM_CANOPY) .* 40 ./ 7)
    "Carbon-based constituents in lma `[g cm⁻²]`"
    cbc::NTuple{DIM_CANOPY,FT} = NTuple{DIM_CANOPY,FT}(zeros(FT,DIM_CANOPY))
    "Leaf chlorophyll contents `[μg cm⁻²]`"
    chl::NTuple{DIM_CANOPY,FT} = NTuple{DIM_CANOPY,FT}(ones(FT,DIM_CANOPY) .* 40)
    "Zeaxanthin fraction in Carotenoid (1=all Zeaxanthin, 0=all Violaxanthin) `[-]`"
    f_zeax::NTuple{DIM_CANOPY,FT} = NTuple{DIM_CANOPY,FT}(zeros(FT,DIM_CANOPY))
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    lma::NTuple{DIM_CANOPY,FT} = NTuple{DIM_CANOPY,FT}(ones(FT,DIM_CANOPY) .* FT(0.012))
    "Leaf mesophyll structural parameter that describes mesophyll reflectance and transmittance"
    mesophyll::NTuple{DIM_CANOPY,FT} = NTuple{DIM_CANOPY,FT}(ones(FT,DIM_CANOPY) .* FT(1.4))
    "Protein content in lma (pro = lma - cbc) `[g cm⁻²]`"
    pro::NTuple{DIM_CANOPY,FT} = NTuple{DIM_CANOPY,FT}(zeros(FT,DIM_CANOPY))
end
