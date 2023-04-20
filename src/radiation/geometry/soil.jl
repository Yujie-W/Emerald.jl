
function soil_albedo(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT}) where {FT}
    @assert 1<= state.soil_color <= 20 "Soil color $(state.soil_color) is not within [1,20]!";

    #=
    # use Yujie's method via relative soil water content
    _rwc = state.θ / LAYERS[1].VC.Θ_SAT;
    albedo.ρ_sw[1] = SOIL_ALBEDOS[COLOR,1] * (1 - _rwc) + _rwc * SOIL_ALBEDOS[COLOR,3];
    albedo.ρ_sw[2] = SOIL_ALBEDOS[COLOR,2] * (1 - _rwc) + _rwc * SOIL_ALBEDOS[COLOR,4];

    @assert !config.HYPER_SOIL "Hyperspectral soil albedo is not supported at the moment!";
    =#

    return nothing
end
