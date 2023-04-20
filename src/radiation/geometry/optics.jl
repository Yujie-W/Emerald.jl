"""

    average_transmittance(α::FT, nr::FT) where {FT}

Return the average transmittance of isotropic radiation across an interface between two dielectrics, given
- `α` angle of incidence
- `nr` Index of refraction

# References
- Stern (1964) Transmission of isotropic radiation across an interface between two dielectrics. Applied Optics 3(1): 111-113.
- Allen (1973) Transmission of isotropic light across a dielectric surface in two and three dimensions. Journal of the Optical Society of America 63(6): 664-666.

"""
function average_transmittance(α::FT, nr::FT) where {FT}
    @assert 0 < α <= 90;

    # some shortcuts to avoid overly comlicated equation
    _a     = (nr + 1) ^ 2 / 2;
    _a³    = _a ^ 3;
    _n²    = nr ^ 2;
    _n⁴    = nr ^ 4;
    _n⁶    = nr ^ 6;
    _n²p   = _n² + 1;
    _n²p²  = _n²p ^ 2;
    _n²p³  = _n²p ^ 3;
    _n²m²  = (_n² - 1) ^ 2;
    _k     = -1 * _n²m² / 4;
    _k²    = _k ^ 2;
    _sin²α = sind(α) ^ 2;
    _b₂    = _sin²α - _n²p/2;
    _b₁    = (α==90 ? 0 : sqrt(_b₂ ^ 2 + _k));
    _b     = _b₁ - _b₂;
    _b³    = _b ^ 3;
    _npanm = 2 * _n²p * _a - _n²m²;
    _npbnm = 2 * _n²p * _b - _n²m²;

    # S polarization
    _ts  = ( _k² / (6*_b³) + _k/_b - _b/2 ) - ( _k² / (6*_a³) + _k/_a - _a/2 );

    # P polarization
    _tp₁ = -2 * _n² * (_b - _a) / _n²p²;
    _tp₂ = -2 * _n² * _n²p * log(_b / _a) / _n²m²;
    _tp₃ = _n² * (1/_b - 1/_a) / 2;
    _tp₄ = 16 * _n⁴ * (_n⁴ + 1) * log(_npbnm / _npanm) / (_n²p³ * _n²m²);
    _tp₅ = 16 * _n⁶ * (1/_npbnm - 1/_npanm) / _n²p³;
    _tp  = _tp₁ + _tp₂ + _tp₃ + _tp₄ + _tp₅;

    return  (_ts + _tp) / (2 * _sin²α)
end


"""

    leaf_optical_properties(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT}) where {FT}

Return a DIM_WL-element vector of DIM_CANOPY-element vector of leaf optics, given
- `config` `EmeraldConfiguration` struct
- `state` All state variables in a multiple layer SPAC

"""
function leaf_optical_properties end

leaf_optical_properties(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT}) where {FT} = (
    _optics = leaf_optical_properties.((config.LHA,), (state,), eachindex(config.LHA.K_CAB));

    return SVector(_optics)
);

leaf_optical_properties(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT}, j::Int) where {FT} = leaf_optical_properties(config.LHA, state, j);

"""
This method return a vector of tuples of leaf optics for all canopy layers for one wavelength.
This method is supposed to parallelize the calculations, and thus is not meant for public use.
"""
leaf_optical_properties(lha::HyperspectralAbsorption{FT}, state::MultipleLayerSPACState{FT}, j::Int) where {FT} = (
    _optics = leaf_optical_properties.(
                lha.K_ANT[j],
                lha.K_BROWN[j],
                lha.K_CAB[j],
                lha.K_CAR_V[j],
                lha.K_CAR_Z[j],
                lha.K_CBC[j],
                lha.K_H₂O[j],
                lha.K_LMA[j],
                lha.K_PRO[j],
                lha.NR[j],
                state.ant,
                state.brown,
                state.car,
                state.cbc,
                state.chl,
                state.f_zeax,
                state.lma,
                state.mesophyll,
                state.pro,
                state.water);

    return SVector(_optics)
);

"""
This method return the shortwave reflectance, transmittance, and PPAR:APAR ratio for given coefficients and pigment contents.
This method is supposed to parallelize the calculations, and thus is not meant for public use.
"""
leaf_optical_properties(
            K_ANT::FT,
            K_BROWN::FT,
            K_CAB::FT,
            K_CAR_V::FT,
            K_CAR_Z::FT,
            K_CBC::FT,
            K_H₂O::FT,
            K_LMA::FT,
            K_PRO::FT,
            NR::FT,
            ant::FT,
            brown::FT,
            car::FT,
            cbc::FT,
            chl::FT,
            f_zeax::FT,
            lma::FT,
            mesophyll::FT,
            pro::FT,
            water::FT;
            α::FT = FT(40)
) where {FT} = (
    # calculate the average absorption feature and relative Cab and Car partitions
    _a_lhc = (K_CAB * chl + K_CAR_V * car * (1 - f_zeax) + K_CAR_Z * car * f_zeax);
    _a_non = (K_ANT * ant + K_BROWN * brown + K_H₂O * water + K_CBC * cbc + K_PRO * pro + K_LMA * (lma - cbc - pro));
    _a_all = _a_lhc + _a_non;
    _k_all = _a_all / mesophyll;
    _k_lhc = _a_lhc / _a_all;

    # calculate the reflectance and transmittance at the interfaces of one layer
    _τ   = max(0, (1 - _k_all) * exp(-_k_all) + _k_all ^ 2 * expint(_k_all + eps(FT)));
    _τ_α = average_transmittance(α, NR);
    _ρ_α = 1 - _τ_α;
    _τ₁₂ = average_transmittance(FT(90), NR);
    _ρ₁₂ = 1 - _τ₁₂;
    _τ₂₁ = _τ₁₂ / (NR ^ 2);
    _ρ₂₁ = 1 - _τ₂₁;

    # top and bottom side transmittance and reflectance
    _denom = 1 - (_τ * _ρ₂₁) ^ 2;
    _τ_top = _τ_α * _τ * _τ₂₁ / _denom;
    _τ_btm = _τ₁₂ * _τ * _τ₂₁ / _denom;
    _ρ_top = _ρ_α + _τ * _ρ₂₁ * _τ_top;
    _ρ_btm = _ρ₁₂ + _τ * _ρ₂₁ * _τ_btm;

    # calculate the reflectance and transmittance at the interfaces of N layer
    _d     = sqrt((1 + _ρ_btm + _τ_btm) * (1 + _ρ_btm - _τ_btm) * (1 - _ρ_btm + _τ_btm) * (1 - _ρ_btm - _τ_btm));
    _ρ²    = _ρ_btm ^ 2;
    _τ²    = _τ_btm ^ 2;
    _a     = (1 + _ρ² - _τ² + _d) / (2 * _ρ_btm);
    _b     = (1 - _ρ² + _τ² + _d) / (2 * _τ_btm);
    _bⁿ⁻¹  = _b ^ (mesophyll - 1);
    _b²ⁿ⁻² = _bⁿ⁻¹ ^ 2;
    _a²    = _a ^ 2;
    _denom = _a² * _b²ⁿ⁻² - 1;
    _ρ_sub = _a * (_b²ⁿ⁻² - 1) / _denom;
    _τ_sub = _bⁿ⁻¹ * (_a² - 1) / _denom;

    # avoid case of zero absorption
    if _ρ_btm + _τ_btm >= 1
        _τ_sub = _τ_btm / (_τ_btm + (1 - _τ_btm) * (mesophyll - 1));
        _ρ_sub = 1 - _τ_sub;
    end;

    # reflectance & transmittance of the leaf: combine top layer with next N-1 layers
    _denom = 1 - _ρ_sub * _ρ_btm;
    _τ_sw  = _τ_top * _τ_sub / _denom;
    _ρ_sw  = _ρ_top + _τ_top * _ρ_sub * _τ_btm / _denom;

    return _ρ_sw, _τ_sw, _k_lhc
);




leaf_optical_properties(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT}, i::Int, j::Int) where {FT} = (
    (; LHA) = config;

    return leaf_optical_properties(LHA, state.ant[i], state.brown[i], state.car[i], state.cbc[i], state.chl[i], state.f_zeax[i], state.lma[i], state.mesophyll[i], state.pro[i], state.water[i], j)
);

leaf_optical_properties(lha::HyperspectralAbsorption{FT}, ant::FT, brown::FT, car::FT, cbc::FT, chl::FT, f_zeax::FT, lma::FT, mesophyll::FT, pro::FT, water::FT, j::Int; α::FT = FT(40)) where {FT} = (
    (; K_ANT, K_BROWN, K_CAB, K_CAR_V, K_CAR_Z, K_CBC, K_H₂O, K_LMA, K_PRO, NR) = lha;

    # calculate the average absorption feature and relative Cab and Car partitions
    _a_lhc = (K_CAB[j] * chl + K_CAR_V[j] * car * (1 - f_zeax) + K_CAR_Z[j] * car * f_zeax);
    _a_non = (K_ANT[j] * ant + K_BROWN[j] * brown + K_H₂O[j] * water + K_CBC[j] * cbc + K_PRO[j] * pro + K_LMA[j] * (lma - cbc - pro));
    _a_all = _a_lhc + _a_non;
    _k_all = _a_all / mesophyll;
    _k_lhc = _a_lhc / _a_all;

    # calculate the reflectance and transmittance at the interfaces of one layer
    _τ   = max(0, (1 - _k_all) * exp(-_k_all) + _k_all ^ 2 * expint(_k_all + eps(FT)));
    _τ_α = average_transmittance(α, NR[j]);
    _ρ_α = 1 - _τ_α;
    _τ₁₂ = average_transmittance(FT(90), NR[j]);
    _ρ₁₂ = 1 - _τ₁₂;
    _τ₂₁ = _τ₁₂ / (NR[j] ^ 2);
    _ρ₂₁ = 1 - _τ₂₁;

    # top and bottom side transmittance and reflectance
    _denom = 1 - (_τ * _ρ₂₁) ^ 2;
    _τ_top = _τ_α * _τ * _τ₂₁ / _denom;
    _τ_btm = _τ₁₂ * _τ * _τ₂₁ / _denom;
    _ρ_top = _ρ_α + _τ * _ρ₂₁ * _τ_top;
    _ρ_btm = _ρ₁₂ + _τ * _ρ₂₁ * _τ_btm;

    # calculate the reflectance and transmittance at the interfaces of N layer
    _d     = sqrt((1 + _ρ_btm + _τ_btm) * (1 + _ρ_btm - _τ_btm) * (1 - _ρ_btm + _τ_btm) * (1 - _ρ_btm - _τ_btm));
    _ρ²    = _ρ_btm ^ 2;
    _τ²    = _τ_btm ^ 2;
    _a     = (1 + _ρ² - _τ² + _d) / (2 * _ρ_btm);
    _b     = (1 - _ρ² + _τ² + _d) / (2 * _τ_btm);
    _bⁿ⁻¹  = _b ^ (mesophyll - 1);
    _b²ⁿ⁻² = _bⁿ⁻¹ ^ 2;
    _a²    = _a ^ 2;
    _denom = _a² * _b²ⁿ⁻² - 1;
    _ρ_sub = _a * (_b²ⁿ⁻² - 1) / _denom;
    _τ_sub = _bⁿ⁻¹ * (_a² - 1) / _denom;

    # avoid case of zero absorption
    if _ρ_btm + _τ_btm >= 1
        _τ_sub = _τ_btm / (_τ_btm + (1 - _τ_btm) * (mesophyll - 1));
        _ρ_sub = 1 - _τ_sub;
    end;

    # reflectance & transmittance of the leaf: combine top layer with next N-1 layers
    _denom = 1 - _ρ_sub * _ρ_btm;
    _τ_sw  = _τ_top * _τ_sub / _denom;
    _ρ_sw  = _ρ_top + _τ_top * _ρ_sub * _τ_btm / _denom;

    return _ρ_sw, _τ_sw, _k_lhc
);
