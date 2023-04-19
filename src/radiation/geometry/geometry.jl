
function shortwave_radiation(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT}, sza::FT, vza::FT, raa::FT, p_incl::SVector{DIM_INCL,FT}) where {FT,DIM_INCL}
    DIM_AZI = dim_azi(config);
    DIM_CANOPY = dim_canopy(state);

    # 0. compute the clumping index
    _ci = clumping_index(state.ci, sza);

    # 1. update the canopy optical properties related to extinction and scattering coefficients
    _coefs = extinction_coefficient.(sza, vza, raa, config.CAN.Θ_INCL);
    _ko  = p_incl' * [_coef[1] for _coef in _coefs];
    _ks  = p_incl' * [_coef[2] for _coef in _coefs];
    _sob = p_incl' * [_coef[3] for _coef in _coefs];
    _sof = p_incl' * [_coef[4] for _coef in _coefs];
    _bf  = p_incl' * (cosd.(config.CAN.Θ_INCL) .^ 2);
    _sdb = (_ks + _bf) / 2;
    _sdf = (_ks - _bf) / 2;
    _dob = (_ko + _bf) / 2;
    _dof = (_ko - _bf) / 2;
    _ddb = (1   + _bf) / 2;
    _ddf = (1   - _bf) / 2;

    # 2. update the matrices fs and fo
    _Co = [_coef[5] for _coef in _coefs];
    _Cs = [_coef[6] for _coef in _coefs];
    _So = [_coef[7] for _coef in _coefs];
    _Ss = [_coef[8] for _coef in _coefs];
    _fo = (_Co * ones(FT,1,DIM_AZI) .+ _So * (cosd.(config.CAN.Θ_AZI .- raa))') ./ cosd(vza);
    _fs = (_Cs * ones(FT,1,DIM_AZI) .+ _Ss * (cosd.(config.CAN.Θ_AZI))') ./ cosd(sza);
    _abs_fo = abs.(_fo);
    _abs_fs = abs.(_fs);

    # 3. update the viewing fraction ps, po, pso, and p_sunlit
    _Σlai = FT[0; [sum(state.δlai[1:_i]) for _i in 1:DIM_CANOPY]];
    _po = exp.(-_ko .* _ci .* _Σlai);
    _p_sunlit = [1 / (_ks * _ci * state.δlai[_i]) * (exp(-_ks * _ci * _Σlai[_i]) - exp(-_ks * _ci * _Σlai[_i+1])) for _i in 1:DIM_CANOPY];

    _dso = sqrt( tand(sza) ^ 2 + tand(vza) ^ 2 - 2 * tand(sza) * tand(vza) * cosd(raa) );
    @inline _pdf(x::FT) where {FT} = (
        _Σk = _ko + _ks;
        _Πk = _ko * _ks;
        _cl = _ci * state.lai;
        _α  = _dso / config.CAN.HOT_SPOT * 2 / _Σk;

        if _dso == 0
            return exp( (_Σk - sqrt(_Πk)) * _cl * x )
        end;

        return exp( _Σk * _cl * x + sqrt(_Πk) * _cl / _α * (1 - exp(_α * x)) )
    );
    _pso = [quadgk(_pdf, -1 * _Σlai[_i+1] / state.lai, -1 * _Σlai[_i] / state.lai; rtol = 1e-2)[1] / state.δlai[_i] for _i in 1:DIM_CANOPY];

    return _abs_fo, _abs_fs, _po, _pso, _p_sunlit
end
