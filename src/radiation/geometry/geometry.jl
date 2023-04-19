
function shortwave_radiation(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT}, sza::FT, vza::FT, raa::FT, p_incl::SVector{FT}) where {FT}
    DIM_AZI = dim_azi(config);
    DIM_CANOPY = dim_canopy(state);

    # 0. compute the clumping index
    _ci = 1;

    # 1. update the canopy optical properties related to extinction and scattering coefficients
    _coefs = extinction_coefficient.(sza, vza, raa, config.CAN.Θ_INCL);
    # _ko, _ks, _sb, _sf, _Co, _Cs, _So, _Ss
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
    _fs_fo = _fs .* _fo;
    _abs_fo = abs.(_fo);
    _abs_fs = abs.(_fs);
    _abs_fs_fo = abs.(_fs_fo);

    # 3. update the viewing fraction ps, po, pso, and p_sunlit
    # _fac_s = (1 - exp(-_ks * _ci * state.lai / DIM_CANOPY)) / (_ks * _ci * state.lai / DIM_CANOPY);
    # _fac_o = (1 - exp(-_ko * _ci * state.lai / DIM_CANOPY)) / (_ko * _ci * state.lai / DIM_CANOPY);
    # _po = exp.(can._x_bnds .* _ko .* _ci .* state.lai) .* _fac_o;
    # _ps = exp.(can._x_bnds .* _ks .* _ci .* state.lai) .* _fac_s;
    # _p_sunlit = (_ps[1:DIM_CANOPY] .+ _ps[2:DIM_CANOPY+1]) ./ 2;

    return nothing
end
