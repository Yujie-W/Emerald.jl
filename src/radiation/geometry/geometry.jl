
function sun_geometry(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT,DIM_CANOPY}, sza::FT, p_incl::SVector{DIM_INCL,FT}) where {FT,DIM_CANOPY,DIM_INCL}
    DIM_AZI = dim_azi(config);

    # 0. compute the clumping index
    _ci = clumping_index(state.ci, sza);

    # 1. update the canopy optical properties related to extinction and scattering coefficients
    _coefs = extinction_coefficient.(sza, config.CAN.Θ_INCL);
    __ks = [_coef[1] for _coef in _coefs];
    __Cs = [_coef[2] for _coef in _coefs];
    __Ss = [_coef[3] for _coef in _coefs];
    _ks = p_incl' * __ks;
    _bf = p_incl' * (cosd.(config.CAN.Θ_INCL) .^ 2);

    # 2. update the matrices fs and fo
    _fs = (__Cs * ones(FT,1,DIM_AZI) .+ __Ss * (cosd.(config.CAN.Θ_AZI))') ./ cosd(sza);
    _abs_fs = abs.(_fs);

    # 3. update the viewing fraction ps, po, pso, and p_sunlit
    _Σlai = FT[0; [sum(state.δlai[1:_i]) for _i in 1:DIM_CANOPY]];
    _f_sunlit = [1 / (_ks * _ci * state.δlai[_i]) * (exp(-_ks * _ci * _Σlai[_i]) - exp(-_ks * _ci * _Σlai[_i+1])) for _i in 1:DIM_CANOPY];


    return (_ks,_bf,_ci), _abs_fs, _f_sunlit
end


function canopy_scatter_matrices(
            config::EmeraldConfiguration{FT},
            state::MultipleLayerSPACState{FT,DIM_CANOPY},
            leaf_optics::SMatrix{DIM_CANOPY,DIM_WL},
            coefs::NTuple{6,FT},
            i::Int,
            j::Int
) where {FT,DIM_CANOPY,DIM_WL}
    (_ko,_ks,_sob,_sof,_bf,_ci) = coefs;
    _sdb = (_ks + _bf) / 2;
    _sdf = (_ks - _bf) / 2;
    _dob = (_ko + _bf) / 2;
    _dof = (_ko - _bf) / 2;
    _ddb = (1   + _bf) / 2;
    _ddf = (1   - _bf) / 2;
    _ilai = state.δlai[i] * _ci;

    # 1. update the scattering coefficients for ith canopy layer and jth wavelength
    _σ_ddb = _ddb * leaf_optics[i,j][1] + _ddf * leaf_optics[i,j][2];
    _σ_ddf = _ddf * leaf_optics[i,j][1] + _ddb * leaf_optics[i,j][2];
    _σ_sdb = _sdb * leaf_optics[i,j][1] + _sdf * leaf_optics[i,j][2];
    _σ_sdf = _sdf * leaf_optics[i,j][1] + _sdb * leaf_optics[i,j][2];
    _σ_dob = _dob * leaf_optics[i,j][1] + _dof * leaf_optics[i,j][2];
    _σ_dof = _dof * leaf_optics[i,j][1] + _dob * leaf_optics[i,j][2];
    _σ_so  = _sob * leaf_optics[i,j][1] + _sof * leaf_optics[i,j][2];

    # 2. update the transmittance and reflectance for single directions per layer (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    _τ_ss = exp(-1 * _ks * _ilai);
    _τ_dd = exp(-1 * (1 - _σ_ddf) * _ilai);
    _τ_sd = 1 - exp(-1 * _σ_sdf * _ilai);
    _ρ_dd = 1 - exp(-1 * _σ_ddb * _ilai);
    _ρ_sd = 1 - exp(-1 * _σ_sdb * _ilai);
end

#=

function canopy_geometry(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT,DIM_CANOPY}, sza::FT, vza::FT, raa::FT, p_incl::SVector{DIM_INCL,FT}) where {FT,DIM_CANOPY,DIM_INCL}
    DIM_AZI = dim_azi(config);

    # 0. compute the clumping index
    _ci = clumping_index(state.ci, sza);

    # 1. update the canopy optical properties related to extinction and scattering coefficients
    _coefs = extinction_coefficient.(sza, vza, raa, config.CAN.Θ_INCL);
    _ko  = p_incl' * [_coef[1] for _coef in _coefs];
    _ks  = p_incl' * [_coef[2] for _coef in _coefs];
    _sob = p_incl' * [_coef[3] for _coef in _coefs];
    _sof = p_incl' * [_coef[4] for _coef in _coefs];
    _bf  = p_incl' * (cosd.(config.CAN.Θ_INCL) .^ 2);

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
    _f_sunlit = [1 / (_ks * _ci * state.δlai[_i]) * (exp(-_ks * _ci * _Σlai[_i]) - exp(-_ks * _ci * _Σlai[_i+1])) for _i in 1:DIM_CANOPY];

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
    _pso = [quadgk(_pdf, -1 * _Σlai[_i] / state.lai - FT(0.01), -1 * _Σlai[_i] / state.lai; rtol = 1e-2)[1] / (state.lai * FT(0.01)) for _i in 1:DIM_CANOPY+1];
    _f_so = [quadgk(_pdf, -1 * _Σlai[_i+1] / state.lai, -1 * _Σlai[_i] / state.lai; rtol = 1e-2)[1] / state.δlai[_i] for _i in 1:DIM_CANOPY];
    for _i in 1:DIM_CANOPY
        @info "Debugging" _pso[_i] _f_so[_i] _pso[_i+1];
    end;

    return (_ko,_ks,_sob,_sof,_bf,_ci), _abs_fo, _abs_fs, _po, _pso, _f_so, _f_sunlit
end


function canopy_scatter_matrices(
            config::EmeraldConfiguration{FT},
            state::MultipleLayerSPACState{FT,DIM_CANOPY},
            leaf_optics::SMatrix{DIM_CANOPY,DIM_WL},
            coefs::NTuple{6,FT},
            i::Int,
            j::Int
) where {FT,DIM_CANOPY,DIM_WL}
    (_ko,_ks,_sob,_sof,_bf,_ci) = coefs;
    _sdb = (_ks + _bf) / 2;
    _sdf = (_ks - _bf) / 2;
    _dob = (_ko + _bf) / 2;
    _dof = (_ko - _bf) / 2;
    _ddb = (1   + _bf) / 2;
    _ddf = (1   - _bf) / 2;
    _ilai = state.δlai[i] * _ci;

    # 1. update the scattering coefficients for ith canopy layer and jth wavelength
    _σ_ddb = _ddb * leaf_optics[i,j][1] + _ddf * leaf_optics[i,j][2];
    _σ_ddf = _ddf * leaf_optics[i,j][1] + _ddb * leaf_optics[i,j][2];
    _σ_sdb = _sdb * leaf_optics[i,j][1] + _sdf * leaf_optics[i,j][2];
    _σ_sdf = _sdf * leaf_optics[i,j][1] + _sdb * leaf_optics[i,j][2];
    _σ_dob = _dob * leaf_optics[i,j][1] + _dof * leaf_optics[i,j][2];
    _σ_dof = _dof * leaf_optics[i,j][1] + _dob * leaf_optics[i,j][2];
    _σ_so  = _sob * leaf_optics[i,j][1] + _sof * leaf_optics[i,j][2];

    # 2. update the transmittance and reflectance for single directions per layer (it was 1 - k*Δx, and we used exp(-k*Δx) as Δx is not infinitesmal)
    _τ_ss = exp(-1 * _ks * _ilai);
    _τ_dd = exp(-1 * (1 - _σ_ddf) * _ilai);
    _τ_sd = 1 - exp(-1 * _σ_sdf * _ilai);
    _ρ_dd = 1 - exp(-1 * _σ_ddb * _ilai);
    _ρ_sd = 1 - exp(-1 * _σ_sdb * _ilai);
end
=#
