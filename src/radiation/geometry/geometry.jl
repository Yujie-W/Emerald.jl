
function sun_geometry(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT}, sza::FT, p_incl::SVector{DIM_INCL,FT}) where {FT,DIM_INCL}
    DIM_AZI = dim_azi(config);
    DIM_CANOPY = dim_canopy(state);

    # compute the clumping index
    _ci = clumping_index(state.ci, sza);

    # compute the canopy optical properties related to extinction coefficients
    _coefs = extinction_coefficient.(sza, config.CAN.Θ_INCL);
    __ks = [_coef[1] for _coef in _coefs];
    __Cs = [_coef[2] for _coef in _coefs];
    __Ss = [_coef[3] for _coef in _coefs];
    _ks = p_incl' * __ks;
    _bf = p_incl' * (cosd.(config.CAN.Θ_INCL) .^ 2);

    # compute the matrix fs
    _fs = (__Cs * ones(FT,1,DIM_AZI) .+ __Ss * (cosd.(config.CAN.Θ_AZI))') ./ cosd(sza);
    _abs_fs = abs.(_fs);

    # compute the sunlit fraction f_sunlit
    _Σlai = FT[0; [sum(state.δlai[1:_i]) for _i in 1:DIM_CANOPY]];
    _f_sunlit = [1 / (_ks * _ci * state.δlai[_i]) * (exp(-_ks * _ci * _Σlai[_i]) - exp(-_ks * _ci * _Σlai[_i+1])) for _i in 1:DIM_CANOPY];

    # long wave reflectance and transmittance
    _ddb = (1 + _bf) / 2;
    _ddf = (1 - _bf) / 2;
    _σ_lwb = _ddb * ρ_LEAF_LW(config.CONST) + _ddf * τ_LEAF_LW(config.CONST);
    _σ_lwf = _ddf * ρ_LEAF_LW(config.CONST) + _ddb * τ_LEAF_LW(config.CONST);
    _ρ_lw = 1 .- exp.(-1 .* _σ_lwb .* _ci .* state.δlai);
    _τ_lw = exp.(-1 .* (1 - _σ_lwf) .* _ci .* state.δlai);

    return (_ks,_bf,_ci), _abs_fs, _f_sunlit, _ρ_lw, _τ_lw
end


function longwave_coefs(config::EmeraldConfiguration{FT}, state::MultipleLayerSPACState{FT}, ext_coefs::NTuple{3,FT}) where {FT}
    (_,_bf,_ci) = ext_coefs;

    # long wave reflectance and transmittance based on leaf area
    _ddb = (1 + _bf) / 2;
    _ddf = (1 - _bf) / 2;
    _σ_lwb = _ddb * ρ_LEAF_LW(config.CONST) + _ddf * τ_LEAF_LW(config.CONST);
    _σ_lwf = _ddf * ρ_LEAF_LW(config.CONST) + _ddb * τ_LEAF_LW(config.CONST);
    _ρ_lw = 1 .- exp.(-1 .* _σ_lwb .* _ci .* state.δlai);
    _τ_lw = exp.(-1 .* (1 - _σ_lwf) .* _ci .* state.δlai);

    return _ρ_lw, _τ_lw
end


function effective_longwave_coefs(config::EmeraldConfiguration{FT}, ρ_lw::VT_CANOPY, τ_lw::VT_CANOPY) where {FT,VT_CANOPY}
    DIM_CANOPY = length(ρ_lw);

    _r_lw = zeros(FT, DIM_CANOPY+1);
    _t_lw = zeros(FT, DIM_CANOPY);

    # set last reflectance to be soil albedo
    _r_lw[end] = ρ_SOIL_LW(config.CONST);

    # iterate from bottom to top
    for _i in DIM_CANOPY:-1:1
        _ρ_lw_i = ρ_lw[_i];
        _τ_lw_i = τ_lw[_i];
        _r_lw_j = _r_lw[_i+1];

        _dnorm_lw = 1 - _ρ_lw_i * _r_lw_j;

        _t_lw[_i] = (_τ_lw_i) / _dnorm_lw;                      # it, then rescale
        _r_lw[_i] = _ρ_lw_i + _τ_lw_i * _r_lw_j * _t_lw[_i];    # ir + it-jr-it
    end;

    return _r_lw, _t_lw
end


"""

    shortwave_coefs(leaf_optics::VT_WL_VT_CANOPY_T3, δlai::VT_CANOPY, ext_coefs::NTuple{3,FT}) where {FT,VT_CANOPY,VT_WL_VT_CANOPY_T3}

Return a DIM_WL-element vector of DIM_CANOPY-element vector of shortwave scattering coefficients, given
- `leaf_optics` Leaf optical properties
- `δlai` Leaf area index per canopy layer
- `ext_coefs` Canopy extinction coefficients

"""
function shortwave_coefs end

shortwave_coefs(leaf_optics::VT_WL_VT_CANOPY_T3, δlai::VT_CANOPY, ext_coefs::NTuple{3,FT}) where {FT,VT_CANOPY,VT_WL_VT_CANOPY_T3} = (
    return shortwave_coefs.(leaf_optics, (δlai,), ext_coefs...)
);

# This method return a vector of tuple of scattering coefficients for shortwave radiation.
# This method is supposed to parallelize the calculations, and thus is not meant for public use.
shortwave_coefs(leaf_optics::VT_CANOPY_T3, δlai::VT_CANOPY, ks::FT, bf::FT, ci::FT) where {FT,VT_CANOPY,VT_CANOPY_T3} = (
    return shortwave_coefs.(leaf_optics, δlai, ks, bf, ci);
);

# This method return a tuple of scattering coefficients for shortwave radiation.
# This method is supposed to parallelize the calculations, and thus is not meant for public use.
shortwave_coefs(leaf_optics::NTuple{3,FT}, δlai::FT, ks::FT, bf::FT, ci::FT) where {FT} = (
    _ilai = δlai * ci;
    _sdb = (ks + bf) / 2;
    _sdf = (ks - bf) / 2;
    _ddb = (1 + bf) / 2;
    _ddf = (1 - bf) / 2;

    # compute the scattering coefficients for ith canopy layer and jth wavelength
    _σ_ddb = _ddb * leaf_optics[1] + _ddf * leaf_optics[2];
    _σ_ddf = _ddf * leaf_optics[1] + _ddb * leaf_optics[2];
    _σ_sdb = _sdb * leaf_optics[1] + _sdf * leaf_optics[2];
    _σ_sdf = _sdf * leaf_optics[1] + _sdb * leaf_optics[2];

    # update the transmittance and reflectance for single directions per layer
    _τ_ss = exp(-1 * ks * _ilai);
    _τ_dd = exp(-1 * (1 - _σ_ddf) * _ilai);
    _τ_sd = 1 - exp(-1 * _σ_sdf * _ilai);
    _ρ_dd = 1 - exp(-1 * _σ_ddb * _ilai);
    _ρ_sd = 1 - exp(-1 * _σ_sdb * _ilai);

    return _τ_ss, _τ_dd, _τ_sd, _ρ_dd, _ρ_sd
);


function effective_shortwave_coefs end

effective_shortwave_coefs(sca_coefs::VT_WL_VT_CANOPY_T5, ρ_soil_sw::VT_WL) where {VT_WL<:Union{SVector,Vector},VT_WL_VT_CANOPY_T5} = (
    return effective_shortwave_coefs.(sca_coefs, ρ_soil_sw);
)

effective_shortwave_coefs(sca_coefs::VT_CANOPY_T5, ρ_soil_sw::FT) where {FT,VT_CANOPY_T5} = (
    DIM_CANOPY = length(sca_coefs[1]);

    _t_dd = zeros(FT, DIM_CANOPY);
    _t_sd = zeros(FT, DIM_CANOPY);
    _r_dd = zeros(FT, DIM_CANOPY+1);
    _r_sd = zeros(FT, DIM_CANOPY+1);

    # set last reflectance to be soil albedo
    _r_dd[end] = ρ_soil_sw;
    _r_sd[end] = ρ_soil_sw;

    # iterate from bottom to top
    for _i in DIM_CANOPY:-1:1
        (_τ_ss_i,_τ_dd_i,_τ_sd_i,_ρ_dd_i,_ρ_sd_i) = sca_coefs[_i];
        _r_dd_j = _r_dd[_i+1];
        _r_sd_j = _r_sd[_i+1];

        _denom_sw = 1 - _ρ_dd_i * _r_dd_j;

        _t_dd[_i] = (_τ_dd_i) / _denom_sw;                                                  # it, then rescale
        _t_sd[_i] = (_τ_sd_i + _τ_ss_i * _r_sd_j * _ρ_dd_i) / _denom_sw;                    # it + it-jr-ir, then rescale
        _r_dd[_i] = _ρ_dd_i + _τ_dd_i * _r_dd_j * _t_dd[_i];                                # ir + it-jr-it
        _r_sd[_i] = _ρ_sd_i + _τ_ss_i * _r_sd_j * _τ_dd_i + _t_sd[_i] * _r_dd_j * _τ_dd_i;  # ir + it-jr-it(v) + it-jr_dd-it
    end;

    if USE_STATIC_ARRAY
        return SVector{DIM_CANOPY+1,FT}(_r_dd), SVector{DIM_CANOPY+1,FT}(_r_sd), SVector{DIM_CANOPY,FT}(_t_dd), SVector{DIM_CANOPY,FT}(_t_sd)
    else
        return _r_dd, _r_sd, _t_dd, _t_sd
    end;
);








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
