"""
---

    extinction_coefficient(sza::FT, lia::FT) where {FT}

Return the extinction coefficient for direct radiation, given
- `sza` Solar zenith angle in `°`
- `lia` Leaf inclination angle in `°`

---

    extinction_coefficient(lia::FT) where {FT}

Return the extinction coefficient for diffuse radiation (unifrom 18 average angles from 2.5° to 87.5°), given
- `sza` Solar zenith angle in `°`

"""
function extinction_coefficient end

extinction_coefficient(sza::FT, lia::FT) where {FT} = (
    _π = FT(pi);

    _Cs = cosd(lia) * cosd(sza);
    _Ss = sind(lia) * sind(sza);
    _βs = (_Cs >= _Ss ? _π : acos(-_Cs/_Ss));

    return 2/_π / cosd(sza) * (_Cs * (_βs - _π/2) + _Ss * sin(_βs))
);

extinction_coefficient(lia::FT) where {FT} = (
    # compute the mean extinction coefficient for diffuse solar radiation from 18 angles
    _kd::FT = 0;
    for _sza in 0:5:89
        _kd += extinction_coefficient(_sza + FT(2.5), lia) * sind(_sza);
    end;

    return _kd / 18
);
