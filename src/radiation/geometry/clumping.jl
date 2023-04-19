"""

    clumping_index(ci::ClumpingIndexPinty{FT}, sza::FT) where {FT}

Return the clumping index, given
- `ci` `ClumpingIndexPinty` model
- `sza` Solar zenith angle

"""
function clumping_index(ci::ClumpingIndexPinty{FT}, sza::FT) where {FT}
    (; Ω_A, Ω_B) = ci;

    return Ω_A + Ω_B * (1 - cosd(sza))
end
