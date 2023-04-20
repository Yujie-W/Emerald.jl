"""

$(TYPEDEF)

Hierarchy of AbstractSoilVC:
- [`BrooksCorey`](@ref)
- [`VanGenuchten`](@ref)

"""
abstract type AbstractSoilVC{FT} end


"""

$(TYPEDEF)

Brooks Corey soil parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct BrooksCorey{FT} <:AbstractSoilVC{FT}
    # General model information
    "Maximum soil hydraulic conductivity at 25 °C `[mol m⁻¹ s⁻¹ MPa⁻¹]`"
    K_MAX::FT
    "Soil b"
    B::FT
    "Soil type"
    TYPE::String
    "Potential at saturation `[MPa]`"
    Ψ_SAT::FT
    "Saturated soil volumetric water content"
    Θ_SAT::FT
    "Residual soil volumetric water content"
    Θ_RES::FT
end


"""

$(TYPEDEF)

van Genuchten soil parameters

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct VanGenuchten{FT} <:AbstractSoilVC{FT}
    # General model information
    "Maximum soil hydraulic conductivity at 25 °C `[mol m⁻¹ s⁻¹ MPa⁻¹]`"
    K_MAX::FT
    "Soil n is Measure of the pore-size distribution"
    N::FT
    "Soil type"
    TYPE::String
    "Soil α is related to the inverse of the air entry suction, α > 0"
    α::FT
    "Residual soil volumetric water content"
    Θ_RES::FT
    "Saturated soil volumetric water content"
    Θ_SAT::FT

    # Parameters based on the ones above
    "Soil m = 1 - 1/n"
    M::FT = 1 - 1 / N
end

VanGenuchten{FT}(name::String) where {FT} = (
    # https://structx.com/Soil_Properties_007.html
    # Parameters from Loam soil
    _p = [ 367.3476, 1.56, 0.43, 0.078, 7.19e-6];

    # switch name
    if name=="Sand"
        _p = [1479.5945, 2.68, 0.43, 0.045, 1.76e-4];
    elseif name=="Loamy Sand"
        _p = [1265.3084, 2.28, 0.41, 0.057, 1.56e-4];
    elseif name=="Sandy Loam"
        _p = [ 765.3075, 1.89, 0.41, 0.065, 3.45e-5];
    elseif name=="Loam"
        _p = [ 367.3476, 1.56, 0.43, 0.078, 6.94e-6];
    elseif name=="Sandy Clay Loam"
        _p = [ 602.0419, 1.48, 0.39, 0.100, 6.31e-6];
    elseif name=="Silt Loam"
        _p = [ 204.0820, 1.41, 0.45, 0.067, 7.19e-6];
    elseif name=="Silt"
        _p = [ 163.2656, 1.37, 0.46, 0.034, 7.19e-6];
    elseif name=="Clay Loam"
        _p = [ 193.8779, 1.31, 0.41, 0.095, 2.45e-6];
    elseif name=="Silty Clay Loam"
        _p = [ 102.0410, 1.23, 0.43, 0.089, 1.70e-6];
    elseif name== "Sandy Clay"
        _p = [ 275.5107, 1.23, 0.38, 0.100, 2.17e-6];
    elseif name=="Silty Clay"
        _p = [  51.0205, 1.09, 0.36, 0.070, 1.02e-6];
    elseif name=="Clay"
        _p = [  81.6328, 1.09, 0.38, 0.068, 1.28e-6];
    else
        @warn "Soil type $(name) not recognized, use Loam instead.";
        name = "Loam";
    end;

    # return a new struct
    return VanGenuchten{FT}(K_MAX = _p[5] * ρ_H₂O() / M_H₂O() / ρg_MPa(), N = _p[2], TYPE = name, α = _p[1], Θ_RES = _p[4], Θ_SAT = _p[3])
);
