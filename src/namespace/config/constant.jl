"""

$(TYPEDEF)

Structure to save universal constants.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef struct UniversalConstants{FT}
    "Avogadro's number `[molecule mol⁻¹]`"
    AVOGADRO::FT = 6.02214076e23
    "Isobaric specific heat of dry air `[J kg⁻¹ K⁻¹]`"
    CP_D::FT = 1003.5
    "Isobaric specific heat of ice water `[J kg⁻¹ K⁻¹]`"
    CP_I::FT = 2108
    "Isobaric specific heat of liquid water `[J kg⁻¹ K⁻¹]`"
    CP_L::FT = 4181
    "Isobaric specific heat of water vapor `[J kg⁻¹ K⁻¹]`"
    CP_V::FT = 1859
    "O₂ fraction in air `[-]`"
    F_O₂::FT = 0.2095
    "Universal gas constant `[J mol⁻¹ K⁻¹]`"
    GAS_R::FT = 8.3144598
    "Gravity of the Earth `[m s⁻²]`"
    GRAVITY::FT = 9.81
    "Planck constant `[m² kg s⁻¹]`"
    H_PLANCK::FT = 6.626e-34
    "Boltzmann constant `[m² kg s⁻² K⁻¹]`"
    K_BOLTZMANN::FT = 1.381e-23
    "Stefan-Boltzmann constant `[W m⁻² K⁻⁴]`"
    K_STEFAN::FT = 5.670e-8
    "Von Karman constant `[-]`"
    K_VON_KARMAN::FT = 0.4
    "Latent heat vaporization at T₀ `[K kg⁻¹]`"
    LH_V₀::FT = 2.5008e6
    "Light speed in vacuum `[m s⁻¹]`"
    LIGHT_SPEED::FT = 2.99792458e8
    "Molar mass of dry air `[kg mol⁻¹]`"
    M_DRYAIR::FT = 28.97e-3
    "Molar mass of water `[kg mol⁻¹]`"
    M_H₂O::FT = 18.01528e-3
    "Mean atmospheric pressure at sea level `[Pa]`"
    P_ATM::FT = 1.01325e5
    "Water vapor pressure at triple temperature `[Pa]`"
    PRESS_TRIPLE::FT = 611.657
    "Freezing temperature of water `[K]`"
    T₀::FT = 273.15
    "Triple temperature of water `[K]`"
    T_TRIPLE::FT = 273.16
    "Mean number of days per year [day]"
    YEAR_D::FT = 365.2422222
    "Thermal conductivity of water `[W m⁻¹ K⁻¹]`"
    Λ_THERMAL_H₂O::FT = 0.57
    "Density of liquid water `[kg m⁻³]`"
    ρ_H₂O::FT = 1000
    "Broadband leaf thermal reflectance, related to blackbody emittance `[-]`"
    ρ_LEAF_LW::FT = 0.01
    "Broadband soil thermal reflectance, related to blackbody emittance `[-]`"
    ρ_SOIL_LW::FT = 0.06
    "Broadband leaf thermal transmission, related to blackbody emittance `[-]`"
    τ_LEAF_LW::FT = 0.01
end


const UNIVERSAL_CONSTANTS = UniversalConstants{Float64}();

""" Avogadro's number `[molecule mol⁻¹]` """
AVOGADRO(FT=Float64) = FT(UNIVERSAL_CONSTANTS.AVOGADRO);
AVOGADRO(uc::UniversalConstants) = uc.AVOGADRO;

""" Isobaric specific heat of dry air `[J kg⁻¹ K⁻¹]` """
CP_D(FT=Float64) = FT(UNIVERSAL_CONSTANTS.CP_D);;
CP_D(uc::UniversalConstants) = uc.CP_D;

""" Isobaric specific heat of dry air per mole `[J mol⁻¹ K⁻¹]` """
CP_D_MOL(FT=Float64) = CP_D(FT) * M_DRYAIR(FT);
CP_D_MOL(uc::UniversalConstants) = uc.CP_D * uc.M_DRYAIR;

""" Isobaric specific heat of ice water `[J kg⁻¹ K⁻¹]` """
CP_I(FT=Float64) = FT(UNIVERSAL_CONSTANTS.CP_I);
CP_I(uc::UniversalConstants) = uc.CP_I;

""" Isobaric specific heat of ice water per mole `[J mol⁻¹ K⁻¹]` """
CP_I_MOL(FT=Float64) = CP_I(FT) * M_H₂O(FT);
CP_I_MOL(uc::UniversalConstants) = uc.CP_I * uc.M_H₂O;

""" Isobaric specific heat of liquid water `[J kg⁻¹ K⁻¹]` """
CP_L(FT=Float64) = FT(UNIVERSAL_CONSTANTS.CP_L);
CP_L(uc::UniversalConstants) = uc.CP_L;

""" Isobaric specific heat of liquid water per mole `[J mol⁻¹ K⁻¹]` """
CP_L_MOL(FT=Float64) = CP_L(FT) * M_H₂O(FT);
CP_L_MOL(uc::UniversalConstants) = uc.CP_L * uc.M_H₂O;

""" Isobaric specific heat of water vapor `[J kg⁻¹ K⁻¹]` """
CP_V(FT=Float64) = FT(UNIVERSAL_CONSTANTS.CP_V);
CP_V(uc::UniversalConstants) = uc.CP_V;

""" Isobaric specific heat of water vapor per mole `[J mol⁻¹ K⁻¹]` """
CP_V_MOL(FT=Float64) = CP_V(FT) * M_H₂O(FT);
CP_V_MOL(uc::UniversalConstants) = uc.CP_V * uc.M_H₂O;

""" O₂ fraction in air `[-]` """
F_O₂(FT=Float64) = FT(UNIVERSAL_CONSTANTS.F_O₂);
F_O₂(uc::UniversalConstants) = uc.F_O₂;

""" Universal gas constant `[J mol⁻¹ K⁻¹]` """
GAS_R(FT=Float64) = FT(UNIVERSAL_CONSTANTS.GAS_R);
GAS_R(uc::UniversalConstants) = uc.GAS_R;

""" Gravity of the Earth `[m s⁻²]` """
GRAVITY(FT=Float64) = FT(UNIVERSAL_CONSTANTS.GRAVITY);
GRAVITY(uc::UniversalConstants) = uc.GRAVITY;

""" Planck constant `[m² kg s⁻¹]` """
H_PLANCK(FT=Float64) = FT(UNIVERSAL_CONSTANTS.H_PLANCK);
H_PLANCK(uc::UniversalConstants) = uc.H_PLANCK;

""" Boltzmann constant `[m² kg s⁻² K⁻¹]` """
K_BOLTZMANN(FT=Float64) = FT(UNIVERSAL_CONSTANTS.K_BOLTZMANN);
K_BOLTZMANN(uc::UniversalConstants) = uc.K_BOLTZMANN;

""" Stefan-Boltzmann constant `[W m⁻² K⁻⁴]` """
K_STEFAN(FT=Float64) = FT(UNIVERSAL_CONSTANTS.K_STEFAN);
K_STEFAN(uc::UniversalConstants) = uc.K_STEFAN;

""" Von Karman constant `[-]` """
K_VON_KARMAN(FT=Float64) = FT(UNIVERSAL_CONSTANTS.K_VON_KARMAN);
K_VON_KARMAN(uc::UniversalConstants) = uc.K_VON_KARMAN;

""" Latent heat vaporization at T₀ `[K kg⁻¹]` """
LH_V₀(FT=Float64) = FT(UNIVERSAL_CONSTANTS.LH_V₀);
LH_V₀(uc::UniversalConstants) = uc.LH_V₀;

""" Light speed in vacuum `[m s⁻¹]` """
LIGHT_SPEED(FT=Float64) = FT(UNIVERSAL_CONSTANTS.LIGHT_SPEED);
LIGHT_SPEED(uc::UniversalConstants) = uc.LIGHT_SPEED;

""" Molar mass of dry air `[kg mol⁻¹]` """
M_DRYAIR(FT=Float64) = FT(UNIVERSAL_CONSTANTS.M_DRYAIR);
M_DRYAIR(uc::UniversalConstants) = uc.M_DRYAIR;

""" Molar mass of water `[kg mol⁻¹]` """
M_H₂O(FT=Float64) = FT(UNIVERSAL_CONSTANTS.M_H₂O);
M_H₂O(uc::UniversalConstants) = uc.M_H₂O;

""" Mean atmospheric pressure at sea level `[Pa]` """
P_ATM(FT=Float64) = FT(UNIVERSAL_CONSTANTS.P_ATM);
P_ATM(uc::UniversalConstants) = uc.P_ATM;

""" Water vapor pressure at triple temperature `[Pa]` """
PRESS_TRIPLE(FT=Float64) = FT(UNIVERSAL_CONSTANTS.PRESS_TRIPLE);
PRESS_TRIPLE(uc::UniversalConstants) = uc.PRESS_TRIPLE;

""" Gas constant water vapor `[J kg⁻¹ K⁻¹]` """
R_V(FT=Float64) = GAS_R(FT) / M_H₂O(FT);
R_V(uc::UniversalConstants) = uc.GAS_R / uc.M_H₂O;

""" Gas constant times 298.15 K `[J mol⁻¹]` """
RT₂₅(FT=Float64) = GAS_R(FT) * T₂₅(FT);
RT₂₅(uc::UniversalConstants) = uc.GAS_R * (uc.T₀ + 25);

""" Freezing temperature of water `[K]` """
T₀(FT=Float64) = FT(UNIVERSAL_CONSTANTS.T₀);
T₀(uc::UniversalConstants) = uc.T₀;

""" Kelvin temperature at 25 Celcius `[K]` """
T₂₅(FT=Float64) = T₀(FT) + 25;
T₂₅(uc::UniversalConstants) = uc.T₀ + 25;

""" Triple temperature of water `[K]` """
T_TRIPLE(FT=Float64) = FT(UNIVERSAL_CONSTANTS.T_TRIPLE);
T_TRIPLE(uc::UniversalConstants) = uc.T_TRIPLE;

""" Molar volume of liqiud water """
V_H₂O(FT=Float64) = M_H₂O(FT) / ρ_H₂O(FT);
V_H₂O(uc::UniversalConstants) = uc.M_H₂O / uc.ρ_H₂O;

""" Mean number of days per year [day] """
YEAR_D(FT=Float64) = FT(UNIVERSAL_CONSTANTS.YEAR_D);
YEAR_D(uc::UniversalConstants) = uc.YEAR_D;

""" Thermal conductivity of water `[W m⁻¹ K⁻¹]` """
Λ_THERMAL_H₂O(FT=Float64) = FT(UNIVERSAL_CONSTANTS.Λ_THERMAL_H₂O);
Λ_THERMAL_H₂O(uc::UniversalConstants) = uc.Λ_THERMAL_H₂O;

""" Density of liquid water `[kg m⁻³]` """
ρ_H₂O(FT=Float64) = FT(UNIVERSAL_CONSTANTS.ρ_H₂O);
ρ_H₂O(uc::UniversalConstants) = uc.ρ_H₂O;

""" Broadband leaf thermal reflectance, related to blackbody emittance `[-]` """
ρ_LEAF_LW(FT=Float64) = FT(UNIVERSAL_CONSTANTS.ρ_LEAF_LW);
ρ_LEAF_LW(uc::UniversalConstants) = uc.ρ_LEAF_LW;

""" Broadband soil thermal reflectance, related to blackbody emittance `[-]` """
ρ_SOIL_LW(FT=Float64) = FT(UNIVERSAL_CONSTANTS.ρ_SOIL_LW);
ρ_SOIL_LW(uc::UniversalConstants) = uc.ρ_SOIL_LW;

""" Broadband leaf thermal transmission, related to blackbody emittance `[-]` """
τ_LEAF_LW(FT=Float64) = FT(UNIVERSAL_CONSTANTS.τ_LEAF_LW);
τ_LEAF_LW(uc::UniversalConstants) = uc.τ_LEAF_LW;

""" Density of water times gravity `[MPa m⁻¹]` """
ρg_MPa(FT=Float64) = ρ_H₂O(FT) * GRAVITY(FT) * FT(1e-6);
ρg_MPa(uc::UniversalConstants{FT}) where {FT} = uc.ρ_H₂O * uc.GRAVITY * FT(1e-6);
