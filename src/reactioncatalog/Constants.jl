module Constants

# Physical constants
const k_CtoK          = 273.15     # convert temperature in Celsius to Kelvin
const k_molVolIdealGas= 22.4136    # l/mol  molar volume of ideal gas at STP (0C, 1 atm)
const k_Rgas          = 8.3144621  # J/K/mol gas constant
const k_Avogadro      = 6.0221409e23 # molecules mol-1

# Earth system present-day constants
const age_present_yr  = 4.5e9       # Age of Earth at present-day
const k_solar_presentday  = 1368.0 # present-day solar insolation W/m^2
const k_secpyr        = 3.15569e7  # present-day seconds per year
const k_secpday       = 24.0*3600.0 # sec per day
const k_daypyr        = k_secpyr/k_secpday  # days in a year
const k_g_earth       = 9.80665     # gravitational field strength for Earth m/s^2
const k_SurfAreaEarth = 5.101e14 # Earth surface area in m

const pCO2atm0        = 280e-6     # ppm pre-industrial pCO2
const k16_PANtoO      = 3.762      #
const k18_oceanmass   = 1.397e21   # kg Ocean total mass
const k_moles1atm     = 1.77e20    # Moles in 1 atm
const k_preindCO2atm  = 280e-6     # pre-industrial pCO2 (atm)

# Atmospheric composition from Sarmiento & Gruber (2006) Table 3.1.1 (which cites Weast & Astle 1982)
const k_atmmixrN2    = 0.78084     # N2 atmospheric mixing ratio (moles / moles dry air)
const k_atmmixrO2    = 0.20946     # O2 atmospheric mixing ratio (moles / moles dry air)

"""
    STANDARD_ATOMIC_WEIGHTS

IUPAC recommended values of relative atomic masses of sources in the local environment of the Earth's crust and atmosphere
(ie with Earth-specific isotope composition)
"""
const STANDARD_ATOMIC_WEIGHTS = (
    H   = 1.0080,
    He  = 4.0026,
    Li  = 6.94,
    C   = 12.011,
    N   = 14.007,
    O   = 15.999,
    F   = 18.998,
    Na  = 22.990,
    P   = 30.974,
    S   = 32.06,
    Cl  = 35.45,
    Ar  = 39.95,
    K   = 39.098,
    Ti  = 47.867,
    V   = 50.942,
    Cr  = 51.996,
    Fe  = 55.845,
    Rb  = 85.468,
    Cs  = 132.91,
)

end