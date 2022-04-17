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

const pCO2atm0        = 280e-6     # ppm pre-industrial pCO2
const k16_PANtoO      = 3.762      #
const k18_oceanmass   = 1.397e21   # kg Ocean total mass
const k_moles1atm     = 1.77e20    # Moles in 1 atm
const k_preindCO2atm  = 280e-6     # pre-industrial pCO2 (atm)

# Atmospheric composition from Sarmiento & Gruber (2006) Table 3.1.1 (which cites Weast & Astle 1982)
const k_atmmixrN2    = 0.78084     # N2 atmospheric mixing ratio (moles / moles dry air)
const k_atmmixrO2    = 0.20946     # O2 atmospheric mixing ratio (moles / moles dry air)


end