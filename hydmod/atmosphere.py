import numpy as np
from hydmod.conversions import *

THERMODYNAMIC_VAPOR_CONSTANT = 0.4615 #kJ/kg/K
VON_KARMON_CONSTANT = 0.41

#wind roughness parameters
D = 0.0 # m - height of zero plane displacement
ZM = 0.001 # m - momentum roughness parameter
ZH = 0.0002 # m - heat and vapor roughness parameter
ZU = 2.0 # m height of wind measurements
ZT = 2.0 # m height of temp measurements

def RelativeHumidity(tc, td):
    """
    Relative humidity
    Args:
        tc: temperature (degrees C)
        td: dewpoint temperature (degrees C)

    Returns:
        relative humidity (proportion 0.0 - 1.0)

    """
    e = VaporPressure(td)
    es = VaporPressure(tc)
    rh = np.divide(e, es)
    return rh

def VaporPressure(tc):
    """
    Vapor pressure (kPa)
    Args:
        tc: temperature (degrees C)

    Returns:
        vapor pressure (kPa)

    """
    es = np.exp(np.divide(np.subtract(np.multiply(16.78, tc), 116.9), np.add(tc, 237.3)))
    return es # kPa

def VaporDensity(e, tc):
    vd = np.divide(np.divide(e, CelciusToKelvin(tc)), THERMODYNAMIC_VAPOR_CONSTANT)
    return vd #kg/m^3

def WindRoughness(windv, pfc):
    """
    Wind roughness (s/m)

    Args:
        windv: wind speed (m/s)
        pfc: percent forest cover (%)

    Returns:
        wind roughness (s/m)

    """
    pt1 = np.log(np.divide(np.add(np.subtract(ZU, D), ZM), ZM))
    pt2 = np.log(np.divide(np.add(np.subtract(ZT, D), ZH), ZH))
    pt3 = np.multiply(np.multiply(VON_KARMON_CONSTANT**2.0, windv), np.subtract(1.0, np.divide(pfc, 101.0)))
    wr = np.divide(np.multiply(pt1, pt2), pt3)
    return wr