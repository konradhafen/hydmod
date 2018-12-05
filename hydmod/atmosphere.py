import numpy as np

#wind roughness parameters
D = 0.0 # m - height of zero plane displacement
ZM = 0.001 # m - momentum roughness parameter
ZH = 0.0002 # m - heat and vapor roughness parameter
ZU = 2.0 # m height of wind measurements
ZT = 2.0 # m height of temp measurements

def ActualVaporPressure(td):
    """
    Actual vapor pressure (kPa)
    Args:
        td: dewpoint temperature (degrees C)

    Returns:
        Actual vapor pressure (kPa)

    """
    e = np.multiply(6.11, np.power(10.0, np.divide(np.multiply(7.5, td), np.add(237.7, td)))) # millibars
    e = np.divide(e, 10.0) # kPa
    return e

def RelativeHumidity(tc, td):
    """
    Relative humidity
    Args:
        tc: temperature (degrees C)
        td: dewpoint temperature (degrees C)

    Returns:
        relative humidity (proportion 0.0 - 1.0)

    """
    e = ActualVaporPressure(td)
    es = SaturationVaporPressure(tc)
    rh = np.divide(e, es)
    return rh

def SaturationVaporPressure(tc):
    """
    Saturation vapor pressure (kPa)
    Args:
        tc: temperature (degrees C)

    Returns:
        Saturation vapor pressure (kPa)

    """
    es = np.exp(np.divide(np.subtract(np.multiply(16.78, tc), 116.9), np.add(tc, 237.3)))
    return es # kPa

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