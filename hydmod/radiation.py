import math
from conversions import *

SOLAR_CONSTANT = 0.0820 #MJ m^-2 min^-1
LATENT_HEAT_VAPORIZATION = 0.408 #MJ/kg

def ExtraterrestrialRadiation(phi, doy, units='radians'):
    """
    Daily extraterrestrial radiation

    Args:
        phi: Latitude (degrees)
        doy: Day of year
        units: Angle units for phi. One of 'radians' (default) or 'degrees'.

    Returns:
        Extraterrestrial solar radiation in mm/day(numpy array)

    """
    if units == 'degrees':
        phi = DegreesToRadians(phi)

    dr = EarthSunIRD(doy)
    delta = SolarDeclination(doy)
    omegas = SunsetHourAngle(phi, delta)

    pt1 = np.multiply(((24.0*60.0)/PI)*SOLAR_CONSTANT, dr)
    pt2 = np.multiply(np.multiply(omegas, np.sin(phi)), np.sin(delta))
    pt3 = np.multiply(np.multiply(np.cos(phi), np.cos(delta)), np.sin(omegas))
    Ra = np.multiply(pt1, np.add(pt2, pt3))

    return Ra*LATENT_HEAT_VAPORIZATION

def ExtraterrestrialRadiation_2d(phi, doy, units='radians'):
    """
    Daily extraterrestrial radiation

    Args:
        phi: Latitude (degrees)
        doy: Day of year
        units: Angle units for phi. One of 'radians' (default) or 'degrees'.

    Returns:
        Extraterrestrial solar radiation in mm/day(numpy array)

    """
    if units == 'degrees':
        phi = DegreesToRadians(phi)

    dr = EarthSunIRD(doy)
    delta = SolarDeclination(doy)
    omegas = SunsetHourAngle(phi[np.newaxis, :, :], delta[:, np.newaxis, np.newaxis])

    pt1 = np.multiply(((24.0*60.0)/PI)*SOLAR_CONSTANT, dr)
    pt2 = np.multiply(np.multiply(omegas, np.sin(phi)[np.newaxis, :, :]), np.sin(delta)[:, np.newaxis, np.newaxis])
    pt3 = np.multiply(np.multiply(np.cos(phi)[np.newaxis, :, :], np.cos(delta)[:, np.newaxis, np.newaxis]), np.sin(omegas))
    Ra = np.multiply(pt1[:, np.newaxis, np.newaxis], np.add(pt2, pt3))

    return Ra*LATENT_HEAT_VAPORIZATION

def EarthSunIRD(doy):
    """
    The inverse relative distance Earth-Sun

    Args:
        doy: Day of year (1-365 or 366)

    Returns:
        Inverse relative distance Earth-Sun

    """
    dr = np.add(1.0, np.multiply(0.033, np.cos(np.multiply(np.divide(2*PI, 365), doy))))
    return dr

def SolarDeclination(doy, units='radians'):
    """
    Solar declination angle

    Args:
        doy: Day of year (1-365 or 366
        units One of 'radians' (default) or 'degrees'

    Returns:
        Solar declination angle

    """
    delta = np.multiply(0.4093, np.sin(np.subtract(np.multiply(np.divide(2*PI, 365), doy), 1.39)))

    if units == 'degrees':
        return RadiansToDegrees(delta)
    else:
        return delta

def SunsetHourAngle(phi, delta, units='radians'):
    """
    Sunset hour angle

    Args:
        phi: Latitude
        delta: Solar declination
        units: Specifies input and output units. One of 'radians' (default) or 'degrees'

    Returns:
        Sunset hour angle
    """
    if units == 'degrees':
        phi = DegreesToRadians(phi)
        delta = DegreesToRadians(delta)
    ntanphi = np.multiply(-1.0, np.tan(phi))
    tandelta = np.tan(delta)
    value = np.arccos(np.multiply(ntanphi, tandelta))

    if units == 'degrees':
        return RadiansToDegrees(value)
    else:
        return value