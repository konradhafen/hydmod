import math
from conversions import *

SOLAR_CONSTANT_MIN = 0.0820 #MJ m^-2 min^-1
SOLAR_CONSTANT_DAY = 118.1 #MJ m^-2 day^-1
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

    pt1 = np.multiply(((24.0*60.0)/PI)*SOLAR_CONSTANT_MIN, dr)
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

    pt1 = np.multiply(((24.0*60.0)/PI)*SOLAR_CONSTANT_MIN, dr)
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

def SolarAzimuthAngle(phi, lat, t=12, t0=12, units="radians"):
    """
    Horizontal angle between due south and the sun

    Args:
        phi: solar elevation or altitude angle
        lat: latitude
        t: time of day (hours; default=12)
        t0: time of solar noon (hours; default=12)
        units: units for phi and lat, "radians" (default) or "degrees"

    Returns: Solar Azimuth Angle (radians)

    """


def SolarDeclination(doy, units='radians'):
    """
    Solar declination angle

    Args:
        doy: Day of year (1-365 or 366)
        units One of 'radians' (default) or 'degrees'

    Returns:
        Solar declination angle

    """
    delta = np.multiply(0.4093, np.sin(np.subtract(np.multiply(np.divide(2*PI, 365), doy), 1.39)))

    if units == 'degrees':
        return RadiansToDegrees(delta)
    else:
        return delta

def SolarIncidenceAngle_2d(slope, aspect, units="radians"):
    """
    The angle the suns' rays make with a horizontal surface

    Args:
        slope: Slope of land surface
        aspect: Aspect of sloping surface
        units: "radians" (default) or "degrees", units must be in radians for calculation

    Returns: Solar incidence angle (radians)

    """

def SolarElevationAngle(lat, delta, tod=12, tsn=12, units='radians'):
    """
    The angle the sun's rays make with a horizontal surface. Defaults for tod and tsn calculate the daily value (i.e. value at solar noon)

    Args:
        lat: latitude
        delta: solar declination angle (radians)
        tod: time of day (hours; default=12)
        tsn: time of solar noon (hours; default=12)
        units: units for latitude, one of 'radians' (default) or 'degrees'

    Returns: Solar elevation angle (radians)

    """

    pt1 = np.multiply(np.sin(lat), np.sin(delta))
    pt2 = np.multiplty(np.cos(lat), np.cos(delta))
    pt3 = np.cos(np.divide(np.multiply(PI, np.subtract(tod-tsn)),12.0))
    sinphi = np.add(pt1, np.multiply(pt2,pt3))
    return np.arcsin(sinphi)

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