import math
from conversions import *

SOLAR_CONSTANT_MIN = 0.0820 #MJ m^-2 min^-1
SOLAR_CONSTANT_DAY = 118.1 #MJ m^-2 day^-1
LATENT_HEAT_VAPORIZATION = 0.408 #MJ/kg

def DirectSolarRadiation(lat, doy, slope, aspect, dls=0.0, cf=0.0, units='radians'):
    """
    Direct solar radiation incident on a sloping surface

    Args:
        lat: latitude
        doy: day of year
        dls: days since last snow (default=0)
        cf: forest canopy factor (default=0)
        units: units for lat, slope and aspect; one of 'radians' (default) or 'degrees'

    Returns: Solar radiation

    """

    lat = ConvertToRadians(lat, units)
    slope = ConvertToRadians(slope, units)
    aspect = ConvertToRadians(aspect, units)
    delta = SolarDeclination(doy) #solar declination angle
    phi = SolarElevationAngle(lat, delta) #solar elevation angle
    azs = SolarAzimuthAngle(phi, lat, delta) #solar azimuth angle
    qd = ExtraterrestrialRadiation(lat, doy, 'radians') #solar radiation incident on a flat surface
    alpha = SnowAlbedo(dls) #snow albedo
    i = SolarIncidenceAngle_2d(slope, aspect, azs, phi) #angle sun's rays make with sloping surface

    pt1 = np.divide(np.sin(i), np.sin(phi))
    pt1[pt1 <= 0] = 0.0
    pt2 = np.multiply(qd, np.multiply(np.subtract(1.0, cf), np.subtract(1.0, alpha)))
    qs = np.multiply(pt1, pt2)
    qs[qs < 0.1] = 0.1
    return qs

def ExtraterrestrialRadiation(phi, doy, units='radians'):
    """
    Daily extraterrestrial radiation

    Args:
        phi: Latitude
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

def SnowAlbedo(dls):
    """
    Albedo of snow by days since last snowfall (dls)

    Args:
        dls: days since last snowfal

    Returns: Albedo

    """

    alpha = np.multiply(0.738, np.power(dls, -0.1908))
    return alpha

def SolarAzimuthAngle(phi, lat, delta, tod=12, tsn=12, units="radians"):
    """
    Horizontal angle between due south and the sun

    Args:
        phi: solar elevation or altitude angle (radians)
        lat: latitude
        delta: solar declination angle (radians)
        tod: time of day (hours; default=12)
        tsn: time of solar noon (hours; default=12)
        units: units for lat, "radians" (default) or "degrees"

    Returns: Solar Azimuth Angle (radians)

    """

    lat = ConvertToDegrees(lat, units)
    top = np.subtract(np.multiply(np.sin(phi), np.sin(lat)), np.sin(delta))
    bottom = np.multiply(np.cos(phi), np.cos(lat))
    azs = np.arccos(np.divide(top, bottom))
    azs = np.where(tod>tsn, np.add(PI, azs), np.subtract(PI, azs))
    return azs

def SolarDeclination(doy):
    """
    Solar declination angle

    Args:
        doy: Day of year (1-365 or 366)

    Returns:
        Solar declination angle (radians)

    """
    delta = np.multiply(0.4093, np.sin(np.subtract(np.multiply(np.divide(2*PI, 365), doy), 1.39)))
    return delta

def SolarElevationAngle(lat, delta, tod=12.0, tsn=12.0, units='radians'):
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

    lat = ConvertToDegrees(lat, units)
    pt1 = np.multiply(np.sin(lat), np.sin(delta))
    pt2 = np.multiplty(np.cos(lat), np.cos(delta))
    pt3 = np.cos(np.divide(np.multiply(PI, np.subtract(tod-tsn)),12.0))
    sinphi = np.add(pt1, np.multiply(pt2,pt3))
    return np.arcsin(sinphi)

def SolarIncidenceAngle_2d(slope, aspect, az, phi, units="radians"):
    """
    The angle the suns' rays make with a horizontal surface

    Args:
        slope: Slope of land surface
        aspect: Aspect of sloping surface
        az: solar azimuth angle (radians)
        phi: solar elevation angle (radians)
        units: units for slope and aspect "radians" (default) or "degrees"

    Returns: Solar incidence angle (radians)

    """

    slope = ConvertToRadians(slope, units)
    aspect = ConvertToRadians(aspect, units)
    i = np.arcsin(np.multiply(np.sin(np.arctan(slope)), np.multiply(np.cos(phi), np.cos(np.subtract(az, slope)))))
    np.where(i>0.0, i, 0.0)
    return i

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