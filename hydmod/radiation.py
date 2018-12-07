from hydmod.conversions import *
from hydmod.atmosphere import *

SOLAR_CONSTANT_MIN = 0.0820 # MJ m^-2 min^-1
SOLAR_CONSTANT_DAY = 118.1 # MJ m^-2 day^-1
LATENT_HEAT_VAPORIZATION_MJ = 0.408 # MJ/kg
LATENT_HEAT_VAPORIZATION = 2500.0 # KJ/kg
LATENT_HEAT_FUSION = 335.0 # KJ/kg
EMISSIVITY_SNOW = 0.97
DENSITY_WATER = 1000.0 # kg/m^3
DENSITY_AIR = 1.29 # kg/m^3
HEAT_CAPACITY_AIR = 1.0 # KJ/kg/C
HEAT_FROM_GROUND = 173 # kJ/m^2/day
STEFAN_BOLTZMANN_CONSTANT = 5.6697 * 10.0 ** (-8.0) # kW m^-2 K^-4

def ClearSkyRadiation(qo, transmissivity=0.75):
    """
    Predicted radiation at earth's surface

    Args:
        qo: potential solar radiation above Earth's atmosphere
        transmissivity: transmissivity factor (default: 0.75)

    Returns:
        estimate clear sky radiation

    """
    return np.multiply(transmissivity, qo)

def CloudCoverFraction(qd, qcs):
    """
    Estimate of the fraction of cloud cover using observed radiation and predicted clear sky radiation
    Args:
        qd: observed solar radiation
        qcs: predicted clear sky solar radiation

    Returns:
        fraction of cloud cover (cf)

    """

    cf = np.subtract(1.0, np.divide(qd, qcs))
    return np.where(cf<0.0, 0.0, cf)

def DirectSolarRadiation_Adjustment(qin, pfc=0.0, snowAlbedo=0.0, cc=0.0):
    adj = np.multiply((qin), np.multiply(np.subtract(1.0, np.divide(pfc, 100.0)), np.subtract(1.0, snowAlbedo)))
    return adj #kJ/m^2

def DirectSolarRadiation_SlopingSurface(lat, doy, qin, slope, aspect, units='radians'):
    """
    Direct solar radiation incident on a sloping surface (MJ/m^2)

    Args:
        lat: latitude
        doy: day of year
        qin: input solar radiation (kJ/m^2)
        slope: slope of land surface
        aspect: aspect of land surface
        units: units for lat and aspect; one of 'radians' (default) or 'degrees'

    Returns:
        direct solar radiation (kJ/m^2)

    """
    lat = ConvertToRadians(lat, units)
    aspect = ConvertToRadians(aspect, units)
    delta = SolarDeclination(doy) #solar declination angle
    phi = SolarElevationAngle(lat, delta) #solar elevation angle
    azs = SolarAzimuthAngle(phi, lat, delta) #solar azimuth angle
    i = SolarIncidenceAngle_2d(slope, aspect, azs, phi) #angle sun's rays make with sloping surface
    pt1 = np.divide(np.sin(i), np.sin(phi))
    pt1 = np.where(pt1 <= 0.0, 0.0, pt1) #[pt1 <= 0.0] = 0.0
    #pt2 = np.multiply((qcs), np.multiply(np.subtract(1.0, np.divide(pfc, 100.0)), np.subtract(1.0, snowAlbedo)))
    #qs = np.multiply(pt1, pt2)
    qs = np.multiply(pt1, qin)
    qs = np.where(qs < 0.1, 0.1, qs) #[qs < 0.1] = 0.1
    return qs #kJ/m^2

def EarthSunIRD(doy):
    """
    The inverse relative distance between the Earth and the Sun

    Args:
        doy: Day of year (1-365 or 366)

    Returns:
        Inverse relative distance Earth-Sun

    """
    dr = np.add(1.0, np.multiply(0.033, np.cos(np.multiply(np.divide(2*PI, 365), doy))))
    return dr

def Emissivity(tavg, qd, qo, fce=0.92):
    qcs = ClearSkyRadiation(qo)
    cf = CloudCoverFraction(qd, qcs)
    ea = Emissivity_CloudCover(tavg, cf)
    ets_pt1 = np.multiply(ea, np.subtract(1.0, cf))
    ets_pt2 = np.multiply(fce, cf)
    ets = np.add(ets_pt1, ets_pt2)
    return ets

def Emissivity_CloudCover(tavg, cf):
    """
    Change in emissivity due to cloud cover
    Args:
        tavg: Average daily air temperature
        cf: Fraction of cloud cover

    Returns:

    """
    pt1 = np.add(0.72, np.multiply(0.005, tavg))
    pt2 = np.subtract(1.0, np.multiply(0.84, cf))
    pt3 = np.multiply(0.84, cf)
    ea = np.add(np.multiply(pt1, pt2), pt3)
    return ea

def ExtraterrestrialRadiation(lat, doy, units='radians'):
    """
    Daily potential solar radiation

    Args:
        lat: Latitude (radians)
        doy: Day of year

    Returns:
        Extraterrestrial solar radiation in kJ/m^2/day (numpy array)

    """

    lat = ConvertToRadians(lat, units)
    dr = EarthSunIRD(doy)
    delta = SolarDeclination(doy)
    omegas = SunsetHourAngle(lat, delta)

    pt1 = np.multiply(((24.0*60.0)/PI)*SOLAR_CONSTANT_MIN, dr)
    pt2 = np.multiply(np.multiply(omegas, np.sin(lat)), np.sin(delta))
    pt3 = np.multiply(np.multiply(np.cos(lat), np.cos(delta)), np.sin(omegas))
    Ra = np.multiply(pt1, np.add(pt2, pt3))

    return Ra*1000 # kJ/m^2/day

def ExtraterrestrialRadiation_2d(lat, doy, units='radians'):
    """
    Daily potential solar radiation

    Args:
        lat: Latitude (degrees)
        doy: Day of year
        units: Angle units for phi. One of 'radians' (default) or 'degrees'.

    Returns:
        Extraterrestrial solar radiation (kJ/m^2/day)

    """
    if units == 'degrees':
        lat = DegreesToRadians(lat)

    dr = EarthSunIRD(doy)
    delta = SolarDeclination(doy)
    omegas = SunsetHourAngle(lat[np.newaxis, :, :], delta[:, np.newaxis, np.newaxis])

    pt1 = np.multiply(((24.0*60.0)/PI)*SOLAR_CONSTANT_MIN, dr)
    pt2 = np.multiply(np.multiply(omegas, np.sin(lat)[np.newaxis, :, :]), np.sin(delta)[:, np.newaxis, np.newaxis])
    pt3 = np.multiply(np.multiply(np.cos(lat)[np.newaxis, :, :], np.cos(delta)[:, np.newaxis, np.newaxis]), np.sin(omegas))
    Ra = np.multiply(pt1[:, np.newaxis, np.newaxis], np.add(pt2, pt3))

    return Ra*1000 # kJ/m^2/day

def HalfDayLength(lat, delta):
    """

    Args:
        lat: latitude (radians)
        delta: solar declination (radians)

    Returns:
        half day length (radians)

    """
    hdl = np.arccos(np.multiply(-1.0, np.multiply(np.tan(delta), np.tan(lat))))
    return hdl

def LatentRadiation(td, tsnow, windv, pfc):
    vda = VaporDensity(VaporPressure(td), td)
    vdsnow = VaporDensity(VaporPressure(tsnow), tsnow)
    wr = WindRoughness(windv, pfc)
    qlatent = np.multiply(LATENT_HEAT_VAPORIZATION, np.divide(np.subtract(vda, vdsnow), np.divide(wr, 3600.0*24.0))) # kJ/m^2
    return qlatent

def LongwaveRadiation(tavg, pfc, cc, fce=0.92):
    """
    Net longwave radiation (KJ m^2)
    Args:
        tavg: average daily temperature
        qd: observed shorwave radiation
        qo: potential shortwave radiation
        es: emissivity of snow (default: 0.98)
        ts: temperature of snow (default: 0.0 C)
        fce: forest cover emissivity (default: 0.92)

    Returns:
        net longwave radiation (KJ/m^2)
    """
    pt1 = (STEFAN_BOLTZMANN_CONSTANT * 3600.0 * 24.0) / 1000.0
    intm = ((0.72+0.005*tavg)*(1-0.84*cc)+0.84*cc)*(1-pfc/100.0)+fce*pfc/100.0
    pt2 = intm*np.power(CelciusToKelvin(tavg), 4.0)-EMISSIVITY_SNOW*np.power(CelciusToKelvin(0.0), 4.0)
    qlw = pt1*pt2
    return qlw # KJ/m^2

def SensibleRadiation(tavg, windr):
    qsens = np.multiply(HEAT_CAPACITY_AIR, np.multiply(DENSITY_AIR, np.divide(tavg, np.divide(windr, 3600.0*24.0))))
    return qsens # KJ/m^2

def SnowAlbedo(snowage, swe=0.0, exponent=-0.1908):
    """
    Albedo value for snow as a function of days since last snowfall (snowage)

    Args:
        snowage: days since last snowfall (days)
        swe: snow water equivalent (default: 0.0)
        exponent: exponential decay parameter (default: -0.1908)

    Returns:
        snow albedo value (alpha)

    """

    alpha = np.where(snowage == 0.0, 0.738, 0.0)
    alpha = np.where((alpha == 0.0) & (swe > 0.0), np.multiply(0.738, np.power(snowage, exponent)), alpha)
    return alpha

def SolarAzimuthAngle(phi, lat, delta, tod=12, tsn=12, units='radians'):
    """
    Horizontal angle between due south and the sun

    Args:
        phi: solar elevation or altitude angle (radians)
        lat: latitude
        delta: solar declination angle (radians)
        tod: time of day (hours; default: 12)
        tsn: time of solar noon (hours; default: 12)
        units: units for lat, "radians" (default) or "degrees"

    Returns: Solar Azimuth Angle (radians)

    """

    lat = ConvertToRadians(lat, units)
    if len(phi.shape) == 3:
        delta = delta[:,None,None]
    top = np.subtract(np.multiply(np.sin(phi), np.sin(lat)), np.sin(delta))
    bottom = np.multiply(np.cos(phi), np.cos(lat))
    azs = np.arccos(np.round(np.divide(top, bottom), 6))
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
    #delta = np.multiply(0.4093, np.sin(np.multiply(np.divide(2*PI, 365), doy-80)))
    delta = np.multiply(0.4093, np.sin(np.subtract(np.multiply(np.divide(2 * PI, 365), doy), 1.39)))
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

    lat = ConvertToRadians(lat, units)
    if len(lat.shape) == 2 and len(delta.shape) == 1:
        delta = delta[:,None, None]
    pt1 = np.multiply(np.sin(lat), np.sin(delta))
    pt2 = np.multiply(np.cos(lat), np.cos(delta))
    pt3 = np.cos(np.divide(np.multiply(PI, np.subtract(tod, tsn)),12.0))
    sinphi = np.add(pt1, np.multiply(pt2,pt3))
    # print("Elevation angle", np.arcsin(sinphi))
    return np.arcsin(sinphi)

def SolarIncidenceAngle_2d(slope, aspect, az, phi):
    """
    The angle the suns' rays make with a horizontal surface

    Args:
        slope: Slope of land surface (percent)
        aspect: Aspect of sloping surface
        az: solar azimuth angle (radians)
        phi: solar elevation angle (radians)

    Returns: Solar incidence angle (radians)

    """

    i_pt1 = np.multiply(np.sin(np.arctan(slope/100.0)), np.multiply(np.cos(phi), np.cos(np.subtract(az, aspect))))
    i_pt2 = np.multiply(np.cos(np.arctan(slope/100.0)), np.sin(phi))
    i = np.arcsin(np.add(i_pt1, i_pt2))
    # print(i)
    i = np.where(i>0.0, i, 0.0)
    # print("Incidence angle", i, "az-aspect", az-aspect, "aspect", aspect, 'az', az, 'elevation angle', phi)
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
