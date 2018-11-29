from hydmod.conversions import *

SOLAR_CONSTANT_MIN = 0.0820 # MJ m^-2 min^-1
SOLAR_CONSTANT_DAY = 118.1 # MJ m^-2 day^-1
LATENT_HEAT_VAPORIZATION_MJ = 0.408 # MJ/kg
LATENT_HEAT_VAPORIZATION = 2500.0 #kJ/kg
LATENT_HEAT_FUSION = 335.0 #KJ/kg
EMISSIVITY_SNOW = 0.97
DENSITY_WATER = 1000.0 #kg/m^3
STEFAN_BOLTZMANN_CONSTANT = 5.6697 * 10.0 ** (-8.0) # kW m^-2 K^-4
VON_KARMON_CONSTANT = 0.41

#wind roughness parameters
D = 0.0 # m - height of zero plane displacement
ZM = 0.001 # m - momentum roughness parameter
ZH = 0.0002 # m - heat and vapor roughness parameter
ZU = 2.0 # m height of wind measurements
ZT = 2.0 # m height of temp measurements

def ClearSkyRadiation(qo, transmissivity=0.75):
    """
    Predicted radiation at earth's surface

    Args:
        qo: potential solar radiation above Earth's atmosphere
        transmissivity: transmissivity factor (default: 0.75)

    Returns:
        estimate clear sky radiation

    """
    return transmissivity*qo

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
    return cf

def DirectSolarRadiation(lat, doy, slope, aspect, swe=0.0, dls=0.1, cf=0.0, units='radians'):
    """
    Direct solar radiation incident on a sloping surface (MJ/m^2)

    Args:
        lat: latitude
        doy: day of year
        slope: slope of land surface
        aspect: aspect of land surface
        swe: snow water equivalent (default: 0.0)
        dls: days since last snow (default: 0.1)
        cf: forest canopy factor (default: 0)
        units: units for lat, slope and aspect; one of 'radians' (default) or 'degrees'

    Returns:
        direct solar radiation (MJ/m^2)

    """

    lat = ConvertToRadians(lat, units)
    aspect = ConvertToRadians(aspect, units)
    delta = SolarDeclination(doy) #solar declination angle
    phi = SolarElevationAngle(lat, delta) #solar elevation angle
    azs = SolarAzimuthAngle(phi, lat, delta) #solar azimuth angle
    qo = ExtraterrestrialRadiation(lat, doy) #solar radiation incident on a flat surface
    qcs = ClearSkyRadiation(qo)
    alpha = SnowAlbedo(dls, swe) #snow albedo
    i = SolarIncidenceAngle_2d(slope, aspect, azs, phi) #angle sun's rays make with sloping surface
    pt1 = np.divide(np.sin(i), np.sin(phi))
    pt1 = np.where(pt1 <= 0.0, 0.0, pt1) #[pt1 <= 0.0] = 0.0
    pt2 = np.multiply(qcs, np.multiply(np.subtract(1.0, cf), np.subtract(1.0, alpha)))
    qs = np.multiply(pt1, pt2)
    qs = np.where(qs < 0.1, 0.1, qs) #[qs < 0.1] = 0.1
    # print('q0', qo, 'qs', qs)
    return qs

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

def ExtraterrestrialRadiation(lat, doy):
    """
    Daily potential solar radiation

    Args:
        lat: Latitude (radians)
        doy: Day of year

    Returns:
        Extraterrestrial solar radiation in mm/day(numpy array)

    """

    dr = EarthSunIRD(doy)
    # print('ird to sun', dr)
    delta = SolarDeclination(doy)
    # print('solar declination', delta)
    omegas = SunsetHourAngle(lat, delta)
    # print('sunset hour angle', omegas)

    pt1 = np.multiply(((24.0*60.0)/PI)*SOLAR_CONSTANT_MIN, dr)
    pt2 = np.multiply(np.multiply(omegas, np.sin(lat)), np.sin(delta))
    pt3 = np.multiply(np.multiply(np.cos(lat), np.cos(delta)), np.sin(omegas))
    Ra = np.multiply(pt1, np.add(pt2, pt3))

    return Ra #*LATENT_HEAT_VAPORIZATION_MJ

def ExtraterrestrialRadiation_2d(lat, doy, units='radians'):
    """
    Daily potential solar radiation

    Args:
        lat: Latitude (degrees)
        doy: Day of year
        units: Angle units for phi. One of 'radians' (default) or 'degrees'.

    Returns:
        Extraterrestrial solar radiation in mm/day(numpy array)

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

    return Ra # *LATENT_HEAT_VAPORIZATION_MJ

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

def LongwaveRadiation(tavg, qd, qo, es=0.98, ts=0.0, fce=0.92):
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
    PercentForestCover = 80.0
    CloudCover = 0.038948
    tavg = 7.8
    pt1 = (STEFAN_BOLTZMANN_CONSTANT * 3600 * 24) / 1000.0
    intm = ((0.72+0.005*tavg)*(1-0.84*CloudCover)+0.84*CloudCover)*(1-PercentForestCover/100.0)+0.92*PercentForestCover/100.0
    pt2 = intm*np.power(CelciusToKelvin(tavg), 4.0)-EMISSIVITY_SNOW*np.power(CelciusToKelvin(0.0), 4.0)
    qlw = pt1*pt2
    # ets = Emissivity(tavg, qd, qo, fce)
    # qlw_pt1 = np.multiply(np.power(CelciusToKelvin(tavg), 4.0), np.multiply(ets, STEFAN_BOLTZMANN_CONSTANT))
    # qlw_pt2 = np.multiply(np.power(ts, 4.0), np.multiply(es, STEFAN_BOLTZMANN_CONSTANT))
    # qlw = np.subtract(qlw_pt1, qlw_pt2)
    return qlw # KJ/m^2

def SnowAlbedo(dsls, swe=0.0, exponent=-0.1908):
    """
    Albedo value for snow as a function of days since last snowfall (dsls)

    Args:
        dsls: days since last snowfall (days)
        swe: snow water equivalent (default: 0.0)
        exponent: exponential decay parameter (default: -0.1908)

    Returns:
        snow albedo value (alpha)

    """

    if swe > 0.0:
        alpha = np.multiply(0.738, np.power(dsls, exponent))
    else:
        alpha = 0.0
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

def WindRoughness(windv):
    pt1 = np.log(np.divide(np.add(np.subtract(ZU, D), ZM)), ZM)
    pt2 = np.log(np.add(np.subtract(ZT, D), ZH))
    pt3 = np.multiply(VON_KARMON_CONSTANT**2.0, np.multiply(windv, np.divide(np.subtract(1, D), 101.0)))
    wr = np.divide(np.multiply(pt1, pt2), pt3)
    return wr