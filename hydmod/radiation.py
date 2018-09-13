import numpy as np
import math
from conversions import *

SOLAR_CONSTANT = 0.0820 #MJ m^-2 min^-1

def DayOfYear(m, d, y):
    """
    Day of year (e.g. 1-365 or 366)

    Args:
        m: Month (1-12)
        d: Day (1-31)
        y: Year (e.g. 2012)

    Returns:
        Day of year (1-365 or 366)

    """
    j = d-32+round(275*(m/9))+2*round(3/(m+1))+round((m/100)-((y%4)/4)+0.975)
    return j

def ExtraterrestrialRadiation(dr, omegas, phi, delta, doy):
    """
    Daily extraterrestrial radiation

    Args:
        dr: Inverse relative distance Earth-Sun
        omegas: Sunset hour angle (degrees)
        phi: Latitude (degrees)
        delta: Solar declination (degrees)
        doy: Day of year

    Returns:
        Extraterrestrial solar radiation (numpy array)

    """
    omegas = DegreesToRadains(omegas)
    phi = DegreesToRadians(phi)
    delta = DegreesToRadians(delta)

    Ra = ((24.0*60.0)/PI)*SOLAR_CONSTANT*dr*\
         (omegas*math.sin(phi)*math.sin(delta)+math.cos(phi)*math.cos(delta)*math.sin(omegas))
    return Ra