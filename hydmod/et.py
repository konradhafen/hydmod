import numpy as np
from radiation import *

def ET_theta(pet, fc, wp, wc):
    """
    Model evapotranspiration, multiply potential ET by factor based on soil water content
    Args:
        pet: Potential evapotranspiration
        fc: Field capacity
        wp: Wilting point
        wc: Water content

    Returns:
        ET estimate

    """
    theta = 1.0
    if wc < 0.8*fc and wc > wp:
        theta = (0.8*(fc - wc))/(0.8*fc-wp)
    elif wc <= wp:
        theta = 0.0
    #print("wc", wc, "fc", fc, "wp", wp, "theta", theta)
    return(pet*theta)

def PET_Hargreaves1985(tmax, tmin, tmean, Ra):
    """
    Estimate reference/potential evapotranspiration with the Hargreaves 1985 method.

    Args:
        tmax: Maximum daily air temperature (deg C)
        tmin: Minimum daily air temperature (deg C)
        tmean: Mean daily air temperature (deg C)
        Ra: Extra terrestrial radiation (mm/day)

    Returns:

    """
    pt1 = np.multiply(0.0023, np.power(np.subtract(tmax, tmin), 0.5))
    pt2 = np.multiply(np.add(tmean, 17.8), Ra)
    ETref = np.multiply(pt1, pt2)
    return ETref