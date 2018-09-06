import numpy as np
import os

def meltDegreeDay_USACE(k, temp, tbase=0.0):
    """
    Degree day snowmelt with US Army Corps of Engineers empirical model

    Args:
        k: Empirically derived constant with units length/degree/day
        temp: Daily mean (or min or max) air temperature (numpy array)
        tbase: Base temperature. USACE recommends 0.0 (C) for forested areas and -4.4 (C) for open areas (default = 0.0)

    Returns:
        numpy array of daily melt values

    """
    melt = k*(temp-tbase)
    melt[melt < 0.0] = 0.0
    return melt

def modelSWE(ppt_snow, swe_melt):
    """
    Change in snow water equivalent when accounting for melt

    Args:
        ppt_snow: Precipitation that falls as snow
        swe_melt: Snowmelt

    Returns:
        Cumulative snow water equivalent (SWE), and the actual amount of snow that melted

    """
    swe_cum = np.zeros(ppt_snow.shape, dtype=np.float32)
    act_melt = np.zeros(swe_melt.shape, dtype=np.float32)
    for i in range(1,ppt_snow.shape[0]):
        swe_inc = swe_cum[i-1] + ppt_snow[i] - swe_melt[i]
        if swe_inc > 0.0:
            swe_cum[i] = swe_inc
            act_melt[i] = swe_melt[i]
        else:
            act_melt[i] = swe_cum[i-1]

    return swe_cum, act_melt

def precipDaily(acprecip):
    """
    Calculate daily precipitation from accumulated precipitation. Assumes accumulated precipitation for day preceding record is 0
    Args:
        acprecip: Accumulated precipitation (numpy array 1D)

    Returns:
        Daily precipitation as a 1D numpy array

    """

    acprecip = np.insert(acprecip, 0, 0)
    return(np.ediff1d(acprecip))

def precipInput(ppt_rain, swe_melt):
    """
    Rain plus melted snow

    Args:
        ppt_rain: Precipitation falling as rain
        swe: Cumulative snow water equivalent
        swe_melt: Melted snow water equivalent

    Returns:
        Melted SWE plus rainfall

    """
    melt = np.where(np.subtract(swe, swe_melt)>0.0, swe_melt, swe)
    ppt_in = np.add(melt, ppt_rain)
    return ppt_in

def precipPhase(precip, temp, train=3.0, tsnow=0.0):
    """
    Determine the amount of precipitation falling as snow

    Args:
        precip: Daily precipitation (numpy array)
        temp: Daily mean (or min or max) air temperature (numpy array)
        train: Temperature at which all precipitation becomes rain (default = 3.0 C)
        tsnow: Temperature at which all precipitation becomes snow (default = 0.0 C)

    Returns:
        Daily precipitation falling as snow

    """
    snow = np.where(temp < train, np.where(temp>tsnow, np.multiply(precip, np.divide(np.subtract(temp,tsnow),
                                                                                     (train-tsnow))), precip), 0.0)
    rain = np.subtract(precip, snow)
    return snow, rain