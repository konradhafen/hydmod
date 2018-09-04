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
    melt[melt < 0] = 0
    return melt

def precipAsSnow(precip, temp, train=3.0, tsnow=0.0):
    """
    Determine the amount of precipitation falling as snow

    Args:
        precip: Daily precipitation (numpy array)
        temp: Daily mean (or min or max) air temperature (numpy array)
        train: Temperature at which all precipitation becomes rain (default = 3.0 C)
        tsnow: Temperature at which all precipitation becomes snow (default = 0.0 C)

    Returns:

    """
    snow = np.where(temp < train, np.where(temp>tsnow, (temp-tsnow)/(train-tsnow), 0.0), 0.0)
    return snow