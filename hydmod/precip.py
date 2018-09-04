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
    return 0