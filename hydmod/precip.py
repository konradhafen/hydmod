import numpy as np
import os

def meltDegreeDay_USACE(k, tbase, tmean):
    melt = k*(tmean-tbase)
    melt[melt < 0] = 0
    return melt