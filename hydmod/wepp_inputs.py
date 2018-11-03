import hydmod.precip as precip
import hydmod.et as et
import hydmod.stats as stats
import hydmod.conversions as conv
import hydmod.groundwater as gw
import hydmod.flow_routing as fr
import hydmod.radiation as rad
from datetime import datetime
import pandas as pd
from osgeo import gdal
import numpy as np
import math
import matplotlib.pyplot as plt
import richdem as rd
import os

os.chdir('C:/temp/smr/')
dempath = 'dem/dem.tif'
demds = gdal.Open(dempath)
geot = demds.GetGeoTransform()
demnp = demds.GetRasterBand(1).ReadAsArray()
slpds = gdal.Open('export/arcmap/gtiffs/FVSLOP.tif')
slpnp = slpds.GetRasterBand(1).ReadAsArray()
aspds = gdal.Open('export/arcmap/gtiffs/TASPEC.tif')
aspnp = aspds.GetRasterBand(1).ReadAsArray()
clipath = 'climate/265191.cli'
clidat = np.genfromtxt(clipath, skip_header=15)
d = clidat[:,0]
m = clidat[:,1]
y = clidat[:,2]
ppt = clidat[:,3]
tmax = clidat[:,7]
tmin = clidat[:,8]

nrow = demds.RasterYSize
ncol = demds.RasterXSize

ndays = 365
doy = conv.DayOfYear(m, d, y)
ppt2d = np.reshape(np.repeat(ppt, nrow*ncol), (ndays-1, nrow, ncol))
tmin2d = np.reshape(np.repeat(tmin, nrow*ncol), (ndays-1, nrow, ncol))
tmax2d = np.reshape(np.repeat(tmax, nrow*ncol), (ndays-1, nrow, ncol))
tavg2d = 0.5 * (tmin2d + tmax2d)

#from optimization in R
k = 1.16
tbase = 0.
train = 3.
tsnow = 0.

#calculate melted snow
swe_melt2d = np.apply_along_axis(precip.MeltDegreeDay_USACE, 0, tavg2d, k=k, tbase=tbase)
ppt_snow2d, ppt_rain2d = precip.PrecipPhase_2d(ppt2d, tavg2d, train, tsnow)

swe_mod2d, act_melt2d = precip.ModelSWE_2d(ppt_snow2d, swe_melt2d)

ppt_in2d = np.add(ppt_rain2d, act_melt2d)

