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

#calculate ET radiation
Ra2d = rad.ExtraterrestrialRadiation_2d(np.full((nrow, ncol),conv.DegreesToRadians(41.97)), doy)

#calculate PET
pet2d = et.PET_Hargreaves1985(tmax2d, tmin2d, tavg2d, Ra2d)/1000.0 # m/day

s = np.zeros(ppt_in2d.shape)
r = np.zeros(ppt_in2d.shape)
ra = np.zeros(ppt_in2d.shape) #accumulated runoff
hwt = np.zeros(ppt_in2d.shape)
qlat_out = np.zeros(ppt_in2d.shape)
qlat_in = np.zeros(ppt_in2d.shape)
perc = np.zeros(ppt_in2d.shape)
sb = np.zeros(ppt_in2d.shape)
bf = np.zeros(ppt_in2d.shape)
q = np.zeros(ppt_in2d.shape)
r[0,:,:] = 0 # set initial runoff m
ra[0,:,:] = 0
s[0,:,:] = 0.3 # set initial storage (i.e water content) m
sb[0,:,:] = 0 # set initial aquifer storage (baseflow source) m
soildepth = np.full((nrow, ncol), 1.0) #depth of soil profile m

aet = np.zeros(ppt_in2d.shape)
ksat = np.full((nrow, ncol), 1.0) # m/day
slpds = gdal.Open(dirpath + "/testslopeper.tif")
geot = slpds.GetGeoTransform()
slope = slpds.GetRasterBand(1).ReadAsArray()
# print("slope", slope)
por = np.full((nrow, ncol), 0.5)
fc = np.full((nrow, ncol), 0.3)
wp = np.full((nrow, ncol), 0.1)
ksub = np.full((nrow, ncol), 0.001) # m/day
alpha = 0.02
porl = np.multiply(por, soildepth)
fcl = np.multiply(fc, soildepth)
wpl = np.multiply(wp, soildepth)
smax = porl
qlat_nr = np.zeros(ndays)

fprop = fr.FlowProportions(demnp)
aet[0,:,:] = et.ET_theta_2d(pet2d[0,:,:], fcl, wpl, s[0,:,:])

for i in range(1, ppt_in2d.shape[0]):
    s[i, :, :] = s[i - 1, :, :] + ppt_in2d[i,:,:]
    hwt[i, :, :] = gw.WaterTableHeight(por, fc, np.divide(s[i, :, :], soildepth), soildepth)
    qlat_out[i, :, :] = gw.LateralFlow_Darcy_2d(ksat, slope, hwt[i, :, :], geot[1], geot[1])
    qlat_in[i, :, :], qlat_nr[i] = fr.RouteFlow(fprop, qlat_out[i, :, :])
    s[i, :, :] = s[i, :, :] + qlat_in[i,:,:] - qlat_out[i,:,:]
    aet[i, :, :] = et.ET_theta_2d(pet2d[i, :, :], fcl, wpl, s[i - 1, :, :])
    s[i,:,:] = s[i,:,:] - aet[i,:,:]
    perc[i, :, :] = gw.Percolation_2d(ksub, gw.WaterTableHeight(por, fc, np.divide(s[i, :, :], 1.0), soildepth))
    s[i, :, :] = s[i, :, :] - perc[i,:,:]
    r[i, :, :] = np.where(s[i, :, :] > (soildepth*por), (s[i, :, :] - (soildepth*por)), 0.0)
    ra[i, :, :] = rd.FlowAccumulation(rd.LoadGDAL(dempath), method='D8', weights=r[i,:,:])
    # ra[i, :, :] = intm
    s[i, :, :] = np.where(s[i, :, :] > (soildepth * por), (soildepth * por), s[i, :, :])

# print("ppt", ppt_in2d)
# print("pet", pet)
# print("pet", pet2d)
# print("et", et1)
# print("et", et)
# rd.FillDepressions(rddem, in_place=True)
# accum = rd.FlowAccumulation(rddem, method="Dinf")
# dinf_fig = rd.rdShow(accum, axes=False, cmap='jet')
# accum = rd.FlowAccumulation(rddem, method="D8")
# d8_fig = rd.rdShow(accum, axes=False, cmap='jet')
# print("ppt", ppt_in2d)
# print("hwt", hwt)
# print("qlat out", qlat_out)
# print("qlat in", qlat_in)
# print("r", r)
# print("s", s)
plt.plot(date, qlat_in[:,4,2], 'g', date, ra[:,4,2], 'c', date, (r[:,4,2]+qlat_in[:,4,2]), 'b')
plt.show()
plt.plot(date, hwt[:,4,2], 'b')
plt.show()
plt.plot(date, s[:,4,2], 'b')
plt.show()

