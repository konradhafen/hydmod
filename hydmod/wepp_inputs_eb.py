import hydmod.precip as precip
import hydmod.et as et
import hydmod.stats as stats
import hydmod.conversions as conv
import hydmod.groundwater as gw
import hydmod.flow_routing as fr
import hydmod.radiation as rad
import datetime
import pandas as pd
from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt
import richdem as rd
from contextlib import contextmanager
import sys, os

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

os.chdir('C:/temp/smr/temple_fork')
dempath = 'dem/dem.tif'
demds = gdal.Open(dempath)
geot = demds.GetGeoTransform()
demnp = demds.GetRasterBand(1).ReadAsArray()
nrow = demds.RasterYSize
ncol = demds.RasterXSize
demds = None
with suppress_stdout():
    demfil = rd.FillDepressions(rd.LoadGDAL(dempath))
    rdaccum = rd.FlowAccumulation(demfil, method='Dinf')

slpds = gdal.Open('export/arcmap/gtiffs/FVSLOP.tif')
slpnp = slpds.GetRasterBand(1).ReadAsArray()*100.0
aspds = gdal.Open('export/arcmap/gtiffs/TASPEC.tif')
aspnp = aspds.GetRasterBand(1).ReadAsArray()
clipath = 'climate/424856.cli'
clidat = np.genfromtxt(clipath, skip_header=15)

lat = 41.82 #latitude
lat2d = np.full((nrow, ncol), lat)
pfc = 0.0 #percent forest cover

#set simulation time frame
daystart = 1
dayend = 730
ndays = dayend-daystart
dayindex = np.arange(ndays)

#read input data
d = clidat[daystart-1:dayend-1, 0]
m = clidat[daystart-1:dayend-1, 1]
y = clidat[daystart-1:dayend-1, 2]+2000

df = pd.DataFrame({'year': y, 'month': m, 'day': d})
dt = pd.DatetimeIndex(pd.to_datetime(df))
#print(dt)

doy = conv.DayOfYear(m, d, y)
ppt = np.reshape(np.repeat(clidat[daystart-1:dayend-1, 3] * 0.0254, nrow*ncol), (ndays, nrow, ncol))
tmax =  np.reshape(np.repeat(clidat[daystart-1:dayend-1, 7], nrow*ncol), (ndays, nrow, ncol))
tmin =  np.reshape(np.repeat(clidat[daystart-1:dayend-1, 8], nrow*ncol), (ndays, nrow, ncol))
tavg = 0.5 * (tmin + tmax)
wind =  np.reshape(np.repeat(clidat[daystart-1:dayend-1, 10], nrow*ncol), (ndays, nrow, ncol))
td =  np.reshape(np.repeat(clidat[daystart-1:dayend-1, 12], nrow*ncol), (ndays, nrow, ncol))
radobs =  np.reshape(np.repeat(conv.LangleyTokJsqm(clidat[daystart-1:dayend-1, 9]), nrow*ncol), (ndays, nrow, ncol)) #langleys -> kJ/m^2
qcs = rad.ClearSkyRadiation(rad.ExtraterrestrialRadiation_2d(np.full((nrow, ncol),conv.DegreesToRadians(lat)), doy))
#radmod_max = rad.DirectSolarRadiation_SlopingSurface(lat, doy, qcs, slpnp, aspnp, 'degrees')
cc = rad.CloudCoverFraction(radobs, qcs)
print(doy.shape, ppt.shape, ndays)
#parameters for precip phase
train = 3.0
tsnow = 0.0

#calc precip phase and snow properties
snow, rain = precip.PrecipPhase_2d(ppt, tavg, train, tsnow)
snowage = precip.SnowAge(snow) #assumes snow never melts
snowalbedo = rad.SnowAlbedo(snowage, 1.0) #assumes snow never melts
radadj = rad.DirectSolarRadiation_SlopingSurface(lat2d, doy, radobs, slpnp, aspnp)
radadj = rad.DirectSolarRadiation_Adjustment(radadj, pfc, snowalbedo)

rl = rad.LongwaveRadiation(tavg, pfc, cc)
wr = rad.WindRoughness(wind, pfc)
qs = rad.SensibleRadiation(tavg, wr)
es = rad.VaporPressure(tavg)
e = rad.VaporPressure(td)
vda = rad.VaporDensity(e, td)
vds = rad.VaporDensity(es, tavg)
vdsnow = rad.VaporDensity(rad.VaporPressure(np.where(tavg < 0.0, tavg, 0.0)), np.where(tavg < 0.0, tavg, 0.0))
ql = rad.LatentRadiation(td, tsnow, wind, pfc)
qtotal = radadj + rl + qs + ql + rad.HEAT_FROM_GROUND
qtotal = np.where(qtotal>0.0, qtotal, 0.0)
maxmelt = qtotal/(rad.LATENT_HEAT_FUSION*rad.DENSITY_WATER)
swe, actmelt = precip.ModelSWE_2d(snow, maxmelt)
pptin = rain+actmelt

#recalculate albedo and radiation after swe calculated
snowalbedo = rad.SnowAlbedo(snowage, swe)
radadj = rad.DirectSolarRadiation_Adjustment(radobs, pfc, snowalbedo)
qtotal = radadj + rl + qs + ql + rad.HEAT_FROM_GROUND
qtotal = np.where(qtotal>0.0, qtotal, 0.0)
pet = et.PET_Hargreaves1985(tmax, tmin, tavg, qtotal/rad.LATENT_HEAT_VAPORIZATION)/1000.0 # m/day

#initiate arrays for variables
s = np.zeros((ndays, nrow, ncol))
r = np.zeros((ndays, nrow, ncol))
ra = np.zeros((ndays, nrow, ncol))  # accumulated runoff
hwt = np.zeros((ndays, nrow, ncol))
qlat_out = np.zeros((ndays, nrow, ncol))
qlat_in = np.zeros((ndays, nrow, ncol))
perc = np.zeros((ndays, nrow, ncol))
sb = np.zeros((ndays, nrow, ncol))
bf = np.zeros((ndays, nrow, ncol))
q = np.zeros((ndays, nrow, ncol))

#set initial values
r[0, :, :] = 0.0  # set initial runoff m
ra[0, :, :] = 0.0
s[0, :, :] = 0.3  # set initial storage (i.e water content) m
sb[0, :, :] = 0.0  # set initial aquifer storage (baseflow source) m
soildepth = np.full((nrow, ncol), 2.0)  # depth of soil profile m

# Set initial and constant/preset values for simulation
lat = np.full((nrow, ncol), 41.97)
pfc = np.full((nrow, ncol), 80.0)

#set initial values
aet = np.zeros((ndays, nrow, ncol))
ksat = np.full((nrow, ncol), 10.0) # m/day
por = np.full((nrow, ncol), 0.5)
fc = np.full((nrow, ncol), 0.3)
wp = np.full((nrow, ncol), 0.1)
ksub = np.full((nrow, ncol), 0.001) # m/day
alpha = 0.02 # baseflow recession coefficient
porl = np.multiply(por, soildepth)
fcl = np.multiply(fc, soildepth)
wpl = np.multiply(wp, soildepth)
smax = porl
qlat_nr = np.zeros(ndays)

fprop = fr.FlowProportions(demfil)

with suppress_stdout():
    for i in range(1, ppt.shape[0]):
        s[i, :, :] = s[i - 1, :, :] + ppt[i, :, :]
        hwt[i, :, :] = gw.WaterTableHeight(por, fc, np.divide(s[i, :, :], soildepth), soildepth)
        qlat_out[i, :, :] = gw.LateralFlow_Darcy_2d(ksat, slpnp / 100.0, hwt[i, :, :], geot[1], geot[1])
        qlat_in[i, :, :], qlat_nr[i] = fr.RouteFlow(fprop, qlat_out[i, :, :])
        s[i, :, :] = s[i, :, :] + qlat_in[i, :, :] - qlat_out[i, :, :]
        aet[i, :, :] = et.ET_theta_2d(pet[i, :, :], fcl, wpl, s[i - 1, :, :])
        s[i, :, :] = s[i, :, :] - aet[i, :, :]
        perc[i, :, :] = gw.Percolation_2d(ksub, gw.WaterTableHeight(por, fc, np.divide(s[i, :, :], 1.0), soildepth))
        s[i, :, :] = s[i, :, :] - perc[i, :, :]
        r[i, :, :] = np.where(s[i, :, :] > (soildepth * por), (s[i, :, :] - (soildepth * por)), 0.0)
        ra[i, :, :] = rd.FlowAccumulation(rd.LoadGDAL(dempath), method='Dinf', weights=r[i, :, :])
        s[i, :, :] = np.where(s[i, :, :] > (soildepth * por), (soildepth * por), s[i, :, :])

#basin outlet
outletrow = 130
outletcol = 34
#another point
outrow = 13
outcol = 17

plt.plot(dayindex, qlat_in[:, outrow, outcol], 'b', dayindex, qlat_out[:, outrow, outcol], 'r')
plt.show()

plt.plot(dayindex, hwt[:, outrow, outcol], 'b')
plt.show()
plt.plot(dayindex, pet[:, outrow, outcol], 'r', dayindex, aet[:, outrow, outcol], 'g')
plt.show()

plt.plot(dayindex, radadj[:, outrow, outcol], 'r',
         dayindex, rl[:, outrow, outcol], 'b',
         dayindex, ql[:, outrow, outcol], 'c',
         dayindex, qs[:, outrow, outcol], 'g',
         dayindex, qtotal[:, outrow, outcol], 'y',
         dayindex, radobs[:, outrow, outcol], 'k')
plt.show()

plt.plot(dayindex, qlat_in[:, outletrow, outletcol]-qlat_out[:, outletrow, outletcol], 'g',
         dayindex, ra[:, outletrow, outletcol], 'c',
         dayindex, (ra[:, outletrow, outletcol]+(qlat_in[:, outletrow, outletcol]-qlat_out[:, outletrow, outletcol])), 'b')
plt.show()

flow = ((ra[:, outletrow, outletcol]+(qlat_in[:, outletrow, outletcol]-qlat_out[:, outletrow, outletcol]))
        *geot[1]*np.abs(geot[5]))/(3600*24) #cms
plt.plot(dt, flow, 'b')
plt.show()
