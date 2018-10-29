from precip import *
from et import *
from stats import *
from conversions import *
from groundwater import *
from datetime import datetime
import pandas as pd
from osgeo import gdal
import matplotlib.pyplot as plt
import richdem as rd

def CreateTestDEM(dempath):
    driver = gdal.GetDriverByName("GTiff")
    dem = np.array([[10.0, 9.0, 10.0],
              [9.0, 8.0, 9.0],
              [8.0, 7.0, 8.0]])
    dem = np.array([[12.0, 11.0, 10.0, 11.0, 12.0],
                    [11.0, 10.0, 9.0, 10.0, 11.0],
                    [10.0, 9.0, 8.0, 9.0, 10.0],
                    [9.0, 8.0, 7.0, 8.0, 9.0],
                    [8.0, 7.0, 6.0, 7.0, 8.0]])
    ds = driver.Create(dempath, xsize=5, ysize=5, bands=1, eType=gdal.GDT_Float32)
    geot = [1000.0, 10.0, 0, 1000.0, 0, -10.0]
    ds.SetGeoTransform(geot)
    ds.GetRasterBand(1).WriteArray(dem)
    ds.GetRasterBand(1).FlushCache()
    ds.GetRasterBand(1).SetNoDataValue(-9999.0)
    ds = None
    return dem


nrow = 5
ncol = 5

dirpath = "C:/Users/khafe/Desktop/Classes/WR_502_EnviroHydroModeling/data"
fn = dirpath + "/snotel_klondike_0918.csv"
dempath = dirpath + "/testdem.tif"
dem = CreateTestDEM(dempath)
rddem = rd.LoadGDAL(dempath, no_data=-9999)
rdprop = rd.FlowProportions(dem=rddem, method='D4')
gdal.DEMProcessing(dirpath + "/testslopedeg.tif",
                   dempath,'slope', computeEdges=True, slopeFormat='degree')
gdal.DEMProcessing(dirpath + "/testslopeper.tif",
                   dempath,'slope', computeEdges=True, slopeFormat='percent')
gdal.DEMProcessing(dirpath + "/testaspect.tif", dempath, 'aspect',
                   computeEdges=True, trigonometric=False, zeroForFlat=True)
#read data
str2date = lambda x: datetime.strptime(x.decode("utf-8"), '%m/%d/%Y')
indate = pd.DatetimeIndex(
    np.genfromtxt(fn, delimiter=",", skip_header=1, converters={0: str2date}, usecols=(0)))
indat = np.genfromtxt(fn, delimiter=",", skip_header=1, usecols=(1,2,3,4,5,6))

doy = DayOfYear(indate.month.values, indate.day.values, indate.year.values)

#convert to mm
ndays = math.ceil(5)#indat.shape[0]
swe = indat[1:ndays,0]*25.4
ppt = indat[:ndays-1,5]*25.4
ppt2d = np.reshape(np.repeat(ppt, nrow*ncol), (ndays-1, nrow, ncol))

doy = doy[:ndays-1]
date = indate[:ndays-1]

#convert to degrees C
tmin = np.multiply(np.subtract(indat[:ndays-1,3],32.0), (5.0/9.0))
tmax = np.multiply(np.subtract(indat[:ndays-1,2],32.0), (5.0/9.0))
tavg = np.multiply(np.subtract(indat[:ndays-1,4],32.0), (5.0/9.0))

tmin2d = np.reshape(np.repeat(tmin, nrow*ncol), (ndays-1, nrow, ncol))
tmax2d = np.reshape(np.repeat(tmax, nrow*ncol), (ndays-1, nrow, ncol))
tavg2d = np.reshape(np.repeat(tavg, nrow*ncol), (ndays-1, nrow, ncol))

#from optimization in R
k = 1.1456
tbase = 4.77
train = 2.22
tsnow = 0.76

#calculate melted snow
swe_melt = MeltDegreeDay_USACE(tavg, k, tbase)
swe_melt2d = np.apply_along_axis(MeltDegreeDay_USACE, 0, tavg2d, k=k, tbase=tbase)
ppt_snow, ppt_rain = PrecipPhase(ppt, tavg, train, tsnow)
ppt_snow2d, ppt_rain2d = PrecipPhase_2d(ppt2d, tavg2d, train, tsnow)

swe_mod, act_melt = ModelSWE(ppt_snow, swe_melt)
swe_mod2d, act_melt2d = ModelSWE_2d(ppt_snow2d, swe_melt2d)

ppt_in = np.add(ppt_rain, act_melt)
ppt_in2d = np.add(ppt_rain2d, act_melt2d)

#calculate ET radiation
Ra = ExtraterrestrialRadiation(DegreesToRadians(41.97), doy)
Ra2d = ExtraterrestrialRadiation_2d(np.full((nrow, ncol),DegreesToRadians(41.97)), doy)

#calculate PET
pet = PET_Hargreaves1985(tmax, tmin, tavg, Ra)
pet2d = PET_Hargreaves1985(tmax2d, tmin2d, tavg2d, Ra2d)

s = np.zeros(ppt_in2d.shape)
r = np.zeros(ppt_in2d.shape)
hwt = np.zeros(ppt_in2d.shape)
qlat_out = np.zeros(ppt_in2d.shape)
qlat_in = np.zeros(ppt_in2d.shape)
perc = np.zeros(ppt_in2d.shape)
sb = np.zeros(ppt_in2d.shape)
bf = np.zeros(ppt_in2d.shape)
q = np.zeros(ppt_in2d.shape)
r[0,:,:] = 0 #set initial runoff
s[0,:,:] = 500.0 #set initial storage (i.e water content)
sb[0,:,:] = 0 #set initial aquifer storage (baseflow source)
soildepth = np.full((nrow, ncol), 1000.0) #mm

et1 = np.zeros(pet.shape)
et = np.zeros(ppt_in2d.shape)
ksat = np.full((nrow, ncol), 1000.0) #mm/day
slpds = gdal.Open(dirpath + "/testslopeper.tif")
geot = slpds.GetGeoTransform()
slope = slpds.GetRasterBand(1).ReadAsArray()
print("slope", slope)
por = np.full((nrow, ncol), 0.5)
fc = np.full((nrow, ncol), 0.3)
wp = np.full((nrow, ncol), 0.1)
ksub = np.full((nrow, ncol), 0.1) #mm/day
alpha = 0.02
porl = np.multiply(por, soildepth)
fcl = np.multiply(fc, soildepth)
wpl = np.multiply(wp, soildepth)
smax = porl #mm

et[0,:,:] = ET_theta_2d(pet2d[0,:,:], fcl, wpl, s[0,:,:])

for i in range(1, ppt_in2d.shape[0]):
    et1[i] = ET_theta(pet[i], fcl[1,1], wpl[1,1], s[i-1,1,1])
    et[i,:,:] = ET_theta_2d(pet2d[i,:,:], fcl, wpl, s[i-1,:,:])
    hwt[i,:,:] = WaterTableHeight_2d(por, fc, np.divide(s[i-1,:,:],1000.0), soildepth)
    perc[i,:,:] = Percolation_2d(ksub, hwt[i,:,:])
    qlat_out[i,:,:] = LateralFlow_Darcy_2d(ksat, slope, hwt[i,:,:], geot[1], geot[1])
    qlat_in[i,:,:] = np.multiply(qlat_out[i,:,:], 1.0)
    s[i,:,:] = s[i-1,:,:] - et[i,:,:] + ppt_in2d[0,:,:] - perc[i,:,:] - qlat_out[i,:,:] + qlat_in[i,:,:]
    print(i)

# print("ppt", ppt_in2d)
# print("pet", pet)
# print("pet", pet2d)
# print("et", et1)
# print("et", et)
rd.FillDepressions(rddem, in_place=True)
accum = rd.FlowAccumulation(rddem, method="Dinf")
dinf_fig = rd.rdShow(accum, axes=False, cmap='jet')
accum = rd.FlowAccumulation(rddem, method="D8")
d8_fig = rd.rdShow(accum, axes=False, cmap='jet')
print("hwt", hwt)
print("qlat out", qlat_out)
print("qlat in", qlat_in)
print("s", s)
print("flowprop", rdprop)
print("rddem", rddem)