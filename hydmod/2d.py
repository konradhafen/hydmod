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

def CreateTestDEM():
    dem = np.array([[10.0, 9.0, 10.0],
              [9.0, 8.0, 9.0],
              [8.0, 7.0, 8.0]])
    return dem

dem = CreateTestDEM()
nrow = 3
ncol = 3

fn = "C:/Users/khafe/Desktop/Classes/WR_502_EnviroHydroModeling/data/snotel_klondike_0918.csv"
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

