from precip import *
from et import *
from stats import *
from conversions import *
from groundwater import *
from datetime import datetime
import pandas as pd
#import richdem

import matplotlib.pyplot as plt

fn = "C:/Users/konrad/Desktop/Classes/WR_502_EnviroHydroModeling/data/snotel_klondike_0918.csv"
#read data
str2date = lambda x: datetime.strptime(x.decode("utf-8"), '%m/%d/%Y')
indate = pd.DatetimeIndex(
    np.genfromtxt(fn, delimiter=",", skip_header=1, converters={0: str2date}, usecols=(0)))
indat = np.genfromtxt(fn, delimiter=",", skip_header=1, usecols=(1,2,3,4,5,6))

doy = DayOfYear(indate.month.values, indate.day.values, indate.year.values)

#convert to mm
ndays = math.ceil(365*8.5)#indat.shape[0]
swe = indat[1:ndays,0]*25.4
ppt = indat[:ndays-1,5]*25.4

#convert to degrees C
tmin = np.multiply(np.subtract(indat[:ndays-1,3],32.0), (5.0/9.0))
tmax = np.multiply(np.subtract(indat[:ndays-1,2],32.0), (5.0/9.0))
tavg = np.multiply(np.subtract(indat[:ndays-1,4],32.0), (5.0/9.0))

doy = doy[:ndays-1]
date = indate[:ndays-1]

def RunModel(swe, ppt, tmin, tmax, tavg, doy):

    #from optimization in R
    k = 1.1456
    tbase = 4.77
    train = 2.22
    tsnow = 0.76

    #calculate melted snow
    swe_melt = MeltDegreeDay_USACE(tavg, k, tbase)
    ppt_snow, ppt_rain = PrecipPhase(ppt, tavg, train, tsnow)
    swe_mod, act_melt = ModelSWE(ppt_snow, swe_melt)
    ppt_in = np.add(ppt_rain, act_melt)

    #calculate ET radiation
    Ra = ExtraterrestrialRadiation(DegreesToRadians(41.97), doy)
    #calculate PET
    pet = PET_Hargreaves1985(tmax, tmin, tavg, Ra)

    print("precip in", np.sum(ppt_in))
    print("precip", np.sum(ppt))
    print("rain + snow", np.sum(ppt_rain) + np.sum(ppt_snow))

    #Storage and runoff
    s = np.zeros(ppt_in.shape)
    r = np.zeros(ppt_in.shape)
    hwt = np.zeros(ppt_in.shape)
    qlat = np.zeros(ppt_in.shape)
    qlatin = np.zeros(ppt_in.shape)
    perc = np.zeros(ppt_in.shape)
    sb = np.zeros(ppt_in.shape)
    bf = np.zeros(ppt_in.shape)
    q = np.zeros(ppt_in.shape)
    r[0] = 0 #set initial runoff
    s[0] = 0 #set initial storage (i.e water content)
    sb[0] = 0 #set initial aquifer storage (baseflow source)
    soildepth = 1000 #mm
    ponddepth = 1000-soildepth #mm

    #Soil info and ET info
    et = np.zeros(ppt_in.shape)
    ksat = 1000.0 #mm/day
    slope = 0.1
    por = 0.5
    fc = 0.3
    wp = 0.1
    ksub = 1.0 #mm/day
    alpha = 0.02
    porl = por*soildepth
    fcl = fc*soildepth
    wpl = wp*soildepth
    smax = porl + ponddepth #mm
    print("smax", smax)

    et[0] = ET_theta(pet[0], fc, wp, s[0])
    s[0] = s[0] + ppt_in[0] - et[0]
    for i in range(1, ppt_in.shape[0], 1):
        et[i]= ET_theta(pet[i], fcl, wpl, s[i-1])

        hwt[i] = WaterTableHeight(por, fc, (s[i-1]/1000.0), soildepth)
        if hwt[i] > 0:
            perc[i] = Percolation(ksub)
        bf[i] = Baseflow(alpha, sb[i-1])
        sb[i] = sb[i-1] + perc[i] - bf[i]
        qlat[i] = LateralFlow_Darcy(ksat, slope, hwt[i], length=10000.0)
        qlatin[i] = 1.0*qlat[i]
        s[i] = s[i-1] + ppt_in[i] - et[i] - qlat[i] + qlatin[i] - perc[i]
        if s[i] > smax:
            r[i] = s[i] - smax
            s[i] = smax
        q[i] = r[i] + qlat[i] + bf[i]

    #pond depth
    p = np.where(s>porl, s-porl, 0.0)
    swc = np.subtract(s, p)

    index = np.arange(0, swe.shape[0], 1)
    plt.subplot(4,1,1)
    plt.plot(index, s, 'r', index, swc, 'k', index, p, 'b')
    plt.subplot(4,1,2)
    plt.plot(index, pet, 'r', index, et, 'b')
    plt.subplot(4, 1, 3)
    plt.plot(index, q, 'b', index, qlat, 'g', index, r, 'c')
    plt.subplot(4,1,4)
    plt.plot(date, bf, 'r', date, sb, 'k')
    # plt.subplot(3,1,3)
    # plt.plot(index,et,'g')
    # plt.subplot(4,1,4)
    # plt.plot(index,swe,'b',index,swe_mod,'r')
    plt.show()

    plt.plot(date, bf, 'r', date, qlat, 'g', date, r, 'c', date, q, 'b')
    plt.show()
    plt.plot(date, ((q/1000.0)*10.0*10.0)/(60.0*60.0*24.0), 'b')
    plt.show()
    print("precip", np.sum(ppt_in))
    print("et", np.sum(et))
    print("r", np.sum(r))
    print("storage", s[-1])
    print("balance", s[-1]+np.sum(et)+np.sum(r)-s[0]+np.sum(np.subtract(qlat, qlatin))+np.sum(perc))

RunModel(swe, ppt, tmin, tmax, tavg, doy)