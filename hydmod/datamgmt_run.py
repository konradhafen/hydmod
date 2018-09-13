from precip import *
from et import *
from stats import *
#import richdem

import matplotlib.pyplot as plt

fn = "C:/Users/khafe/Desktop/Classes/WR_502_EnviroHydroModeling/data/snotel_klondike_0918.csv"
#read data
indat = np.genfromtxt(fn, delimiter=",", skip_header=1)

#convert to cm
ndays = indat.shape[0]
swe = indat[1:ndays,1]*2.54
ppt = indat[:ndays-1,6]*2.54

#convert to degrees C
tmin = np.multiply(np.subtract(indat[:ndays-1,4],32.0), (5.0/9.0))
tmax = np.multiply(np.subtract(indat[:ndays-1,3],32.0), (5.0/9.0))
tavg = np.multiply(np.subtract(indat[:ndays-1,5],32.0), (5.0/9.0))

#from optimization in R
k = 1.1456
tbase = 4.77
train = 2.22
tsnow = 0.76

#calculate melted snow
swe_melt = MeltDegreeDay_USACE(k, tavg, tbase)
ppt_snow, ppt_rain = PrecipPhase(ppt, tavg, train, tsnow)
swe_mod, act_melt = ModelSWE(ppt_snow, swe_melt)
ppt_in = np.add(ppt_rain, act_melt)

print("NSE", NSE(swe_mod, swe))
print("RMSE", RMSE(swe_mod, swe))
print("MD", MeanDifference(swe_mod, swe))

print("precip in", np.sum(ppt_in))
print("precip", np.sum(ppt))
print("rain + snow", np.sum(ppt_rain) + np.sum(ppt_snow))

#Storage and runoff
s = np.zeros(ppt_in.shape)
r = np.zeros(ppt_in.shape)
r[0] = 0 #set initial runoff
s[0] = 40 #set initial storage (i.e water content)
soildepth = 80 #cm
ponddepth = 100-soildepth #cm

#Soil info and ET info
et = np.zeros(ppt_in.shape)
pet = 0.4 #cm
fc = 0.4
wp = 0.1
fcl = 0.4*soildepth
wpl = 0.1*soildepth
smax = fcl + ponddepth #cm
print("smax", smax)

et[0] = ModelET(pet, fc, wp, s[0])
s[0] = s[0] + ppt_in[0] - et[0]
for i in range(1, ppt_in.shape[0], 1):
    et[i]= ModelET(pet, fcl, wpl, s[i-1])
    s[i] = s[i-1] + ppt_in[i] - et[i]
    if s[i] > smax:
        r[i] = s[i] - smax
        s[i] = smax

#runoff/excess

#pond depth
p = np.where(s>fcl, s-fcl, 0.0)
swc = np.subtract(s, p)

index = np.arange(0, swe.shape[0], 1)
plt.subplot(4,1,1)
plt.plot(index, s, 'r', index, swc, 'k', index, p, 'b')
plt.subplot(4,1,2)
plt.plot(index, ppt_in, 'b',index, r, 'c')
plt.subplot(4,1,3)
plt.plot(index,et,'g')
plt.subplot(4,1,4)
plt.plot(index,swe,'b',index,swe_mod,'r')
plt.show()
