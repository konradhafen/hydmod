from precip import *
from et import *
from stats import *
import matplotlib.pyplot as plt

fn = "C:/Users/khafe/Desktop/Classes/WR_502_EnviroHydroModeling/data/snotel_klondike_0918.csv"
#read data
indat = np.genfromtxt(fn, delimiter=",", skip_header=1)

#convert to cm
swe = indat[1:,1]*2.54
ppt = indat[:-1,6]*2.54

#convert to degrees C
tmin = np.multiply(np.subtract(indat[:-1,4],32.0), (5.0/9.0))
tmax = np.multiply(np.subtract(indat[:-1,3],32.0), (5.0/9.0))
tavg = np.multiply(np.subtract(indat[:-1,5],32.0), (5.0/9.0))

#from optimization in R
k = 1.1456
tbase = 4.77
train = 2.22
tsnow = 0.76

#calculate melted snow
swe_melt = meltDegreeDay_USACE(k, tavg, tbase)
ppt_snow, ppt_rain = precipPhase(ppt, tavg, train, tsnow)
swe_mod, act_melt = modelSWE(ppt_snow, swe_melt)
ppt_in = np.add(ppt_rain, act_melt)

print("NSE", nse(swe_mod, swe))
print("RMSE", rmse(swe_mod, swe))
print("MDSE", md(swe_mod, swe))

print("precip in", np.sum(ppt_in))
print("precip", np.sum(ppt))
print("rain + snow", np.sum(ppt_rain) + np.sum(ppt_snow))

#Soil info and ET info (mm)
et = np.zeros(ppt_in.shape)
pet = 1
fc = 20
wp = 5

#Storage (mm)
s = np.zeros(ppt_in.shape)
s[0] = 50 #set initial storage (i.e water content)
smax = 100

et[0] = modelET(pet, fc, wp, s[0])
s[0] = s[0] + ppt_in[0] - et[0]
for i in range(1, ppt_in.shape[0], 1):
    et[i] = modelET(pet, fc, wp, s[i-1])
    s[i] = s[i-1] + ppt_in[i] - et[i]

#runoff/excess
r = np.where(s>smax, s-smax, 0.0)

index = np.arange(0, swe.shape[0], 1)
plt.plot(index, s, 'r', index, r, 'b')
plt.show()
