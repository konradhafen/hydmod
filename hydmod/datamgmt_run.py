from precip import *
from stats import *

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

print(ppt[4:10])
print(ppt_rain[4:10])
print(ppt_snow[4:10])
print(swe_melt[4:10])
