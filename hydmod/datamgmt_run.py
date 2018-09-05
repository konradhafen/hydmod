from precip import *
from stats import *

fn = "C:/Users/khafe/Desktop/Classes/WR_502_EnviroHydroModeling/data/snotel_klondike_0918.csv"
indat = np.genfromtxt(fn, delimiter=",", skip_header=1)

swe = indat[1:,1]*2.54
ppt = indat[:-1,6]*2.54
tmin = (indat[:-1,4]-32)*(5/9)
tmax = (indat[:-1,3]-32)*(5/9)
tavg = (indat[:-1,2]-32)*(5/9)
