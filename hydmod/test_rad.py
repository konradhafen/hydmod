from hydmod.radiation import *
from hydmod.atmosphere import *

slope = 0
aspect = 0
tavg = 7.8
td = -4.0322
doy = 274
lat = 46.75
windv = 1.622
pfc = 80.0
cc = 0.04


ra = DirectSolarRadiation(lat, doy, slope, aspect, units='degrees')
qd = ra
rl = LongwaveRadiation(tavg, pfc, cc, qd, ra)
wr = WindRoughness(windv, pfc)
qs = SensibleRadiation(tavg, wr)
es = SaturationVaporPressure(tavg)
e = ActualVaporPressure(td)

print("es", es, "e", e)
print('shortwave', ra, 'longwave', rl, 'sensible', qs)