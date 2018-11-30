from hydmod.radiation import *

slope = 0
aspect = 0
tavg = 7.8
doy = 274
lat = 46.75
windv = 1.622
pfc = 80.0


ra = DirectSolarRadiation(lat, doy, slope, aspect, units='degrees')
qd = ra
rl = LongwaveRadiation(tavg, qd, ra)
wr = WindRoughness(windv, pfc)
qs = SensibleRadiation(tavg, wr)

print('shortwave', ra, 'longwave', rl, 'sensible', qs)