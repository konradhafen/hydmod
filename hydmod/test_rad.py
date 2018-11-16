from hydmod.radiation import *

slope = 100.0
aspect = 180.0
tavg = 7.8
doy = 274
lat = 46.75


ra = DirectSolarRadiation(lat, doy, slope, aspect, units='degrees')
qd = ra
rl = LongwaveRadiation(tavg, qd, ra)

print('shortwave', ra, 'longwave', rl)