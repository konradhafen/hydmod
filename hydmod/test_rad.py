from hydmod.radiation import *

slope = 0.0
aspect = 0.0
tavg = 8.0
doy = 274
lat = 46.75


ra = DirectSolarRadiation(lat, doy, slope, aspect, units='degrees')
qd = ra
rl = LongwaveRadiation(tavg, qd, ra)

print('shortwave', ra, 'longwave', rl)