from hydmod.radiation import *
from hydmod.atmosphere import *

slope = 100
aspect = 180
tavg = 7.8
td = -4.0322
doy = 274
lat = 46.75
windv = 1.622
pfc = 80.0
cc = 0.04
tsnow = 0.0


ra = DirectSolarRadiation(lat, doy, slope, aspect, units='degrees')
qd = ra
rl = LongwaveRadiation(tavg, pfc, cc, qd, ra)
wr = WindRoughness(windv, pfc)
qs = SensibleRadiation(tavg, wr)
es = VaporPressure(tavg)
e = VaporPressure(td)
vda = VaporDensity(e, td)
vds = VaporDensity(es, tavg)
vdsnow = VaporDensity(VaporPressure(0.0), 0.0)
ql = LatentRadiation(td, tsnow, windv, pfc)

print("es", es, "e", e, 'vds', vds, 'vda', vda, 'vdsnow', vdsnow)
print('shortwave', ra, 'longwave', rl, 'sensible', qs, 'latent', ql)