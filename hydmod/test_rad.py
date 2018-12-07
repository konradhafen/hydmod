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
tsnow = 0.0
dsls = 1
swe = 1


sa = SnowAlbedo(dsls, swe)

qo = ExtraterrestrialRadiation(lat, doy, units='degrees') #solar radiation incident on a flat surface
qcs = ClearSkyRadiation(qo) #kJ/m^2
cc = CloudCoverFraction(186.1 / 1000.0 * 3600.0 * 24.0, qcs)
print('qcs', qo, 'cc', cc, lat, doy)
ra = DirectSolarRadiation_SlopingSurface(lat, doy, qcs, slope, aspect, units='degrees')
ra = DirectSolarRadiation_Adjustment(ra, pfc, sa, cc)
# ra = DirectSolarRadiation_SlopingSurface(lat, doy, 186.1 / 1000.0 * 3600.0 * 24.0, slope, aspect, pfc,
#                                          swe=1, snowage=1, units='degrees')
qd = ra
rl = LongwaveRadiation(tavg, pfc, cc)
wr = WindRoughness(windv, pfc)
qs = SensibleRadiation(tavg, wr)
es = VaporPressure(tavg)
e = VaporPressure(td)
vda = VaporDensity(e, td)
vds = VaporDensity(es, tavg)
vdsnow = VaporDensity(VaporPressure(0.0), 0.0)
ql = LatentRadiation(td, tsnow, windv, pfc)
# radj = (1-(pfc/100.0))*(1.0-sa)*ra#*1000.0
rtotal = ra + rl + qs + ql + HEAT_FROM_GROUND
melt = rtotal/(LATENT_HEAT_FUSION*DENSITY_WATER)

print('snow albedo', sa)
print("es", es, "e", e, 'vds', vds, 'vda', vda, 'vdsnow', vdsnow)
print('shortwave', ra, 'longwave', rl, 'sensible', qs, 'latent', ql, 'total', rtotal)
print('melt (m)', melt)