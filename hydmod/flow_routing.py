import numpy as np
import richdem as rd
from osgeo import gdal

def FlowProportions(dem, nodata=-9999):
    #calculate proportion of flow from a cell to neighbors
    eprop = np.zeros((8, dem.shape[0], dem.shape[1]))
    eprop[0, :, 0:-1] = (dem[:, 0:-1] - dem[:, 1:]) # to east
    eprop[1, 0:-1, 0:-1] = (dem[0:-1, 0:-1] - dem[1:, 1:])/np.sqrt(2.0)# to southeast
    eprop[2, 0:-1, :] = (dem[0:-1, :] - dem[1:, :]) # to south
    eprop[3, 0:-1, 1:] = (dem[0:-1, 1:] - dem[1:, 0:-1]) / np.sqrt(2.0) # to southwest
    eprop[4, :, 1:] = (dem[:, 1:] - dem[:, 0:-1]) # to west
    eprop[5, 1:, 1:] = (dem[1:, 1:] - dem[0:-1, 0:-1]) / np.sqrt(2.0) # to northwest
    eprop[6, 1:, :] = (dem[1:, :] - dem[0:-1, :]) # to north
    eprop[7, 1:, 0:-1] = (dem[1:,0:-1] - dem[0:-1,1:])/np.sqrt(2.0) # to northeast

    eprop[eprop < 0.0] = 0.0
    esum = np.sum(eprop, axis=0)
    fprop = np.where(esum > 0.0, np.divide(eprop, esum), 0.0)

    return fprop

def RouteFlow(flowprop, flow):
    rflow = np.zeros((flowprop.shape[1], flowprop.shape[2]))
    #from west
    flowdec = flow.copy()
    rflow[:, 1:] = rflow[:, 1:] + flowprop[0, :, 0:-1] * flow[:, 0:-1]  # flow from west
    flowdec[:, 0:-1] = flowdec[:, 0:-1] - flowprop[0, :, 0:-1] * flow[:, 0:-1]
    rflow[1:, 1:] = rflow[1:, 1:] + flowprop[1, 0:-1, 0:-1] * flow[0:-1, 0:-1]  # flow from northwest
    flowdec[0:-1, 0:-1] = flowdec[:-1, 0:-1] - flowprop[1, 0:-1, 0:-1] * flow[0:-1, 0:-1]
    rflow[1:, :] = rflow[1:, :] + flowprop[2, 0:-1, :] * flow[0:-1, :]  # flow from north
    flowdec[0:-1, :] = flowdec[0:-1, :] - flowprop[2, 0:-1, :] * flow[0:-1, :]
    rflow[1:, 0:-1] = rflow[1:, 0:-1] + flowprop[3, 0:-1, 1:] * flow[0:-1, 1:]  # flow from northeast
    flowdec[0:-1, 1:] = flowdec[0:-1, 1:] - flowprop[3, 0:-1, 1:] * flow[0:-1, 1:]
    rflow[:, 0:-1] = rflow[:, 0:-1] + flowprop[4, :, 1:] * flow[:, 1:]  # flow from east
    flowdec[:, 1:] = flowdec[:, 1:] - flowprop[4, :, 1:] * flow[:, 1:]
    rflow[0:-1, 0:-1] = rflow[0:-1, 0:-1] + flowprop[5, 1:, 1:] * flow[1:, 1:]  # flow from southeast
    flowdec[1:, 1:] = flowdec[1:, 1:] - flowprop[5, 1:, 1:] * flow[1:, 1:]
    rflow[0:-1, :] = rflow[0:-1, :] + flowprop[6, 1:, :] * flow[1:, :]  # flow from south
    flowdec[1:, :] = flowdec[1:, :] - flowprop[6, 1:, :] * flow[1:, :]
    rflow[0:-1, 1:] = rflow[0:-1, 1:] + flowprop[7, 1:, 0:-1] * flow[1:, 0:-1]  # flow from southwest
    flowdec[1:, 0:-1] = flowdec[1:, 0:-1] - flowprop[7, 1:, 0:-1] * flow[1:, 0:-1]
    print("flow not routed", np.sum(flowdec))

    return rflow

# dirpath = "C:/Users/khafe/Desktop/Classes/WR_502_EnviroHydroModeling/data"
# dempath = dirpath + "/testdem.tif"
# demds = gdal.Open(dempath)
# demnp = demds.GetRasterBand(1).ReadAsArray()
#
# print(demnp)
# print("fprops", FlowProportions(demnp))
# print("routed", RouteFlow(FlowProportions(demnp), np.ones((5,5))))
# print("sum routed", np.sum(RouteFlow(FlowProportions(demnp), np.ones((5,5)))))


