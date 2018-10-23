import numpy as np

def Baseflow(alpha, storage, beta=1):
    return alpha*storage**beta

def LateralFlow_Darcy(ksat, slope, hwt, length=1, width=1):
    qlat = ((ksat*slope*hwt*width)/(width*length))
    return qlat

def Percolation(ksub):
    return ksub

def WaterTableHeight(por, fc, theta, depth):
    if (theta < fc):
        hwt = 0
    elif (theta >= por):
        hwt = depth
    else:
        hwt = depth*((theta-fc)/(por-fc))
    return hwt

def WaterTableHeight_2d(por, fc, theta, depth):
    # print("por", por)
    # print("fc", fc)
    # print("theta", theta)
    # print("depth", depth)
    hwt = np.where(theta < fc, 0.0, depth)
    hwt = np.where((theta < por) & (theta >= fc), depth*((theta-fc)/(por-fc)), hwt)
    return hwt