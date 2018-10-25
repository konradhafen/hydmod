import numpy as np

def Baseflow(alpha, storage, beta=1):
    return alpha*storage**beta

def LateralFlow_Darcy_2d(ksat, slope, hwt, length=1, width=1):
    qlat = np.divide(np.multiply(ksat, np.multiply(np.divide(slope, 100.0), np.multiply(hwt, width))),
                     np.multiply(width, length))
    return qlat

def LateralFlow_Darcy(ksat, slope, hwt, length=1, width=1):
    qlat = ((ksat*slope*hwt*width)/(width*length))
    return qlat

def Percolation(ksub):
    return ksub

def Percolation_2d(ksub, hwt):
    perc = np.where(hwt > 0.0, ksub, 0.0)
    return perc

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