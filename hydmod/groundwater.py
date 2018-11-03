import numpy as np

def Baseflow(alpha, storage, beta=1):
    return alpha*storage**beta

def LateralFlow_Darcy_2d(ksat, slope, hwt, length=1, width=1, convf=1.0):
    """
    2-dimensional darcy flow

    Args:
        ksat: saturated hydraulic conductivity
        slope: slope of land surface or water table
        hwt: height of the water table
        length: distance from one cell to another
        width: cell width
        convf: unit conversion factor to convert width and length to units of hwt and ksat (default: 1.0)

    Returns:

    """
    qlat = np.divide(np.multiply(ksat, np.multiply(slope, np.multiply(hwt, np.multiply(width, convf)))),
                     np.multiply(np.multiply(width, convf), np.multiply(length, convf)))
    return qlat

def LateralFlow_Darcy(ksat, slope, hwt, length=1, width=1):
    qlat = ((ksat*slope*hwt*width)/(width*length))
    return qlat

def Percolation(ksub):
    return ksub

def Percolation_2d(ksub, hwt):
    perc = np.where(hwt > 0.0, ksub, 0.0)
    return perc

# def WaterTableHeight(por, fc, theta, depth):
#     if (theta < fc):
#         hwt = 0
#     elif (theta >= por):
#         hwt = depth
#     else:
#         hwt = depth*((theta-fc)/(por-fc))
#     return hwt

def WaterTableHeight(por, fc, theta, depth):
    """

    Args:
        por: soil porosity
        fc: soil field capacity
        theta: soil moisture content
        depth: depth of soil profile

    Returns:

    """
    # where water content is less than field capacity no water table, set water table height to depth of the
    # soil profile everywhere else
    hwt = np.where(theta < fc, 0.0, depth)
    # where water content is between field capacity and porosity, water table height is determined linearly
    # the previous line accounts for anywhere water content is >= porosity
    hwt = np.where((theta < por) & (theta >= fc), depth*((theta-fc)/(por-fc)), hwt)
    return hwt