
def WaterTableHeight(por, fc, theta, depth):
    if (theta < fc):
        hwt = 0
    elif (theta >= por):
        hwt = depth
    else:
        hwt = depth*((theta-fc)/(por-fc))
    return hwt

def LateralFlow_Darcy(ksat, slope, hwt, width=1):
    qlat = -1.0*ksat*slope*hwt*width
    return qlat