
def Baseflow(alpha, storage):
    return alpha*storage

def LateralFlow_Darcy(ksat, slope, hwt, length=1, width=1):
    qlat = ((ksat*slope*hwt*width)/length)
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