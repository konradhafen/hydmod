import numpy as np

PI = 3.14159

def ConvertToDegrees(var, units):
    if units == "radians":
        return RadiansToDegrees(var)
    else:
        return var

def ConvertToRadians(var, units):
    if units == "degrees":
        return DegreesToRadians(var)
    else:
        return var

def DayOfYear(m, d, y):
    """
    Day of year (e.g. 1-365 or 366)

    Args:
        m: Month (1-12)
        d: Day (1-31)
        y: Year (e.g. 2012)

    Returns:
        Day of year (1-365 or 366)

    """
    #j2 = d-32.0+math.floor(275.0*(m/9.0))+2.0*math.floor(3.0/(m+1.0))+math.floor((m/100.0)-((y%4.0)/4.0)+0.975)

    pt1 = np.subtract(d, 32.0)
    pt2 = np.trunc(np.multiply(275.0, np.divide(m, 9.0)))
    pt3 = np.multiply(2.0, np.trunc(np.divide(3.0, np.add(m, 1.0))))
    pt4 = np.trunc(np.add(np.subtract(np.divide(m, 100.0), np.divide(np.remainder(y, 4.0), 4.0)), 0.975))

    j = np.add(np.add(pt1, pt2), np.add(pt3, pt4))

    return j.astype(int)

def DegreesToRadians(deg):
    """
    Convert degrees to radians

    Args:
        deg: Value in degrees

    Returns:
        Value in radians

    """
    rad = deg * (PI/180.0)
    return rad

def RadiansToDegrees(rad):
    """
    Convert radians to degrees

    Args:
        rad: Value in radians

    Returns:
        Value in degrees

    """
    deg = rad * (180.0/PI)