PI = 3.14159

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