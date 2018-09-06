import numpy as np

def modelET(pet, fc, wp, wc):
    """
    Model evapotranspiration, multiply potential ET by factor based on soil water content
    Args:
        pet: Potential evapotranspiration
        fc: Field capacity
        wp: Wilting point
        wc: Water content

    Returns:
        ET estimate

    """
    theta = 1.0
    if wc < 0.8*fc and wc > wp:
        theta = (0.8*fc - wc)/(0.8*fc-wp)
    elif wc <= wp:
        theta = 0.0
    return(pet*theta)