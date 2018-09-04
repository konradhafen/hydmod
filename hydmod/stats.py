import numpy as np

def rmse(modeled, observed):
    """
    Root mean square error

    Args:
        modeled: Modeled/predicted values (numpy array)
        observed: Observed values (numpy array)

    Returns:

    """
    return(np.nanmean(np.sqrt(modeled - observed)^2))