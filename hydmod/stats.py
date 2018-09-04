import numpy as np

def md(modeled, observed):
    """
    Mean difference between modeled and observed values

    Args:
        modeled: Modeled values (numpy array 1D)
        observed: Observed values (numpy array 1D)

    Returns:
        Mean Difference value

    """
    return(np.mean(np.subtract(modeled, observed)))

def nse(modeled, observed):
    """
    Nash-Sutcliffe Efficiency

    Args:
        modeled: Modeled values (numpy array 1D)
        observed: Observed values (numpy array 1D)

    Returns:
        NSE value

    """

    return(1-(np.sum(np.square(np.subtract(observed, modeled)))/np.sum(np.square(np.subtract(observed, np.sum(observed))))))

def rmse(modeled, observed):
    """
    Root mean square error

    Args:
        modeled: Modeled/predicted values (numpy array)
        observed: Observed values (numpy array)

    Returns:

    """
    return(np.nanmean(np.sqrt(np.square(np.subtract(modeled - observed^2)))))