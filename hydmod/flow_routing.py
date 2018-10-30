import numpy as np
import richdem as rd

def FlowProportions(dem, nodata=-9999):
    #calculate proportion of flow from a cell to neighbors
    eprop = np.zeros((8, dem.shape[0], dem.shape[1]))
    eprop[0, :, 0:-1] = (dem[:, 0:-1] - dem[:, 1:]) #east
    eprop[1, 0:-1, 0:-1] = (dem[0:-1, 0:-1] - dem[1:, 1:])/np.sqrt(2.0)#southeast
    eprop[2, 0:-1, :] = (dem[0:-1, :] - dem[1:, :]) #south
    eprop[3, 0:-1, 1:] = (dem[0:-1, 1:] - dem[1:, 0:-1]) / np.sqrt(2.0) #southwest
    eprop[4, :, 1:] = (dem[:, 1:] - dem[:, 0:-1]) #west
    eprop[5, 1:, 1:] = (dem[1:, 1:] - dem[0:-1, 0:-1]) / np.sqrt(2.0) #northwest
    eprop[6, 1:, :] = (dem[1:, :] - dem[0:-1, :]) #north
    eprop[7, 1:, 0:-1] = (dem[1:,0:-1] - dem[0:-1,1:])/np.sqrt(2.0) #northeast


