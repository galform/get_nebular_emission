#! /usr/bin/env python

"""
Some useful functions
  percentiles(xbins,xarray,yarray,per): obtains percentiles of yarray in xbins
  convert_to_stdev(grid): normalised a grid to cumulative standard deviations.
  n_gt_x(x,array): returns the number of elements in the array larger than each of the values in x.
  chi2(obs,model,err): returns the chi^2 for a model
NOTE: this module requires the numpy and scipy libraries to be
      available for import!
"""
import sys
import numpy as np


def percentiles(val, data, weights=None):
    """
    Examples
    --------
    >>> import numpy as np
    >>> import stats
    >>> data = np.array(np.arange(0.,100.,10.))
    >>> stats.percentiles(0.5,data)
    >>> 45.0
    """

    if (val < 0 or val > 1):
        sys.exit('STOP percentiles: 0<val<1')

    if (weights is None):
        ws = np.zeros(shape=(len(data)));
        ws.fill(1.)
    else:
        ws = weights

    data = np.array(data);
    ws = np.array(ws)
    ind_sorted = np.argsort(data)  # Median calculation from wquantiles
    sorted_data = data[ind_sorted];
    sorted_weights = ws[ind_sorted]

    num = np.cumsum(sorted_weights) - 0.5 * sorted_weights
    den = np.sum(sorted_weights)
    if (den != 0):
        pn = num / den
        percentiles = np.interp(val, pn, sorted_data)
    else:
        sys.exit('STOP percentiles: problem with weights')
    return percentiles


def perc_2arrays(xedges, xarray, yarray, val, weights=None, nmin=None):
    """
    Returns percentiles of yarray over xbins
    Parameters
    ----------
    xedges : array of floats
        Bin edges on the x-axis
    xarray : array of floats
        Values for the x-axis
    yarray : array of floats
        Values for the y-axis
    val : float from 0 to 1
        Value used to determine the percentile to be calculated
    weights : array of floats
        Weights for the yarray values
    nmin : integer
        Minimal number of points to be considered in a bin
    Returns
    -------
    apercentile : string of floats
       Percentiles of the yarray within xbins
    Examples
    --------
    >>> import numpy as np
    >>> import stats
    >>> xedges = np.array([0.,1.,2.])
    >>> xarray = np.array(np.arange(0.,2.,0.1))
    >>> yarray = np.append(np.array(np.arange(1.,11.,1.)),np.array(np.arange(1.,11.,1.)))
    >>> stats.perc_2arrays(xedges,xarray,yarray,0.5)
    >>> array([5.5, 5.5])
    """
    xlen = len(xedges) - 1
    apercentile = np.zeros(shape=(xlen));
    apercentile.fill(-999.)

    if (len(xarray) != len(yarray)):
        sys.exit('ERROR @ perc_2arrays: The lenght of the input arrays should be equal.')

    if (nmin is None):
        nmin = 1

    for i in range(xlen):
        ind = np.where((xarray >= xedges[i]) & (xarray < xedges[i + 1]))
        # We require at least nmin points per bin
        if (np.shape(ind)[1] >= nmin):
            data = yarray[ind]

            if (weights is None):
                apercentile[i] = percentiles(val, data)
            else:
                if (len(weights) != len(yarray)):
                    sys.exit(
                        'ERROR @ perc_2arrays: The lenght of the weights array should be equal to the input array.')

                ws = weights[ind]
                apercentile[i] = percentiles(val, data, weights=ws)

    return apercentile


def av_2arrays(xbins, xarray, yarray, weights, nmin):
    """ Returns average of yarray over xbins"""
    xlen = len(xbins) - 1
    av_2arrays = np.zeros(shape=(xlen));
    av_2arrays.fill(-999.)

    if len(xarray) != len(yarray):
        sys.exit('ERROR @ perc_2arrays: The lenght of the input arrays should be equal.')

    for i in range(xlen):
        ind = np.where((xarray >= xbins[i]) & (xarray < xbins[i + 1]))
        # We require at least nmin points per bin
        num = np.shape(ind)[1]
        if (num > nmin):
            data = yarray[ind];
            ws = weights[ind]
            dw = ws * data
            av_2arrays[i] = np.sum(dw) / np.sum(ws)

    return av_2arrays


def convert_to_stdev(grid):
    """
    From Coleman Krawczyk
    Based on https://pypi.python.org/simple/astroml/
    Given a grid of values, convert them to cumulative standard deviation.
    This is useful for drawing contours with standard deviations as the levels.
    """
    shape = grid.shape
    grid = grid.ravel()
    # Obtain indices to sort and unsort the flattened array
    i_sort = np.argsort(grid)[::-1]
    i_unsort = np.argsort(i_sort)
    grid_cumsum = grid[i_sort].cumsum()
    grid_cumsum /= grid_cumsum[-1]

    return grid_cumsum[i_unsort].reshape(shape)


def n_gt_x(xedges, array):
    y = np.zeros(len(xedges))

    for i, xedge in enumerate(xedges):
        ind = np.where(array > xedge)
        y[i] = np.shape(ind)[1]

    return y


def get_interval(val, low, high):
    '''
    Get the index, i, of the interval where low[i]<val<high[i]

    Parameters:
    val : int or float
        Value
    low : array of int or floats
        Array of the low edges of the intervals
    high : array of int or floats
        Array of the high edges of the intervals
    Returns:
    ind : integer
       Index of the interval the value belongs to or -999 if outside
    '''

    ind = -999

    if (len(low) != len(high)):
        return ind

    if (val == high[-1] and val >= low[-1]):
        ind = len(high) - 1
    else:
        linds = np.where(val >= low)
        hinds = np.where(val < high)

        if (np.shape(linds)[1] > 0 and np.shape(hinds)[1] > 0):
            lind = linds[0]
            hind = hinds[0]
            common = list(set(lind).intersection(hind))
            if (len(common) == 1):
                ind = common[0]

    return ind


def chi2(obs, model, err2):
    '''
    Get the chi^2 for a given model

    Parameters:
    obs : array of floats
        The observatioins or target values
    model : array of floats
        The model values (should be the same length as obs)
    err2 : array of floats
        The error**2 of the observations
    Returns:
    val : float
       chi^2 values
    '''
    val = 0.
    for i, iobs in enumerate(obs):
        val = val + (iobs - model[i]) ** 2 / err2[i]
    return val


def get_err2Pk(k, Pk, dk, N, vol):
    '''
    Get the error of the Power Spectrum
    Parameters:
    -----------
    k : numpy array of floats
       Wavenumber of the modes
    Pk : numpy array of floats
       Power spectrum at each k
    dk : float
      Size of the step for k
    N : float
      Number of elements
    vol : float
      Considered volume
    Returns:
    --------
    err2Pk : float
      Square of the power spectrum error
    '''

    err2Pk = None

    if (len(k) != len(Pk)):
        print('STOP (stats.get_err2Pk): k and Pk are different lengths')
        return err2Pk

    norm = (2 * np.pi) ** 2 / (k * k * dk * vol)
    err2Pk = norm * (Pk + vol / N) ** 2

    return err2Pk