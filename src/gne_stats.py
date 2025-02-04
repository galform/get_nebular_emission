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
import src.gne_const as c


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

            if len(data)<5:
                apercentile[i] = -999.
            elif (weights is None):
                apercentile[i] = percentiles(val, data)
                # print(val, apercentile[i])
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


def locate_interval(val, edges):
    '''
    Get the index, i, of the interval, [), to which val belongs.
    If outside the limits, using values -1 or the number of bins+1.

    Parameters
    ----------
    val : int or float or array of ints or floats
        Value to evaluate
    edges : array of int or floats
        Array of the n edges for the (n-1) intervals
        
    Returns
    -------
    jl : integer
        Index of the interval, [edges(jl),edges(jl+1)), where val is place.
    '''

    n = edges.size
    jl = np.searchsorted(edges, val, side='right') - 1
    jl = np.clip(jl, -1, n - 1)
    return jl


def interpl_weights(xx,edges):
    '''
    Get linear interpolation weights: xd=(x-x1)/(x2-x1)
    Values outside the edges limits are given the weights
    corresponding to the minimum and maximum edge values.
    
    Parameters
    ----------
    xx : float (or int) or array of floats (or int)
        Values to be evaluated
    edges : array of floats (or int)
        Array of the n edges for the (n-1) intervals
        
    Returns
    -------
    xd : float or list of float (or int)
        Weights for linear interpolation
    ix : int or list of ints
        Lower index of the interval the value belongs to
    '''
    # Size of the 1D grid
    n = edges.size

    # If scalar, turn it into array
    scalar = False
    if isinstance(xx, (float, int)): # Floats
        scalar = True
        xx = np.array([xx])

    xd = []; ix = []
    for x in xx[:]:
        jl = locate_interval(x,edges)
        if jl<0: # Use first value in the grid
            xd.append(0.0)
            jl = 0
        elif jl > n - 2: # Use last value in the grid
            xd.append(1.0)
            jl = n - 2
        else:
            d = (x - edges[jl]) / (edges[jl + 1] - edges[jl])
            xd.append(d)
        ix.append(jl)

    outxd = np.asarray(xd)
    outix = np.asarray(ix,dtype=int)
    if scalar:
        outxd = outxd[0]
        outix = outix[0]
    return outxd, outix


def bilinear_interpl_matrix(vals,grid):
    '''
    Get linear interpolation eights: xd=(x-x1)/(x2-x1)
    Values outside the edges limits are given the weights
    corresponding to the minimum and maximum edge values.
    
    Parameters
    ----------
    xx : float (or int) or array of floats (or int)
        Values to be evaluated
    edges : array of floats (or int)
        Array of the n edges for the (n-1) intervals
        
    Returns
    -------
    xd : float or list of float (or int)
        Weights for linear interpolation
    ix : int or list of ints
        Lower index of the interval the value belongs to
    '''
    n = edges.size
    
    if isinstance(xx, (float, int)): # Floats
        x = xx
        jl = locate_interval(x,edges)
        if jl<0: # Use first value in the grid
            xd = 0.0
            jl = 0
        elif jl > n - 2: # Use last value in the grid
            xd = 1.0
            jl = n - 2
        else:
            xd = (x - edges[jl]) / (edges[jl + 1] - edges[jl])
        ix = jl
        
    else: # Arrays
        xd = []; ix = []
        for x in xx[:]:
            jl = locate_interval(x,edges)
            if jl<0: # Use first value in the grid
                xd.append(0.0)
                jl = 0
            elif jl > n - 2: # Use last value in the grid
                xd.append(1.0)
                jl = n - 2
            else:
                d = (x - edges[jl]) / (edges[jl + 1] - edges[jl])
                xd.append(d)
            ix.append(jl)

    return np.asarray(xd), np.asarray(ix,dtype=int)


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
