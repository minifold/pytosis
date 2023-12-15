"""Various RFI reduction techniques that can be applied as a quick
"""
import os
from itertools import islice
import numpy as np

from astropy.io import fits
from astropy.convolution import convolve, Gaussian1DKernel
from astropy.stats import SigmaClip, sigma_clipped_stats

from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

from scipy.ndimage import median_filter


def gauss_conv(data, std=2, showplot=False):
    """
    First take the median of the data along the time axis, then run a smooth
    convolution along the full spectrum

    Parameters
    ----------
    data : 1D array-like
        Data to be convolved. Will be the median of the 2D data.
    std : int, optional
        Standard deviation of the gaussian kernel. Default is 2.
    showplot : boolean, optional
        Determines whether the plot of the graph for gaussian convolution is shown.
    """

    gauss_kernel = Gaussian1DKernel(stddev=std)

    conv = np.zeros(data.shape)
    conv = convolve(data, gauss_kernel)

    if showplot is True:
        plt.plot(data)
        plt.plot(conv)
        plt.show()

        # Plot the data divided by the convolved data.
        plt.plot(data / conv)
        plt.show()

    return conv


def robust_avg(data, sigma: int = 3, maxiters: int = 15, showplot: bool = False):
    """
    Takes an average using the astropy SigmaClip function.

    Parameters
    ----------
    data : 1D array-like
        Array of the median
    sigma : int, optional
        Number of standard deviations away from the median the algorithm is
        allowed to stray. Keep this number small.
    maxiters : int, optional
        Number of iterations the averaging algorithm will perform.
    showplot : boolean, optional
        Determines whether the plot of the graph for gaussian convolution is shown.
    """

    sigclip = SigmaClip(sigma=sigma, maxiters=maxiters)
    filtered_data = sigclip(data)

    if showplot is True:
        plt.plot(data)
        plt.plot(filtered_data)
        plt.show()

    return filtered_data


def med_filter(data, boxes: int = 6, showplot: bool = False):
    """
    Conducts a running median filter across the data.

    Essentially a convolution function, but where instead of a gaussian filter,
    the median of each set of points is taken and applied as a kernel.

    Parameters
    ----------
    data : array-like
       The median data of the 2D array passed in from `read_fits()`. For
       practicality, I also check if the array is 2D then take the median for
       you, since the idea of taking the median and then taking the median
       again through this filter is kind of confusing.
    boxes : int, optional
       Number of boxes that are taken to perform one "averaging" calculation
       over. Default is 6.
    showplot : boolean, optional
        Determines whether the plot of the graph for gaussian convolution is shown.

    Returns
    -------
    med : 1D array
    """

    if len(data.shape) == 2:
        data = np.median(data, axis=0)

    med = median_filter(data, 6)

    if showplot is True:
        plt.plot(med)
        plt.show()

    return med


def plot_averages(meandata, line=0):
    """
    """
    (rows, cols) = meandata.shape
    pols = cols // 2


def plot_polarizations(data, time=0, y_lim=[1e6, 5e6]) -> None:
    """
    Plot the polarizations of a single line of the given FITS file.
    Can be called multiple times.

    Parameters
    ----------
    data : 3d np array
    time : int
    y_lim : list (2,)
        * lower bound y limit
        * upper bound y limit
    """
    (rows, pnum, cols) = data.shape
    pols = cols // 2

    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.set_ylim(lims)
    ax2.set_ylim(lims)
    ax1.plot(d[time][:pols + 1])
    ax2.plot(d[time][pols:-1])


if __name__ == "__main__":
    # data = read_fits(path + "a4002.20230628.b0s1g0.00200.fits")
    # med_data = np.median(data, axis=0)

    # gauss_conv(med_data)
    # robust_avg(med_data)
