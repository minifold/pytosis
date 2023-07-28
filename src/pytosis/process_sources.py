"""

Authors
-------
Alex Cingoranelli, University of Central Florida
Dr. Emmanuel J. Morales Butler, University of Puerto Rico,
"""

import os
import numpy as np
from astropy.io import fits
from astropy.stats import histogram
from scipy.stats import pearson3

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from kurtosis import SpectralKurtosis as SK


def read_fits(fname, showplot=False, path="/share/pdata1/pdev/"):
    """
    Reads in a FITS file and creates a variety of plots of the data.

    Parameters
    ----------
    fname : string
        File name that will be opened.

    Returns
    -------
    data : 2d nparray
        Data of size (time, channels) to be analyzed.
    showplot : boolean, optional
        Determines whether the 2-D plotting image is shown.
    """

    # First check if file exists.
    # We will assume that data is being stored in the `pdev` directory:
    # We want to accept the string whether its an absolute path or a relative path.
    if not os.path.exists(os.path.realpath(fname)) or \
       not os.path.exists(os.path.realpath(path + fname)):
        raise FileNotFoundError("The given FITS file does not exist.")

    with fits.open(fname, memmap=True, lazy_load_hdus=True,
                   output_verify="silentfix") as f:

        rows = len(f[1].data['DATA'])
        # cols is the number of channels.
        cols = len(f[1].data['DATA'][0])
        # Find the polarizations and save that to pols.
        # Easier than converting to a 2D datacube and then reconverting it back
        pols = cols // 2

        data = np.zeros((rows, cols))
        print(rows, cols)

        # This must be done because when working with the astropy Table the
        # values will not be written and astropy will throw an error because
        # the DATA column is classified as an object, and fits.writeto() will
        # not allow the `object` datatype to be written.
        for i in range(rows):
            data[i, :] = f[1].data['DATA'][i]

            f[1].data['DATA'][i] = data[i, :]

    if showplot is True:
        plt.imshow(data, aspect="auto", norm=LogNorm(
            np.min(data), np.max(data)))
        plt.show()

    # f.writeto("temp" + ext, overwrite=True)
    return data


def histbins(M, sk):
    """
    """
    # When M = 300 seconds
    if (M == 300):
        sample = sk[np.where(sk < 1.5)]
    else:
        sample = sk[np.where(sk < 5.0)]

    ym, bin_edges = histogram(sample, bins="scott", density=True)
    return ym, bin_edges


def process_source_on_off(data, M, source_nums, channels=16384, wlen=1450,
                          bandlen=400):
    """
    Parameters
    ----------
    data : array-like
    M :
    source_nums :
        number of sources.
    channels : int, optional
        The number of channels in the FITS file.
    wlen : int, optional
    """

    # TODO This variable represents the length of the band artifact that is
    # present on all 12 meter data at Arecibo.
    # If a telescope does not have this specific band shape, this may be set to
    # `None`. There should probably be a better default here.
    if bandlen is None:
        bandlen = 0

    # We split the data into its two separate polarizations and save them as A
    # and B respectively.
    A = np.zeroes(M, channels // 2 - wlen * 2 - bandlen * 2)
    B = np.zeroes(M, channels // 2 - wlen * 2 - bandlen * 2)

    # I have already converted the data into a 2d array through the `read_fits`
    # method. Now we need to remove the artifacts that are known issues.
    # Clip the data to remove the central peak and "wings" of the spectra.
    A = data[wlen:(channels // 4) - bandlen, :]
    B = data[(channels // 4) + bandlen:(channels // 2) - wlen + 1, :]

    # The following are both taken care of when the kurtosis class is initialized.
    k = SK()
    sk_red = k.sk_estimator(A)
    sk_grn = k.sk_estimator(B)

    ym_red, edges_red = histbins(M, sk_red)
    ym_grn, edges_grn = histbins(M, sk_grn)


if __name__ == "__main__":
    pass

# K = Nd
# theta = a / N
# mu = Nd * a/N
# mu = da
# var = k * theta ^ 2
# var = da^2 / N
# d = mu ^ 2 / N var
