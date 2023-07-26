"""

Authors
-------
Alex Cingoranelli, University of Central Florida
Dr. Emmanuel J. Morales Butler, University of Puerto Rico,
"""

import numpy as np
import mat2py as mp
from mat2py.core import *
from scipy.stats import pearson3
import matplotlib.pyplot as plt
from kurtosis import Kurtosis as K


def process_source_on_off(tabledata, M, source_nums, channels=16384, wlen=1450,
                          bandlen=None):
    """
    Parameters
    ----------
    tabledata :
    M :
    source_nums :

    channels : int, optional
        The number of channels in the FITS file.
    wlen : int, optional
    """

    # TODO This variable represents the length of the band artifact that is
    # present on all 12 meter data at Arecibo.
    # This will be added to the function definition such that if a telescope
    # does not have this specific band shape, this may be set to None.
    bandlen = 400

    # We split the data into its two separate polarizations and save them as A
    # and B respectively.
    A = np.zeroes(M, channels // 2 - wlen * 2 - bandlen * 2)
    B = np.zeroes(M, channels // 2 - wlen * 2 - bandlen * 2)
    x = M[
        [
            M[(windlength + 1):((16384 / 4) - middleBandlength)],
            M[(((16384 / 4) + 1) + middleBandlength):((16384 / 2) - windlength)],
        ]
    ]

    for i in range(M):
        aux_vec = tabledata[1, 1, source_nums][i, 1]
        A[i, :] = aux_vec[]
        B[i, :] = aux_vec[]

    s1_red = np.sum(A, axis=1)
    s1_grn = np.sum(B, axis=1)

    s2_red = np.sum(A ** 2, axis=1)
    s2_grn = np.sum(B ** 2, axis=1)

    mu_red = np.mean(A, axis=1)
    mu_grn = np.mean(B, axis=1)
    var_red = np.var(A)
    var_grn = np.var(B)

    d_red = mu_red ** 2 // np.median(var_red)
    d_grn = mu_grn ** 2 // np.median(var_grn)
