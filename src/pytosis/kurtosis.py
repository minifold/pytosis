"""
Kurtosis class for conducting spectral kurtosis to reduce RFI from an input of PSD data.

Author
------
Alex Cingoranelli, University of Central Florida
"""

import numpy as np
from scipy.stats import pearson3, skew, kurtosis


class Kurtosis:
    def __init__(self, data=None, M=300, n: int = 1, d: [float, None] = None):
        """
        Initializes the Kurtosis class.

        Parameters
        ----------
        data : array-like, optional
            The spectral data to be passed in by the user. Optional in case one
            desires to call the class and then initialize the S1 and S2 parameters
            independently.
        M : int, optional
            The number of spectral estimates. Defaults to 5 minutes, or 300 seconds.
        n : int, optional
            The number of spectral channels.
        d : int, optional

        """
        self.data = data
        self.n = n  # the number of channels
        self.d = d
        self.M = M  # the number of seconds

        if self.data is not None:
            return self.sk_estimator(self.data, self.M)

    def sk_estimator(self, data, M):
        """
        Parameters
        ----------
        data : array-like
            The given power-sepctral density estimate $\hat{P_k}$
        M : int
            The number of spectral estimates

        Returns
        -------
        sk_est : array-like
            Spectral kurtosis estimator.

        Notes
        -----
        The $SK$ estimator for a set of M samples in the same freq. channel is defined as
        .. math: $\hat{SK} = \hat{V}_k^2 = \frac{\hat{\sigma}^2}{\hat{\mu}^2}$
        The summations $S_1$ and $S_2$ are the summations of the
        .. math: $S_1 = \Sum{\Sum{\hat{P_ki}}}$
        .. math: $S_2 = \Sum{\Sum{\left(\hat{P_ki}\right)}^{2}}$
        ..
        If N spectra are averaged together before passing through the $\hat{SK}$
        detection, then the previous equation becomes
        ... math: $\hat{SK} = \frac{MNd + 1}{M-1}\left(\frac{MS_1}{S_2^2 - 1}\right)$

        See
        ---
        Section 3.1 Spectral Kurtosis Estimator, eqns 19-21 of Nita et. al 2007
        https://www.jstor.org/stable/10.1086/520938
        """

        a = data[:, :M] * self.n
        s1 = np.sum(np.sum(a, axis=1), axis=0)
        s2 = np.sum(np.sum(a, axis=1) ** 2, axis=0)

        if self.d is None:
            self.d = np.mean(a, axis=1) ** 2 / np.var(a, ddof=1)
            # Based on previous variance calculation where we take the median
            # of the variance
            # self.d = np.mean(a, axis=1) ** 2 / np.median(np.var(a, ddof=1))

        return ((M * self.n * self.d + 1) / (M - 1)) * ((M * (s2 / np.square(s1))) - 1.)

    def single_scale(self, data, M):
        """
        Calculates the single-scale Spectral Kurtosis routine as specified in
        Smith 2022 and Nita 2010.

        Parameters
        ----------
        data : 2D array-like
            Input data of the format (channels, spectra)
        M : int, optional
            spectral kurtosis M-value

        Returns
        -------
        out_sk : array-like
        out_s : array-like

        See
        ---
        Section 1.6 Spectral Kurtosis, eqns 1.23-1.32 of Smith 2020
        https://researchrepository.wvu.edu/etd/11467
        """

        num_spectra = data.shape[1] // M

        for i in range(num_spectra):
            s_i = data[:, i * M:(i + 1) * M]

            if i == 0:
                out_sk = np.expand_dims(self.sk_estimator(s_i, M), axis=1)
                out_s = np.expand_dims(np.mean(s_i, axis=1), axis=1)
            else:
                out_sk = np.c_[out_sk, self.sk_estimator(s_i, M)]
                out_s = np.c_[out_s, np.mean(s_i, axis=1)]

        return out_sk, out_s

    def pearson_criterion():
        mu_2 = np.var(sample)
        mu_3 = skew(sample)

    def sk_edit(self, data, w):
        """
        data : 2D array-like
        w :
        """

        kurtosis = np.array([self.sk_estimator(i) for i in data.T])
        thresh = (1. + 100/data.T.shape[1]) * np.nanmedian(kurtosis)
        w[:, kurtosis >= thresh] = 0
        return kurtosis, w

    # Threshold
