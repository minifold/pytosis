"""
Kurtosis class for conducting spectral kurtosis to reduce RFI from an input of PSD data.

Author
------
Alex Cingoranelli, University of Central Florida
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.stats import kurtosis, pearson3, skew


class SpectralKurtosis:
    """Functions for computing Spectral Kurtosis.

    Methods
    -------
    * pearson3_ccdf(x, a, b, d)
    * sk_estimator(data)
    * single_scale(data, M=None)
    * sk_edit(data, w)
    * estimate_threshold(data, rho=0.0013499)
    """

    def __init__(self, data=None, M: int = 300, N: int = 2929, d: [float, None] = None):
        """Initialize the Kurtosis class.

        Parameters
        ----------
        data : array-like, optional
            The spectral data to be passed in by the user. Optional in case one
            desires to call the class and then initialize the S1 and S2 parameters
            independently.
        M : int, optional
            The number of spectral estimates. Defaults to 5 minutes, or 300 seconds.
        N : int, optional
            The number of instantaneous FFT PSD spectra taken per second. Defaults to 2929 PSDs,
            which is the amount of spectra produced at a bandwidth of 24 MHz.
        d : int, optional

        Returns
        -------
        See
        ---
        """
        self.data = data
        # the number of instantaneous FFT spectra taken by the mock spectrograph per second.
        self.N = N
        self.d = d  #
        self.M = M  # the number of seconds

        if self.data is not None:
            return self.sk_estimator(self.data)

    def pearson3_ccdf(x, a, b, d):
        r"""Perform the complementary Pearson III CDF.

        Parameters
        ----------
        x : array-like
        a : float
            $\alpha$ shape parameter
        b : float
            $\beta$ scale parameter
        d : float
            $\delta$ location parameter

        Returns
        -------
        _ : array-like
        """
        # TODO check that the values passed in here are in the correct order.
        return 1 - pearson3.cdf(x - d, a, b)

    def sk_estimator(self, data, dtype="averaged"):
        r"""Perform the spectral kurtosis estimation function $\textbb{\hat{SK}}$.

        The `sk_estimator` function makes several assumptions about the data:
        It is defined as:
        .. math: $\langle P \rangle = \Sum{P_i}_{i=1} / N$
        For the Arecibo 12m telescope, the values are taken at
        Which is the PSD averaged over 24 * 2 MHz. One PSD will take 0.34 ms to complete.

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
        # <P> for M spectra
        P = data[:, : self.M]
        s1 = np.sum(P, axis=1)
        s2 = np.sum(np.square(P), axis=1)

        if self.d is None:
            self.d = np.mean(P, axis=1) ** 2 / (np.var(P, ddof=1) * self.N)
            self.d = np.ones(np.shape(d)) * np.median(d)
            # In previous variance calculation, we divided by the median
            # of the variance, which produced inconsistent approximations.

        return ((self.M * self.N * self.d + 1) / (self.M - 1)) * (
            (self.M * (s2 / np.square(s1))) - 1.0)

    def single_scale(self, data, M=None):
        """Calculate the single-scale Spectral Kurtosis routine.

        Perform single-scale spectral kurtosis as specified in the papers
        Smith 2022 and Nita 2010.

        Parameters
        ----------
        data : 2D array-like
            Input data of the format (channels, spectra)
        M : int or None, optional
            spectral kurtosis M-value. If None, the class M attribute will be used.

        Returns
        -------
        out_sk : array-like
        out_s : array-like

        See
        ---
        Section 1.6 Spectral Kurtosis, eqns 1.23-1.32 of Smith 2020
        https://researchrepository.wvu.edu/etd/11467
        """
        if M is None:
            M = self.M

        num_spectra = data.shape[1] // M

        for i in range(num_spectra):
            s_i = data[:, i * M : (i + 1) * M]

            if i == 0:
                out_sk = np.expand_dims(self.sk_estimator(s_i, M), axis=1)
                out_s = np.expand_dims(np.mean(s_i, axis=1), axis=1)
            else:
                out_sk = np.c_[out_sk, self.sk_estimator(s_i, M)]
                out_s = np.c_[out_s, np.mean(s_i, axis=1)]

        return out_sk, out_s

    def sk_edit(self, data, w):
        """
        data : 2D array-like
        w :
        """
        kurtosis = np.array([self.sk_estimator(i) for i in data.T])
        thresh = (1.0 + 100 / data.T.shape[1]) * np.nanmedian(kurtosis)
        w[:, kurtosis >= thresh] = 0
        return kurtosis, w

    def estimate_threshold(self, data, rho: float = 0.0013499):
        r"""Estimate the upper and lower thresholds of the spectral curve.

        Finds the threshold values of the Pearson III CDF and Pearson III CCDF
        using the `scipy.optimize` root solver function.

        Parameters
        ----------
        data : array-like
        rho : float, optional
             1 - 3$\sigma$ = 0.0013499

        Returns
        -------
        upper : int
            Upper threshold as given by the Pearson III CCDF
        lower : int
            Lower threshold as given by the Pearson III CDF

        Notes
        -----
        .. math: $\beta_1 = \frac{\mu_3^2}{\mu_2^2}$
        .. math: $\beta_2 = \frac{\mu_4}{\mu_2^2}$
        The scale parameter is defined as:
        .. math: $\alpha = \frac{\mu_3}{2 \mu_2}$
        The shape parameter is defined as:
        .. math: $\beta = \frac{4 \mu_2^3}{mu_3 ** 2}$
        The location parameter is defined as:
        .. math: $\delta = 1 - \frac{2 \mu_2^2}{\mu_3}$
        The Pearson Criterion is defined as (Pearson 1985):
        .. math: $\kappa = \frac{\beta_1 (\beta_2 + 3)^2}{4 (4 \beta_2 - 3 \beta_1)(2 \beta_2 - 3\beta_1 - 6)}
        When $\kappa$ is > 1, the Pearson distribution is type III.
        """
        mu_2 = np.var(data)
        mu_3 = skew(data)
        _mu_4 = kurtosis(data)

        upper = fsolve(self.pearson3_ccdf(), 1, args=(mu_2, mu_3, rho))
        lower = fsolve(pearson3.cdf(), 1, args=(mu_2, mu_3, rho))
        return upper, lower
