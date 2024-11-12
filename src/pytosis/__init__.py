#!/usr/bin/env python3
"""Package for excising unwanted signals using the spectral kurtosis method.

Radio-frequency interference (RFI) is any source of transmission that is
observed in the frequency band other than the desired celestial sources.
Because transmitters on and around the Earth can be many times stronger than
the signal of interest, RFI is a major concern when analyzing radio astronomy data.

Modules exported by this package:

- `kurtosis`: Provides the Spectral Kurtosis class.
- `reduction`: Provides several noise reduction methods.
- `cli`: The command-line interface for working with files.
- `process sources`: Provides methods for reading and manipulating Arecibo FITS files.
"""
# read version from installed package
from importlib.metadata import version

__version__ = version("pytosis")

