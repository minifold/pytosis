#! /usr/bin/env python3
# pytosis/cli.py
"""Provides a command-line interface for the user to interact with the package.

This module includes an argument parser for various inputs that may be passed
by the user.

Contains the following functions:
- `parse_cli()` - Takes string inputs from STDIN and calls functions in the pytosis package.
"""

from __future__ import annotations
from argparse import ArgumentParser
from typing import Any


def parse_cli() -> dict[str, Any]:
    """
    Takes the cli inputs and calls functions based on the given arguments.

    Returns
    -------
    args : dict
        The arguments passed through the command-line.
    """

    parser = ArgumentParser(
        prog="readfits",
        description="Reads in a fits file from the pdata directory and then "
        + "performs basic averaging calculations and optionally "
        + "plots them.",
    )

    # Arguments passed in through the CLI.
    parser.add_argument(
        "file(s)", nargs="+"
    )  # There has to be at least one FITS file passed.
    parser.add_argument(
        "-p", "--plot", action="store_true", help="Shows the plots for every function."
    )
    parser.add_argument(
        "-f",
        "--filter-method",
        choices=["conv", "median", "sigmaclip", "kurtosis"],
        default="median",
        help="Algorithm for reducing RFI. (default: %(default)s)",
    )
    # Version information.
    parser.add_argument("--version", action="version", version="%(prog)s 0.1.0")

    options, args = parser.parse_args()
    return vars(options)
