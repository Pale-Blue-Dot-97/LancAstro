# =================================================== HISTOPLOT =======================================================
#
# Script to plot histograms of data stored in FITS files, ascii tables
#
# Author: Harry Baker
# Date: 31.07.2019
# Version: 1.4
# ==VERSION 1==
#
# Part of the LancAstro.py project for XGAL
#
# Edited from HISTO_Environment.py by Dr David Sobral 2008
#
# Reworked from HistoPlot.py to utilise Plot2D's functions
#

# =====================================================================================================================
#                                                     IMPORTS
# =====================================================================================================================

import numpy as np
import Dependencies.progress_bar as pb

# =====================================================================================================================
#                                                    FUNCTIONS
# =====================================================================================================================


def _get_outer_edges(a):
    """
    Determine the outer bin edges to use, from either the data or the range
    argument
    
    Edited from Numpy source code:
    https://github.com/numpy/numpy/blob/master/numpy/lib/histograms.py

    """

    RANGE = np.min(a), np.max(a)
    if RANGE is not None:
        first_edge, last_edge = RANGE
        if first_edge > last_edge:
            raise ValueError(
                'Maximum value must be larger than minimum in range parameter.')
        if not (np.isfinite(first_edge) and np.isfinite(last_edge)):
            raise ValueError(
                "Supplied range of [{}, {}] is not finite".format(
                    first_edge, last_edge))
    elif a.size == 0:
        # handle empty arrays. Can't determine range, so use 0-1.
        first_edge, last_edge = 0, 1
    else:
        first_edge, last_edge = np.min(a), np.max(a)
        if not (np.isfinite(first_edge) and np.isfinite(last_edge)):
            raise ValueError(
                "Autodetected range of [{}, {}] is not finite".format(
                    first_edge, last_edge))

    # expand empty range to avoid divide by zero
    if first_edge == last_edge:
        first_edge = first_edge - 0.5
        last_edge = last_edge + 0.5

    return first_edge, last_edge


def make_bins(a, bin_min=None, bin_max=None, bin_wid=1, n=10):
    """
    """

    # Calculates the outer edges of the bins
    if auto_bins:
        bin_min, bin_max = _get_outer_edges(a)

    if not use_width:
        bin_wid = (bin_max - bin_min) / float(n)

    # Constructs the centre of the bins from the defined mid, max and width
    bins = np.arange(bin_min + bin_wid / 2.0, bin_max - bin_wid / 2.0, bin_wid)

    # Creates an array of the same length as bins to hold the counts per bin
    counts = bins - bins

    # Loops through all data points, for each one cycling through the bins to
    # see if it resides within it. If so, it is added to counts for that bin
    for j in range(len(a)):
        for k in range(len(bins)):
            if a[j] < bins[k] + bin_wid / 2.0:
                if a[j] > bins[k] - bin_wid / 2.0:
                    counts[k] = counts[k] + 1.0

    # Normalises the data
    norm = sum(counts)

    if normalise:
        counts = counts / norm

    return counts, bins, bin_min, bin_max


def output_bins_to_plot(a, bin_min=None, bin_max=None, bin_wid=1, n=10):
    """
    """

    counts, bins, bin_min, bin_max = make_bins(a, bin_min, bin_min, bin_wid, n)

    # Constructs the extra points either side of the mid-point of the bin to
    # plot
    x = [bin_min]  # Array of x-axis plot points (bin edges)
    y = [np.max(counts)]  # Array of y-axis plot points (count values for each bin)

    for j in range(len(bins)):
        # Left edge
        x.append(bins[j] - bin_wid / 2.0)
        y.append(counts[j])

        # Central value
        x.append(bins[j])
        y.append(counts[j])

        # Right edge
        x.append(bins[j] + bin_wid / 2.0)
        y.append(counts[j])

    x.append(bin_max)
    y.append(np.min(counts))

    return x, y, bin_min, bin_max


# =====================================================================================================================
#                                                   PARAMETERS
# =====================================================================================================================
# BIN PARAMETERS ======================================================================================================

use_width = False  # Defines if WID is used (True) to construct bins
# or n is used (False)

auto_bins = True  # Auto create MIN, MAX, WID
normalise = False  # Normalises counts if true

# END PROGRAM =========================================================================================================
