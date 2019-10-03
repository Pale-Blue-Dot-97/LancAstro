# =================================================== ASCII_READ ======================================================
"""Set of definitions to read columns from a whitespace delimited file with '#' comments (e.g. Sextractor file)

By Jim Geach
Modified by Philip Best (August 2007)
Modified by David Sobral (November 2007)
Modified by Harry Baker (October 2019) to comply with PEP
"""
# =====================================================================================================================
#                                                       IMPORTS
# =====================================================================================================================

import numpy as np


# =====================================================================================================================
#                                                      FUNCTIONS
# =====================================================================================================================
def count(f):
    """Counts number of lines in a file

    Args:
        f (str): Name of file

    Returns:
        Number of lines in file

    """

    file_object = open(f)
    alldata_lines = file_object.readlines()
    file_object.close()

    return len(alldata_lines)


def skip(f):
    """Counts commented lines at start of file

    Args:
        f (str): Name of file

    Returns:
        Number of commented lines in file

    """

    file_object = open(f)
    alldata_lines = file_object.readlines()
    file_object.close()

    n = 0
    while alldata_lines[n][0] == '#':
        n += 1

    return n


def read_col(c, f):
    """Reads a column "c" from a file, skipping the commented lines at start of file

    Args:
        c (str): Name of column
        f (str): Name of file

    Returns:
        Array representing column 'c' in file 'f'

    """

    L = count(f)
    sk = skip(f)
    L -= sk

    file_object = open(f)
    alldata_lines = file_object.readlines()
    file_object.close()

    D = np.array(range(L), np.float32)

    for i in range(L):
        line = alldata_lines[i + sk]
        sp = line.split()
        D[i] = float(sp[c - 1])

    return D


def mread_col(cs, f):
    """Reads multiple columns "cs" from a file, skipping the commented lines at start of file

    Args:
        cs ([str]): List of names of columns to read
        f (str): Name of file to read from

    Returns:
        2D array of columns 'cs' of data read from file 'f'

    """

    L = count(f)
    sk = skip(f)
    L -= sk

    file_object = open(f)
    alldata_lines = file_object.readlines()
    file_object.close()

    O = []

    for j in cs:
        D = np.array(range(L), np.float32)
        for i in range(L):
            line = alldata_lines[i + sk]
            sp = line.split()
            D[i] = float(sp[j - 1])
        O.append(D)
    return O
