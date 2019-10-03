# =================================================== DATALOAD ========================================================
"""
Module to load data from FITS files into a 2D array.
Part of LancAstro.py project

Author: Harry Baker

Example:

Todo:
    * Upgrade data_load() to handle ASCII, csv and excel files

"""

# =====================================================================================================================
#                                                     IMPORTS
# =====================================================================================================================

from astropy.io import fits as pyfits       # To handle FITS files
import Plot2D as laplt                      # Uses Plot2D's settings as master
import Dependencies.progress_bar as pb      # To display a progress bar in output if print_on==True in Plot2D


# =====================================================================================================================
#                                                    FUNCTIONS
# =====================================================================================================================
# DATA LOADING ========================================================================================================

def data_load(DATANAME, COLUMNNAME, PATH=['']):
    """Load data in from columns in a FITS file
    
    Args:
        DATANAME ([str]): Array of strings with names of the FITS file containing the data to be loaded

        COLUMNNAME ([[str]]): 2D array of strings with first axis being of equal length to DATANAME, and the second axis
            being the names of the columns to be loaded from each FITS file in DATANAME

        PATH: Optional variable that is a list of strings for paths to each FITS file. Should be same length as DATANAME

    Returns:
        A 2D array of data requested from FITS files

    """
    # Defines the array that will hold all the variables to plot
    DATASET = []

    # Checks that PATH is the same length as DATANAME and corrects if not
    PATH = laplt.default_lengths(PATH, len(DATANAME))

    if laplt.print_on:
        # Initialises the progress bar for loading the data
        pb.printProgressBar(0, 2 * len(DATANAME) * sum(len(x) for x in COLUMNNAME), prefix='LOADING DATA:',
                            suffix='COMPLETE', length=40)

    # Loops through the FITS files
    for i in range(len(DATANAME)):
        # Exception handling in case file cannot be opened
        try:
            # Opens file into memory
            FILE = pyfits.open("%s%s" % (PATH[i], DATANAME[i]))

            if laplt.print_on:
                # Updates progress bar
                pb.printProgressBar(i, 2 * len(DATANAME) * sum(len(x) for x in COLUMNNAME),
                                    prefix='OPENING DATA FILE %s%s:' % (PATH[i], DATANAME[i]),
                                    suffix='COMPLETE', length=40)

            # Loads Extension 1 of the FITS file into memory
            TABLE = FILE[1].data

            # Loops through and loads all the columns in the file to be loaded
            for j in range(len(COLUMNNAME[i])):

                if laplt.print_on:
                    # Progress bar updated to reflect this
                    pb.printProgressBar(i + j, 2 * len(DATANAME) * sum(len(x) for x in COLUMNNAME),
                                        prefix='LOADING COLUMN %s:' % COLUMNNAME[i][j], suffix='COMPLETE', length=40)

                # Loads column into memory
                VAR = TABLE.field("%s" % COLUMNNAME[i][j])

                # Removes any 'nan' (Not an Actual Number) strings that are
                # placeholders if there are blank entries as these would throw
                # errors
                CleanVAR = [x for x in VAR if str(x) != 'nan']

                # Adds 'cleaned' variable into DATASET
                DATASET.append(CleanVAR)

                # If file cannot be opened, throws an exception
        except:
            print("FILE COULD NOT BE OPENED")

    if laplt.print_on:
        # Completes progress bar once all variables are loaded into DATASET
        pb.printProgressBar(2 * len(DATANAME) * sum(len(x) for x in COLUMNNAME),
                            2 * len(DATANAME) * sum(len(x) for x in COLUMNNAME),
                            prefix='DATASET LOADED READY FOR PLOTTING!',
                            suffix='COMPLETE!', length=40)

    return DATASET  # Returns 2D array of variables
