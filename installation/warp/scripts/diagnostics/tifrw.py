"""Set of functions to read and write tif photos from numpy arrays in Python
# Module tifrw.py
#
# by: Rami Kishek and A. Valfells
# created: April 20, 2005
#
#        Last Modified: 4/20/2005
#
# Set of functions to read and write tif photos from numpy arrays in Python
# save_tif ...  Saves photo array into 'tif' file
# read_tif ...  Reads photo array from 'tif' file
# ==================
"""
# GLOBALS + INITIALIZATION

yes = 1
no = 0
import numpy


def tifrwdoc():
    import tifrw
    print tifrw.__doc__


################## Save Array to Tif ##########################################

def save_tif(matrix, filename=None):
    """ save_tif(matrix, filename = "temp.tif")
    Saves a 2-D array "matrix" into a tif picture file.
    """
    S = numpy.ravel(matrix)
    M, N = numpy.shape(matrix)

    min_val = float(numpy.minimum.reduce(S))
    max_val = float(numpy.maximum.reduce(S))
    if max_val != min_val:
        matrix = (matrix - min_val)/(max_val - min_val) * 255  # Preprocessor
    matrix = matrix.astype(numpy.ubyte)                        # Convert to binary
    matrix = numpy.transpose(matrix)                           # Preprocess for tif-ization

    if filename is None:
        filename = "temp.tif"

    tif = "P5\n#TIF version of array\n%d %d\n255\n%s" % (M, N,
                                                         numpy.ravel(matrix).tostring())
    with open(filename, 'wb') as f_tif:
        f_tif.write(tif)


################## Read Array from Tif ########################################

def read_tif(phpath):
    """ read_tif(phpath): read tif photo speicified by phpath and
        return as a 2-D array
    """
    with open(phpath, 'rb') as f_tif:
        tif = f_tif.read()

    header, matrix = tif.split('\n255\n')
    dims = tuple([int(s) for s in header.split('\n')[-1].split()])

    matrix = numpy.array(list(matrix))
    matrix = numpy.reshape(matrix, dims)
    dummy = matrix.astype('l')
    matrix = numpy.where(dummy < 0, dummy+255, dummy)
    return numpy.transpose(matrix), dims
