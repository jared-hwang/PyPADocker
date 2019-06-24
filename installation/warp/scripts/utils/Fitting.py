""" Module Fitting.py

by:      Rami A. Kishek
Created: July 17, 2001

Last Modified: 7/17/01

This module contains the following fitting functions of general use:

lsqfit ... Least Squares fit of a straight line to approximate a set of pts.
"""
import numpy


def Fittingdoc():
    import Fitting
    print Fitting.__doc__


def lsqfit(x, y):
    """ (a, b, L) = lsqfit(x, y)
        Given arrays 'x' and 'y', fits a straight line of the form
        y = ax + b in order to minimize the least squares parameter
        L = (1/N)sum[(yn - a*xn - b)^2, n=1..N]

        Note that x and y are converted to arrays of float before operation.
    """
    N = len(x)
    x = numpy.array(x, 'd')
    y = numpy.array(y, 'd')
    xm = numpy.sum(x)/N
    ym = numpy.sum(y)/N
    a = ((numpy.sum(x*y)/N)-xm*ym)/((numpy.sum(x*x)/N) - xm**2)
    b = ym - a*xm
    L = numpy.sum((y - a*x - b)**2)/N
    return (a, b, L)
