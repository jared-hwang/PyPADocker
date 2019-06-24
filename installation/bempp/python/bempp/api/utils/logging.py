"""This module contains functions to initialize the Python logger"""
#pylint: disable=import-self
#pylint: disable=no-member
#pylint: disable=redefined-builtin
#pylint: disable=invalid-name
from __future__ import absolute_import
import logging as _logging
import time as _time
from mpi4py import MPI as _MPI

# Logging levels

DEBUG = _logging.DEBUG
INFO = _logging.INFO
WARNING = _logging.WARNING
ERROR = _logging.ERROR
CRITICAL = _logging.CRITICAL


DEFAULT_FORMAT = ("%(asctime)s:Host {0}:Rank {1}:".format(
    _MPI.Get_processor_name(), _MPI.COMM_WORLD.Get_rank()) +
                  "%(name)s:%(levelname)s: %(message)s")

def _init_logger():
    """Initialize the BEM++ logger."""

    logger = _logging.getLogger()
    logger.setLevel(_logging.DEBUG)
    logger.addHandler(_logging.NullHandler())
    return logger

def log(message, level=INFO, flush=True):
    """Log including default flushing for IPython."""
    from bempp.api import _LOGGER
    _LOGGER.log(level, message)
    if flush:
        flush_log()

def flush_log():
    """Flush all handlers. Necessary for Jupyter."""
    from bempp.api import _LOGGER
    for handler in _LOGGER.handlers:
        handler.flush()


def enable_console_logging(level=DEBUG):
    """Enable console logging and return the console handler."""

    import bempp.api
    if not bempp.api._console_logging_handler:
        ch = _logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(_logging.Formatter(DEFAULT_FORMAT, "%H:%M:%S"))
        bempp.api._LOGGER.addHandler(ch)
        bempp.api._console_logging_handler = ch
    return bempp.api._console_logging_handler


def enable_file_logging(file_name, level=DEBUG, format=DEFAULT_FORMAT):
    """Enable logging to a specific file."""

    from bempp.api import _LOGGER
    fh = _logging.FileHandler(file_name)
    fh.setLevel(level)
    fh.setFormatter(_logging.Formatter(format, "%H:%M:%S"))
    _LOGGER.addHandler(fh)
    return fh


def set_logging_level(level):
    """Set the logging level."""

    from bempp.api import _LOGGER
    _LOGGER.setLevel(level)


def timeit(message):
    """Decorator to time a method in BEM++"""
    from bempp.api import global_parameters

    def timeit_impl(fun):
        """Implementation of timeit."""
        def timed_fun(*args, **kwargs):
            """The actual timer function."""
            if not global_parameters.verbosity.extended_verbosity:
                return fun(*args, **kwargs)

            st = _time.time()
            res = fun(*args, **kwargs)
            et = _time.time()
            bempp.api.log(message + " : {0:.3e}s".format(et-st))
            return res

        return timed_fun

    return timeit_impl

#pylint: disable=too-few-public-methods
class Timer:
    """Context manager to measure time in BEM++."""

    def __init__(self):
        """Constructor."""
        self.start = 0
        self.end = 0
        self.interval = 0

    def __enter__(self):
        self.start = _time.time()
        return self


    def __exit__(self, *args):
        self.end = _time.time()
        self.interval = self.end - self.start
