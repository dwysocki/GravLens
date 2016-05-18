"""
Generic utilities.
"""
import numpy as np

from os import makedirs
from os.path import isdir


def mirror(X, negate=False):
    """
    Given an array X = [x0, x1, x2, ...], returns a "mirrored" array,
    [..., x2, x1, x0, x1, x2, ...].
    """
    sign = -1 if negate else +1
    X_rev = np.flipud(X[1:])
    return np.concatenate((sign*X_rev, X))


def make_sure_path_exists(path):
    """
    Creates the supplied *path* if it does not exist.
    Raises *OSError* if the *path* cannot be created.

    **Parameters**
    path : str
        Path to create.

    **Returns**
    None
    """
    try:
        makedirs(path)
    except OSError:
        if not isdir(path):
            raise
