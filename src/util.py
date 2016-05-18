"""
Generic utilities.
"""
import numpy as np


def mirror(X):
    """
    Given an array X = [x0, x1, x2, ...], returns a "mirrored" array,
    [..., x2, x1, x0, x1, x2, ...].
    """
    X_rev = np.flipud(X[1:])
    return np.concatenate((X_rev, X))
