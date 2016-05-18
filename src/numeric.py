"""
Numeric routines.
"""

def arr_roots(x):
    """
    Find the roots of an array.
    """
    l = x[ :-1]
    r = x[1:  ]

    return ((l > 0) & (r < 0)) | ((l < 0) & (r > 0)) | (l == 0)
