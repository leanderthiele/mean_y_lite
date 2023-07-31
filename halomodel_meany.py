# Wrapper around the C code to be called from python

from sys import stdout

import ctypes as ct

import numpy as np


_libhalomodel_meany = ct.CDLL('./libhalomodel_meany.so')

# the only function to be used
meany = _libhalomodel_meany.meany

# returns negative if error encountered
meany.restype = ct.c_double

meany.argtypes = [ct.c_char_p, # class.ini file
                  ct.c_char, # either 'y' or 'T', whether we want <y> or <yT>/<y>
                  ct.c_int, # flag whether custom hmf to be used
                  ct.c_int, # number of z sample points
                  ct.c_int, # number of mass sample points

                  # the redshift sample points
                  np.ctypeslib.ndpointer(dtype=ct.c_double, ndim=1, flags='C_CONTIGUOUS'),

                  # the log(mass) sample points
                  np.ctypeslib.ndpointer(dtype=ct.c_double, ndim=1, flags='C_CONTIGUOUS'),

                  # the data, of shape [z x logM]
                  np.ctypeslib.ndpointer(dtype=ct.c_double, ndim=1, flags='C_CONTIGUOUS'),
                 ]

# utility function to get a C string from a python byte/str instance
def c_str(some_str) :
    
    assert isinstance(some_str, (str, bytes))

    if isinstance(some_str, str) :
        some_str = some_str.encode(stdout.encoding)
    
    return ct.c_char_p(some_str)

# utility function to get a C char from a python 1-character string
def c_char(some_str) :

    assert isinstance(some_str, (str, bytes)) 
    assert len(some_str) == 1

    if isinstance(some_str, str) :
        some_str = some_str.encode(stdout.encoding)

    return ct.c_char(some_str)
