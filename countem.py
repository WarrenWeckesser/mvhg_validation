import numpy as np
from cffi import FFI
from _countem_cffi import lib

ffi = FFI()

def countem(indexmap, nsample, samples):
    num_samples, ncolors = samples.shape
    num_support = len(indexmap)
    nbins = indexmap[:, -1].max() + 1
    counts = np.zeros(nbins, dtype=int)

    lib.countem(num_support, ncolors, nsample, nbins,
                ffi.cast("int *", ffi.from_buffer(indexmap)),
                num_samples,
                ffi.cast("int *", ffi.from_buffer(samples)),
                ffi.cast("int *", ffi.from_buffer(counts)))
    return counts


def countem2(indexmap, nsample, samples, counts):
    """
    counts must be a contiguous 1D integer array with length (at least)
    `indexmap[:, -1].max() + 1`.
    """
    num_samples, ncolors = samples.shape
    num_support = len(indexmap)
    nbins = indexmap[:, -1].max() + 1

    lib.countem(num_support, ncolors, nsample, nbins,
                ffi.cast("int *", ffi.from_buffer(indexmap)),
                num_samples,
                ffi.cast("int *", ffi.from_buffer(samples)),
                ffi.cast("int *", ffi.from_buffer(counts)))
