
from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef("""
    int countem(int num_support, int ncolors, int nsample, int nbins, int *indexmap,
                int num_samples, int *samples, int *counts);
""")

ffibuilder.set_source("_countem_cffi",
    """
    #include "countem.h"
    """,
    sources=['countem.c'],
)


if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
