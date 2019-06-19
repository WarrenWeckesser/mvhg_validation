# mvhg_validation
Scripts for validating the multivariate hypergeometric distribution in NumPy.

The scripts "check1.py" and "check2.py" are used to evaluate the multivariate
hypergeometric distribution that is currently in my `new-mvhg` branch of NumPy.

They require `mpmath` and `mpsci` (https://github.com/WarrenWeckesser/mpsci).
"check2.py" also requires matplotlib.

A small C function is used to count the frequencies of the variates in
a random sample.  The Python interface to this function uses cffi.
To build it, run

    $ python countem_build.py
