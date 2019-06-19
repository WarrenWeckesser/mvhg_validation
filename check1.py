
import numpy as np
from numpy.random import Generator, PCG64
import mpmath
from mpsci.distributions import multivariate_hypergeometric
from mpsci.stats import gtest


mpmath.mp.dps = 20


# How many p-values to compute for each set of parameters.
ntrials = 3

# Minimum expected count of any element in the support for the G-test.
min_count = 100

gen = Generator(PCG64(11223344556677))

for method in ['marginals', 'count']:
    print()
    print("=== method:", method, "===")
    # Be sure any parameters for which this code is run have a PMF whose
    # smallest nonzero value is large enough such min_count/min(pmf.values())
    # gives a reasonable sample size.
    for colors, nsample in [([10, 90], 40), ([5, 10, 10], 12),
                            ([40, 6, 5], 19), ([40, 6, 5], 30),
                            ([3, 4, 5, 6], 8), ([2, 3, 3, 3, 3, 4, 5], 13)]:

        pmf = multivariate_hypergeometric.pmf_dict(colors, nsample)

        pmin = float(min(pmf.values()))
        size = int(min_count/pmin + 0.5)
        expected = {key:size*p for key, p in pmf.items()}

        print()
        print("colors:", colors, "  nsample:", nsample, "size:", size, end=' ')
        print("  min expected count:", float(min(expected.values())))

        for i in range(ntrials):
            rvs = gen.multivariate_hypergeometric(colors, nsample,
                                                  size=size, method=method)
            count_array = np.zeros(np.array(colors)+1, dtype=int)
            np.add.at(count_array, tuple(rvs.T), 1)
            indices = np.nonzero(count_array)
            observed_samples = np.column_stack(indices)
            counts = count_array[indices]
            support = np.array(list(pmf.keys()))
            assert np.all(support == observed_samples)
            stat, pvalue = gtest(counts, expected.values())  
            print("p-value:", pvalue)
