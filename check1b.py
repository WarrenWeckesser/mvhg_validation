import sys
import numpy as np
from numpy.random import Generator, PCG64
import mpmath
from mpsci.distributions import multivariate_hypergeometric
from mpsci.stats import gtest


if len(sys.argv) != 4:
    print("Required command line arguments: ntrials seed pvalue_filename")
    sys.exit(-1)

# How many p-values to compute for each set of parameters.
ntrials = int(sys.argv[1])

seed = int(sys.argv[2])
pvalue_filename = sys.argv[3]

mpmath.mp.dps = 20

# Minimum expected count of any element in the support for the G-test.
min_count = 100

gen = Generator(PCG64(seed))

colors = [1, 2, 1, 2, 1, 2, 3, 4, 3, 3, 2]
nsample = 12
method = 'count'

pmf = multivariate_hypergeometric.pmf_dict(colors, nsample)

pmin = float(min(pmf.values()))
size = int(min_count/pmin + 0.5)
expected = {key:size*p for key, p in pmf.items()}

print("colors:", colors, "  nsample:", nsample, "method:", method, "size:", size, end=' ')
print("  min expected count:", float(min(expected.values())))

batch_size = 20000000
out = open(pvalue_filename, 'w')

for i in range(ntrials):
    count_array = np.zeros(np.array(colors)+1, dtype=int)

    num_batches, final_batch_size = divmod(size, batch_size)
    sizes = [batch_size] * num_batches + ([final_batch_size] if final_batch_size else [])
    for k, current_batch_size in enumerate(sizes):
        print(' %04d/%04d' % (k, len(sizes)), end='\r', flush=True)
        rvs = gen.multivariate_hypergeometric(colors, nsample,
                                              size=current_batch_size,
                                              method=method)
        np.add.at(count_array, tuple(rvs.T), 1)

    indices = np.nonzero(count_array)
    observed_samples = np.column_stack(indices)
    counts = count_array[indices]
    support = np.array(list(pmf.keys()))
    assert np.all(support == observed_samples)
    stat, pvalue = gtest(counts, expected.values())  
    print("p-value:", pvalue)
    print(pvalue, file=out)

out.close()
