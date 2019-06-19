
import mpmath
import numpy as np
from numpy.random import Generator, PCG64
from mpsci.stats import gtest

from create_mvhg_sample_groups import create_mvhg_sample_groups
from countem import countem, countem2
import matplotlib.pyplot as plt


def run_tests(seed, ntrials, colors, nsample, method, size, indexmap,
              expected):
    gen = Generator(PCG64(seed))

    pvals = []
    for i in range(ntrials):
        s = gen.multivariate_hypergeometric(colors, nsample, size=size,
                                            method=method)
        counts = countem(indexmap, nsample, s)
        stat, p = gtest(counts, expected)
        pvals.append(float(p))
        if (i + 1) % 10 == 0:
            print("%15r %d" % (seed, i+1))

    return pvals


def run_tests2(seed, ntrials, colors, nsample, method, size, indexmap,
               expected):
    nbins = indexmap[:, -1].max() + 1
    counts = np.empty(nbins, dtype=int)

    batchsize = 1000000
    n, r = divmod(size, batchsize)
    batches = [batchsize]*n + [r]*(r > 0)

    gen = Generator(PCG64(seed))

    pvals = []
    for i in range(ntrials):
        counts[...] = 0
        for b in batches:
            s = gen.multivariate_hypergeometric(colors, nsample, size=b,
                                                method=method)
            countem2(indexmap, nsample, s, counts)
        stat, p = gtest(counts, expected)
        pvals.append(float(p))
        if (i + 1) % 10 == 0:
            print("%15r %d" % (seed, i+1))

    return pvals



mpmath.mp.dps = 80

#colors = np.array([5, 5, 25, 25, 15])
#nsample = 25
colors = np.array([1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3])
nsample = 11
size = 20000000

threshold = 100

method = 'count'

print("colors =", colors)
print("sum(colors) =", sum(colors))
print("nsample =", nsample)
print("method =", method)
print("size =", size)
print("threshold =", threshold)

indexmap, binned_expected = create_mvhg_sample_groups(colors, nsample, size,
                                                      threshold, dir='pmf')
indexmap_array = np.column_stack((list(indexmap.keys()),
                                  list(indexmap.values())))

ntrials = 750


def func(seed):
    return run_tests2(seed, ntrials, colors, nsample, method, size,
                      indexmap_array, binned_expected)


if __name__ == "__main__":
    import multiprocessing
    #seeds = [129055329838793, 487126922734978, 342291325629100, 912315119763521]
    seeds = [219055379838793, 487926922034978, 342291325329130, 912395119263521]
    pool = multiprocessing.Pool(processes=len(seeds))
    all_pvals = pool.map(func, seeds)

    pvals = np.array(all_pvals).ravel()

    nsamples = len(pvals)
    nbins = min(int(len(pvals)/10 + 0.5), 20)
    plt.hist(pvals, range=[0, 1], bins=nbins, alpha=0.75)
    plt.grid()
    plt.axhline(nsamples/nbins, color='k', linestyle='--', alpha=0.5)
    plt.title(("Histogram of p-values of multivariate hypergeometric samples\n"
               "colors=%s, nsample=%s, size=%s, method=%r\n"
               "threshold=%s, nsamples=%s") % (colors, nsample, size, method,
                                               threshold, nsamples),
              fontsize=8)
    plt.savefig('histogram.png')
    plt.close('all')
