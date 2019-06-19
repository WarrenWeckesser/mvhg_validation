import os
import mpmath
import numpy as np
from mpsci.stats import gtest
from mpsci.distributions import multivariate_hypergeometric


mpmath.mp.dps = 80


def create_mvhg_sample_groups(colors, nsample, size, threshold, dir=''):
    colors = np.asarray(colors)
    filename = "pmf" + "-".join(str(c) for c in colors) + "ns" + str(nsample) + ".txt"
    if dir != '' and not os.path.exists(dir):
        os.mkdir(dir)
    fullname = os.path.join(dir, filename)
    if os.path.exists(fullname):
        # We already computed the PMF for these parameters.
        print("Reading PMF from '%s'" % (fullname,))
        with open(fullname, 'r') as f:
            expected = {}
            for line in f:
                fields = line.split()
                sample = tuple(int(field) for field in fields[:-1])
                pval = mpmath.mp.mpf(fields[-1])
                expected[sample] = pval
    else: 
        print("Generating PMF as a dictionary")
        expected = multivariate_hypergeometric.pmf_dict(colors.tolist(), nsample)
        print("Saving PMF to '%s'" % (fullname,))
        items = sorted(expected.items())
        with open(fullname, 'w') as f:
            for sample, pval in items:
                f.write(" ".join(str(c) for c in sample))
                f.write(" %s\n" % pval)

    #print("generating pmf as a dictionary")
    #expected = mvhg_pmf(colors.tolist(), nsample)
    print("Scaling values to expected frequencies")
    for key in expected:
        expected[key] = float(expected[key]*size)

    print("Sorting")
    vals = list(expected.values())
    sortorder = sorted(range(len(expected)), key=vals.__getitem__, reverse=True)
    #keys = sorted(expected.keys(), key=vals.__getitem__)
    tmpkeys = list(expected.keys())
    keys = [tmpkeys[k] for k in sortorder]
    #values = sorted(vals)
    values = np.array([vals[k] for k in sortorder])

    print("Size of support is", len(expected))
    print("Minimum expected value is", min(values))
    print("Maximum expected value is", max(values))
    print()

    # Form groups of samples whose individual expected value
    # is less than the threshold.

    lowmask = values < threshold

    if np.any(lowmask):
        firstlow = np.argmax(values < threshold)
        lowvals = values[firstlow:]

        print("number of expected values below threshold:", len(lowvals))

        if sum(lowvals) < threshold:
            # The sum of all the expected frequencies that are less than
            # the threshold is still less than the threshold, so we have
            # to include a point whose expected frequency is above the
            # threshold.
            firstlow -= 1
            lowvals = values[firstlow:]
            threshold = values[firstlow]

        start = 0
        end = len(lowvals)

        groups = []
        while start < end:
            v0 = lowvals[start]
            group = [start]
            start += 1
            s = v0
            while s < threshold:
                end -= 1
                v1 = lowvals[end]
                group.append(end)
                s += v1
            if lowvals[start:end].sum() < threshold:
                group.extend(range(start, end))
                end = start
            groups.append(group)

        print("number of groups:", len(groups))
        print("number of bins:", len(expected) - len(lowvals) + len(groups))

        m = firstlow + len(groups)

        indexmap = {}
        for index, key in enumerate(keys[:firstlow]):
            indexmap[key] = index

        for k, g in enumerate(groups):
            for lowindex in g:
                keysindex = lowindex + firstlow
                key = keys[keysindex]
                indexmap[key] = k + firstlow

    else:
        # There are no expected frequencies below the threshold
        m = len(keys)

        indexmap = {}
        for index, key in enumerate(keys):
            indexmap[key] = index

    binned_expected = np.zeros(m)
    for t, v in zip(keys, values):
        k = indexmap[t]
        binned_expected[k] += v

    return indexmap, binned_expected
