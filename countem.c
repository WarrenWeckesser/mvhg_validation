
#include <stdio.h>
#include <inttypes.h>
#include "khash.h"


KHASH_MAP_INIT_INT64(intmap, uint64_t)


/*
 *  indexmap is (in effect) a 2-d array with shape (num_support, ncolors + 1).
 *  The first ncolors columns are elements of the support of the distribution.
 *  The last column is the index into counts to which that element of the
 *  support contributes.  The maximum value in the last column must be
 *  less than nbins.
 *
 *  samples is the 2-d array of samples from the distribution.  It has shape
 *  (num_samples, ncolors).
 *
 *  counts is a 1-d array of the counts of the occurrences of the elements
 *  in samples.  It has length nbins.  counts is *not* initialized to zero
 *  in the function, so that it can be called more than once to accumulate
 *  the results for several samples.  Be sure to initialize the values in
 *  counts to 0 before making the first call to countem().
 */

int countem(int num_support, int ncolors, int nsample, int nbins, int64_t *indexmap,
            int num_samples, int64_t *samples, int64_t *counts)
{
    int i, j;
    khiter_t k;
    int status;
    khash_t(intmap) *h = kh_init(intmap);

    /*
    printf("countem: num_support = %d\n", num_support);
    printf("countem: ncolors     = %d\n", ncolors);
    printf("countem: nsample     = %d\n", nsample);
    printf("countem: num_samples = %d\n", num_samples);
    */

    /*
    printf("indexmap[:2,:] =");
    for (i = 0; i < ncolors + 1; ++i) {
        printf("  %lld", indexmap[i]);
    }
    printf("\n                ");
    for (i = 0; i < ncolors + 1; ++i) {
        printf("  %lld", indexmap[i + ncolors + 1]);
    }
    printf("\n");
    */

    for (i = 0; i < num_support; ++i) {
        int64_t key = 0;
        for (j = 0; j < (ncolors - 1); ++j) {
            key = key*(nsample + 1) + *indexmap;
            ++indexmap;
        }
        ++indexmap;  // Skip the last component of the current sample.

        k = kh_put(intmap, h, key, &status);
        kh_value(h, k) = *indexmap;
        //printf("key = %lld   value = %lld\n", key, *indexmap);
        ++indexmap;
    }

    for (i = 0; i < num_samples; ++i) {
        int index;
        int64_t key = 0;
        for (j = 0; j < (ncolors - 1); ++j) {
            //printf("%lld ", *samples);
            key = key*(nsample + 1) + *samples;
            ++samples;
        }
        ++samples;  // Skip the last component of the current sample.
        //printf("   key = %lld", key);

        k = kh_get(intmap, h, key);
        index = kh_value(h, k);
        //printf(" value = %d\n", index);
        counts[index] += 1;

    }

    kh_destroy(intmap, h);

    return 0; 
}


#ifdef MAIN

int test_countem()
{
    int i;

    int num_support = 10;
    int ncolors = 3;
    int nsample = 3;
    int nbins = 8;
    int indexmap[] = {3, 0, 0,  0,
                      2, 1, 0,  1,
                      2, 0, 1,  2,
                      1, 2, 0,  3,
                      1, 1, 1,  4,
                      1, 0, 2,  5,
                      0, 3, 0,  0,
                      0, 2, 1,  6,
                      0, 1, 2,  7,
                      0, 0, 3,  0};

    int num_samples = 10;
    int samples[] = {3, 0, 0,
                     0, 0, 3,
                     0, 2, 1,
                     2, 1, 0,
                     0, 1, 2,
                     0, 2, 1,
                     2, 0, 1,
                     0, 0, 3,
                     2, 0, 1,
                     2, 1, 0};
    int counts[] = {0, 0, 0, 0, 0, 0, 0, 0};
    int expected_counts[] = {3, 2, 2, 0, 0, 0, 2, 1};

    countem(num_support, ncolors, nsample, nbins, indexmap, num_samples, samples, counts);

    for (i = 0; i < nbins; ++i) {
        printf("counts[%d] = %d   expected: %d\n", i, counts[i], expected_counts[i]);
    }
    return 0;
}

int main(int argc, char *argv[])
{
    test_countem();

    return 0;
}
#endif
