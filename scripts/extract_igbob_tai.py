#!/usr/bin/python3
import json
import numpy as np
import scipy.stats

from tai_makespan_arpd_ref import tai_makespan_arpd_ref


res_dirname = "../experiments/igbob_tai__e90dd172980af75bfb58030929e31760add9cba0/"

TAI_TYPES = [
    (20, 5),
    (20, 10),
    (20, 20),
    (50, 5),
    (50, 10),
    (50, 20),
    (100, 5),
    (100, 10),
    (100, 20),
    (200, 10),
    (200, 20),
    (500, 20)
]

I = range(0,10)

TIMES = [30, 60, 90]


def extract_igbob_tai():
    res = {}
    for (n,m) in TAI_TYPES:
        for i in I:
            for t in TIMES:
                with open("{}perfprofile_igbob_{}_{}_{}_{}".format(res_dirname, t, n, m, i)) as f:
                    data = json.load(f)
                    res[(n,m,i,t)] = float(data["avgmakespan"])
                # res[(n,m,i,t)] = tmp
    return res

def compute_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

def compute_igbob_tai_arpd():
    res = extract_igbob_tai()
    tai_ref = tai_makespan_arpd_ref()
    aggregated = {}
    arpd = {}
    for n,m in TAI_TYPES:
        for t in TIMES:
            aggregated[(n,m,t)] = compute_confidence_interval([res[(n,m,i,t)] for i in I])
            arpd[(n,m,t)] = [
                (e-tai_ref[(n,m)])/(tai_ref[(n,m)])*100
                for e in aggregated[(n,m,t)]
            ]
    return arpd

if __name__ == "__main__":

    arpd = compute_igbob_tai_arpd()
    for e in arpd:
        print("{}\t{}".format(e, arpd[e]))