#!/usr/bin/python3
import json
import numpy as np
import scipy.stats

from vfr_arpd_ref import vfr_arpd_ref, vfr_arpd_refs
from helper_stats import compute_confidence_interval, apply_arpd_tab


res_dirname = "../experiments/igbob_vrf__2838421f1bb08a8464e6118b7c7dcb848f4827d4/"

N = [100,200,300,400,500,600,700,800]
M = [20,40,60]
I = range(1,11)

TIMES = [30, 60, 90]


def extract_igbob_vfr():
    res = {}
    for n in N:
        for m in M:
            for i in I:
                for t in TIMES:
                    with open("{}perfprofile_igbob_{}_{}_{}_{}".format(res_dirname, t, n, m, i)) as f:
                        data = json.load(f)
                        res[(n,m,i,t)] = float(data["avgmakespan"])
                    # res[(n,m,i,t)] = tmp
    return res


def compute_igbob_vfr_arpd():
    res = extract_igbob_vfr()
    ref = vfr_arpd_ref()
    aggregated = {}
    arpd = {}
    for n in N:
        for m in M:
            for t in TIMES:
                aggregated[(n,m,t)] = compute_confidence_interval([res[(n,m,i,t)] for i in I])
                arpd[(n,m,t)] = [
                    (e-ref[(n,m)])/(ref[(n,m)])*100
                    for e in aggregated[(n,m,t)]
                ]
    return arpd


def igbob_confidence_interval(n, m, t):
    res = extract_igbob_vfr()
    ref_tab = vfr_arpd_refs(n,m)
    tab = [res[(n,m,i,t)] for i in I]
    print(tab)
    arpd_tab = apply_arpd_tab(tab, ref_tab)
    print(arpd_tab)
    print(compute_confidence_interval(arpd_tab))
    return compute_confidence_interval(arpd_tab)


if __name__ == "__main__":

    arpd = compute_igbob_vfr_arpd()
    for e in arpd:
        print("{}\t{}".format(e, arpd[e][0]))