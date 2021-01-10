#!/usr/bin/python3
from compare_perf_makespan import p_value_tab, get_ibs_results_vfr, get_ibs_results_tai, humanize_arpd
from tai_makespan_arpd_ref import tai_ub
from vfr_arpd_ref import vfr_arpd_refs


def is_first_better(arpd_obj, pvalue_limit=0.05):
    return arpd_obj["arpd_1"] < arpd_obj["arpd_2"] and arpd_obj["p"] < pvalue_limit

if __name__ == "__main__":
    pvalue_limit = 0.05
    T = [30,60,90]
    # taillard results
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
    N = [100,200,300,400,500,600,700,800]
    M = [20,40,60]
    branchings = ["forward", "bi_min"]
    guides = ["bound","idle","alpha","wfrontalpha","gap"]
    variants = []
    for b in branchings:
        for g in guides:
            s = "{}_{}".format(b,g)
            variants.append(s)
    # TAILLARD
    for (n,m) in TAI_TYPES:
        ref = tai_ub[(n,m)]
        ibs_lists = [get_ibs_results_tai(n,m,90*n*m/2000,v) for v in variants]
        ibs_arpds = [p_value_tab(ref,e,e)['arpd_1'] for e in ibs_lists]
        print("{}\t& {} \\\\".format(
            "TAI{}\_{}".format(n,m),
            " & ".join([ humanize_arpd(
                e,
                sum(
                    [int(is_first_better(p_value_tab(ref,ibs_lists[i],ibs_lists[j]))) for j in range(len(ibs_lists))]
                ) == len(ibs_lists)-1
            ) for i,e in enumerate(ibs_arpds) ])
        ))
    # VFR
    print("\\hline")
    for n in N:
        for m in M:
            ref = vfr_arpd_refs(n,m)
            ibs_lists = [get_ibs_results_vfr(n,m,90*n*m/2000,v) for v in variants]
            ibs_arpds = [p_value_tab(ref,e,e)['arpd_1'] for e in ibs_lists]
            print("{}\t& {} \\\\".format(
            "VFR{}\_{}".format(n,m),
            " & ".join([ humanize_arpd(
                e,
                sum(
                    [int(is_first_better(p_value_tab(ref,ibs_lists[i],ibs_lists[j]))) for j in range(len(ibs_lists))]
                ) == len(ibs_lists)-1
            ) for i,e in enumerate(ibs_arpds) ])
        ))