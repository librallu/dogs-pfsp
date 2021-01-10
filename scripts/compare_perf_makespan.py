#!/usr/bin/python3
from extract_igbob_vfr import extract_igbob_vfr, igbob_confidence_interval
from extract_igbob_tai import extract_igbob_tai
from helper_stats import read_perfprofile_file, best_known_at_time, apply_arpd_tab
from vfr_arpd_ref import vfr_arpd_refs
from tai_makespan_arpd_ref import tai_ub
from extract_igirms_arpd_ref import average_igirms_arpd_tai, average_igirms_arpd_vfr
from vbih_extract_sols import extract_vbih_sols

from scipy.stats import wilcoxon

ratio_cpu_ibs_igrms = 2584/737 # ~ 3.5
ratio_cpu_ibs_vbih = 2584/2047 # ~ 1.26

VFR_I = list(range(1,11))
TAI_I = list(range(0,10))

def get_igbob_results_vfr(n,m,t):
    return [ extract_igbob_vfr()[(n,m,i,t)] for i in VFR_I ]

def get_igbob_results_tai(n,m,t):
    return [ extract_igbob_tai()[(n,m,i,t)] for i in TAI_I ]



def get_ibs_results_vfr(n,m,t,algoname):
    dirname = "../experiments/EXPERIMENTS_IBS_VFR/"
    filenames = map(
        (lambda i: "{}IBS_VFR_{}_{}_{}_{}.perfprofile.json".format(dirname,algoname,n,m,i)), 
        VFR_I
    )
    points_list = map(lambda e:read_perfprofile_file(e), filenames)
    return list(map(lambda e:best_known_at_time(e,t), points_list))

def get_ibs_results_tai(n,m,t,algoname):
    dirname = "../experiments/IBS_makespan_tai__e1ee51bac12047452ba609a589ef85a64d235cab/"
    filenames = map(
        (lambda i: "{}IBS_{}_{}_{}_{}.perfprofile.json".format(dirname,algoname,n,m,i)), 
        TAI_I
    )
    points_list = map(lambda e:read_perfprofile_file(e), filenames)
    return list(map(lambda e:best_known_at_time(e,t), points_list))


def p_value_tab(ref,res1,res2):
    arpd_1 = sum(apply_arpd_tab(res1, ref))/10.
    arpd_2 =   sum(apply_arpd_tab(res2, ref))/10.
    res = {'arpd_1': arpd_1,'arpd_2': arpd_2}
    to_continue = False
    for i in range(len(res1)):
        if res1[i] != res2[i]:
            to_continue = True
            break
    if not to_continue:
        res['p'] = 1
    else:
        _,p = wilcoxon([res2[i]-res1[i] for i in range(10)])
        res['p'] = p
    return res




def p_value_igbob_vfr(n,m,t,algoname):
    igbob = get_igbob_results_vfr(n,m,t)
    ibs_t = n*m*t/2000
    ibs = get_ibs_results_vfr(n, m, ibs_t, algoname)
    ref = vfr_arpd_refs(n,m)
    tmp = p_value_tab(ref,igbob,ibs)
    return {
        'p': tmp["p"], 
        'arpd_igbob': tmp["arpd_1"],
        'arpd_ibs': tmp["arpd_2"],
    }

def p_value_igbob_tai(n,m,t,algoname):
    igbob = get_igbob_results_tai(n,m,t)
    ibs_t = n*m*t/2000
    ibs = get_ibs_results_tai(n, m, ibs_t, algoname)
    ref = tai_ub[(n,m)]
    tmp = p_value_tab(ref,igbob,ibs)
    return {
        'p': tmp["p"], 
        'arpd_igbob': tmp["arpd_1"],
        'arpd_ibs': tmp["arpd_2"],
    }


def humanize_arpd(v, isbold=False):
    v_text = "{:0.2f}".format(v)
    if isbold:
        return "\\textbf{"+v_text+"}"
    else:
        return v_text

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


    print("##### VBIH COMPARISON")
    vbih30 = extract_vbih_sols("./input/Permutations_VBIH_15nm.txt")
    vbih60 = extract_vbih_sols("./input/Permutations_VBIH_30nm.txt")
    vbih90 = extract_vbih_sols("./input/Permutations_VBIH_45nm.txt")
    for n in N:
        for m in M:
            ref = vfr_arpd_refs(n,m)
            vbih30_l = [ vbih30[(n,m,i)] for i in range(1,11) ]
            vbih60_l = [ vbih60[(n,m,i)] for i in range(1,11) ]
            vbih90_l = [ vbih90[(n,m,i)] for i in range(1,11) ]
            ibsw_30_l = get_ibs_results_vfr(n, m, n*m*30/2000/ratio_cpu_ibs_vbih, "bi_min_wfrontalpha")
            ibsw_60_l = get_ibs_results_vfr(n, m, n*m*60/2000/ratio_cpu_ibs_vbih, "bi_min_wfrontalpha")
            ibsw_90_l = get_ibs_results_vfr(n, m, n*m*90/2000/ratio_cpu_ibs_vbih, "bi_min_wfrontalpha")
            ibsa_30_l = get_ibs_results_vfr(n, m, n*m*30/2000/ratio_cpu_ibs_vbih, "bi_min_gap")
            ibsa_60_l = get_ibs_results_vfr(n, m, n*m*60/2000/ratio_cpu_ibs_vbih, "bi_min_gap")
            ibsa_90_l = get_ibs_results_vfr(n, m, n*m*90/2000/ratio_cpu_ibs_vbih, "bi_min_gap")
            pw_30 = p_value_tab(ref, vbih30_l, ibsw_30_l)
            pw_60 = p_value_tab(ref, vbih60_l, ibsw_60_l)
            pw_90 = p_value_tab(ref, vbih90_l, ibsw_90_l)
            pa_30 = p_value_tab(ref, vbih30_l, ibsa_30_l)
            pa_60 = p_value_tab(ref, vbih60_l, ibsa_60_l)
            pa_90 = p_value_tab(ref, vbih90_l, ibsa_90_l)
            vbih_best_30 = pw_30["p"]<pvalue_limit and (pa_30["p"]) and pw_30["arpd_1"] < pw_30["arpd_2"] and pa_30["arpd_1"] < pa_30["arpd_2"]
            vbih_best_60 = pw_60["p"]<pvalue_limit and (pa_60["p"]) and pw_60["arpd_1"] < pw_60["arpd_2"] and pa_60["arpd_1"] < pa_60["arpd_2"]
            vbih_best_90 = pw_90["p"]<pvalue_limit and (pa_90["p"]) and pw_90["arpd_1"] < pw_90["arpd_2"] and pa_90["arpd_1"] < pa_90["arpd_2"]
            print("{}\t & {} & {} & {}  & {} & {} & {}  & {} & {} & {} \\\\".format(
                "VFR{}\_{}".format(n,m),
                humanize_arpd(pw_30["arpd_1"],isbold=vbih_best_30), 
                humanize_arpd(pw_30["arpd_2"],isbold=pw_30["arpd_2"]<pw_30["arpd_1"] and pw_30["p"]<pvalue_limit),
                humanize_arpd(pa_30["arpd_2"],isbold=pa_30["arpd_2"]<pa_30["arpd_1"] and pa_30["p"]<pvalue_limit),

                humanize_arpd(pw_60["arpd_1"],isbold=vbih_best_60), 
                humanize_arpd(pw_60["arpd_2"],isbold=pw_60["arpd_2"]<pw_60["arpd_1"] and pw_60["p"]<pvalue_limit),
                humanize_arpd(pa_60["arpd_2"],isbold=pa_60["arpd_2"]<pa_60["arpd_1"] and pa_60["p"]<pvalue_limit),

                humanize_arpd(pw_90["arpd_1"],isbold=vbih_best_90), 
                humanize_arpd(pw_90["arpd_2"],isbold=pw_90["arpd_2"]<pw_90["arpd_1"] and pw_90["p"]<pvalue_limit),
                humanize_arpd(pa_90["arpd_2"],isbold=pa_90["arpd_2"]<pa_90["arpd_1"] and pa_90["p"]<pvalue_limit),
            ))



    print("##### IGRMS COMPARISON")
    avg_igirms = average_igirms_arpd_tai()
    for n,m in TAI_TYPES:
        igirms_60 = avg_igirms[(n,m)][0]
        igirms_120 = avg_igirms[(n,m)][1]
        igirms_240 = avg_igirms[(n,m)][2]
        ibsw_60 = sum(apply_arpd_tab(
            get_ibs_results_tai(n,m,60*n*m/ratio_cpu_ibs_igrms/2000,"bi_min_wfrontalpha"),
            tai_ub[(n,m)]
        ))/10
        ibsw_120 = sum(apply_arpd_tab(
            get_ibs_results_tai(n,m,120*n*m/ratio_cpu_ibs_igrms/2000,"bi_min_wfrontalpha"),
            tai_ub[(n,m)]
        ))/10
        ibsw_240 = sum(apply_arpd_tab(
            get_ibs_results_tai(n,m,240*n*m/ratio_cpu_ibs_igrms/2000,"bi_min_wfrontalpha"),
            tai_ub[(n,m)]
        ))/10
        ibsg_60 = sum(apply_arpd_tab(
            get_ibs_results_tai(n,m,60*n*m/ratio_cpu_ibs_igrms/2000,"bi_min_gap"),
            tai_ub[(n,m)]
        ))/10
        ibsg_120 = sum(apply_arpd_tab(
            get_ibs_results_tai(n,m,120*n*m/ratio_cpu_ibs_igrms/2000,"bi_min_gap"),
            tai_ub[(n,m)]
        ))/10
        ibsg_240 = sum(apply_arpd_tab(
            get_ibs_results_tai(n,m,240*n*m/ratio_cpu_ibs_igrms/2000,"bi_min_gap"),
            tai_ub[(n,m)]
        ))/10
        print("{}\t & {:0.2f} & {:0.2f} & {:0.2f}  & {:0.2f} & {:0.2f} & {:0.2f}  & {:0.2f} & {:0.2f} & {:0.2f} \\\\".format(
            "TAI{}\_{}".format(n,m),
            igirms_60, ibsw_60, ibsg_60,
            igirms_120, ibsw_120, ibsg_120,
            igirms_240, ibsw_240, ibsg_240,
        ))
    print("\\hline")
    avg_igirms_vfr = average_igirms_arpd_vfr()
    for n in N:
        for m in M:
            ref = vfr_arpd_refs(n,m)
            igirms_60 =  avg_igirms_vfr[(n,m)][0]
            igirms_120 = avg_igirms_vfr[(n,m)][1]
            igirms_240 = avg_igirms_vfr[(n,m)][2]
            ibsw_60 = sum(apply_arpd_tab(
                get_ibs_results_vfr(n,m,60*n*m/ratio_cpu_ibs_igrms/2000,"bi_min_wfrontalpha"),
                ref
            ))/10
            ibsw_120 = sum(apply_arpd_tab(
                get_ibs_results_vfr(n,m,120*n*m/ratio_cpu_ibs_igrms/2000,"bi_min_wfrontalpha"),
                ref
            ))/10
            ibsw_240 = sum(apply_arpd_tab(
                get_ibs_results_vfr(n,m,240*n*m/ratio_cpu_ibs_igrms/2000,"bi_min_wfrontalpha"),
                ref
            ))/10
            ibsg_60 = sum(apply_arpd_tab(
                get_ibs_results_vfr(n,m,60*n*m/ratio_cpu_ibs_igrms/2000,"bi_min_gap"),
                ref
            ))/10
            ibsg_120 = sum(apply_arpd_tab(
                get_ibs_results_vfr(n,m,120*n*m/ratio_cpu_ibs_igrms/2000,"bi_min_gap"),
                ref
            ))/10
            ibsg_240 = sum(apply_arpd_tab(
                get_ibs_results_vfr(n,m,240*n*m/ratio_cpu_ibs_igrms/2000,"bi_min_gap"),
                ref
            ))/10
            print("{}\t & {:0.2f} & {:0.2f} & {:0.2f}  & {:0.2f} & {:0.2f} & {:0.2f}  & {:0.2f} & {:0.2f} & {:0.2f} \\\\".format(
            "VFR{}\_{}".format(n,m),
            igirms_60, ibsw_60, ibsg_60,
            igirms_120, ibsw_120, ibsg_120,
            igirms_240, ibsw_240, ibsg_240,
        ))

    print("##### IGBOB COMPARISON")
    for n,m in TAI_TYPES:
        p30w = p_value_igbob_tai(n,m,30,"bi_min_wfrontalpha")
        p30g = p_value_igbob_tai(n,m,30,"bi_min_gap")
        p60w = p_value_igbob_tai(n,m,60,"bi_min_wfrontalpha")
        p60g = p_value_igbob_tai(n,m,60,"bi_min_gap")
        p90w = p_value_igbob_tai(n,m,90,"bi_min_wfrontalpha")
        p90g = p_value_igbob_tai(n,m,90,"bi_min_gap")
        igbob_is_best30 = (p30w["p"]<pvalue_limit) and (p30g["p"]<pvalue_limit) and (p30w["arpd_igbob"]<p30w["arpd_ibs"]) and (p30g["arpd_igbob"]<p30g["arpd_ibs"])
        igbob_is_best60 = (p60w["p"]<pvalue_limit) and (p60g["p"]<pvalue_limit) and (p60w["arpd_igbob"]<p60w["arpd_ibs"]) and (p60g["arpd_igbob"]<p60g["arpd_ibs"])
        igbob_is_best90 = (p90w["p"]<pvalue_limit) and (p90g["p"]<pvalue_limit) and (p90w["arpd_igbob"]<p90w["arpd_ibs"]) and (p90g["arpd_igbob"]<p90g["arpd_ibs"])
        print("{}\t & {} & {} & {}  & {} & {} & {}  & {} & {} & {} \\\\".format(
            "TAI{}\_{}".format(n,m),
            humanize_arpd(p30w["arpd_igbob"],isbold=igbob_is_best30), 
            humanize_arpd(p30w["arpd_ibs"],isbold=p30w["arpd_ibs"]<p30w["arpd_igbob"] and p30w["p"]<pvalue_limit),
            humanize_arpd(p30g["arpd_ibs"],isbold=p30g["arpd_ibs"]<p30g["arpd_igbob"] and p30g["p"]<pvalue_limit),
            humanize_arpd(p60w["arpd_igbob"],isbold=igbob_is_best60), 
            humanize_arpd(p60w["arpd_ibs"],isbold=p60w["arpd_ibs"]<p60w["arpd_igbob"] and p60w["p"]<pvalue_limit),
            humanize_arpd(p60g["arpd_ibs"],isbold=p60g["arpd_ibs"]<p60g["arpd_igbob"] and p60g["p"]<pvalue_limit),
            humanize_arpd(p90w["arpd_igbob"],isbold=igbob_is_best90), 
            humanize_arpd(p90w["arpd_ibs"],isbold=p90w["arpd_ibs"]<p90w["arpd_igbob"] and p90w["p"]<pvalue_limit),
            humanize_arpd(p90g["arpd_ibs"],isbold=p90g["arpd_ibs"]<p90g["arpd_igbob"] and p90g["p"]<pvalue_limit),
        ))
    print("\\hline")
    # VFR results
    for n in N:
        for m in M:
            p30w = p_value_igbob_vfr(n,m,30,"bi_min_wfrontalpha")
            p30g = p_value_igbob_vfr(n,m,30,"bi_min_gap")
            p60w = p_value_igbob_vfr(n,m,60,"bi_min_wfrontalpha")
            p60g = p_value_igbob_vfr(n,m,60,"bi_min_gap")
            p90w = p_value_igbob_vfr(n,m,90,"bi_min_wfrontalpha")
            p90g = p_value_igbob_vfr(n,m,90,"bi_min_gap")
            igbob_is_best30 = (p30w["p"]<pvalue_limit) and (p30g["p"]<pvalue_limit) and (p30w["arpd_igbob"]<p30w["arpd_ibs"]) and (p30g["arpd_igbob"]<p30g["arpd_ibs"])
            igbob_is_best60 = (p60w["p"]<pvalue_limit) and (p60g["p"]<pvalue_limit) and (p60w["arpd_igbob"]<p60w["arpd_ibs"]) and (p60g["arpd_igbob"]<p60g["arpd_ibs"])
            igbob_is_best90 = (p90w["p"]<pvalue_limit) and (p90g["p"]<pvalue_limit) and (p90w["arpd_igbob"]<p90w["arpd_ibs"]) and (p90g["arpd_igbob"]<p90g["arpd_ibs"])
            print("{}\t & {} & {} & {}  & {} & {} & {}  & {} & {} & {} \\\\".format(
                "VFR{}\_{}".format(n,m),
                humanize_arpd(p30w["arpd_igbob"],isbold=igbob_is_best30), 
                humanize_arpd(p30w["arpd_ibs"],isbold=p30w["arpd_ibs"]<p30w["arpd_igbob"] and p30w["p"]<pvalue_limit),
                humanize_arpd(p30g["arpd_ibs"],isbold=p30g["arpd_ibs"]<p30g["arpd_igbob"] and p30g["p"]<pvalue_limit),

                humanize_arpd(p60w["arpd_igbob"],isbold=igbob_is_best60), 
                humanize_arpd(p60w["arpd_ibs"],isbold=p60w["arpd_ibs"]<p60w["arpd_igbob"] and p60w["p"]<pvalue_limit),
                humanize_arpd(p60g["arpd_ibs"],isbold=p60g["arpd_ibs"]<p60g["arpd_igbob"] and p60g["p"]<pvalue_limit),

                humanize_arpd(p90w["arpd_igbob"],isbold=igbob_is_best90), 
                humanize_arpd(p90w["arpd_ibs"],isbold=p90w["arpd_ibs"]<p90w["arpd_igbob"] and p90w["p"]<pvalue_limit),
                humanize_arpd(p90g["arpd_ibs"],isbold=p90g["arpd_ibs"]<p90g["arpd_igbob"] and p90g["p"]<pvalue_limit),
            ))
