#!/usr/bin/python3
import csv
import matplotlib.pyplot as plt

from helper_stats import points_from_perfprofile, aggregate_points, apply_arpd_tab, get_confidence_curves
from vfr_arpd_ref import vfr_arpd_refs
from extract_igbob_vfr import compute_igbob_vfr_arpd, igbob_confidence_interval
from vbih_extract_sols import vbih_arpd, vbih_confidence_interval

def get_ibs_vfr_confidence(algoname,n,m,t_end=0):
    I = list(range(1,11))
    filenames = map(
        (lambda i: "../experiments/EXPERIMENTS_IBS_VFR/IBS_VFR_{}_{}_{}_{}.perfprofile.json".format(algoname,n,m,i)), 
        I
    )
    points_list = map(lambda e:points_from_perfprofile(e,t_end), filenames)
    aggregated_points = aggregate_points(list(points_list))
    ref = vfr_arpd_refs(n,m)
    aggregated_arpd = {}
    for k in aggregated_points:
        aggregated_arpd[k] = apply_arpd_tab(aggregated_points[k], ref)
    return get_confidence_curves(aggregated_arpd)

def read_igrms_arpd():
    with open("input/igirms_vfr.csv") as csvfile:
        res = {}
        reader = csv.reader(csvfile)
        for row in reader:
            if row[0] == "instance_class": # skip header
                pass
            else:
                res[tuple([ int(e) for e in row[0][3:].split("_")])] = [float(row[1]), float(row[2]), float(row[3])]
        return res


def read_igall_arpd():
    with open("input/iga_vfr.csv") as csvfile:
        res = {}
        reader = csv.reader(csvfile)
        for row in reader:
            if row[0] == "instance_class": # skip header
                pass
            else:
                res[tuple([ int(e) for e in row[0][3:].split("_")])] = [float(row[1]), float(row[2]), float(row[3])]
        return res


def gen_pareto_tai_makespan(n,m,g):
    ratio_cpu_ibs_igrms = 2584/737 # ~ 3.5
    
    # PRINT IGBOB
    igbob_30 = igbob_confidence_interval(n,m,30)
    igbob_60 = igbob_confidence_interval(n,m,60)
    igbob_90 = igbob_confidence_interval(n,m,90)

    data = [
        [30/2*n*m/1000, float(igbob_30[1]), "IGbob", igbob_30[2]-igbob_30[1]],
        [60/2*n*m/1000, float(igbob_60[1]), "IGbob", igbob_60[2]-igbob_60[1]],
        [90/2*n*m/1000, float(igbob_90[1]), "IGbob", igbob_90[2]-igbob_90[1]],
    ]
    x,y,t,e = list(zip(*data))
    # plt.errorbar(x, y, yerr=e, fmt="x", capsize=3, c="g")
    plt.plot(x,y,"+")
    for i,text in enumerate(t):
        plt.annotate(text, (x[i],y[i]))

    # PRINT IGRMS
    igrms_arpd = read_igrms_arpd()
    igrms_60 = igrms_arpd[(n,m)][0]
    igrms_120 = igrms_arpd[(n,m)][1]
    igrms_240 = igrms_arpd[(n,m)][2]
    # print(igrms_60, igrms_120, igrms_240)
    data = [
        [60/2/ratio_cpu_ibs_igrms*n*m/1000,  float(igrms_60), "IGrms"],
        [120/2/ratio_cpu_ibs_igrms*n*m/1000, float(igrms_120), "IGrms"],
        [240/2/ratio_cpu_ibs_igrms*n*m/1000, float(igrms_240), "IGrms"],
    ]
    x,y,t = list(zip(*data))
    plt.plot(x,y,"+")
    for i,text in enumerate(t):
        plt.annotate(text, (x[i],y[i]-0.01))

    # PRINT VBIH
    ratio_cpu_ibs_vbih = 2584/2047 # ~ 1.26
    vbih_15 = vbih_confidence_interval("./input/Permutations_VBIH_15nm.txt", n, m)
    vbih_30 = vbih_confidence_interval("./input/Permutations_VBIH_30nm.txt", n, m)
    vbih_45 = vbih_confidence_interval("./input/Permutations_VBIH_45nm.txt", n, m)
    print(vbih_15)
    data = [
        [15/ratio_cpu_ibs_vbih*n*m/1000, vbih_15[1], "VBIH", vbih_15[2]-vbih_15[1]],
        [30/ratio_cpu_ibs_vbih*n*m/1000, vbih_30[1], "VBIH", vbih_30[2]-vbih_30[1]],
        [45/ratio_cpu_ibs_vbih*n*m/1000, vbih_45[1], "VBIH", vbih_45[2]-vbih_45[1]],
    ]
    x,y,t,e = list(zip(*data))
    # plt.errorbar(x, y, yerr=e, fmt="x", capsize=3, c="r")
    plt.plot(x,y,"+")
    for i,text in enumerate(t):
        plt.annotate(text, (x[i]+0.01,y[i]-.01))

    # PRINT IBS
    if g == 'gap' or g == 'all':
        conf_curve = get_ibs_vfr_confidence("bi_min_gap", n, m, t_end=int(n*m*0.045))
        ibs_a_t = []
        ibs_a_vmin = []
        ibs_a_v = []
        ibs_a_vmax = []
        for t in conf_curve:
            if t > 10:
                ibs_a_t.append(t)
                ibs_a_vmin.append(conf_curve[t][0])
                ibs_a_v.append(conf_curve[t][1])
                ibs_a_vmax.append(conf_curve[t][2])
        plt.plot(ibs_a_t, ibs_a_v, "b--", label="IBS bi_gap", drawstyle='steps-post')
        # plt.fill_between(ibs_a_t, ibs_a_vmin, ibs_a_vmax, alpha=0.2, color="b")
    if g == 'wfrontalpha' or g == 'all':
        conf_curve = get_ibs_vfr_confidence("bi_min_wfrontalpha", n, m, t_end=int(n*m*0.045))
        ibs_a_t = []
        ibs_a_vmin = []
        ibs_a_v = []
        ibs_a_vmax = []
        for t in conf_curve:
            if t > 10:
                ibs_a_t.append(t)
                ibs_a_vmin.append(conf_curve[t][0])
                ibs_a_v.append(conf_curve[t][1])
                ibs_a_vmax.append(conf_curve[t][2])
        plt.plot(ibs_a_t, ibs_a_v, "b-", label="IBS bi_wfrontalpha", drawstyle='steps-post')
        # plt.fill_between(ibs_a_t, ibs_a_vmin, ibs_a_vmax, alpha=0.2, color="b")

    # draw 0 line
    # plt.plot(0,0,max(ibs_a_t),0, "-")
    plt.grid()


    # DISPLAY PLOT
    plt.xlabel("CPU-regularized running time (s)")
    plt.ylabel("Average relative percentage deviation (ARPD, smaller is better)")
    # plt.yscale("symlog")
    plt.legend()
    plt.show()

from sys import argv
if __name__ == "__main__":
    print("USAGE: {} N M IBS_guide".format(argv[0]))
    print("\tIBS_guide (gap, wfrontalpha, all)")
    gen_pareto_tai_makespan(int(argv[1]),int(argv[2]),argv[3])
