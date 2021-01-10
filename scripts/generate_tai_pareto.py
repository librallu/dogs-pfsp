#!/usr/bin/python3
import csv
import matplotlib.pyplot as plt

from extract_igbob_tai import compute_igbob_tai_arpd
from extract_arpd_perfprofile_makespan_tai import gen_arpd_perfprofile_makespan

def read_igrms_arpd():
    with open("input/igirms_tai.csv") as csvfile:
        res = {}
        reader = csv.reader(csvfile)
        for row in reader:
            if row[0] == "instance_class": # skip header
                pass
            else:
                res[tuple([ int(e) for e in row[0].split("_")])] = [float(row[1]), float(row[2]), float(row[3])]
        return res

def gen_pareto_tai_makespan(n,m):
    ratio_cpu_ibs_igrms = 2584/737 # ~ 3.5
    arpd_igbob = compute_igbob_tai_arpd()
    igbob_30 = arpd_igbob[(n,m,30)]
    igbob_60 = arpd_igbob[(n,m,60)]
    igbob_90 = arpd_igbob[(n,m,90)]
    igrms_arpd = read_igrms_arpd()
    igrms_60 = igrms_arpd[(n,m)][0]
    igrms_120 = igrms_arpd[(n,m)][1]
    igrms_240 = igrms_arpd[(n,m)][2]
    # print(igrms_60, igrms_120, igrms_240)
    
    # PRINT IGBOB
    data = [
        [30/2*n*m/1000, float(igbob_30[0]), "IGbob"],
        [60/2*n*m/1000, float(igbob_60[0]), "IGbob"],
        [90/2*n*m/1000, float(igbob_90[0]), "IGbob"],
    ]
    x,y,t = list(zip(*data))
    plt.plot(x,y,"+")
    for i,text in enumerate(t):
        plt.annotate(text, (x[i],y[i]))
    print(data)

    # PRINT IGRMS
    data = [
        [60/2/ratio_cpu_ibs_igrms*n*m/1000,  float(igrms_60), "IGrms"],
        [120/2/ratio_cpu_ibs_igrms*n*m/1000, float(igrms_120), "IGrms"],
        [240/2/ratio_cpu_ibs_igrms*n*m/1000, float(igrms_240), "IGrms"],
    ]
    x,y,t = list(zip(*data))
    plt.plot(x,y,"x")
    for i,text in enumerate(t):
        plt.annotate(text, (x[i],y[i]))

    # PRINT IBS
    ibs_a_perf = gen_arpd_perfprofile_makespan(n,m,"bi_min_gap")
    ibs_a_t = []
    ibs_a_v = []
    for i in range(1,len(ibs_a_perf)):
        ibs_a_t.append(i)
        ibs_a_v.append(ibs_a_perf[i])
    plt.plot(ibs_a_t, ibs_a_v, "g--", label="bi-min_g5")

    ibs_a_perf = gen_arpd_perfprofile_makespan(n,m,"bi_min_bound")
    ibs_a_t = []
    ibs_a_v = []
    for i in range(1,len(ibs_a_perf)):
        ibs_a_t.append(i)
        ibs_a_v.append(ibs_a_perf[i])
    plt.plot(ibs_a_t, ibs_a_v, "r-.", label="bi-min_g1")

    ibs_a_perf = gen_arpd_perfprofile_makespan(n,m,"bi_min_wfrontalpha")
    ibs_a_t = []
    ibs_a_v = []
    for i in range(1,len(ibs_a_perf)):
        ibs_a_t.append(i)
        ibs_a_v.append(ibs_a_perf[i])
    plt.plot(ibs_a_t, ibs_a_v, "b-", label="bi-min_g4")


    # DISPLAY PLOT
    plt.xlabel("regularized running time (s)")
    plt.ylabel("Average relative percentage deviation (ARPD, smaller is better)")
    plt.yscale("log")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    gen_pareto_tai_makespan(500,20)
