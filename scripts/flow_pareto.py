#!/usr/bin/python3
import csv
import matplotlib.pyplot as plt

from helper_stats import points_from_perfprofile, aggregate_points, apply_arpd_tab, get_confidence_curves
from flow_arpd_ref import flow_arpd_ref
from gen_bks_tabs_flow import read_bks_file

def get_ibs_flow_confidence(algoname,n,m,t_end=0):
    I = list(range(0,10))
    filenames = map(
        (lambda i: "../experiments/IBS_flowtime_tai__74aed6e9568574f41e722b3e3b98329846fd2f8f/IBS_{}_{}_{}_{}.perfprofile.json".format(algoname,n,m,i)), 
        I
    )
    points_list = map(lambda e:points_from_perfprofile(e,t_end), filenames)
    aggregated_points = aggregate_points(list(points_list))
    ref = flow_arpd_ref[(n,m)]
    aggregated_arpd = {}
    for k in aggregated_points:
        aggregated_arpd[k] = apply_arpd_tab(aggregated_points[k], ref)
    return get_confidence_curves(aggregated_arpd)

def read_igrms_arpd():
    with open("input/flow_algirtct.csv") as csvfile:
        res = {}
        reader = csv.reader(csvfile)
        for row in reader:
            if row[0] == "instance_class": # skip header
                pass
            else:
                res[tuple([ int(e) for e in row[0].split("_")])] = [float(row[1]), float(row[2]), float(row[3])]
        return res


def read_iga_arpd():
    with open("input/flow_iga.csv") as csvfile:
        res = {}
        reader = csv.reader(csvfile)
        for row in reader:
            if row[0] == "instance_class": # skip header
                pass
            else:
                res[tuple([ int(e) for e in row[0].split("_")])] = [float(row[1]), float(row[2]), float(row[3])]
        return res


def pareto_flowtime(n,m,g):
    ratio_cpu_ibs_igrms = 2584/737 # ~ 3.5

    # PRINT IGRMS
    igrms_arpd = read_igrms_arpd()
    igrms_60 = igrms_arpd[(n,m)][0]
    igrms_120 = igrms_arpd[(n,m)][1]
    igrms_240 = igrms_arpd[(n,m)][2]
    data = [
        [60/2/ratio_cpu_ibs_igrms*n*m/1000,  float(igrms_60),  "ALGirtct"],
        [120/2/ratio_cpu_ibs_igrms*n*m/1000, float(igrms_120), "ALGirtct"],
        [240/2/ratio_cpu_ibs_igrms*n*m/1000, float(igrms_240), "ALGirtct"],
    ]
    x,y,t = list(zip(*data))
    plt.plot(x,y,"x")
    for i,text in enumerate(t):
        plt.annotate(text, (x[i],y[i]-0.01))

    # PRINT IGA
    iga_arpd = read_iga_arpd()
    iga_60 =   iga_arpd[(n,m)][0]
    iga_120 =  iga_arpd[(n,m)][1]
    iga_240 =  iga_arpd[(n,m)][2]
    data = [
        [60/2/ratio_cpu_ibs_igrms*n*m/1000,  float(iga_60),  "IGA"],
        [120/2/ratio_cpu_ibs_igrms*n*m/1000, float(iga_120), "IGA"],
        [240/2/ratio_cpu_ibs_igrms*n*m/1000, float(iga_240), "IGA"],
    ]
    x,y,t = list(zip(*data))
    plt.plot(x,y,"x")
    for i,text in enumerate(t):
        plt.annotate(text, (x[i],y[i]))

    # PRINT BSCH(n)
    ratio_cpu_ibs_bsch = 2584/2070 # ~ 1.24
    bsch_arpd = {
        (20, 5):    1.25,
        (20, 10):   0.75,
        (20, 20):   0.75,
        (50, 5):    0.75,
        (50, 10):   1.04,
        (50, 20):   1.48,
        (100, 5):   0.30,
        (100, 10):  0.57,
        (100, 20):  1.14,
        (200, 10):  -0.61,
        (200, 20):  -0.76,
        (500, 20):  -1.87
    }

    bsch_time = {
        (20, 5):    0.00,
        (20, 10):   0.00,
        (20, 20):   0.00,
        (50, 5):    0.03,
        (50, 10):   0.03,
        (50, 20):   0.06,
        (100, 5):   0.31,
        (100, 10):  0.40,
        (100, 20):  0.68,
        (200, 10):  7.25,
        (200, 20):  8.57,
        (500, 20):  215.44
    }

    mrsils_arpd = {
        (20, 5):    0.01,
        (20, 10):   0.00,
        (20, 20):   0.00,
        (50, 5):    0.28,
        (50, 10):   0.47,
        (50, 20):   0.63,
        (100, 5):   0.22,
        (100, 10):  0.27,
        (100, 20):  0.83,
        (200, 10):  -0.71,
        (200, 20):  -0.83,
        (500, 20):  -1.90
    }

    mrsils_time = {
        (20, 5):    3.00,
        (20, 10):   6.00,
        (20, 20):   12.00,
        (50, 5):    7.50,
        (50, 10):   15.00,
        (50, 20):   30.00,
        (100, 5):   15.00,
        (100, 10):  30.00,
        (100, 20):  60.00,
        (200, 10):  60.00,
        (200, 20):  120.00,
        (500, 20):  300.00
    }
    plt.plot(bsch_time[(n,m)]/ratio_cpu_ibs_bsch,bsch_arpd[(n,m)],"x", color="r")
    plt.annotate("BSCH(n)", (bsch_time[(n,m)]/ratio_cpu_ibs_bsch,bsch_arpd[(n,m)]))
    
    plt.plot(mrsils_time[(n,m)]/ratio_cpu_ibs_bsch,mrsils_arpd[(n,m)],"x", color="r")
    plt.annotate("MRSILS", (mrsils_time[(n,m)]/ratio_cpu_ibs_bsch,mrsils_arpd[(n,m)]))

    # PRINT shake_ls
    # bks_flow = read_bks_file()
    # shake_ls = sum([bks_flow["TA{}_{}_{}".format(n,m,i)]["shake_ls"] for i in range(0,10)])/10
    # ref = sum(flow_arpd_ref[(n,m)])/10
    # shake_ls_arpd = (shake_ls-ref)/ref*100
    # print(shake_ls_arpd)
    # plt.plot(n*m*0.045,shake_ls_arpd,"x", color="r")
    # plt.annotate("shake_LS (1h)", (n*m*0.045,shake_ls_arpd))

    # PRINT IBS
    if g == 'alpha' or g == 'all':
        conf_curve = get_ibs_flow_confidence("forward_alpha", n, m, t_end=int(n*m*0.045))
        ibs_a_t = []
        ibs_a_vmin = []
        ibs_a_v = []
        ibs_a_vmax = []
        for t in conf_curve:
            if t > -1:
                ibs_a_t.append(t)
                ibs_a_vmin.append(conf_curve[t][0])
                ibs_a_v.append(conf_curve[t][1])
                ibs_a_vmax.append(conf_curve[t][2])
        plt.plot(ibs_a_t, ibs_a_v, "b--", label="IBS alpha", drawstyle='steps-post')
        # plt.fill_between(ibs_a_t, ibs_a_vmin, ibs_a_vmax, alpha=0.2, color="b")
    if g == 'walpha' or g == 'all':
        conf_curve = get_ibs_flow_confidence("forward_walpha", n, m, t_end=int(n*m*0.045))
        ibs_a_t = []
        ibs_a_vmin = []
        ibs_a_v = []
        ibs_a_vmax = []
        for t in conf_curve:
            if t > -1:
                ibs_a_t.append(t)
                ibs_a_vmin.append(conf_curve[t][0])
                ibs_a_v.append(conf_curve[t][1])
                ibs_a_vmax.append(conf_curve[t][2])
        plt.plot(ibs_a_t, ibs_a_v, "g-", label="IBS walpha", drawstyle='steps-post')
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
    print("\tIBS_guide (alpha, walpha, all)")
    pareto_flowtime(int(argv[1]),int(argv[2]),argv[3])
