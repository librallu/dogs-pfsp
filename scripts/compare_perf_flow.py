#!/usr/bin/python3
from flow_arpd_ref import flow_arpd_ref
from flow_pareto import read_igrms_arpd, read_iga_arpd
from helper_stats import points_from_perfprofile, read_perfprofile_file, best_known_at_time, apply_arpd_tab

ratio_cpu_ibs_bsch = 2584/2070 # ~ 1.24
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

TAI_I = range(10)

def get_ibs_results_tai(n,m,t,algoname):
    dirname = "../experiments/IBS_flowtime_tai__74aed6e9568574f41e722b3e3b98329846fd2f8f/"
    filenames = map(
        (lambda i: "{}IBS_{}_{}_{}_{}.perfprofile.json".format(dirname,algoname,n,m,i)), 
        TAI_I
    )
    points_list = map(lambda e:read_perfprofile_file(e), filenames)
    return list(map(lambda e:best_known_at_time(e,t), points_list))

if __name__ == "__main__":
    print("##### ALGIRTCT COMPARISON")
    ratio_cpu_ibs_igrms = 2584/737 # ~ 3.5
    igrms_arpd = read_igrms_arpd()
    iga_arpd = read_iga_arpd()
    for n,m in TAI_TYPES:
        ref = flow_arpd_ref[(n,m)]
        igrms_60 = igrms_arpd[(n,m)][0]
        igrms_120 = igrms_arpd[(n,m)][1]
        igrms_240 = igrms_arpd[(n,m)][2]
        iga_60 =   iga_arpd[(n,m)][0]
        iga_120 =  iga_arpd[(n,m)][1]
        iga_240 =  iga_arpd[(n,m)][2]
        alpha_60 =   sum(apply_arpd_tab(get_ibs_results_tai(n,m,n*m*60/2000/ratio_cpu_ibs_igrms,"forward_alpha")  , ref))/10
        walpha_60 =  sum(apply_arpd_tab(get_ibs_results_tai(n,m,n*m*60/2000/ratio_cpu_ibs_igrms,"forward_walpha") , ref))/10
        alpha_120 =  sum(apply_arpd_tab(get_ibs_results_tai(n,m,n*m*120/2000/ratio_cpu_ibs_igrms,"forward_alpha") , ref))/10
        walpha_120 = sum(apply_arpd_tab(get_ibs_results_tai(n,m,n*m*120/2000/ratio_cpu_ibs_igrms,"forward_walpha"), ref))/10
        alpha_240 =  sum(apply_arpd_tab(get_ibs_results_tai(n,m,n*m*240/2000/ratio_cpu_ibs_igrms,"forward_alpha") , ref))/10
        walpha_240 = sum(apply_arpd_tab(get_ibs_results_tai(n,m,n*m*240/2000/ratio_cpu_ibs_igrms,"forward_walpha"), ref))/10
        print("{} \t & {:0.2f} & {:0.2f} & {:0.2f} & {:0.2f}  & {:0.2f} & {:0.2f} & {:0.2f} & {:0.2f}  & {:0.2f} & {:0.2f} & {:0.2f} & {:0.2f} \\\\".format(
            "TAI{}\_{}".format(n,m),
            igrms_60, iga_60, alpha_60, walpha_60,
            igrms_120, iga_120, alpha_120, walpha_120,
            igrms_240, iga_240, alpha_240, walpha_240,
        ))
    print("##### CBSH COMPARISON")
    ratio_cpu_ibs_bsch = 2584/2070 # ~ 1.24
    for n,m in TAI_TYPES:
        ref = flow_arpd_ref[(n,m)]
        mrsils = mrsils_arpd[(n,m)]
        alpha = sum(apply_arpd_tab(get_ibs_results_tai(n,m,n*m*60/2000/ratio_cpu_ibs_bsch,"forward_alpha")  , ref))/10
        walpha = sum(apply_arpd_tab(get_ibs_results_tai(n,m,n*m*60/2000/ratio_cpu_ibs_bsch,"forward_walpha")  , ref))/10
        print("{} \t & {:0.2f} & {:0.2f} & {:0.2f} \\\\".format(
            "TAI{}\_{}".format(n,m),
            mrsils, alpha, walpha
        ))