#!/usr/bin/python3
import csv
import json

N = [100,200,300,400,500,600,700,800]
M = [20,40,60]
I = range(1,11)

def read_bks_file(filename="./input/bks_vfr.csv"):
    with open(filename) as csvfile:
        res = {}
        reader = csv.reader(csvfile)
        for row in reader:
            if row[0] == "instance": # skip header
                pass
            else:
                inst = (row[0].split("_Gap")[0]).replace(" ",'')
                res_igrms = int(row[1])
                res_vbih = int(row[2])
                res[inst] = {'igrms':res_igrms, 'vbih':res_vbih}
        return res

def extract_igbob_vrf(bks_res, dir_results="../experiments/igbob_vrf__2838421f1bb08a8464e6118b7c7dcb848f4827d4/"):
    for n in N:
        for m in M:
            for i in I:
                t = 90
                with open("{}perfprofile_igbob_{}_{}_{}_{}".format(dir_results, t, n, m, i)) as f:
                    data = json.load(f)
                    bks_res['VFR{}_{}_{}'.format(n,m,i)]["igbob"] = min(data["objlist"])


def extract_ibs_vrf(bks_res, dir_results="../experiments/EXPERIMENTS_IBS_VFR/"):
    for n in N:
        for m in M:
            for i in I:
                with open("{}IBS_VFR_bi_min_gap_{}_{}_{}.perfprofile.json".format(dir_results, n, m, i)) as f:
                    data = json.load(f)
                    bks_res['VFR{}_{}_{}'.format(n,m,i)]["IBS_gap"] = data["points"][-1]["v"]
                with open("{}IBS_VFR_bi_min_wfrontalpha_{}_{}_{}.perfprofile.json".format(dir_results, n, m, i)) as f:
                    data = json.load(f)
                    bks_res['VFR{}_{}_{}'.format(n,m,i)]["IBS_wfrontalpha"] = data["points"][-1]["v"]


def is_min_array(array):
    maxi = min(array)
    return list(map(lambda e: e==maxi, array))

def humanize_int(v):
    if v < 1000:
        return str(v)
    else:
        a = int(v/1000)
        b = int(v%1000)
        # res = humanize_int(a)
        return humanize_int(a)+"."+"{:03}".format(b)

def format_values(array):
    res = []
    for v in array:
        if v >= 1e9:
            res.append(" - ")
        else:
            res.append(humanize_int(v))
    return res

def apply_bold(values, flags):
    res = []
    for i,v in enumerate(values):
        if flags[i]:
            res.append("\\textbf{"+str(v)+"}")
        else:
            res.append(str(v))
    return res

if __name__ == "__main__":
    bks_res = read_bks_file()
    extract_igbob_vrf(bks_res)
    extract_ibs_vrf(bks_res)
    tab = []
    for n in N:
        for m in M:
            for i in I:
                k = "VFR{}_{}_{}".format(n,m,i)
                print("{}\t{}\t{}\t{}\t{}\t{}".format(
                    k,
                    bks_res[k]['igrms'],
                    bks_res[k]['vbih'],
                    bks_res[k]['igbob'],
                    bks_res[k]['IBS_gap'],
                    bks_res[k]['IBS_wfrontalpha']
                ))
                # tab.append([
                #     k,
                #     bks_res[k]['igrms'],
                #     bks_res[k]['vbih'],
                #     bks_res[k]['igbob'],
                #     bks_res[k]['IBS_gap'],
                #     bks_res[k]['IBS_wfrontalpha']
                # ])

    gmys_bks = {
        (100,20,1): 6121,
        (100,20,2): 6224,
        (100,20,3): 6157,
        (100,20,4): 6173,
        (100,20,5): 6221,
        (100,20,6): 6247,
        (100,20,7): 6358,
        (100,20,8): 6023,
        (100,20,9): 6286,
        (100,20,10): 6048,

        (200,20,1): 11181,
        (200,20,2): 11254,
        (200,20,3): 11233,
        (200,20,4): 11090,
        (200,20,5): 11076,
        (200,20,6): 11208,
        (200,20,7): 11266,
        (200,20,8): 11041,
        (200,20,9): 11008,
        (200,20,10): 11193,
        
        (300,20,1): 15996,
        (300,20,2): 16409,
        (300,20,3): 16010,
        (300,20,4): 16052,
        (300,20,5): 21399,
        (300,20,6): 16021,
        (300,20,7): 16188,
        (300,20,8): 16287,
        (300,20,9): 16203,
        (300,20,10): 16780,

        (400,20,1): 20952,
        (400,20,2): 21346,
        (400,20,3): 21379,
        (400,20,4): 21125,
        (400,20,5): 16245,
        (400,20,6): 21075,
        (400,20,7): 21507,
        (400,20,8): 21198,
        (400,20,9): 21236,
        (400,20,10): 21456,

        (500,20,1): 26253,
        (500,20,2): 26555,
        (500,20,3): 26268,
        (500,20,4): 25994,
        (500,20,5): 26703,
        (500,20,6): 26325,
        (500,20,7): 26313,
        (500,20,8): 26217,
        (500,20,9): 26345,
        (500,20,10): 26345,

        (600,20,1): 31303,
        (600,20,2): 31281,
        (600,20,3): 31374,
        (600,20,4): 31417,
        (600,20,5): 31323,
        (600,20,6): 31613,
        (600,20,7): 31461,
        (600,20,8): 31414,
        (600,20,9): 31473,
        (600,20,10): 31021,

        (700,20,1): 36285,
        (700,20,2): 36220,
        (700,20,3): 36419,
        (700,20,4): 36361,
        (700,20,5): 36496,
        (700,20,6): 36556,
        (700,20,7): 36540,
        (700,20,8): 36418,
        (700,20,9): 36212,
        (700,20,10): 36362,

        (800,20,1): 41413,
        (800,20,2): 41282,
        (800,20,3): 41319,
        (800,20,4): 41375,
        (800,20,5): 41626,
        (800,20,6): 41919,
        (800,20,7): 41342,
        (800,20,8): 41390,
        (800,20,9): 41697,
        (800,20,10): 41489,

    }

    N = [100,200]
    # N = [300,400]
    # N = [500,600]
    # N = [700,800]
    M = [20,40,60]

    print("\\begin{tabular}{c|cccc|cc}")
    print("instance & Gmys\_{B\&B} & IGrms & VBIH & IGbob & IBS\_gap & IBS\_wfrontalpha \\\\")
    print("\\hline")
    for n in N:
        for m in M:
            for i in I:
                k = "VFR{}_{}_{}".format(n,m,i)
                res_gmys = 1e10
                if (n,m,i) in gmys_bks:
                    res_gmys = gmys_bks[(n,m,i)]
                tab = [
                    res_gmys,
                    bks_res[k]['igrms'],
                    bks_res[k]['vbih'],
                    bks_res[k]['igbob'],
                    bks_res[k]['IBS_gap'],
                    bks_res[k]['IBS_wfrontalpha']
                ]
                bold_flag = is_min_array(tab)
                print("{} & {} \\\\".format(k.replace("_","\\_"), " & ".join(apply_bold(format_values(tab), bold_flag))))
    print("\\end{tabular}")
    # print(tab)
        