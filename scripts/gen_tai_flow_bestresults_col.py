#!/usr/bin/python3
import csv
import json

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

def read_bks_file(filename="./input/bks_flowtime.csv"):
    tai_name_to_standard = {}
    id_inst = 1
    for (n,m) in TAI_TYPES:
        for i in I:
            tai_name_to_standard["TA{}".format(id_inst)] = "TA{}_{}_{}".format(n,m,i)
            id_inst += 1
    with open(filename) as csvfile:
        res = {}
        reader = csv.reader(csvfile)
        for row in reader:
            if row[0] == "instance": # skip header
                pass
            else:
                inst = tai_name_to_standard[row[0]]
                res[inst] = {'algirtct':int(row[1]), 'shake_ls':int(row[2])}
        return res


def extract_ibs_flowtime(bks_res, dir_results="../experiments/IBS_flowtime_tai__74aed6e9568574f41e722b3e3b98329846fd2f8f/"):
    for (n,m) in TAI_TYPES:
        for i in I:
            with open("{}IBS_forward_alpha_{}_{}_{}.perfprofile.json".format(dir_results, n, m, i)) as f:
                data = json.load(f)
                bks_res['TA{}_{}_{}'.format(n,m,i)]["IBS_alpha"] = data["points"][-1]["v"]
            with open("{}IBS_forward_walpha_{}_{}_{}.perfprofile.json".format(dir_results, n, m, i)) as f:
                data = json.load(f)
                bks_res['TA{}_{}_{}'.format(n,m,i)]["IBS_walpha"] = data["points"][-1]["v"]

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
    extract_ibs_flowtime(bks_res)
    tab = []

    types = [
        [
            (20, 5),
            (20, 10),
            (20, 20),
            (50, 5)
        ],
        [
            (50, 10),
            (50, 20),
            (100, 5),
            (100, 10)
        ],
        [
            (100, 20),
            (200, 10),
            (200, 20),
            (500, 20)
        ]
    ]

    id_inst = 1
    for t in types:
        for (n,m) in t:
            for i in I:
                k = "TA{}_{}_{}".format(n,m,i)
                v = min([
                    bks_res[k]['algirtct'],
                    bks_res[k]['shake_ls'],
                    bks_res[k]['IBS_alpha'],
                    bks_res[k]['IBS_walpha']
                ])
                print("{},{}".format(id_inst, k, v))
                id_inst += 1
