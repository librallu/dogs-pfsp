#!/usr/bin/python3

N = [100,200,300,400,500,600,700,800]
M = [20,40,60]

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

def extract_igirms_arpd_ref(filename):
    with open(filename, "r") as f:
        l = f.readlines()[1:]
        l = [ [ a.replace(" ", "").replace("\n","").replace("_Gap","") for a in e.split(",") ][1:] for e in l ]
        l = [ (e[0],int(e[1])) for e in l ]
        res = {}
        for e in l:
            res[e[0]] = e[1] 
        return res

def average_igirms_arpd(igirms_arpd_ref):
    I = range(1,10)
    res = {}
    for n in N:
        for m in M:
            avg = 0
            for i in I:
                avg += igirms_arpd_ref["VFR{}_{}_{}".format(n,m,i)]
            res["VFR{}_{}".format(n,m)] = avg/len(I)
    return res

def average_igirms_arpd_tai():
    with open("input/igirms_tai.csv", "r") as f:
        res = {}
        l = f.readlines()[1:]
        l = [ [ a.replace("\n","") for a in e.split(",") ] for e in l ]
        # print(l)
        for e in l:
            inst = [int(i) for i in e[0].split('_')]
            vals = [float(i) for i in e[1:]]
            res[tuple(inst)] = vals
        return res


def average_igirms_arpd_vfr():
    with open("input/igirms_vfr.csv", "r") as f:
        res = {}
        l = f.readlines()[1:]
        l = [ [ a.replace("\n","").replace("VFR","") for a in e.split(",") ] for e in l ]
        # print(l)
        for e in l:
            inst = [int(i) for i in e[0].split('_')]
            vals = [float(i) for i in e[1:]]
            res[tuple(inst)] = vals
        return res


if __name__ == "__main__":
    igirms_arpd_ref = extract_igirms_arpd_ref("input/vfr_igirms.csv")
    avg_igirms_arpd = average_igirms_arpd(igirms_arpd_ref)
    for k in avg_igirms_arpd.keys():
        print("{}\t{}".format(k, avg_igirms_arpd[k]))