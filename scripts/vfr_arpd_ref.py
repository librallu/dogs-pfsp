#!/usr/bin/python3
from generate_bks_tabs_vfr import read_bks_file

N = [100,200,300,400,500,600,700,800]
M = [20,40,60]
I = range(1,11)

def vfr_arpd_ref():
    f = read_bks_file()
    res = {}
    for n in N:
        for m in M:
            v = sum([f["VFR{}_{}_{}".format(n,m,i)]["igrms"] for i in I])/len(I)
            res[(n,m)] = v
    return res


def vfr_arpd_refs(n,m):
    f = read_bks_file()
    return [f["VFR{}_{}_{}".format(n,m,i)]["igrms"] for i in I]

if __name__ == "__main__":
    refs = vfr_arpd_ref()
    for k in refs:
        print("{}\t{}".format(k,refs[k]))