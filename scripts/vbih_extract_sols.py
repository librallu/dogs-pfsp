#!/usr/bin/python3
from vfr_arpd_ref import vfr_arpd_ref, vfr_arpd_refs
from helper_stats import compute_confidence_interval, apply_arpd_tab

N = [100,200,300,400,500,600,700,800]
M = [20,40,60]
I = range(1,11)

def extract_vbih_sols(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        i = 0
        lines = lines[1:] # skip header
        res = {}
        for l in lines:
            a,b,_ = l.split("\t")
            if a == "" or b == "":
                break
            n,m = a.split("x")
            v = int(b)
            res[(int(n),int(m),i+1)] = v
            i = (i+1) % 10
        return res

def vbih_arpd(filename, n, m):
    tmp = extract_vbih_sols(filename)
    ref = vfr_arpd_ref()[(n,m)]
    v = sum([tmp[(n,m,i)] for i in I])/len(I)
    return (v-ref)/ref*100


def vbih_confidence_interval(filename, n, m):
    tmp = extract_vbih_sols(filename)
    tmp = [ tmp[(n,m,i)] for i in I ]
    ref_tab = vfr_arpd_refs(n,m)
    arpd_tab = apply_arpd_tab(tmp, ref_tab)
    return compute_confidence_interval(arpd_tab)


if __name__ == "__main__":
    res = extract_vbih_sols("./input/Permutations_VBIH_15nm.txt")
    for k in res:
        print("{}\t{}".format(k,res[k]))
    print(vbih_arpd("./input/Permutations_VBIH_15nm.txt", 800, 60))