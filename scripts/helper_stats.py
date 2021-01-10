#!/usr/bin/python3
import json
from math import ceil
import numpy as np
import scipy.stats


def compute_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m-h, m, m+h


def get_confidence_curves(data, confidence=0.95):
    """
    @param data[t:int -> Vec<f64>]
    @return[t:int -> {min_interval:f64, avg:f64, max_interval:f64}]
    """
    res = {}
    for t in data:
        res[t] = compute_confidence_interval(data[t], confidence)
    return res


def apply_arpd(v, ref):
    return (v-ref)/ref * 100

def apply_arpd_tab(v_list, ref_list):
    return [ apply_arpd(v_list[i], ref_list[i]) for i in range(len(v_list)) ]

def points_from_perfprofile(filename, t_end=0):
    with open(filename) as f:
        res = {}
        perfprofile = json.load(f)
        current_point = perfprofile["points"][0]
        current_t = int(ceil(current_point["t"]))
        current_v = current_point["v"]
        res[current_t] = current_v
        t_max = current_t
        for e in perfprofile["points"][1:]:
            t = int(ceil(e["t"]))
            v = e["v"]
            for t_before in range(current_t, t):
                res[t_before] = current_v
            res[t] = v
            current_t, current_v = t,v
            t_max = max(t_max, t)
        # complete perfprofile between t_max and t_end
        for t in range(t_max+1,t_end):
            res[t] = current_v
        return res

def read_perfprofile_file(filename):
    with open(filename) as f:
        perfprofile = json.load(f)
        return perfprofile["points"]

def aggregate_points(points_list):
    """
    @param points_list[Vec<t:int -> f64>]
    @return [t:int -> Vec<f64>]
    """
    res = {}
    for k in points_list[0]:
        tmp_vec = []
        for e in points_list:
            if k in e:
                tmp_vec.append(e[k])
            else:
                tmp_vec = []
                break
        if len(tmp_vec) == len(points_list):
            res[k] = tmp_vec
    return res


def best_known_at_time(points, t):
    """
    @param points: performance profile points list
    @param t: time at which we want the best-known value
    """
    if points[0]["t"] > t:
        return None
    else:
        prev_v = points[0]["v"]
        i = 1
        while i < len(points) and points[i]["t"] < t:
            prev_v = points[i]["v"]
            i += 1
        return prev_v


