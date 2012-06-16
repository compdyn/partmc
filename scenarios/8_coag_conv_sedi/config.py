#!/usr/bin/env python

import os

ref_dirname = "ref"
run_dirname = "runs"
fig_dirname = "figs"

n_nodes = 1
n_procs_per_node = 1

n_part_list = [
    #("1000", "1k", r"$10^3$"),
    #("10000", "10k", r"$10^4$"),
    ("100000", "100k", r"$10^5$"),
    ]

weight_list = [
    #("power", "0"),
    #("power", "-1"),
    #("power", "-2"),
    #("power", "-3"),
    ("nummass", None),
    ]

plot_weight_list = [
    ("power", "0"),
    ("power", "-3"),
    ("nummass", None),
    ]

def all_runs():
    for (n_part, n_part_name, n_part_tex) in n_part_list:
        for (weight_type, exponent) in weight_list:
            if weight_type == "power":
                name = "%s_%s%s" % (n_part_name, weight_type, exponent)
            else:
                name = "%s_%s" % (n_part_name, weight_type)
            run = {"name": name,
                   "n_part": n_part,
                   "n_part_name": n_part_name,
                   "n_part_tex": n_part_tex,
                   "weight_type": weight_type,
                   "exponent": exponent,
                   }
            yield run

if __name__ == "__main__":
    for run in all_runs():
        print run
