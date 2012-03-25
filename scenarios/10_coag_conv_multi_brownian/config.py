#!/usr/bin/env python

import os

ref_dirname = "ref"
run_dirname = "runs"
fig_dirname = "figs"

n_repeat = 100

n_part_list = [
    ("1000", "1k", r"$10^3$"),
    #("10000", "10k", r"$10^4$"),
    #("100000", "100k", r"$10^5$"),
    ]

weight_list = [
    ("flat", None),
    ("flat_source", None)
    ]

ratio_list = [
    ("conc", "0.999"),
    ("equal", "0.5"),
    ]

def all_runs():
    for (i_part, (n_part, n_part_name, n_part_tex)) in enumerate(n_part_list):
        for (i_weight, (weight_type, exponent)) in enumerate(weight_list):
            for (i_ratio, (ratio_type, ratio)) in enumerate(ratio_list):
                if weight_type == "power":
                    name = "%s_%s_%s%s" % (n_part_name, ratio_type, weight_type, exponent)
                else:
                    name = "%s_%s_%s" % (n_part_name, ratio_type, weight_type)
                run = {"name": name,
                       "i_part": i_part,
                       "n_part": n_part,
                       "n_part_name": n_part_name,
                       "n_part_tex": n_part_tex,
                       "i_weight": i_weight,
                       "weight_type": weight_type,
                       "exponent": exponent,
                       "i_ratio": i_ratio,
                       "ratio_type": ratio_type,
                       "ratio": ratio,
                       }
                yield run

if __name__ == "__main__":
    for run in all_runs():
        print run
