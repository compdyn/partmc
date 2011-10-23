#!/usr/bin/env python

import os

ref_dirname = "ref"
run_dirname = "runs"

n_nodes = 10
n_procs_per_node = 12

def all_runs():
    for (n_part, n_part_name) in [
        ("1000", "1k"),
        ("10000", "10k"),
#        ("100000", "100k"),
        ]:
        for (weight_type, exponent) in [
            ("power", "0"),
            ("power", "-1"),
            ("power", "-2"),
            ("power", "-3"),
            ("nummass", None),
            ]:
            if weight_type == "power":
                name = "%s_%s%s" % (n_part_name, weight_type, exponent)
            else:
                name = "%s_%s" % (n_part_name, weight_type)
            run = {"name": name,
                   "n_part": n_part,
                   "n_part_name": n_part_name,
                   "weight_type": weight_type,
                   "exponent": exponent,
                   }
            yield run

if __name__ == "__main__":
    for run in all_runs():
        print run
