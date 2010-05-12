#!/usr/bin/env python

import os, sys, re

graph_width = 8
color_bar_offset = 0.5

diameter_axis_min = 0.01
diameter_axis_max = 1.0
num_diameter_bins = 70
diameter_axis_label = r'dry diameter $D\ (\rm\mu m)$'

bc_axis_min = 0
bc_axis_max = 80
num_bc_bins = 40

main_dir = os.path.join(os.environ["HOME"], "parallel_runs")
extra_sub_path = "6_urban_plume_parallel/out"
run_re = re.compile("^job\.([^.]+)\.[0-9]+$")
datafile_parallel_re = re.compile("^(.+)_([0-9]{4})_([0-9]{4})_([0-9]{8})\.nc")
datafile_serial_re = re.compile("^(.+)_([0-9]{4})_([0-9]{8})\.nc")

runs = []
main_dir_ls = os.listdir(main_dir)
main_dir_ls.sort()
for main_entry in main_dir_ls:
    main_entry_path = os.path.join(main_dir, main_entry)
    main_entry_path = os.path.join(main_entry_path, extra_sub_path)
    main_match = run_re.search(main_entry)
    if main_match:
        run_name = main_match.group(1)
        run_dir_ls = os.listdir(main_entry_path)
        run_dir_ls.sort()
        run = {"name": run_name,
               "loops": {}}
        runs.append(run)
        for run_entry in run_dir_ls:
            run_entry_path = os.path.join(main_entry_path, run_entry)
            run_match = datafile_parallel_re.search(run_entry)
            if run_match:
                data_base = run_match.group(1)
                data_loop = int(run_match.group(2))
                data_proc = int(run_match.group(3))
                data_index = int(run_match.group(4))
            else:
                run_match = datafile_serial_re.search(run_entry)
                if run_match:
                    data_base = run_match.group(1)
                    data_loop = int(run_match.group(2))
                    data_proc = 1
                    data_index = int(run_match.group(3))
                else:
                    continue
            loop = run["loops"].setdefault(data_loop, {"num": data_loop,
                                                       "indices": {}})
            index = loop["indices"].setdefault(data_index, {"num": data_index,
                                                            "procs": []})
            index["procs"].append({"num": data_proc,
                                   "filename": run_entry_path})
for run in runs:
    run["loops"] = run["loops"].values()
    run["loops"].sort(key=lambda loop: loop["num"])
    for loop in run["loops"]:
        loop["indices"] = loop["indices"].values()
        loop["indices"].sort(key=lambda index: index["num"])
        for index in loop["indices"]:
            index["procs"].sort(key=lambda proc: proc["num"])

if __name__ == "__main__":
    for run in runs:
        print "run name: %s" % run["name"]
        print "run loops: %d" % len(run["loops"])
        for loop in run["loops"]:
            print "loop num: %d" % loop["num"]
            print "loop indices: %d" % len(loop["indices"])
            for index in loop["indices"]:
                print "index num: %d" % index["num"]
                print "index procs: %d" % len(index["procs"])
                for proc in index["procs"]:
                    print "proc num: %d" % proc["num"]
                    print "proc filename: %s" % proc["filename"]

        
