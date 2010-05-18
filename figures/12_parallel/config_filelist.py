#!/usr/bin/env python

import os, sys, re

main_dir = os.path.join(os.environ["HOME"], "parallel_runs_32")
extra_sub_path = "6_urban_plume_parallel/out"
#run_re = re.compile("^job\.([^.]+)\.[0-9]+$")
run_re = re.compile("^job\.(urban_plume_serial_big)\.[0-9]+$")
true_run_name = "urban_plume_serial"
datafile_parallel_re = re.compile("^(.+)_([0-9]{4})_([0-9]{4})_([0-9]{8})\.nc")
datafile_serial_re = re.compile("^(.+)_([0-9]{4})_([0-9]{8})\.nc")

runs = []
true_run = None
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
        if run_name == true_run_name:
            true_run = run
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

runs_by_base = {}
for run in runs:
    match = re.search("^(.+)_([0-9]+)$", run["name"])
    if match:
        run_base = match.group(1)
        run_size = int(match.group(2))
        run_group = runs_by_base.setdefault(run_base, {"name": run_base,
                                                       "run_items": []})
        run_group["run_items"].append({"size": run_size,
                                   "run": run})
runs_by_base = runs_by_base.values()
for run_group in runs_by_base:
    run_group["run_items"].sort(key=lambda run_item: run_item["size"])

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
