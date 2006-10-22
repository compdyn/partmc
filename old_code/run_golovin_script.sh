#!/bin/sh
./run_golovin_adapt
./run_golovin_exact
./run_golovin_fix
./run_golovin_var

./process_out out_golovin_adapt.d
./process_out out_golovin_exact.d
./process_out out_golovin_fix.d
./process_out out_golovin_var.d

