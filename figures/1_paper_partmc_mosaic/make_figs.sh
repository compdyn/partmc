#!/bin/sh

for f in fig_aero_2d_all_no_coag.py \
    fig_aero_2d_all.py \
    fig_aero_2d_n_orig.py \
    fig_aero_2d_oc.py \
    fig_aero_2d_proj.py \
    fig_aero_2d_water.py \
    fig_aero_bc_mixing.py \
    fig_aero_dist_mixing.py \
    fig_aero_dist_size.py \
    fig_aero_particles.py \
    fig_aero_time_species.py \
    fig_aero_time_totals.py \
    fig_aero_water_dist.py \
    fig_env.py \
    fig_gas.py \
    fig_test_brownian.py \
    fig_test_emission.py \
    ; do
    echo
    echo "\vspace{1em}"
    echo "\begin{verbatim}"
    echo $f
    ./$f
    echo "\end{verbatim}"
done

