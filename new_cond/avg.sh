#!/bin/sh

../build/bin_average_size 1e-10 1e-5 25 dry average start/urban_plume_wc_0001_00000002.nc start_size/urban_plume_wc
../build/bin_average_size 1e-10 1e-5 25 dry average start/urban_plume_wc_0001_00000008.nc start_size/urban_plume_wc
../build/bin_average_size 1e-10 1e-5 25 dry average start/urban_plume_wc_0001_00000016.nc start_size/urban_plume_wc
../build/bin_average_size 1e-10 1e-5 25 dry average start/urban_plume_wc_0001_00000025.nc start_size/urban_plume_wc

../build/bin_average_comp 1e-10 1e-5 25 dry start/urban_plume_wc_0001_00000002.nc start_comp/urban_plume_wc
../build/bin_average_comp 1e-10 1e-5 25 dry start/urban_plume_wc_0001_00000008.nc start_comp/urban_plume_wc
../build/bin_average_comp 1e-10 1e-5 25 dry start/urban_plume_wc_0001_00000016.nc start_comp/urban_plume_wc
../build/bin_average_comp 1e-10 1e-5 25 dry start/urban_plume_wc_0001_00000025.nc start_comp/urban_plume_wc

../build/bin_average_size 1e-10 1e-5 25 dry average start_comp/urban_plume_wc_0001_00000002.nc start_both/urban_plume_wc
../build/bin_average_size 1e-10 1e-5 25 dry average start_comp/urban_plume_wc_0001_00000008.nc start_both/urban_plume_wc
../build/bin_average_size 1e-10 1e-5 25 dry average start_comp/urban_plume_wc_0001_00000016.nc start_both/urban_plume_wc
../build/bin_average_size 1e-10 1e-5 25 dry average start_comp/urban_plume_wc_0001_00000025.nc start_both/urban_plume_wc
