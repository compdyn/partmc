#!/bin/sh

cat <<ENDINFO

Urban Plume Test-case
---------------------

This simulates an idealized urban plume with gas and aerosol
chemistry.

ENDINFO
sleep 1

echo "../build/partmc urban_plume_with_coag.spec"
../build/partmc urban_plume_with_coag.spec

echo "../build/partmc urban_plume_no_coag.spec"
../build/partmc urban_plume_no_coag.spec

echo "../build/extract_env out/urban_plume_wc_0001_ out/urban_plume_wc_env.txt"
../build/extract_env out/urban_plume_wc_0001_ out/urban_plume_wc_env.txt
echo "../build/extract_gas out/urban_plume_wc_0001_ out/urban_plume_wc_gas.txt"
../build/extract_gas out/urban_plume_wc_0001_ out/urban_plume_wc_gas.txt
echo "../build/extract_aero_total out/urban_plume_wc_0001_ out/urban_plume_wc_aero_total.txt"
../build/extract_aero_total out/urban_plume_wc_0001_ out/urban_plume_wc_aero_total.txt
echo "../build/extract_aero_species out/urban_plume_wc_0001_ out/urban_plume_wc_aero_species.txt"
../build/extract_aero_species out/urban_plume_wc_0001_ out/urban_plume_wc_aero_species.txt
echo "../build/extract_aero_size_num 1e-9 1e-6 100 out/urban_plume_wc_0001_ out/urban_plume_wc_aero_size_num.txt"
../build/extract_aero_size_num 1e-9 1e-6 100 out/urban_plume_wc_0001_ out/urban_plume_wc_aero_size_num.txt
echo "../build/extract_aero_size_mass 1e-9 1e-6 100 out/urban_plume_wc_0001_ out/urban_plume_wc_aero_size_mass.txt"
../build/extract_aero_size_mass 1e-9 1e-6 100 out/urban_plume_wc_0001_ out/urban_plume_wc_aero_size_mass.txt

echo "../build/extract_aero_particle_mass out/urban_plume_wc_0001_00000001.nc out/urban_plume_wc_aero_particle_mass_00000001.txt"
../build/extract_aero_particle_mass out/urban_plume_wc_0001_00000001.nc out/urban_plume_wc_aero_particle_mass_00000001.txt
echo "../build/extract_aero_particle_mass out/urban_plume_wc_0001_00000301.nc out/urban_plume_wc_aero_particle_mass_00000301.txt"
../build/extract_aero_particle_mass out/urban_plume_wc_0001_00000301.nc out/urban_plume_wc_aero_particle_mass_00000301.txt
echo "../build/extract_aero_particle_mass out/urban_plume_wc_0001_00000421.nc out/urban_plume_wc_aero_particle_mass_00000421.txt"
../build/extract_aero_particle_mass out/urban_plume_wc_0001_00000421.nc out/urban_plume_wc_aero_particle_mass_00000421.txt
echo "../build/extract_aero_particle_mass out/urban_plume_wc_0001_00001441.nc out/urban_plume_wc_aero_particle_mass_00001441.txt"
../build/extract_aero_particle_mass out/urban_plume_wc_0001_00001441.nc out/urban_plume_wc_aero_particle_mass_00001441.txt

echo "../build/extract_env out/urban_plume_nc_0001_ out/urban_plume_nc_env.txt"
../build/extract_env out/urban_plume_nc_0001_ out/urban_plume_nc_env.txt
echo "../build/extract_gas out/urban_plume_nc_0001_ out/urban_plume_nc_gas.txt"
../build/extract_gas out/urban_plume_nc_0001_ out/urban_plume_nc_gas.txt
echo "../build/extract_aero_total out/urban_plume_nc_0001_ out/urban_plume_nc_aero_total.txt"
../build/extract_aero_total out/urban_plume_nc_0001_ out/urban_plume_nc_aero_total.txt
echo "../build/extract_aero_species out/urban_plume_nc_0001_ out/urban_plume_nc_aero_species.txt"
../build/extract_aero_species out/urban_plume_nc_0001_ out/urban_plume_nc_aero_species.txt
echo "../build/extract_aero_size_num 1e-9 1e-6 100 out/urban_plume_nc_0001_ out/urban_plume_nc_aero_size_num.txt"
../build/extract_aero_size_num 1e-9 1e-6 100 out/urban_plume_nc_0001_ out/urban_plume_nc_aero_size_num.txt
echo "../build/extract_aero_size_mass 1e-9 1e-6 100 out/urban_plume_nc_0001_ out/urban_plume_nc_aero_size_mass.txt"
../build/extract_aero_size_mass 1e-9 1e-6 100 out/urban_plume_nc_0001_ out/urban_plume_nc_aero_size_mass.txt

echo "../build/extract_aero_particle_mass out/urban_plume_nc_0001_00000001.nc out/urban_plume_nc_aero_particle_mass_00000001.txt"
../build/extract_aero_particle_mass out/urban_plume_nc_0001_00000001.nc out/urban_plume_nc_aero_particle_mass_00000001.txt
echo "../build/extract_aero_particle_mass out/urban_plume_nc_0001_00000301.nc out/urban_plume_nc_aero_particle_mass_00000301.txt"
../build/extract_aero_particle_mass out/urban_plume_nc_0001_00000301.nc out/urban_plume_nc_aero_particle_mass_00000301.txt
echo "../build/extract_aero_particle_mass out/urban_plume_nc_0001_00000421.nc out/urban_plume_nc_aero_particle_mass_00000421.txt"
../build/extract_aero_particle_mass out/urban_plume_nc_0001_00000421.nc out/urban_plume_nc_aero_particle_mass_00000421.txt
echo "../build/extract_aero_particle_mass out/urban_plume_nc_0001_00001441.nc out/urban_plume_nc_aero_particle_mass_00001441.txt"
../build/extract_aero_particle_mass out/urban_plume_nc_0001_00001441.nc out/urban_plume_nc_aero_particle_mass_00001441.txt
