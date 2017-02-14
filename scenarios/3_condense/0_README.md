
# Condensation scenario

J. Ching, N. Riemer, and M. West (2012) Impacts of black carbon mixing state on black carbon nucleation scavenging: Insights from a particle-resolved model, _J. Geophys. Res. 117(D23209), DOI: <http://dx.doi.org/10.1029/2012JD018269>.

This scenario demonstrates the use of PartMC in a "cloud parcel" mode. It also shows how a Python script (`3_make_specs.py`) can be used to generate spec files from a template to perform many runs.

There are 49 * 4 = 196 cloud parcel runs. Each one starts from one of the 49 hourly outputs of the 2_urban_plume2 scenario, with one of four different averaging techniques applied to the hourly output aerosol population ("ref" = no averaging, "comp" = composition averaging, "size" = size averaging, and "both" = composition and size averaging).

1. Run 1_copy_start.sh. This copies scenarios/urban_plume2/out files to the "start" directory.
2. Run 2_average.py. This produces composition and size-averaged start files (4 cases: ref, comp, size, both).
3. Run 3_make_specs.py. This produces spec files and temperature profiles.
4. Run 4_run.py. This runs the cloud parcel simulations for all 49 starting times, and for all 4 cases.
5. Run 5_clean.sh. This deletes all files created by the processes above.
