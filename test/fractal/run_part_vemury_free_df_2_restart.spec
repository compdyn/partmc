run_type particle               # particle-resolved run
output_prefix out_dimless_t/restart/part_vemury_free_df_2     # prefix of output files
n_repeat 1                      # number of Monte Carlo repeats
n_part 10000                     # total number of particles
restart yes                      # whether to restart from saved state (yes/no)
restart_file out_dimless_t/part_vemury_free_df_2_0001_00000021.nc # saved state file to restart from

t_max 6.5e+4                       # total simulation time (s)
del_t 10                         # timestep (s)
t_output 5e+3                    # output interval (0 disables) (s)
t_progress 1e+3                   # progress printing interval (0 disables) (s)

temp_profile temp_free.dat           # temperature profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas mixing ratios file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file

rel_humidity 0.47               # initial relative humidity (1)
pressure 1e5                    # initial pressure (Pa)
latitude 40                     # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 0                    # start time (s since 00:00 UTC)
start_day 1                     # start day of year (UTC)

do_coagulation yes              # whether to do coagulation (yes/no)
coag_kernel vemury_free         # coagulation kernel
do_condensation no              # whether to do condensation (yes/no)
do_mosaic no                    # whether to do MOSAIC (yes/no)
do_nucleation no                # whether to do nucleation (yes/no)

rand_init 7                     # random initialization (0 to auto-generate)
allow_doubling yes              # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
record_removals yes             # whether to record particle removals (yes/no)
do_parallel no                  # whether to run in parallel (yes/no)

do_fractal yes                   # whether to do fractal treatment
frac_dim 2                     # fractal dimension
prime_radius 5e-10              # radius of monomer
vol_fill_factor 1               # volume filling factor
