run_type particle               # particle-resolved run
output_prefix out/mosaic_restarted # prefix of output files
n_repeat 1                      # number of Monte Carlo repeats
n_part 20                       # total number of particles
restart yes                     # whether to restart from saved state (yes/no)
restart_file out/mosaic_0001_00000013.nc # saved state file to restart from

t_max 43200                     # total simulation time (s)
del_t 300                       # timestep (s)
t_output 3600                   # output interval (0 disables) (s)
t_progress 3600                 # progress printing interval (0 disables) (s)

n_bin 160                       # number of bins
d_min 1e-8                      # minimum diameter (m)
d_max 1e-3                      # maximum diameter (m)

temp_profile temp.dat           # temperature profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas mixing ratios file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file

rel_humidity 0.85               # initial relative humidity (1)
pressure 1e5                    # initial pressure (Pa)
latitude 0                      # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 43200                # start time (s since 00:00 UTC)
start_day 200                   # start day of year (UTC)

do_coagulation no               # whether to do coagulation (yes/no)
do_condensation no              # whether to do condensation (yes/no)
do_mosaic yes                   # whether to do MOSAIC (yes/no)
do_optical yes                  # whether to compute optical props (yes/no)
do_nucleation no                # whether to do nucleation (yes/no)

rand_init 0                     # random initialization (0 to auto-generate)
allow_doubling yes              # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
record_removals no              # whether to record particle removals (yes/no)
do_parallel no                  # whether to run in parallel (yes/no)
