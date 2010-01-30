run_type particle               # particle-resolved run
output_prefix out/mosaic_restarted # prefix of output files
n_loop 1                        # number of Monte Carlo loops
n_part 20                       # total number of particles
kernel golovin                  # coagulation kernel
nucleate none                   # nucleation parameterization
restart yes                     # whether to restart from saved state (yes/no)
restart_file out/mosaic_0001_00000013.nc # saved state file to restart from

t_max 43200                     # total simulation time (s)
del_t 300                       # timestep (s)
t_output 3600                   # output interval (0 disables) (s)
t_progress 600                  # progress printing interval (0 disables) (s)

n_bin 160                       # number of bins
r_min 1e-8                      # minimum radius (m)
r_max 1e-3                      # maximum radius (m)

weight power                    # weighting function
ref_radius 1e-7                 # radius at which weight is 1
exponent -1                     # weighting exponent

gas_data gas_data.dat           # file containing gas data

aerosol_data aero_data.dat      # file containing aerosol data

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

rand_init 0                     # random initialization (0 to auto-generate)
do_coagulation no               # whether to do coagulation (yes/no)
allow_doubling yes              # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
do_condensation no              # whether to do condensation (yes/no)
do_mosaic yes                   # whether to do MOSAIC (yes/no)
record_removals no              # whether to record particle removals (yes/no)
do_parallel no                  # whether to run in parallel (yes/no)
