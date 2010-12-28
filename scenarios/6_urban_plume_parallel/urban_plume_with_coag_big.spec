run_type particle               # particle-resolved run
output_prefix out/urban_plume_wc # prefix of output files
n_repeat 1                      # number of Monte Carlo repeats
n_part 1048576                     # total number of particles
kernel brown                    # coagulation kernel
nucleate none                   # nucleation parameterization
restart no                      # whether to restart from saved state (yes/no)

t_max 86400                     # total simulation time (s)
del_t 60                        # timestep (s)
t_output 3600                    # output interval (0 disables) (s)
t_progress 600                  # progress printing interval (0 disables) (s)

n_bin 160                       # number of bins
d_min 1e-10                     # minimum diameter (m)
d_max 1e-5                      # maximum diameter (m)

weight none                     # weighting function

gas_data gas_data.dat           # file containing gas data
gas_init gas_init.dat           # initial gas concentrations

aerosol_data aero_data.dat      # file containing aerosol data
aerosol_init aero_init_dist.dat # aerosol initial condition file

temp_profile temp.dat           # temperature profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas concentrations file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file

rel_humidity 0.95               # initial relative humidity (1)
pressure 1e5                    # initial pressure (Pa)
latitude 0                      # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 21600                # start time (s since 00:00 UTC)
start_day 200                   # start day of year (UTC)

rand_init 0                     # random initialization (0 to use time)
do_coagulation yes              # whether to do coagulation (yes/no)
allow_doubling yes              # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
do_condensation no              # whether to do condensation (yes/no)
do_mosaic yes                   # whether to do MOSAIC (yes/no)
do_optical no                   # whether to compute optical props (yes/no)
record_removals no              # whether to record particle removals (yes/no)
do_parallel no                  # whether to run in parallel (yes/no)
