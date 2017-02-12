run_type particle               # particle-resolved run
output_prefix out/loss_part_constant # prefix of output files
n_repeat 1                      # number of Monte Carlo repeats
n_part 10000                    # number of Monte Carlo particles
restart no                      # whether to restart from saved state (yes/no)

t_max 3600                      # total simulation time (s)
del_t 60                        # timestep (s)
t_output 600                    # output interval (0 disables) (s)
t_progress 600                  # progress printing interval (0 disables) (s)

gas_data gas_data.dat           # file containing gas data
gas_init gas_init.dat           # initial gas mixing ratios

aerosol_data aero_data.dat      # file containing aerosol data
do_fractal no                   # whether to do fractal treatment
aerosol_init aero_init_dist.dat # aerosol initial condition file

temp_profile temp.dat           # temperature profile file
pressure_profile pressure.dat   # pressure profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas mixing ratios file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file
loss_function constant          # particle loss function

rel_humidity 0.999              # initial relative humidity (1)
latitude 40                     # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 0                    # start time (s since 00:00 UTC)
start_day 1                     # start day of year (UTC)

do_coagulation no               # whether to do coagulation (yes/no)
do_condensation no              # whether to do condensation (yes/no)
do_mosaic no                    # whether to do MOSAIC (yes/no)
do_nucleation no                # whether to do nucleation (yes/no)

rand_init 0                     # random initialization (0 to auto-generate)
allow_doubling yes              # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
do_select_weighting no          # whether to select weighting explicitly (yes/no)
record_removals no              # whether to record particle removals (yes/no)
do_parallel no                  # whether to run in parallel (yes/no)
