run_type particle               # particle-resolved run
output_prefix output/pap1_exp8_acc/freezing_part # prefix of output files
n_repeat 20                      # number of Monte Carlo repeats
n_part 10000                     # total number of particles
restart no                      # whether to restart from saved state (yes/no)
#restart_file output/chiexp_indsize_0.8/restart.nc

t_max 600                            # total simulation time (s)
del_t 1                           # timestep (s)
t_output 10                   # output interval (0 disables) (s)
t_progress 10                  # progress printing interval (0 disables) (s)

do_camp_chem no                 # whether to run the campible chemistry module

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
loss_function none              # particle loss function

rel_humidity 1.000              # initial relative humidity (1)
latitude 40                     # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 0                    # start time (s since 00:00 UTC)
start_day 1                     # start day of year (UTC)

do_coagulation no               # whether to do coagulation (yes/no)
#coag_kernel brown
do_condensation no             # whether to do condensation (yes/no)
#do_init_equilibrate yes         # whether to initially equilibrate water (yes/no) 
do_mosaic no                    # whether to do MOSAIC (yes/no)
do_nucleation no                # whether to do nucleation (yes/no)
do_immersion_freezing yes                 # whether to do freezing (yes/no)
#immersion_freezing_scheme singular
immersion_freezing_scheme ABIFM
#immersion_freezing_scheme const
#freezing_rate -.01123456789 
do_coating no

rand_init 1                     # random initialization (0 to auto-generate)
allow_doubling yes              # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
do_select_weighting yes          # whether to select weighting explicitly (yes/no)
weight_type flat
record_removals no              # whether to record particle removals (yes/no)
do_parallel no                  # whether to run in parallel (yes/no)
