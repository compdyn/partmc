run_type particle               # Monte Carlo run
output_prefix out_we/out_20/condense      # prefix of output files
n_repeat 1                      # number of Monte Carlo repeats
n_part 10000                    # total number of particles
restart yes                      # whether to restart from saved state (yes/no)
restart_file start/urban_plume2_wc_0001_00000020.nc

t_max 600                       # total simulation time (s)
del_t 1                         # timestep (s)
t_output 10                     # output interval (0 disables) (s)
t_progress 60                   # progress printing interval (0 disables) (s)

#gas_data gas_data.dat           # file containing gas data
#gas_init gas_init.dat           # initial gas mixing ratios

#aerosol_data aero_data.dat      # file containing aerosol data
#aerosol_init aero_init_dist.dat # aerosol initial condition file

temp_profile temp_001_00000020.dat           # temperature profile file
pressure_profile pres_001_00000020.dat       # pressure profile file
do_entrain yes                              # whether to do entrainment (yes/no)
q_tot_profile qtot_001_00000020.dat         # total water mixing ratio profile file
q_background q_back.dat         # background spec. humidity profile file
height_profile height.dat       # height profile file

gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas concentrations file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file

rel_humidity 0.95               # initial relative humidity (1)
latitude 0                      # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 81000                # start time (s since 00:00 UTC)
start_day 200                   # start day of year (UTC)

do_coagulation no               # whether to do coagulation (yes/no)
do_condensation yes             # whether to do condensation (yes/no)
do_init_equilibriate yes        # whether to initially equilibriate water (yes/no)
do_mosaic no                    # whether to do MOSAIC (yes/no)
do_nucleation no                # whether to do nucleation (yes/no)

rand_init 0                     # random initialization (0 to use time)
allow_doubling yes              # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
record_removals no              # whether to record particle removals (yes/no)
do_parallel no                  # whether to run in parallel (yes/no)
