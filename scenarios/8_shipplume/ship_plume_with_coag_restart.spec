run_type particle                # particle-resolved run
output_prefix out_redist_HC_from2pm/restart_from600s/ship_plume_wc # prefix of output files
n_repeat 1                        # number of Monte Carlo loops
n_part 100000                    # total number of particles
restart yes                      # whether to restart from saved state (yes/no)
restart_file out_redist_HC_from2pm/ship_plume_wc_0001_00000061.nc   # saved state file to restart from

t_max 49800                     # total simulation time (s)
del_t 60                        # timestep (s)
t_output 600                   # output interval (0 disables) (s)
t_progress 60                   # progress printing interval (0 disables) (s)

temp_profile temp.dat           # temperature profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas concentrations file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file

rel_humidity 0.9               # initial relative humidity (1)
pressure 9.9252e4                    # initial pressure (Pa)
latitude 50.179                     # latitude (degrees, -90 to 90)
longitude -6.3298                     # longitude (degrees, -180 to 180)
altitude 174                      # altitude (m)
start_time 51396                # start time (s since 00:00 UTC)
start_day 165                   # start day of year (UTC)

do_coagulation yes              # whether to do coagulation (yes/no)
coag_kernel brown               # coagulation kernel
do_condensation no              # whether to do condensation (yes/no)
do_mosaic yes                   # whether to do MOSAIC (yes/no)
do_optical yes                  # whether to compute optical props (yes/no)
do_nucleation no                # whether to do nucleation (yes/no)

rand_init 7                     # random initialization (0 to use time)
allow_doubling yes              # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
record_removals yes             # whether to record particle removals (yes/no)
do_parallel no                  # whether to run in parallel (yes/no)

