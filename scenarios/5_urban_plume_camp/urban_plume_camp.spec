run_type particle               # particle-resolved run
output_prefix out/urban_plume   # prefix of output files
n_repeat 1                      # number of Monte Carlo repeats
n_part 1000                     # total number of particles
restart no                      # whether to restart from saved state (yes/no)

t_max 3600                     # total simulation time (s)
del_t 60                        # timestep (s)
t_output 60                   # output interval (0 disables) (s)
t_progress 60                  # progress printing interval (0 disables) (s)

do_camp_chem yes
camp_config config_1.json
gas_init gas_init_camp.dat           # initial gas concentrations

do_fractal no                   # whether to do fractal treatment
aerosol_init aero_init_dist_camp.dat # aerosol initial condition file

temp_profile temp_camp.dat           # temperature profile file
pressure_profile pres.dat            # pressure profile file
height_profile height_camp.dat       # height profile file
gas_emissions gas_emit_camp.dat      # gas emissions file
gas_background gas_back_camp.dat     # background gas concentrations file
aero_emissions aero_emit_camp.dat    # aerosol emissions file
aero_background aero_back_camp.dat   # aerosol background file
loss_function none                   # loss function specification

rel_humidity 0.95               # initial relative humidity (1)
latitude 0                      # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 21600                # start time (s since 00:00 UTC)
start_day 200                   # start day of year (UTC)

do_coagulation yes              # whether to do coagulation (yes/no)
coag_kernel brown               # coagulation kernel
do_condensation no              # whether to do condensation (yes/no)
do_mosaic no                   # whether to do MOSAIC (yes/no)
do_nucleation no                # whether to do nucleation (yes/no)

rand_init 0                     # random initialization (0 to use time)
allow_doubling yes              # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
do_select_weighting no          # whether to select weighting explicitly (yes/no)
record_removals yes             # whether to record particle removals (yes/no)
do_parallel no                  # whether to run in parallel (yes/no)
