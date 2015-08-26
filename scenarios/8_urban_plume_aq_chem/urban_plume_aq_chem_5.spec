run_type particle               # particle-resolved run
output_prefix out/urban_plume_aq_chem_e   # prefix of output files
n_repeat 1                      # number of Monte Carlo repeats
n_part 10000                    # total number of particles
restart yes                     # whether to restart from saved state (yes/no)
restart_file out/urban_plume_aq_chem_0001_00000241.nc

t_max 1200                      # total simulation time (s)
del_t 1                         # timestep (s)
t_output 60                     # output interval (0 disables) (s)
t_progress 600                  # progress printing interval (0 disables) (s)

temp_profile temp_2.dat         # temperature profile file
pressure_profile pres_2.dat     # pressure profile file
height_profile height_2.dat     # height profile file
gas_emissions gas_emit_2.dat    # gas emissions file
gas_background gas_back_2.dat   # background gas concentrations file
aero_emissions aero_emit_2.dat  # aerosol emissions file
aero_background aero_back_2.dat # aerosol background file

rel_humidity 0.95               # initial relative humidity (1)
latitude 0                      # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 108000               # start time (s since 00:00 UTC)
start_day 200                   # start day of year (UTC)

do_coagulation no               # whether to do coagulation (yes/no)
# coag_kernel brown               # coagulation kernel
do_loss no                      # whether to do particle loss (yes/no)
do_condensation yes             # whether to do condensation (yes/no)
do_init_equilibriate yes
do_mosaic no                    # whether to do MOSAIC (yes/no)
do_nucleation no                # whether to do nucleation (yes/no)

rand_init 0                     # random initialization (0 to use time)
allow_doubling yes              # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
record_removals yes             # whether to record particle removals (yes/no)

do_parallel no                  # whether to run in parallel (yes/no)
# output_type single              # parallel output type (central/dist/single)
# mix_timescale 0                 # mixing timescale between processors (s)
# gas_average yes                 # whether to average gases each timestep
# env_average yes                 # whether to average environment each timestep
# parallel_coag local             # parallel coagulation method (local/dist)

do_aq_chem yes                  # whether to do aqueous chemistry (yes/no)
aq_mech capram24_red+_mod.txt       # aqueous chemical mechanism file
aq_spec aq_spec_data_with_abstol.dat    # file containing aqueous species data
aq_map aq_spec_map.dat          # file containing map between species in aq. mechanism and PartMC
aq_init aq_spec_init.dat        # file containing initial and constant species concentrations

