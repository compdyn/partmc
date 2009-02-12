run_type sect                   # sectional code
output_prefix out/brownian_sect # prefix of output files
process_spec process.dat        # processing specification
kernel brown                    # coagulation kernel

t_max 86400                     # total simulation time (s)
del_t 60                        # timestep (s)
t_output 3600                   # output interval (0 disables) (s)
t_progress 600                  # progress printing interval (0 disables) (s)

n_bin 220                       # number of bins
r_min 1e-10                     # minimum radius (m)
r_max 1e-4                      # maximum radius (m)

gas_data gas_data.dat           # file containing gas data
aerosol_data aero_data.dat      # file containing aerosol data
aerosol_init aero_init_dist.dat # initial aerosol distribution

temp_profile temp.dat           # temperature profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas concentrations file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file

rel_humidity 0.999              # initial relative humidity (1)
pressure 1e5                    # initial pressure (Pa)
latitude 40                     # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 0                    # start time (s since 00:00 UTC)
start_day 1                     # start day of year (UTC)

do_coagulation yes              # whether to do coagulation (yes/no)
