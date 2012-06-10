run_type sectional              # sectional code run
output_prefix out/sedi_sect     # prefix of output files

t_max 600                       # total simulation time (s)
del_t 1                         # timestep (s)
t_output 300                    # output interval (0 disables) (s)
t_progress 60                   # progress printing interval (0 disables) (s)

n_bin 100                       # number of bins
d_min 1e-7                      # minimum diameter (m)
d_max 1                         # maximum diameter (m)

gas_data gas_data.dat           # file containing gas data
aerosol_data aero_data.dat      # file containing aerosol data
aerosol_init aero_init_dist.dat # initial aerosol distribution

temp_profile temp.dat           # temperature profile file
pres_profile pres.dat           # pressure profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas mixing ratios file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file

rel_humidity 0.999              # initial relative humidity (1)
latitude 40                     # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 0                    # start time (s since 00:00 UTC)
start_day 1                     # start day of year (UTC)

do_coagulation yes              # whether to do coagulation (yes/no)
coag_kernel sedi                # coagulation kernel
