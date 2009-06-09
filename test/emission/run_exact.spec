run_type exact                  # Monte Carlo
output_prefix out/emission_exact # prefix of output files
num_den 1e9                     # particle number density (#/m^3)

t_max 86400                     # total simulation time (s)
t_output 1200                   # output interval (0 disables) (s)

n_bin 160                       # number of bins
r_min 1e-8                      # minimum radius (m)
r_max 1e-3                      # maximum radius (m)

gas_data gas_data.dat           # file containing gas data
aerosol_data aero_data.dat      # file containing aerosol data

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

soln zero                       # solution type
aerosol_init aero_init_dist.dat # aerosol initial condition file
