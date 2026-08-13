run_type sectional              # sectional code run
output_prefix out/mosaic_sect   # prefix of output files

# timing matches the particle-resolved run_part.spec so the two can be compared
t_max 86400                     # total simulation time (s)
del_t 300                       # timestep (s)
t_output 3600                   # output interval (0 disables) (s)
t_progress 3600                 # progress printing interval (0 disables) (s)

do_tchem no                     # whether to use TChem for chemistry
do_mosaic yes                   # whether to do MOSAIC (yes/no)

n_bin 40                        # number of bins
d_min 1e-8                      # minimum diameter (m)
d_max 1e-5                      # maximum diameter (m)

gas_data gas_data.dat           # file containing gas data
gas_init gas_init.dat           # initial gas mixing ratios

aerosol_data aero_data.dat      # file containing aerosol data
do_fractal no                   # whether to do fractal treatment
aerosol_init aero_init_dist.dat # aerosol initial condition file (monodisperse)

temp_profile temp.dat           # temperature profile file
pressure_profile pressure.dat   # pressure profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas mixing ratios file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file
loss_function none              # particle loss function

rel_humidity 0.85               # initial relative humidity (1)
latitude 0                      # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 43200                # start time (s since 00:00 UTC)
start_day 200                   # start day of year (UTC)

do_coagulation no               # whether to do coagulation (yes/no)
