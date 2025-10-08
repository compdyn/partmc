run_type modal                       # modal run
output_prefix data/modal_dpg_00001000000000000000_sig_2_5_emerson_grass

t_max 28800                          # total simulation time (s)
del_t 60                             # timestep (s)
t_output 3600                        # output interval (0 disables) (s)
t_progress 0                         # progress printing interval (0 disables) (s)

n_bin 1000                           # number of bins (for processing purposes)
d_min 4e-8                           # minimum diameter (m)
d_max 2.5e-3                         # maximum diameter (m)

gas_data gas_data.dat                # file containing gas data
aerosol_data aero_data.dat           # file containing aerosol data
do_fractal no                        # whether to do fractal treatment 
aerosol_init aero_init_dist.dat      # aerosol initial condition file

temp_profile temp.dat                # temperature profile file
pressure_profile pres.dat            # pressure profile file
height_profile height.dat            # height profile file
gas_emissions gas_emit.dat           # gas emissions file
gas_background gas_back.dat          # background gas concentrations file
aero_emissions aero_emit.dat         # aerosol emissions file
aero_background aero_back.dat        # aerosol background file
loss_function drydep                 # loss function specification
drydep_params drydep_grass_emerson.dat

rel_humidity 0.95                    # initial relative humidity (1)
latitude 0                           # latitude (degrees, -90 to 90)
longitude 0                          # longitude (degrees, -180 to 180)
altitude 0                           # altitude (m)
start_time 21600                     # start time (s since 00:00 UTC)
start_day 200                        # start day of year (UTC)
