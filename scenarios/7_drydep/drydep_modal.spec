run_type modal                         # modal run
output_prefix out/modal_dpg_10_0_sig_2_5_emerson_broad

t_max 28800                            # total simulation time (s)
del_t 60                               # timestep (s)
t_output 3600                          # output interval (0 disables) (s)
t_progress 0                           # progress printing interval (0 disables) (s)

n_bin 1000                             # number of bins (for processing purposes)
d_min 4e-8                             # minimum diameter (m)
d_max 2.5e-3                           # maximum diameter (m)

do_camp_chem no                        # whether to use CAMP for chemistry (yes/no)
do_tchem no                            # whether to use TChem for chemistry (yes/no)

gas_data gas_data.dat                  # file containing gas data
aerosol_data aero_data.dat             # file containing aerosol data
do_fractal no                          # whether to do fractal treatment (yes/no)
aerosol_init aero_init_run.dat

temp_profile temp.dat                  # temperature profile file
pressure_profile pres.dat              # pressure profile file
height_profile height.dat              # height profile file
gas_emissions gas_emit.dat             # gas emissions file
gas_background gas_back.dat            # background gas concentrations file
aero_emissions aero_emit.dat           # aerosol emissions file
aero_background aero_back.dat          # aerosol background file
loss_function drydep                   # loss function specification
drydep_params drydep_broad_emerson.dat

rel_humidity 0.95                      # initial relative humidity (1)
latitude 0                             # latitude (degrees, -90 to 90)
longitude 0                            # longitude (degrees, -180 to 180)
altitude 0                             # altitude (m)
start_time 21600                       # start time (s since 00:00 UTC)
start_day 200                          # start day of year (UTC)

do_coagulation no                      # whether to do coagulation (yes/no)
do_condensation no                     # whether to do condensation (yes/no)
do_mosaic no                           # whether to do MOSAIC (yes/no)
do_optical no                          # whether to compute optical props (yes/no)
do_nucleation no                       # whether to do nucleation (yes/no)
do_immersion_freezing no               # whether to do freezing (yes/no)
