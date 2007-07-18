run_type mc                     # Monte Carlo
output_file out/sedi_bidisperse_mc_summary.d #  output filename
state_prefix out/sedi_bidisperse_mc_state # prefix of state files
n_loop 1                        # number of Monte Carlo loops
n_part 10001                    # number of Monte Carlo particles
kernel sedi                     # coagulation kernel

t_max 600                       # total simulation time (s)
del_t 1                         # timestep (s)
t_output 0                      # output interval (0 disables) (s)
t_state 10                      # state output interval (0 disables) (s)
t_progress 10                   # progress printing interval (0 disables) (s)

n_bin 250                       # number of bins
r_min 1e-8                      # minimum radius (m)
r_max 1e0                       # maximum radius (m)

gas_data gas_data_simple.dat    # file containing gas data
gas_init gas_init_simple.dat    # initial gas concentrations

aerosol_data aerosol_data_water.dat # file containing aerosol data
aerosol_init sedi_bidisperse_mc_init.dat # aerosol initial condition file

temp_profile temp_constant_15C.dat # temperature profile file
RH 0.999                        # initial relative humidity (1)
pressure 1e5                    # initial pressure (Pa)
rho_a 1.25                      # initial air density (kg/m^3)
latitude 40                     # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 0                    # start time (s since 00:00 UTC)
start_day 1                     # start day of year (UTC)

gas_emissions gas_none.dat      # gas emissions file
gas_emission_rate 0             # gas emission rate (s^{-1})
gas_background gas_none.dat     # background gas concentrations file
gas_dilution_rate 0             # gas dilution rate with background (s^{-1})

aerosol_emissions aerosol_none.dat # aerosol emissions file
aerosol_emission_rate 0         # aerosol emission rate (s^{-1})
aerosol_background aerosol_none.dat # aerosol background file
aerosol_dilution_rate 0         # aerosol dilution rate with background (s^{-1})

rand_init 17                    # random initialization (0 to auto-generate)
do_coagulation yes              # whether to do coagulation (yes/no)
allow_double no                 # whether to allow doubling (yes/no)
do_condensation no              # whether to do condensation (yes/no)
do_mosaic no                    # whether to do MOSAIC (yes/no)
do_restart no                   # whether to restart from stored state (yes/no)
restart_name XXXX.d             # filename to restart from
