run_type mc                     # Monte Carlo run
output_file out/mosaic_summary.d # name of output file
state_prefix out/mosaic_state   # prefix of state files
n_loop 1                        # number of Monte Carlo loops
n_part 3                        # total number of particles
kernel golovin                  # coagulation kernel

t_max 259200                    # total simulation time (s)
del_t 10                        # timestep (s)
t_output 600                    # output interval (0 disables) (s)
t_state 0                       # state output interval (0 disables) (s)
t_progress 600                  # progress printing interval (0 disables) (s)

n_bin 160                       # number of bins
r_min 1e-8                      # minimum radius (m)
r_max 1e-3                      # maximum radius (m)

gas_data mosaic_gas_data.dat    # file containing gas data
gas_init mosaic_gas_init.dat    # initial gas concentrations

aerosol_data mosaic_aerosol_data.dat # file containing aerosol data
aerosol_init mosaic_aerosol_init.dat # aerosol initial condition file

temp_profile temp_constant_15C.dat # temperature profile file
RH 0.85                         # initial relative humidity (1)
pressure 1e5                    # initial pressure (Pa)
rho_a 1.25                      # initial air density (kg/m^3)
latitude 40                     # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 43200                # start time (s since 00:00 UTC)
start_day 200                   # start day of year (UTC)

gas_emissions mosaic_gas_emissions.dat # gas emissions file
gas_emission_rate 10            # gas emission rate (s^{-1})
gas_background mosaic_gas_background.dat # background gas concentrations file
gas_dilution_rate 0             # gas dilution rate with background (s^{-1})

aerosol_emissions mosaic_aerosol_emissions.dat # aerosol emissions file
aerosol_emission_rate 0         # aerosol emission rate (s^{-1})
aerosol_background mosaic_aerosol_background.dat # aerosol background file
aerosol_dilution_rate 0         # aerosol dilution rate with background (s^{-1})

rand_init 17                    # random initialization (0 to auto-generate)
do_coagulation no               # whether to do coagulation (yes/no)
allow_double yes                # whether to allow doubling (yes/no)
do_condensation no              # whether to do condensation (yes/no)
do_mosaic yes                   # whether to do MOSAIC (yes/no)
do_restart no                   # whether to restart from stored state (yes/no)
restart_name XXXX.d             # filename to restart from
