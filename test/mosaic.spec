run_type mc                     # Monte Carlo run
output_file out/mosaic_out.d    # name of output file
state_prefix out/mosaic_state   # prefix of state files
n_loop 1                        # number of Monte Carlo loops
num_conc 1d9                    # particle concentration (#/m^3)
kernel golovin                  # coagulation kernel

t_max 600                       # total simulation time (s)
del_t 60                        # timestep (s)
t_output 60                     # output interval (0 disables) (s)
t_state 0                       # state output interval (0 disables) (s)
t_progress 60                   # progress printing interval (0 disables) (s)

temp_profile temp_constant_15C.dat # temperature profile file
RH 0.999                        # initial relative humidity (1)
pressure 1d5                    # initial pressure (Pa)
rho_a 1.25                      # initial air density (kg/m^3)
latitude 40                     # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 0                    # start time (s since 00:00 UTC)
start_day 1                     # start day of year (UTC)

gas_init_conc gas_init_mosaic.dat # initial gas concentrations
aerosol_data aerosol_mosaic.dat # file containing aerosol data

n_init_dist 1                   # number of initial aerosol distributions

n_p 10                          # number of particles
vol_frac comp_mosaic.dat        # composition proportions of species
dist_type exp                   # type of distribution
dist_mean_vol 4.1886d-15        # mean volume (m^3)

n_bin 160                       # number of bins
v_min 1d-24                     # volume of smallest bin (m^3)
scal 3                          # scale factor (integer)

rand_init 17                    # random initialization (0 to auto-generate)
do_coagulation no               # whether to do coagulation (yes/no)
allow_double yes                # whether to allow doubling (yes/no)
do_condensation no              # whether to do condensation (yes/no)
do_mosaic yes                   # whether to do MOSAIC (yes/no)
do_restart no                   # whether to restart from stored state (yes/no)
restart_name XXXX.d             # filename to restart from
