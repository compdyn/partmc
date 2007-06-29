run_type mc                     # Monte Carlo run
output_file out/dust_salt_part2_summary.d # name of output file
state_prefix out/dust_salt_part2_state # prefix of state files
n_loop 1                        # number of Monte Carlo loops
num_conc 1d9                    # particle concentration (#/m^3)
kernel sedi                     # coagulation kernel

t_max 800                       # total simulation time (s)
del_t 1                         # timestep (s)
t_output 100                    # output interval (0 disables) (s)
t_state 800                     # state output interval (0 disables) (s)
t_progress 1                    # progress printing interval (0 disables) (s)

temp_profile temp_cooling.dat   # temperature profile file
RH 0.999                        # initial relative humidity (1)
pressure 1d5                    # initial pressure (Pa)
rho_a 1.25                      # initial air density (kg/m^3)
latitude 40                     # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 0                    # start time (s since 00:00 UTC)
start_day 1                     # start day of year (UTC)

gas_init_conc gas_init_standard.dat # initial gas concentrations
aerosol_data aerosol_dust_salt.dat # file containing aerosol data

n_init_dist 2                   # number of initial distributions

n_p 5000                        # number of particles
vol_frac comp_salt.dat          # composition proportions of species
dist_type log_normal            # type of distribution
dist_mean_diam 0.266d-6         # mean diameter (m)
dist_std_dev 0.21               # standard deviation (m)

n_p 5000                        # number of particles
vol_frac comp_dust.dat          # composition proportions of species
dist_type log_normal            # type of distribution
dist_mean_diam 0.05d-6          # mean diameter (m)
dist_std_dev 0.6                # standard deviation (m)

n_bin 160                       # number of bins
v_min 1d-24                     # volume of smallest bin (m^3)
scal 3                          # scale factor (integer)

rand_init 17                    # random initialization (0 to 
do_coagulation yes              # whether to do coagulation (yes/no)
allow_double yes                # double when particle number is small (yes/no)
do_condensation yes             # whether to do condensation (yes/no)
do_mosaic no                    # whether to do MOSAIC (yes/no)
do_restart yes                  # whether to restart from stored state (yes/no)
restart_name dust_salt_part1_state_0001_00000800.d  # filename to restart from
