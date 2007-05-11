run_type mc                     # Monte Carlo
output_file out/out_sedi_bidisperse_mc.d # name of output file
n_loop 1                        # number of Monte Carlo loops
num_conc 1d9                    # particle concentration (#/m^3)
kernel sedi                     # coagulation kernel

t_max 600                       # total simulation time (s)
del_t 1                         # timestep (s)
t_output 60                     # output interval (0 disables) (s)
t_state 10                      # state output interval (0 disables) (s)
t_progress 10                   # progress printing interval (0 disables) (s)

n_spec 1                        # number of species
i_water 1                       # species number that is water
rho 1000                        # density of species (kg/m^3)
nu 0                            # number of ions in solution of each species (1)
eps 0                           # solubility of species (1)
M_w 18d-3                       # molecular weight of species (kg/mole)

n_temps 1                       # number of temperature set-points
temp_times 0                    # times of temperature set-points
temps 288                       # temperatures at temperature set-points (K)
RH 0.999                        # initial relative humidity (1)
pressure 1d5                    # initial pressure (Pa)
rho_a 1.25                      # initial air density (kg/m^3)

n_init_dist 1                   # number of initial distributions
n_p 10001                       # number of particles
vol_frac 1                      # composition proportions of species
dist_type bidisperse            # type of distribution
dist_small_vol 4d-15            # volume of each small particle
dist_big_vol 4d-12              # initial volume of each big particle
dist_big_num 1                  # initial number of big particles

n_bin 250                       # number of bins
v_min 1d-24                     # volume of smallest bin (m^3)
scal 3                          # scale factor (integer)

rand_init 17                    # random initialization (0 to auto-generate)
do_coagulation yes              # whether to do coagulation (yes/no)
allow_double no                 # whether to allow doubling (yes/no)
do_condensation no              # whether to do condensation (yes/no)
do_restart no                   # whether to restart from stored state (yes/no)
restart_name XXXX.d             # filename to restart from
