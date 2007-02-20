run_type sect                   # sectional code
output_name sedi_exp_sect       # name of output files
num_conc 1d9                    # particle concentration (#/m^3)
kernel sedi                     # coagulation kernel

dist_type exp                   # type of initial distribution
dist_mean_vol 4.1886d-15        # mean diameter (m)

t_max 600                       # total simulation time (s)
del_t 1                         # timestep (s)
t_output 60                     # output interval (0 disables) (s)
t_progress 60                   # progress printing interval (0 disables) (s)

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

n_bin 220                       # number of bins
v_min 1d-24                     # volume of smallest bin (m^3)
scal 4                          # scale factor (integer)
