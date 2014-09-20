run_type particle               # particle-resolved run
output_prefix out_aerodyne_0828/barrel_wc # prefix of output files
n_repeat 1                      # number of Monte Carlo repeats
n_part 10000                    # total number of particles
restart no                      # whether to restart from saved state (yes/no)

t_max 7800                      # total simulation time (s)
del_t 60                        # timestep (s)
t_output 120                    # output interval (0 disables) (s)
t_progress 60                   # progress printing interval (0 disables) (s)

gas_data gas_data.dat           # file containing gas data
gas_init gas_init.dat           # initial gas mixing ratios

aerosol_data aero_data.dat      # file containing aerosol data
aerosol_init aero_init_dist_sampled_aerodyne_0828.dat # aerosol initial condition file

temp_profile temp_aerodyne_0828.dat           # temperature profile file
pressure_profile pressure.dat   # pressure profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas mixing ratios file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back_aerodyne_0828.dat   # aerosol background file

rel_humidity 0.1                # initial relative humidity (1)
latitude 0                      # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 57660              # start time (s since 00:00 UTC)
start_day 195                    # start day of year (UTC)

do_coagulation yes              # whether to do coagulation (yes/no)
coag_kernel brown               # coagulation kernel
do_condensation no              # whether to do condensation (yes/no)
do_mosaic no                    # whether to do MOSAIC (yes/no)
do_nucleation no                # whether to do nucleation (yes/no)

rand_init 0                     # random initialization (0 to auto-generate)
allow_doubling yes              # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
record_removals no              # whether to record particle removals (yes/no)
do_parallel no                  # whether to run in parallel (yes/no)

do_chamber yes                 # whether to do chamber loss
V_chamber 0.2297                  # aerosol chamber volume (m^3)
A_diffuse 2.1206                   # diffusional deposition area (m^2)
A_sedi 0.2565                     # sedimentational deposition area (m^2)
prefactor_BL 0.05             # prefactor in diffusive boundary layer thickness (m)
exponent_BL 0.274              # exponent in diffusive boundary layer thickness

do_fractal yes                  # whether to do fractal treatment
frac_dim 2.3                   # fractal dimension
prime_radius 1e-8             # radius of monomer
vol_fill_factor 1.43            # volume filling factor
