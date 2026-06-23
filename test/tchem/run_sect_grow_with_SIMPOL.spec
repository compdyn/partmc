run_type sectional              # sectional code run
output_prefix out/tchem_grow_sect          # prefix of output files

t_max 600                       # total simulation time (s)
del_t 60                        # timestep (s)
t_output 60                     # output interval (0 disables) (s)
t_progress 600                  # progress printing interval (0 disables) (s)

do_tchem yes                    # whether to use TChem for chemistry
tchem_gas_config config_gas_cb05cl_ae5_with_SIMPOL.yaml
tchem_aero_config config_aero_cb05cl_ae5_with_SIMPOL.yaml
tchem_numerics_config solver_cb05cl_ae5_with_SIMPOL.yaml

n_bin 40                        # number of bins (== config "maximum computational particles")
d_min 1e-8                      # minimum diameter (m)
d_max 1e-6                      # maximum diameter (m)

gas_init gas_grow_init.dat                 # initial gas mixing ratios (SOA semivolatiles)

do_fractal no                   # whether to do fractal treatment
aerosol_init aero_grow_dist.dat # initial aerosol distribution (100 nm POA seed)

temp_profile temp_cb05cl_ae5.dat           # temperature profile file
pressure_profile pressure_cb05cl_ae5.dat   # pressure profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit_empty.dat      # gas emissions file
gas_background gas_back.dat     # background gas mixing ratios file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file
loss_function none              # particle loss function

rel_humidity 0.13916579011880265 # initial relative humidity (1)
latitude 40                     # latitude (degrees, -90 to 90)
longitude 0                     # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time 0                    # start time (s since 00:00 UTC)
start_day 1                     # start day of year (UTC)

do_coagulation no               # whether to do coagulation (yes/no)
