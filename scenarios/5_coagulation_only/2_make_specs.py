import numpy as np
import os
import shutil

if not os.path.exists("spec"):
    os.mkdir("spec")
if not os.path.exists("temp"):
    os.mkdir("temp")

values = np.loadtxt('lhs.txt')
n_scenarios = (values.shape[0])
n_modes = int((values.shape[1]) / 3)
print(n_scenarios,n_modes)
for counter in range(n_scenarios):
    print("counter: %d" % counter)
    filename_in = "aero_init_dist_template.dat"
    directory = "./inputs/scenario_%04i" % counter
    if not os.path.exists(directory):
        os.makedirs(directory)
    filename_out = "%s/aero_init_dist.dat" % (directory)
    print("filename_out: %s" % filename_out)
    f_in = open(filename_in, 'r')
    f_out = open(filename_out, 'w')
    for line in f_in:
        for i_mode in range(n_modes):
            line = line.replace('NUM_CONC_%1i' % (i_mode+1), 
                 "%e" %(values[counter,2*n_modes + i_mode]))
            line = line.replace('GEOM_MEAN_DIAM_%1i' % (i_mode+1),
                 "%e" %(values[counter,i_mode]))
            line = line.replace('LOG10_GEO_STD_DEV_%1i' % (i_mode+1),
                 "%f" %(np.log10(values[counter,n_modes+i_mode])))
        f_out.write(line)


    f_in.close()
    f_out.close()
    # Copy files
    input_files = ["gas_back.dat", "temp.dat", "aero_back.dat",
                   "aero_emit.dat", "height.dat", "aero_back_comp.dat",
                   "aero_emit_comp_diesel.dat", "aero_data.dat",
                   "aero_emit_comp_pavedrd.dat", "aero_back_dist.dat",
                   "gas_emit.dat", "aero_emit_comp_gasol.dat", "pres.dat",
                   "gas_init.dat", "aero_emit_dist.dat", "aero_emit_comp_cooking.dat",
                   "aero_init_comp.dat", "gas_data.dat", "coag_brownian.spec", "run.sh"]
    for i_file, filename in enumerate(input_files):
        dest = shutil.copy("%s" %filename, "%s" % directory)
