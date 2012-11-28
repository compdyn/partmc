#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import numpy as np
import mpl_helper
import matplotlib.pyplot as plt

prefactor = 0.0
exponent = 0.0
plt_opt = raw_input("Enter plot option (1 for 1d, 2 for 2d):")
if plt_opt == '1':
   para_opt = raw_input("Enter fixed parameter to plot (p for prefactor, e for exponent):")
   if para_opt == 'p':
      prefactor = float(raw_input("Enter prefactor to plot:"))
   if para_opt == 'e':
      exponent = float(raw_input("Enter exponent to plot:"))

data1 = np.loadtxt("rel_err_num.dat")
data2 = np.loadtxt("rel_err_mass.dat")
data3 = np.loadtxt("rel_err_total.dat")

filename_out_num = "max_rel_err_num.dat"
f_out_num = open(filename_out_num, 'w')
f_out_num.write("# Colume 1: prefactor\n")
f_out_num.write("# Colume 2: exponent\n")
f_out_num.write("# Colume 3: maximum relative error of number distribution\n")
filename_out_mass = "max_rel_err_mass.dat"
f_out_mass = open(filename_out_mass, 'w')
f_out_mass.write("# Colume 1: prefactor\n")
f_out_mass.write("# Colume 2: exponent\n")
f_out_mass.write("# Colume 3: maximum relative error of mass distribution\n")
filename_out_total = "max_rel_err_total.dat"
f_out_total = open(filename_out_total, 'w')
f_out_total.write("# Colume 1: prefactor\n")
f_out_total.write("# Colume 2: exponent\n")
f_out_total.write("# Colume 3: total maximum relative error\n")

for row in range(0, data1.shape[0]):
    data_1d = data1[row, :]
    f_out_num.write("%.3f     " % data_1d[0])
    f_out_num.write("%.3f     " % data_1d[1])
    f_out_num.write("%.3f\n" % max(data_1d[2:]))
    data_1d = data2[row, :]
    f_out_mass.write("%.3f     " % data_1d[0])
    f_out_mass.write("%.3f     " % data_1d[1])
    f_out_mass.write("%.3f\n" % max(data_1d[2:]))
    data_1d = data3[row, :]
    f_out_total.write("%.3f     " % data_1d[0])
    f_out_total.write("%.3f     " % data_1d[1])
    f_out_total.write("%.3f\n" % max(data_1d[2:]))

# make plots
list_prefactor = []
list_exponent = []
list_max_err_num = []
list_max_err_mass = []
list_max_err_total = []

# make plot for max error of number distribution
for row in range(0, data1.shape[0]):
    data_1d = data1[row, :]
    if data_1d[0] == prefactor:
       list_exponent.append(data_1d[1])
       list_max_err_num.append(max(data_1d[2:]))
    if data_1d[1] == exponent:
       list_prefactor.append(data_1d[0])
       list_max_err_num.append(max(data_1d[2:]))

(figure, axes) = mpl_helper.make_fig(colorbar=False)
if para_opt == 'p':
   plt.plot(list_exponent, list_max_err_num, 'bo')
   plt.title = "prefactor = %.3f" % (prefactor)
   plt.xlabel = "exponent"
   plt.ylabel = "maximum relative error of number distribution"
   filename_out = "plot_max_err_num_vs_exponent.pdf"
   plt.savefig(filename_out)
   str_command = "open "+filename_out
   os.system(str_command)
if para_opt == 'e':
   plt.plot(list_prefactor, list_max_err_num, 'bo')
   plt.title = "exponent = %.3f" % (exponent)
   plt.xlabel = "prefactor"
   plt.ylabel = "maximum relative error of number distribution"
   filename_out = "plot_max_err_num_vs_prefactor.pdf"
   plt.savefig(filename_out)
   str_command = "open "+filename_out
   os.system(str_command)


f_out_num.close()
f_out_mass.close()
f_out_total.close()
