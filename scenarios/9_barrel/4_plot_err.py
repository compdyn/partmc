#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import numpy as np
import mpl_helper
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt_opt = raw_input("Enter plot option (1 for 1d, 3 for 3d):")
prefactor = float(raw_input("Enter prefactor:"))
exponent = float(raw_input("Enter exponent:"))

data1 = np.genfromtxt("rel_err_num.dat", skip_header=4)
data2 = np.genfromtxt("rel_err_mass.dat", skip_header=4)
data3 = np.genfromtxt("rel_err_total.dat", skip_header=4)
data4 = np.genfromtxt("rmse_num.dat", skip_header=4)
data5 = np.genfromtxt("rmse_mass.dat", skip_header=4)
data6 = np.genfromtxt("rmse_total.dat", skip_header=4)

list_prefactor = []
list_exponent = []
list_rmse_num_vs_pref = []
list_rmse_mass_vs_pref = []
list_rmse_total_vs_pref = []
list_rmse_num_vs_expo = []
list_rmse_mass_vs_expo = []
list_rmse_total_vs_expo = []
delta_t = 7

for row in range(0, data1.shape[0]):
    data1_1d = data1[row, :]
    data2_1d = data2[row, :]
    data3_1d = data3[row, :]
    if data1_1d[1] == prefactor and data1_1d[2] == exponent:
       array_rmse_num_vs_time = data1_1d[3:]
       array_rmse_mass_vs_time = data2_1d[3:]
       array_rmse_total_vs_time = data3_1d[3:]
       time = range(0,len(data1_1d[3:]))
       time = np.array(time) * delta_t

for row in range(0, data4.shape[0]):
    data4_1d = data4[row, :]
    data5_1d = data5[row, :]
    data6_1d = data6[row, :]
    if data4_1d[1] == prefactor:
       list_exponent.append(data4_1d[2])
       list_rmse_num_vs_expo.append(data4_1d[3])
       list_rmse_mass_vs_expo.append(data5_1d[3])
       list_rmse_total_vs_expo.append(data6_1d[3])
    if data4_1d[2] == exponent:
       list_prefactor.append(data4_1d[1])
       list_rmse_num_vs_pref.append(data4_1d[3])
       list_rmse_mass_vs_pref.append(data5_1d[3])
       list_rmse_total_vs_pref.append(data6_1d[3])

# sort total rmse data
data6_sorted = data6[np.argsort(data6[:,3])]
print "Mininum total error = %.4f, at prefactor = %.3f and exponent = %.3f "\
      "with case %04d" %(data6_sorted[0,3], data6_sorted[0,1], \
      data6_sorted[0,2], data6_sorted[0,0])
filename_out = "rmse_total_sorted.dat"
f_out = open(filename_out, 'w')
f_out.write("# Colume 1: caseID\n")
f_out.write("# Colume 2: prefactor\n")
f_out.write("# Colume 3: exponent\n")
f_out.write("# Colume 4: total root mean square error (num + mass)\n")
for row in range(0,data6_sorted.shape[0]):
    f_out.write("%04d   %.3f   %.3f    %.4f\n" % (data6_sorted[row,0], \
    data6_sorted[row,1], data6_sorted[row,2], data6_sorted[row,3]))
f_out.close()

if plt_opt == '1':
   # plot rmse vs. exponent
   (figure, axes) = mpl_helper.make_fig(colorbar=False)
   axes.plot(list_exponent, list_rmse_num_vs_expo, marker='^')
   axes.plot(list_exponent, list_rmse_mass_vs_expo, marker='D')
   axes.plot(list_exponent, list_rmse_total_vs_expo, marker='o')
   axes.set_title("prefactor = %.3f" % (prefactor))
   axes.set_xlabel("exponent")
   axes.set_ylabel("root mean square error")
   axes.grid()
   box = axes.get_position()
   axes.set_position([box.x0, box.y0, box.width * 0.78, box.height])
   axes.legend(('Num', 'Mass', 'Total'), loc='center left', bbox_to_anchor=(1, 0.5))
   filename_out = "plot_rmse_vs_exponent.pdf"
   figure.savefig(filename_out)
   # plot rmse vs. prefactor
   (figure, axes) = mpl_helper.make_fig(colorbar=False)
   axes.plot(list_prefactor, list_rmse_num_vs_pref, marker='^')
   axes.plot(list_prefactor, list_rmse_mass_vs_pref, marker='D')
   axes.plot(list_prefactor, list_rmse_total_vs_pref, marker='o')
   axes.set_title("exponent = %.3f" % (exponent))
   axes.set_xlabel("prefactor")
   axes.set_ylabel("root mean square error")
   axes.grid()
   box = axes.get_position()
   axes.set_position([box.x0, box.y0, box.width * 0.78, box.height])
   axes.legend(('Num', 'Mass', 'Total'), loc='center left', bbox_to_anchor=(1, 0.5))
   filename_out = "plot_rmse_vs_prefactor.pdf"
   figure.savefig(filename_out)
   # plot rmse vs. time
   (figure, axes) = mpl_helper.make_fig(colorbar=False)
   axes.plot(time, array_rmse_num_vs_time, marker='^')
   axes.plot(time, array_rmse_mass_vs_time, marker='D')
   axes.plot(time, array_rmse_total_vs_time, marker='o')
   axes.set_title("prefactor = %.3f, exponent = %.3f" % (prefactor, exponent))
   axes.set_xlabel("time (min)")
   axes.set_ylabel("relative error")
   axes.grid()
   box = axes.get_position()
   axes.set_position([box.x0, box.y0, box.width * 0.78, box.height])
   axes.legend(('Num', 'Mass', 'Total'), loc='center left', bbox_to_anchor=(1, 0.5))
   filename_out = "plot_rmse_vs_time.pdf"
   figure.savefig(filename_out)
if plt_opt == '3':
   fig = plt.figure()
   ax = fig.add_subplot(111, projection='3d')
   ax.scatter(data6[:,1], data6[:,2], data6[:,3], marker='o')
   ax.set_xlabel('prefactor')
   ax.set_ylabel('exponent')
   ax.set_zlabel('total root mean square error')
   fig.savefig("plot_rmse_3d.pdf")
