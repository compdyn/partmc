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
   prefactor = float(raw_input("Enter prefactor:"))
   exponent = float(raw_input("Enter exponent:"))

data1 = np.genfromtxt("rel_err_num.dat", skip_header=4)
data2 = np.genfromtxt("rel_err_mass.dat", skip_header=4)
data3 = np.genfromtxt("rel_err_total.dat", skip_header=4)
data4 = np.genfromtxt("rmse_num.dat", skip_header=4)
data5 = np.genfromtxt("rmse_mass.dat", skip_header=4)
data6 = np.genfromtxt("rmse_total.dat", skip_header=4)

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
if plt_opt == '1':
   axes.plot(list_exponent, list_max_err_num)
   axes.set_title("prefactor = %.3f" % (prefactor))
   axes.set_xlabel("exponent")
   axes.set_ylabel("maximum relative error of number distribution")
   filename_out = "plot_max_err_num_vs_exponent.pdf"
   figure.savefig(filename_out)
if para_opt == 'e':
   plt.plot(list_prefactor, list_max_err_num, 'bo')
   plt.title = "exponent = %.3f" % (exponent)
   plt.xlabel = "prefactor"
   plt.ylabel = "maximum relative error of number distribution"
   filename_out = "plot_max_err_num_vs_prefactor.pdf"
   plt.savefig(filename_out)
