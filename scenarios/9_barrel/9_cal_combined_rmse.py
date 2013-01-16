#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import shutil
from numpy import *

datasets = ['0212','0909','0925']
data1 = genfromtxt("rel_err_total_"+datasets[0]+".dat", skip_header=5)
data2 = genfromtxt("rel_err_total_"+datasets[1]+".dat", skip_header=5)
data3 = genfromtxt("rel_err_total_"+datasets[2]+".dat", skip_header=5)
filename_out_rmse_total = "rmse_total.dat"
f_out_rmse_total = open(filename_out_rmse_total, 'w')
f_out_rmse_total.write("# Colume 1: caseID\n")
f_out_rmse_total.write("# Colume 2: prefactor\n")
f_out_rmse_total.write("# Colume 3: exponent\n")
f_out_rmse_total.write("# Colume 4: fractal dimension\n")
f_out_rmse_total.write("# Colume 5: total root mean square error\n")

for row in range(0,data1.shape[0]):
    f_out_rmse_total.write("%04d   %.3f   %.2f   %.1f   " % (data1[row,0], data1[row,1], data1[row,2], data1[row,3]))
    list_total_err = []
    for col in range(4,data1.shape[1]):
        list_total_err.append(data1[row,col])
    for col in range(4,data2.shape[1]):
        list_total_err.append(data2[row,col])
    for col in range(4,data3.shape[1]):
        list_total_err.append(data3[row,col])
    array_total_err = array(list_total_err)
    rmse_total = sqrt(sum(array_total_err**2)/len(array_total_err))
    f_out_rmse_total.write("%.4f\n" % (rmse_total))

f_out_rmse_total.close()

# sort total rmse data
data6 = genfromtxt("rmse_total.dat", skip_header=5)
data6_sorted = data6[argsort(data6[:,4])]
print "Mininum total error = %.4f, at prefactor = %.3f and exponent = %.2f "\
      "and frac_dim = %.1f with case %04d" %(data6_sorted[0,4], data6_sorted[0,1], \
      data6_sorted[0,2], data6_sorted[0,3], data6_sorted[0,0])
filename_out = "rmse_total_sorted.dat"
f_out_sort = open(filename_out, 'w')
f_out_sort.write("# Colume 1: caseID\n")
f_out_sort.write("# Colume 2: prefactor\n")
f_out_sort.write("# Colume 3: exponent\n")
f_out_sort.write("# Colume 4: fractal dimension\n")
f_out_sort.write("# Colume 5: total root mean square error (num + mass)\n")
for row in range(0,data6_sorted.shape[0]):
    f_out_sort.write("%04d   %.3f   %.2f    %.1f    %.4f\n" % (data6_sorted[row,0], \
    data6_sorted[row,1], data6_sorted[row,2], data6_sorted[row,3], data6_sorted[row,4]))
f_out_sort.close()
