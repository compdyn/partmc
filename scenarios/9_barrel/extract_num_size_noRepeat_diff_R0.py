import sys
sys.path.append('../../tool/')
import partmc
import scipy.io
import os
import numpy
import math
import mpl_helper
import matplotlib.pyplot as plt

col = 41
dataset_name = '0925'

partmc_data_1 = numpy.loadtxt("out/barrel_wc_0001_aero_size_num_R0_15.txt")
partmc_data_2 = numpy.loadtxt("out/barrel_wc_0001_aero_size_num_R0_45.txt")
partmc_data_3 = numpy.loadtxt("out/barrel_wc_0001_aero_size_num_R0_80.txt")

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.semilogx(partmc_data_1[:,0], partmc_data_1[:,col] * math.log(10), color='g')
axes.semilogx(partmc_data_2[:,0], partmc_data_2[:,col] * math.log(10), color='b')
axes.semilogx(partmc_data_3[:,0], partmc_data_3[:,col] * math.log(10), color='r')
axes.set_title("")
axes.set_xlabel("Dry diameter (m)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.set_ylim(0,)
axes.legend(('$R_0$ = 15 nm','$R_0$ = 45 nm','$R_0$ = 80 nm'),loc='upper left')
filename_out = "aero_num_size.pdf"
figure.savefig(filename_out)
