import sys
sys.path.append('../../tool/')
import partmc
import scipy.io
import os
import numpy
import math
import mpl_helper
import matplotlib.pyplot as plt

col = 21
dataset_name = '0925'

ref_data = numpy.loadtxt("ref_"+dataset_name+"/ref_aero_size_num.txt")
partmc_data = numpy.loadtxt("out/barrel_wc_0001_aero_size_num_R0_15.txt")

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.semilogx(partmc_data[:,0], partmc_data[:,col] * math.log(10), color='k')
axes.semilogx(ref_data[:,0],ref_data[:,col], color='#CC4F1B')
axes.set_title("")
axes.set_xlabel("Dry diameter (m)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.set_ylim(0,)
axes.legend(('PartMC','Barrel'),loc='upper left')
filename_out = "aero_num_size.pdf"
figure.savefig(filename_out)
