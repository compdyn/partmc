import sys
sys.path.append('../../tool/')
import partmc
import scipy.io
import os
import numpy
import math
import mpl_helper
import matplotlib.pyplot as plt

"""
dataset_name = '0925'
prefactor = 0.06
frac_dim = 2.2

dataset_name = '0909'
prefactor = 0.04
frac_dim = 2.4

dataset_name = '0908'
prefactor = 0.07
frac_dim = 2.2

dataset_name = '0322'
prefactor = 0.05
frac_dim = 3.0
"""

prefactor = 0.06
frac_dim = 2.3
#ref_data = numpy.loadtxt("rmse_num_"+dataset_name+".dat")
ref_data = numpy.loadtxt("combined_rmse.txt")
list_a = []
list_rmse = []
for i in range(0,ref_data.shape[0]):
    if (ref_data[i,1] == prefactor and ref_data[i,3] == frac_dim): 
       list_a.append(ref_data[i,2])
       list_rmse.append(ref_data[i,4])

list_a, list_rmse = zip(*sorted(zip(list_a, list_rmse)))
print list_a
print list_rmse

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(list_a, list_rmse, '-ko')
axes.set_xlabel(r"Exponent $a$")
axes.set_ylabel(r"Root mean square error")
axes.set_title(r"prefactor $k_{\rm d}$ = %.3f, fractal dimension $d_{\rm f}$ = %.1f" %(prefactor, frac_dim))
axes.grid()
bbox_props_1 = dict(boxstyle="square,pad=0.3", fc="white", ec="r", lw=1)
#axes.annotate('Exp.1', xy=(0.85, 0.05), xycoords='axes fraction',weight='extra bold', size=14, bbox=bbox_props_1)
filename_out = "rmse_vs_a.pdf"
figure.savefig(filename_out)
