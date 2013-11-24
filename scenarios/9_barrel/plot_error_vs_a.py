import sys
sys.path.append('../../tool/')
import partmc
import scipy.io
import os
import numpy
import math
import mpl_helper
import matplotlib.pyplot as plt

dataset_name = '0925'
prefactor = 0.06
frac_dim = 2.2

ref_data = numpy.loadtxt("rmse_num_"+dataset_name+".dat")
list_a = []
list_rmse = []
for i in range(0,ref_data.shape[0]):
    if (ref_data[i,1] == prefactor and ref_data[i,3] == frac_dim): 
       list_a.append(ref_data[i,2])
       list_rmse.append(ref_data[i,4])

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(list_a, list_rmse, '-ko')
axes.set_xlabel("Exponent")
axes.set_ylabel(r"Root mean square error")
axes.set_title("prefactor = %.3f, fractal dimension = %.1f" %(prefactor, frac_dim))
axes.grid()
filename_out = "rmse_vs_a.pdf"
figure.savefig(filename_out)
