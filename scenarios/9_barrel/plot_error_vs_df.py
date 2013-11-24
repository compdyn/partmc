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
exponent = 0.26

ref_data = numpy.loadtxt("rmse_num_"+dataset_name+".dat")
list_df = []
list_rmse = []
for i in range(0,ref_data.shape[0]):
    if (ref_data[i,1] == prefactor and ref_data[i,2] == exponent): 
       list_df.append(ref_data[i,3])
       list_rmse.append(ref_data[i,4])

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(list_df, list_rmse, '-ro')
axes.set_xlabel("Fractal dimension")
axes.set_ylabel(r"Root mean square error")
axes.set_title("prefactor = %.3f, exponent = %.2f" %(prefactor, exponent))
axes.grid()
filename_out = "rmse_vs_df.pdf"
figure.savefig(filename_out)
