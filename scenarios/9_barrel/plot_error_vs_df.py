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
exponent = 0.26

dataset_name = '0909'
prefactor = 0.04
exponent = 0.25

dataset_name = '0908'
prefactor = 0.07
exponent = 0.22

dataset_name = '0322'
prefactor = 0.05
exponent = 0.22
"""

prefactor = 0.06
exponent = 0.26
#ref_data = numpy.loadtxt("rmse_num_"+dataset_name+".dat")
ref_data = numpy.loadtxt("combined_rmse.txt")
list_df = []
list_rmse = []
for i in range(0,ref_data.shape[0]):
    if (ref_data[i,1] == prefactor and ref_data[i,2] == exponent): 
       list_df.append(ref_data[i,3])
       list_rmse.append(ref_data[i,4])

list_df, list_rmse = zip(*sorted(zip(list_df, list_rmse)))

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(list_df, list_rmse, '-ro')
axes.set_xlabel(r"Fractal dimension $d_{\rm f}$")
axes.set_ylabel(r"Root mean square error")
axes.set_title(r"prefactor $k_{\rm d}$ = %.3f, exponent $a$ = %.2f" %(prefactor, exponent))
axes.grid()
bbox_props_1 = dict(boxstyle="square,pad=0.3", fc="white", ec="r", lw=1)
#axes.annotate('Exp.1', xy=(0.85, 0.05), xycoords='axes fraction',weight='extra bold', size=14, bbox=bbox_props_1)
filename_out = "rmse_vs_df.pdf"
figure.savefig(filename_out)
