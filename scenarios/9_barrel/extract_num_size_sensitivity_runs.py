import sys
sys.path.append('../../tool/')
import partmc
import scipy.io
import os
import numpy
import math
import mpl_helper
import matplotlib.pyplot as plt

col = 11
dataset_name = '0925'
case = 415

ref_data = numpy.loadtxt("ref_"+dataset_name+"/ref_aero_size_num_regrid.txt")

list_mean = []
list_std = []
for row in range(0,ref_data.shape[0]):
    list_data_temp_1d = []
    for i in numpy.arange(1,11,1):
        file_load = "out_"+dataset_name+"/case_%04d_wc_%04d_aero_size_num.txt" % (case, i)
        data_temp = numpy.loadtxt(file_load)
        list_data_temp_1d.append(data_temp[row,col])
    data_temp_1d = numpy.array(list_data_temp_1d)
    data_temp_1d *= math.log(10) # convert to d*_dlogDp
    list_mean.append(numpy.mean(data_temp_1d))
    diff = data_temp_1d - numpy.mean(data_temp_1d)
    list_std.append(numpy.sqrt(1. / float(len(data_temp_1d)-1) * sum(diff**2)))
data1_1d = numpy.array(list_mean)
partmc_std = numpy.array(list_std)

diams_list = []
for i in range(0,ref_data.shape[0]):
    diams_list.append(ref_data[i,0])
diams = numpy.array(diams_list)

data_only_coag = numpy.loadtxt("out_0925/only_coag_0001_aero_size_num.txt")
data_only_coag_1d = data_only_coag[:,col] * math.log(10)
data_only_wall = numpy.loadtxt("out_0925/only_wall_0001_aero_size_num.txt")
data_only_wall_1d = data_only_wall[:,col] * math.log(10)

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.semilogx(diams, data_only_coag_1d, color='r')
axes.semilogx(diams, data_only_wall_1d, color='b')
axes.semilogx(diams, data1_1d, color='k')
axes.set_title("")
axes.set_xlabel("Dry diameter (m)", fontsize=14)
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)", fontsize=14)
axes.grid()
axes.set_ylim(0,)
axes.legend(('Only coag.','Only wall','Coag.+wall'),loc='upper left')
bbox_props = dict(boxstyle="square,pad=0.3", fc="cyan", ec="b", lw=1)
axes.annotate('70 min', xy=(0.82, 0.85), xycoords='axes fraction',weight='extra bold', size=14, bbox=bbox_props)
bbox_props_1 = dict(boxstyle="square,pad=0.3", fc="white", ec="r", lw=1)
axes.annotate('Exp.4', xy=(0.02, 0.15), xycoords='axes fraction',weight='extra bold', size=14, bbox=bbox_props_1)
filename_out = "aero_num_size.pdf"
figure.savefig(filename_out)
