from __future__ import division
import sys
sys.path.append('../../tool/')
import partmc
import scipy.io
import os
import numpy
import math
import mpl_helper
import matplotlib.pyplot as plt

dataset_name = '0909'
case = 1

ref_data = numpy.loadtxt("ref_"+dataset_name+"/ref_aero_size_num_regrid.txt")

list_mean = []
list_std = []
list_time = []
file_load = "out_barrel_sep_fit/out_"+dataset_name+"/case_%04d_wc_%04d_aero_size_num.txt" % (case, 1)
data_temp = numpy.loadtxt(file_load)
diams = data_temp[:,0]
for col in range(1,data_temp.shape[1]):
    list_time.append(7*(col-1))
    list_data_temp_1d = []
    for row in range(0,data_temp.shape[0]):
        list_data_temp_1d.append(data_temp[row,col]*diams[row])
    list_mean.append(sum(list_data_temp_1d) / sum(data_temp[:,col]))
for col in range(1,data_temp.shape[1]):
    list_data_temp_1d = []
    for row in range(0,data_temp.shape[0]):
        list_data_temp_1d.append(data_temp[row,col]*(diams[row]-list_mean[col-1])**2)
    list_std.append(numpy.sqrt(sum(list_data_temp_1d) / sum(data_temp[:,col])))
list_skewness = []
for col in range(1,data_temp.shape[1]):
    list_data_temp_1d = []
    temp_2 = []
    for row in range(0,data_temp.shape[0]):
        list_data_temp_1d.append((diams[row]-list_mean[col-1])**3)
        temp_2.append((diams[row]-list_mean[col-1])**2)
    list_skewness.append(sum(list_data_temp_1d) / data_temp.shape[0] / (sum(temp_2) / (data_temp.shape[0]-1))**(3/2))

list_mean_barrel = []
list_std_barrel = []
diams_barrel = ref_data[:,0]
for col in range(1,ref_data.shape[1]):
    list_data_temp_1d = []
    for row in range(0,ref_data.shape[0]):
        list_data_temp_1d.append(ref_data[row,col]*diams_barrel[row])
    list_mean_barrel.append(sum(list_data_temp_1d) / sum(ref_data[:,col]))
for col in range(1,ref_data.shape[1]):
    list_data_temp_1d = []
    for row in range(0,ref_data.shape[0]):
        list_data_temp_1d.append(ref_data[row,col]*(diams_barrel[row]-list_mean_barrel[col-1])**2)
    list_std_barrel.append(numpy.sqrt(sum(list_data_temp_1d) / sum(ref_data[:,col])))
list_skewness_barrel = []
for col in range(1,ref_data.shape[1]):
    list_data_temp_1d = []
    temp_2 = []
    for row in range(0,ref_data.shape[0]):
        list_data_temp_1d.append((diams_barrel[row]-list_mean_barrel[col-1])**3)
        temp_2.append((diams_barrel[row]-list_mean_barrel[col-1])**2)
    list_skewness_barrel.append(sum(list_data_temp_1d) / ref_data.shape[0] / (sum(temp_2) / (ref_data.shape[0]-1))**(3/2))

print "Initial D_mean = ", list_mean_barrel[0]
print "Initial sigma = ", list_std_barrel[0]

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(list_time, list_mean, color='k')
axes.plot(list_time, list_mean_barrel, color='#CC4F1B')
axes.set_title("")
axes.set_xlabel("Time (min)", fontsize=14)
axes.set_ylabel(r"Mean particle diameter $\bar{D}_{\rm p}$ (m)", fontsize=14)
axes.grid()
axes.legend(('PartMC','Barrel'),loc='upper left')
bbox_props_1 = dict(boxstyle="square,pad=0.3", fc="white", ec="r", lw=1)
axes.annotate('Exp.3', xy=(0.1, 0.1), xycoords='axes fraction',weight='extra bold', size=14, bbox=bbox_props_1)
filename_out = "aero_dmean_time.pdf"
figure.savefig(filename_out)

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(list_time, list_std, color='k')
axes.plot(list_time, list_std_barrel, color='#CC4F1B')
axes.set_title("")
axes.set_xlabel("Time (min)", fontsize=14)
axes.set_ylabel(r"Standard deviation of size $\sigma$ (m)", fontsize=14)
axes.grid()
axes.legend(('PartMC','Barrel'),loc='upper left')
bbox_props_1 = dict(boxstyle="square,pad=0.3", fc="white", ec="r", lw=1)
axes.annotate('Exp.3', xy=(0.1, 0.1), xycoords='axes fraction',weight='extra bold', size=14, bbox=bbox_props_1)
filename_out = "aero_sigma_time.pdf"
figure.savefig(filename_out)

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(list_time, list_skewness, color='k')
axes.plot(list_time, list_skewness_barrel, color='#CC4F1B')
axes.set_title("")
axes.set_xlabel("Time (min)", fontsize=14)
axes.set_ylabel(r"Skewness of size distribution", fontsize=14)
axes.grid()
axes.legend(('PartMC','Barrel'),loc='upper right')
bbox_props_1 = dict(boxstyle="square,pad=0.3", fc="white", ec="r", lw=1)
axes.annotate('Exp.3', xy=(0.1, 0.1), xycoords='axes fraction',weight='extra bold', size=14, bbox=bbox_props_1)
filename_out = "aero_skewness_time.pdf"
figure.savefig(filename_out)
