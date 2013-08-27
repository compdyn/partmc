import sys
sys.path.append('../../tool/')
import partmc
import scipy.io
import os
import numpy
import math
import mpl_helper
import matplotlib.pyplot as plt

col = 66
nbin = 100
dataset_name = 'aerodyne_0828'

ref_data = numpy.loadtxt("ref_"+dataset_name+"/ref_aero_size_num.txt")
partmc_data = numpy.loadtxt("out_"+dataset_name+"/barrel_wc_0001_aero_size_num.txt")

list_mean = []
list_std = []
for row in range(0,nbin):
    list_data_temp_1d = []
    for i in numpy.arange(1,11,1):
        file_load = "out_"+dataset_name+"/barrel_wc_%04d_aero_size_num.txt" % (i)
        data_temp = numpy.loadtxt(file_load)
        list_data_temp_1d.append(data_temp[row,col])
    data_temp_1d = numpy.array(list_data_temp_1d)
    data_temp_1d *= math.log(10) # convert to d*_dlogDp
    list_mean.append(numpy.mean(data_temp_1d))
    diff = data_temp_1d - numpy.mean(data_temp_1d)
    list_std.append(numpy.sqrt(1. / float(len(data_temp_1d)-1) * sum(diff**2)))
data1_1d = numpy.array(list_mean)
partmc_std = numpy.array(list_std)

# calculate the relative error
ref_data_err_list = []
for i in range(0,nbin):
    ref_data_err_list.append(3. * numpy.sqrt(partmc_std[i]**2/10.)) # use 3*sigma
ref_data_err = numpy.array(ref_data_err_list)

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.semilogx(partmc_data[:,0], data1_1d, color='k')
#axes.errorbar(ref_data[:,0],ref_data[:,col],yerr=ref_data_err,color='r')
axes.semilogx(ref_data[:,0],ref_data[:,col], color='#CC4F1B')
axes.fill_between(partmc_data[:,0], data1_1d-ref_data_err, data1_1d+ref_data_err,
    alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
axes.set_title("")
axes.set_xlabel("Dry diameter (m)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.set_ylim(0,)
axes.legend(('PartMC','Barrel'),loc='upper left')
filename_out = "aero_num_size.pdf"
figure.savefig(filename_out)
