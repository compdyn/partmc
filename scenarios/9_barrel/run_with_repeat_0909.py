#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import shutil
from numpy import *

def read_values_from_spec_file(filename_in, str_wanted):
    A = []
    f_in = open(filename_in, 'r')
    for line in f_in:
        s =  line.rsplit(' ')
        if (s[0] == str_wanted):
           A.append(s[1]) # add the value to the list
    f_in.close()
    # Convert strings to an array of values (in this case floats)
    val = zeros(len(A))
    for ii in range(0,len(A)):
        val[ii] =  float(A[ii])
    return val

dataset = ['0909']
t_max = [19740]
t_output = [420]
rh = [0.1275]
index = 0

for dataset_name in dataset:
    if not os.path.exists("spec_"+dataset_name):
        os.mkdir("spec_"+dataset_name)
    else:
        shutil.rmtree("spec_"+dataset_name)
        os.mkdir("spec_"+dataset_name)

    if not os.path.exists("out_"+dataset_name):
        os.mkdir("out_"+dataset_name)

    #prefactor = float(raw_input("Enter prefactor:"))
    #exponent = float(raw_input("Enter exponent:"))
    prime_radius = 1e-8
    vol_fill_factor = 1.43
    rho = 1760
    case = 0
    for prefactor in arange(0.025,0.065,0.005):
        for exponent in arange(0.22,0.27,0.01):
            for frac_dim in arange(2.2,3.1,0.1):
                filename_in = "barrel_template_with_repeat.spec"
                case += 1
                filename_out = "spec_"+dataset_name+"/barrel_wc_case_%04d.spec" % (case)
                f_in = open(filename_in, 'r')
                f_out = open(filename_out, 'w')

                for line in f_in:
                    line = line.replace('%%OUTPUT_PREFIX%%', 'out_'+dataset_name+'/case_%04d_wc' % (case))
                    line = line.replace('%%t_max%%', str(t_max[index]))
                    line = line.replace('%%t_output%%', str(t_output[index]))
                    line = line.replace('%%aero_init_dist_sampled%%', 'aero_init_dist_sampled_'+dataset_name+'.dat')
                    line = line.replace('%%temp_profile%%', 'temp_'+dataset_name+'.dat')
                    line = line.replace('%%aero_background%%', 'aero_back_'+dataset_name+'.dat')
                    line = line.replace('%%rel_humidity%%', str(rh[index]))
                    line = line.replace('%%prefactor%%', str(prefactor))
                    line = line.replace('%%exponent%%', str(exponent))
                    line = line.replace('%%frac_dim%%', str(frac_dim))
                    line = line.replace('%%prime_radius%%', str(prime_radius))
                    line = line.replace('%%vol_fill_factor%%', str(vol_fill_factor))
                    f_out.write(line)

                f_in.close()
                f_out.close()

    str_exec = "../../build/partmc"
    str_extr_size = "../../build/extract_aero_mob_size"
    str_extr_time = "../../build/extract_aero_time"

    # Open file to store the errors
    filename_out_num = "rel_err_num_"+dataset_name+".dat"
    f_out_num = open(filename_out_num, 'w')
    f_out_num.write("# Colume 1: caseID\n")
    f_out_num.write("# Colume 2: prefactor\n")
    f_out_num.write("# Colume 3: exponent\n")
    f_out_num.write("# Colume 4: fractal dimension\n")
    f_out_num.write("# Colume j+4: relative error of number distribution at time(j)\n")
    filename_out_rmse_num = "rmse_num_"+dataset_name+".dat"
    f_out_rmse_num = open(filename_out_rmse_num, 'w')
    f_out_rmse_num.write("# Colume 1: caseID\n")
    f_out_rmse_num.write("# Colume 2: prefactor\n")
    f_out_rmse_num.write("# Colume 3: exponent\n")
    f_out_rmse_num.write("# Colume 4: fractal dimension\n")
    f_out_rmse_num.write("# Colume 5: root mean square error of number distribution\n")

    # run simulations, store the errors
    path="spec_"+dataset_name+"/"
    dirList=os.listdir(path)
    case = 0
    for fname in dirList:

        case += 1
        spec_file = "spec_"+dataset_name+"/barrel_wc_case_%04d.spec" % (case)

        prefactor = read_values_from_spec_file(spec_file, 'prefactor_BL')
        exponent = read_values_from_spec_file(spec_file, 'exponent_BL')
        frac_dim = read_values_from_spec_file(spec_file, 'frac_dim')

        command_1 = [str_exec,spec_file]
        print command_1
        subprocess.check_call(command_1)

        for i in arange(1,11,1):
            ncfile_prefix = "out_"+dataset_name+"/case_%04d_wc_%04d" % (case, i)

            command_2 = [str_extr_size,"--num","--dmin","1e-8","--dmax","1e-6","--nbin","100",ncfile_prefix]
            command_3 = [str_extr_size,"--mass","--dmin","1e-8","--dmax","1e-6","--nbin","100",ncfile_prefix]
            command_4 = [str_extr_time,ncfile_prefix]

            print command_2
            subprocess.check_call(command_2)

            print command_3
            subprocess.check_call(command_3)

            print command_4
            subprocess.check_call(command_4)

        filelist = [ f for f in os.listdir("out_"+dataset_name) if f.endswith(".nc") ]
        for f in filelist: 
            os.remove("out_"+dataset_name+"/"+f)

        data2 = loadtxt("ref_"+dataset_name+"/ref_aero_size_num_regrid.txt")
        raw_counts = loadtxt("ref_"+dataset_name+"/ref_aero_raw_counts_regrid.txt")

        f_out_num.write("%04d   %.3f   %.2f   %.1f    " % (case, prefactor, exponent, frac_dim))
        f_out_rmse_num.write("%04d   %.3f   %.2f   %.1f    " % (case, prefactor, exponent, frac_dim))

        list_num_err = []
        for col in range(1,data2.shape[1]):
            data2_1d = data2[:,col]

            # Calculate raw counts error
            ref_data_err_ratio_counts = []
            for i in raw_counts[:,col]:
                if (i == 0):
                   ref_data_err_ratio_counts.append(0)
                else:
                   ref_data_err_ratio_counts.append(1 / sqrt(i))

            # Calculate size error
            ref_data_err_ratio_size = []
            for i in range(0,data2.shape[0]):
                if (i == 0 and data2[i,col]!=0):
                   ref_data_err_ratio_size.append(0.05 * data2[i,0] * abs(data2[i+1,col] - data2[i,col]) \
                        / (data2[i+1,0] - data2[i,0]) / data2[i,col])
                elif (i == data2.shape[0]-1 and data2[i,col]!=0):
                   ref_data_err_ratio_size.append(0.05 * data2[i,0] * abs(data2[i,col] - data2[i-1,col]) \
                        / (data2[i,0] - data2[i-1,0]) / data2[i,col])
                elif (i > 0 and i < data2.shape[0]-1 and data2[i,col]!=0):
                   ref_data_err_ratio_size.append(0.05 * data2[i,0] * abs(data2[i+1,col] - data2[i-1,col]) \
                        / (data2[i+1,0] - data2[i-1,0]) / data2[i,col])
                else:
                   ref_data_err_ratio_size.append(0)

            # Calculate flow rate error (constant of 1.5%)
            ref_data_err_ratio_flow = [0.015]*data2.shape[0]

            ratio_counts = array(ref_data_err_ratio_counts)
            ratio_size = array(ref_data_err_ratio_size)
            ratio_flow = array(ref_data_err_ratio_flow)

            list_mean = []
            list_std = []
            for row in range(0,data2.shape[0]):
                list_data_temp_1d = []
                for i in arange(1,11,1):
                    file_load = "out_"+dataset_name+"/case_%04d_wc_%04d_aero_size_num.txt" % (case, i)
                    data_temp = loadtxt(file_load)
                    list_data_temp_1d.append(data_temp[row,col])
                data_temp_1d = array(list_data_temp_1d)
                data_temp_1d *= math.log(10) # convert to d*_dlogDp
                list_mean.append(mean(data_temp_1d))
                diff = data_temp_1d - mean(data_temp_1d)
                list_std.append(sqrt(1. / float(len(data_temp_1d)-1) * sum(diff**2)))
            data1_1d = array(list_mean)
            partmc_std = array(list_std)

            # calculate the relative error
            diff_num_list = []
            ref_data_err_list = []
            diams_list = []
            for i in range(0,data2.shape[0]):
                if (data2[i,col] > 0):
                   diams_list.append(data2[i,0])
                   ref_data_err_list.append(sqrt((ratio_counts[i]*data2[i,col])**2 + (ratio_size[i]*data2[i,col])**2 \
                        + (ratio_flow[i]*data2[i,col])**2 + partmc_std[i]**2/10.))
                   diff_num_list.append(data2_1d[i] - data1_1d[i])
            ref_data_err = array(ref_data_err_list)
            diams = array(diams_list)
            diff_num = array(diff_num_list)
            rel_err_num = sqrt(sum(diff_num**2 / ref_data_err**2))

            if col != data2.shape[1]-1:
               f_out_num.write("%.4f   " % (rel_err_num))
            else:
               f_out_num.write("%.4f\n" % (rel_err_num))

            list_num_err.append(rel_err_num)

        array_num_err = array(list_num_err)
        rmse_num = sqrt(sum(array_num_err**2)/float(len(array_num_err)))
        f_out_rmse_num.write("%.4f\n" % (rmse_num))

    f_out_num.close()
    f_out_rmse_num.close()

    # sort rmse data
    data6 = genfromtxt("rmse_num_"+dataset_name+".dat")
    data6_sorted = data6[argsort(data6[:,4])]
    print "Mininum error = %.4f, at prefactor = %.3f and exponent = %.2f "\
          "and frac_dim = %.1f with case %04d" %(data6_sorted[0,4], data6_sorted[0,1], \
          data6_sorted[0,2], data6_sorted[0,3], data6_sorted[0,0])
    filename_out = "rmse_num_sorted_"+dataset_name+".dat"
    f_out_sort = open(filename_out, 'w')
    f_out_sort.write("# Colume 1: caseID\n")
    f_out_sort.write("# Colume 2: prefactor\n")
    f_out_sort.write("# Colume 3: exponent\n")
    f_out_sort.write("# Colume 4: fractal dimension\n")
    f_out_sort.write("# Colume 5: total root mean square error\n")
    for row in range(0,data6_sorted.shape[0]):
        f_out_sort.write("%04d   %.3f   %.2f    %.1f    %.4f\n" % (data6_sorted[row,0], \
        data6_sorted[row,1], data6_sorted[row,2], data6_sorted[row,3], data6_sorted[row,4]))
    f_out_sort.close()

    index = index + 1
