#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import shutil
from numpy import *
#import mpl_helper
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

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

dataset = ['0212','0909','0925']
t_max = [15600, 19740, 20160]
t_output = [600, 420, 420]
rh = [0.0385, 0.1275, 0.0986]
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
    command_0 = ["cp", "aero_init_size_dist_"+dataset_name+".dat", "aero_init_size_dist.dat"]
    subprocess.check_call(command_0)
    for prefactor in arange(0.025,0.055,0.005):
        for exponent in arange(0.23,0.27,0.01):
            for frac_dim in arange(2.2,3.1,0.1):
                filename_in = "barrel_template_with_fractal.spec"
                case += 1
                filename_out = "spec_"+dataset_name+"/barrel_wc_case_%04d.spec" % (case)
                f_in = open(filename_in, 'r')
                f_out = open(filename_out, 'w')

                for line in f_in:
                    line = line.replace('%%OUTPUT_PREFIX%%', 'out_'+dataset_name+'/case_%04d_wc' % (case))
                    line = line.replace('%%t_max%%', str(t_max[index]))
                    line = line.replace('%%t_output%%', str(t_output[index]))
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
    filename_out_mass = "rel_err_mass_"+dataset_name+".dat"
    f_out_mass = open(filename_out_mass, 'w')
    f_out_mass.write("# Colume 1: caseID\n")
    f_out_mass.write("# Colume 2: prefactor\n")
    f_out_mass.write("# Colume 3: exponent\n")
    f_out_mass.write("# Colume 4: fractal dimension\n")
    f_out_mass.write("# Colume j+4: relative error of mass distribution at time(j)\n")
    filename_out_total = "rel_err_total_"+dataset_name+".dat"
    f_out_total = open(filename_out_total, 'w')
    f_out_total.write("# Colume 1: caseID\n")
    f_out_total.write("# Colume 2: prefactor\n")
    f_out_total.write("# Colume 3: exponent\n")
    f_out_total.write("# Colume 4: fractal dimension\n")
    f_out_total.write("# Colume j+4: total relative error (num + mass) at time(j)\n")
    filename_out_rmse_num = "rmse_num_"+dataset_name+".dat"
    f_out_rmse_num = open(filename_out_rmse_num, 'w')
    f_out_rmse_num.write("# Colume 1: caseID\n")
    f_out_rmse_num.write("# Colume 2: prefactor\n")
    f_out_rmse_num.write("# Colume 3: exponent\n")
    f_out_rmse_num.write("# Colume 4: fractal dimension\n")
    f_out_rmse_num.write("# Colume 5: root mean square error of number distribution\n")
    filename_out_rmse_mass = "rmse_mass_"+dataset_name+".dat"
    f_out_rmse_mass = open(filename_out_rmse_mass, 'w')
    f_out_rmse_mass.write("# Colume 1: caseID\n")
    f_out_rmse_mass.write("# Colume 2: prefactor\n")
    f_out_rmse_mass.write("# Colume 3: exponent\n")
    f_out_rmse_mass.write("# Colume 4: fractal dimension\n")
    f_out_rmse_mass.write("# Colume 5: root mean square error of mass distribution\n")
    filename_out_rmse_total = "rmse_total_"+dataset_name+".dat"
    f_out_rmse_total = open(filename_out_rmse_total, 'w')
    f_out_rmse_total.write("# Colume 1: caseID\n")
    f_out_rmse_total.write("# Colume 2: prefactor\n")
    f_out_rmse_total.write("# Colume 3: exponent\n")
    f_out_rmse_total.write("# Colume 4: fractal dimension\n")
    f_out_rmse_total.write("# Colume 5: total root mean square error (num + mass)\n")

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

        ncfile_prefix = "out_"+dataset_name+"/case_%04d_wc_0001" % (case)

        command_1 = [str_exec,spec_file]
        command_2 = [str_extr_size,"--num","--dmin","1e-8","--dmax","1e-6","--nbin","100",ncfile_prefix]
        command_3 = [str_extr_size,"--mass","--dmin","1e-8","--dmax","1e-6","--nbin","100",ncfile_prefix]
        command_4 = [str_extr_time,ncfile_prefix]

        print command_1
        subprocess.check_call(command_1)

        print command_2
        subprocess.check_call(command_2)
        data1 = loadtxt(ncfile_prefix+"_aero_size_num.txt")
        data2 = loadtxt("ref_"+dataset_name+"/ref_aero_size_num_regrid.txt")
        raw_counts = loadtxt("ref_"+dataset_name+"/ref_aero_raw_counts_regrid.txt")

        print command_3
        subprocess.check_call(command_3)
        #data3 = loadtxt(ncfile_prefix+"_aero_size_mass.txt")
        #data4 = loadtxt("ref_"+dataset_name+"/ref_aero_size_mass_regrid.txt")

        print command_4
        subprocess.check_call(command_4)

        f_out_num.write("%04d   %.3f   %.2f   %.1f    " % (case, prefactor, exponent, frac_dim))
        f_out_mass.write("%04d   %.3f   %.2f   %.1f    " % (case, prefactor, exponent, frac_dim))
        f_out_total.write("%04d   %.3f   %.2f   %.1f   " % (case, prefactor, exponent, frac_dim))
        f_out_rmse_num.write("%04d   %.3f   %.2f   %.1f    " % (case, prefactor, exponent, frac_dim))
        f_out_rmse_mass.write("%04d   %.3f   %.2f   %.1f    " % (case, prefactor, exponent, frac_dim))
        f_out_rmse_total.write("%04d   %.3f   %.2f   %.1f   " % (case, prefactor, exponent, frac_dim))

        list_num_err = []
        list_mass_err = []
        list_total_err = []
        for col in range(1,data1.shape[1]):
            diameters = data1[:,0]
            data1_1d = data1[:,col]
            data2_1d = data2[:,col]
            #data3_1d = data3[:,col]
            #data4_1d = data4[:,col]

            # extract_aero_size gives d*_dlnDp, neet to convert to d*_dlogDp
            data1_1d *= math.log(10)
            #data3_1d *= math.log(10)
            data3_1d = data1_1d * math.pi / 6. * rho * diameters**3
            data4_1d = data2_1d * math.pi / 6. * rho * diameters**3

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

            # calculate the relative error
            diff_num_list = []
            diff_mass_list = []
            ref_data_err_list = []
            for i in range(0,data2.shape[0]):
                if (data2[i,col] > 0):
                   ref_data_err_list.append(sqrt(ratio_counts[i]**2 + ratio_size[i]**2 + ratio_flow[i]**2) * data2[i,col])
                   diff_num_list.append(data2_1d[i] - data1_1d[i])
                   diff_mass_list.append(data4_1d[i] - data3_1d[i])
            ref_data_err = array(ref_data_err_list)
            ref_data_err_mass = math.pi / 6. * rho * diameters**3 * ref_data_err
            diff_num = array(diff_num_list)
            rel_err_num = sqrt(sum(diff_num**2 / ref_data_err**2))
            diff_mass = array(diff_mass_list)
            rel_err_mass = sqrt(sum(diff_mass**2 / ref_data_err_mass**2))
            # combine num and mass errors by calculating root mean square error
            rel_err_total = sqrt(0.5*(rel_err_num**2 + rel_err_mass**2))

            if col != data1.shape[1]-1:
               f_out_num.write("%.4f   " % (rel_err_num))
               f_out_mass.write("%.4f   " % (rel_err_mass))
               f_out_total.write("%.4f   " % (rel_err_total))
            else:
               f_out_num.write("%.4f\n" % (rel_err_num))
               f_out_mass.write("%.4f\n" % (rel_err_mass))
               f_out_total.write("%.4f\n" % (rel_err_total))

            list_num_err.append(rel_err_num)
            list_mass_err.append(rel_err_mass)
            list_total_err.append(rel_err_total)

        array_num_err = array(list_num_err)
        array_mass_err  = array(list_mass_err)
        array_total_err = array(list_total_err)
        rmse_num = sqrt(sum(array_num_err**2)/len(array_num_err))
        rmse_mass = sqrt(sum(array_mass_err**2)/len(array_mass_err))
        rmse_total = sqrt(sum(array_total_err**2)/len(array_total_err))
        f_out_rmse_num.write("%.4f\n" % (rmse_num))
        f_out_rmse_mass.write("%.4f\n" % (rmse_mass))
        f_out_rmse_total.write("%.4f\n" % (rmse_total))

    f_out_num.close()
    f_out_mass.close()
    f_out_total.close()
    f_out_rmse_num.close()
    f_out_rmse_mass.close()
    f_out_rmse_total.close()

    # sort total rmse data
    data6 = genfromtxt("rmse_total_"+dataset_name+".dat")
    data6_sorted = data6[argsort(data6[:,4])]
    print "Mininum total error = %.4f, at prefactor = %.3f and exponent = %.2f "\
          "and frac_dim = %.1f with case %04d" %(data6_sorted[0,4], data6_sorted[0,1], \
          data6_sorted[0,2], data6_sorted[0,3], data6_sorted[0,0])
    filename_out = "rmse_total_sorted_"+dataset_name+".dat"
    f_out_sort = open(filename_out, 'w')
    f_out_sort.write("# Colume 1: caseID\n")
    f_out_sort.write("# Colume 2: prefactor\n")
    f_out_sort.write("# Colume 3: exponent\n")
    f_out_sort.write("# Colume 4: fractal dimension\n")
    f_out_sort.write("# Colume 5: total root mean square error (num + mass)\n")
    for row in range(0,data6_sorted.shape[0]):
        f_out_sort.write("%04d   %.3f   %.2f    %.1f    %.4f\n" % (data6_sorted[row,0], \
        data6_sorted[row,1], data6_sorted[row,2], data6_sorted[row,3], data6_sorted[row,4]))
    f_out_sort.close()

    index = index + 1

    #data5 = genfromtxt("rmse_num.dat", skip_header=5)
    #data6 = genfromtxt("rmse_mass.dat", skip_header=5)
    #data7 = genfromtxt("rmse_total.dat", skip_header=5)
    #(figure, axes) = mpl_helper.make_fig(colorbar=False)
    #axes.plot(data5[:,1], data5[:,2], marker='^')
    #axes.plot(data6[:,1], data6[:,2], marker='D')
    #axes.plot(data7[:,1], data7[:,2], marker='o')
    #axes.set_title(r"$\mathrm{k}_{\mathrm{D}}$ = %.3f, a = %.2f, "\
    #               "$\mathrm{R}_0$ = %.1e, f = %.2f" \
    #               % (prefactor, exponent, prime_radius, vol_fill_factor))
    #axes.set_xlabel("fractal dimension")
    #axes.set_ylabel("root mean square error")
    #axes.grid()
    #box = axes.get_position()
    #axes.set_position([box.x0, box.y0, box.width * 0.78, box.height])
    #axes.legend(('Num', 'Mass', 'Total'), loc='center left', bbox_to_anchor=(1, 0.5))
    #filename_out = "plot_rmse_vs_frac_dim.pdf"
    #figure.savefig(filename_out)
