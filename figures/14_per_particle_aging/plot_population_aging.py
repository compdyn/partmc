#!/usr/bin/env python

import os, sys
import scipy.io
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("/Users/nriemer/subversion/partmc/trunk/tool")
import partmc
import mpl_helper

def make_plot(input_file_data1, input_file_data2, input_file_data3, input_file_data4, input_file_data5, input_file_data6, input_file_data7, input_file_data8, input_file_data9, output_figure):

	infile1 = []
	infile1.append(input_file_data1)
	infile1.append(input_file_data2)
	infile1.append(input_file_data3)

	infile2 = []
	infile2.append(input_file_data4)
	infile2.append(input_file_data5)
	infile2.append(input_file_data6)

	infile3 = []
	infile3.append(input_file_data7)
	infile3.append(input_file_data8)
	infile3.append(input_file_data9)

	tau_array_sm_1 = np.zeros([3,1441])
	tau_array_sm_2 = np.zeros([3,145])
	tau_array_sm_3 = np.zeros([3,25])

	smooth_window_width_1 = 60 # in timesteps
	smooth_window_width_2 = 6
	smooth_window_width_3 = 2

	for counter_min in range(0,3):
		if counter_min == 0:
			infile = infile1
			smooth_window_width = smooth_window_width_1
		if counter_min == 1:
			infile = infile2
			smooth_window_width = smooth_window_width_2
		if counter_min == 2:
			infile = infile3
			smooth_window_width = smooth_window_width_3
			
		for counter_ss in range(0,3):	
			
			print "file ", infile[counter_ss]

			data_array = np.loadtxt(infile[counter_ss])
		
			delta_t = data_array[1,0] - data_array[0,0]
			print "delta_t ", delta_t, len(data_array[:,0])

			n_f = data_array[:,2]
			n_f_a_cond = data_array[:,8]
			n_f_a_coag = data_array[:,10]
	
			n_f_sm = partmc.smooth(data_array[:,2], window_len=smooth_window_width)
			n_f_a_cond_sm = partmc.smooth(data_array[:,8], window_len=smooth_window_width)
			n_f_a_coag_sm = partmc.smooth(data_array[:,10], window_len=smooth_window_width)

			n_aging = n_f_a_cond + n_f_a_coag
			n_aging_sm = n_f_a_cond_sm + n_f_a_coag_sm

# time scale calculation #

			tau_sm = delta_t * n_f_sm / n_aging_sm / 3600. # aging time scale in hours
			print "tau ", tau_sm

			if counter_min == 0:
				time_1 = data_array[:,0] / 3600.
				tau_array_sm_1[counter_ss,:] = tau_sm
			if counter_min == 1:
				time_2 = data_array[:,0] / 3600.
				tau_array_sm_2[counter_ss,:] = tau_sm
			if counter_min == 2:
				time_3 = data_array[:,0] / 3600.
				tau_array_sm_3[counter_ss,:] = tau_sm

	fig = plt.figure()
	plt.semilogy(time_1, tau_array_sm_1[0,:], 'r', label = 'SS = 0.6\%, 1 min, 1K')
	plt.semilogy(time_1, tau_array_sm_1[1,:], 'g', label = 'SS = 0.6\%, 1 min, 10K')
	plt.semilogy(time_1, tau_array_sm_1[2,:], 'b', label = 'SS = 0.6\%, 1 min, 100K')
	plt.semilogy(time_2, tau_array_sm_2[0,:], '--r', label = 'SS = 0.6\%, 10 min, 1K')
	plt.semilogy(time_2, tau_array_sm_2[1,:], '--g', label = 'SS = 0.6\%, 10 min, 10K')
	plt.semilogy(time_2, tau_array_sm_2[2,:], '--b', label = 'SS = 0.6\%, 10 min, 100K')
	plt.semilogy(time_3-0.5, tau_array_sm_3[0,:], ':r', label = 'SS = 0.6\%, 60 min, 1K')
	plt.semilogy(time_3-0.5, tau_array_sm_3[1,:], ':g', label = 'SS = 0.6\%, 60 min, 10K')
	plt.semilogy(time_3-0.5, tau_array_sm_3[2,:], ':b', label = 'SS = 0.6\%, 60 min, 100K')
	plt.xlim(0, 24)
	plt.legend(loc = 'lower right')
	plt.xticks([0, 6, 12, 18, 24])
	plt.ylim(1e-1, 1e3)
	plt.grid(True)
	plt.xlabel("time in hours")
	plt.ylabel(r"$\tau$ in hours")
#	fig = plt.gcf()
	fig.savefig(output_figure)

make_plot("aging_data_wc_06_1K_1min.txt", "aging_data_wc_06_10K_1min.txt", "aging_data_wc_06_10K_1min.txt", 
	  "aging_data_wc_06_1K_10min.txt", "aging_data_wc_06_10K_10min.txt", "aging_data_wc_06_100K_10min.txt", 
	  "aging_data_wc_06_1K_60min.txt", "aging_data_wc_06_10K_60min.txt", "aging_data_wc_06_100K_60min.txt", "figs/tau_06.pdf")




