#!/usr/bin/env python

import os, sys
import scipy.io
import numpy as np

sys.path.append("/Users/nriemer/subversion/partmc/trunk/tool")
import partmc

error_detect = False
def set_error():
	global error_detect
	error_detect = True
	print "***********************************************************"
	print "***********************************************************"
	print "ERROR DETECTED"
	print "***********************************************************"
	print "***********************************************************"

def make_plot(netcdf_pattern, aging_ss, output_file_data, output_file_count):
	global error_detect
	netcdf_dir = "/Users/nriemer/subversion/partmc/trunk/local_scenarios/aging_comp/run_100K_60min/out"
	time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)

	data_array = np.zeros([len(time_filename_list),18])
	count_array = np.zeros([len(time_filename_list),18],dtype=int)

	fresh_bc_ids_prev = set()
	aged_bc_ids_prev = set()
	all_bc_ids_prev = set()
	all_ids_prev = set()
	small_ids_prev = set()
	
	comp_vol_prev = 0 

	for (i_counter, [time, filename, key]) in enumerate(time_filename_list):
		print time, filename, key
		ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
		particles = partmc.aero_particle_array_t(ncf)
		removed_info = partmc.aero_removed_info_t(ncf)
		env_state = partmc.env_state_t(ncf)
		ncf.close()

		dry_diameters = particles.dry_diameters()
		s_crit = (particles.critical_rel_humids(env_state) - 1)*100
		comp_vol = particles.comp_vols[0]

		bc = particles.masses(include = ["BC"])
		is_bc = (bc > 0)
		total_number_bc = sum(1/particles.comp_vols[is_bc])

		is_fresh_bc = ((bc > 0) & (s_crit > aging_ss))
		is_aged_bc = ((bc > 0) & (s_crit <= aging_ss))

		is_small = (dry_diameters < 5e-8)
		
		number_fresh_bc = sum(1/particles.comp_vols[is_fresh_bc])
		number_aged_bc = sum(1/particles.comp_vols[is_aged_bc])
		
		fresh_bc_ids = set(particles.ids[is_fresh_bc])
		aged_bc_ids = set(particles.ids[is_aged_bc])
		all_bc_ids = set(particles.ids[is_bc])
		all_ids = set(particles.ids)
		small_ids = set(particles.ids[is_small])

		# calculate things here
		n_f_emit_ids = fresh_bc_ids - all_ids_prev # get the emitted that are fresh
		n_a_emit_ids = aged_bc_ids - all_ids_prev # get the emitted that are aged
		
		n_f_dilute_ids = set()
		n_f_f_cond_ids = set()
		n_f_a_cond_ids = set()
		n_f_f_coag_ids = set()
		n_f_a_coag_ids = set()
		n_f_coag_ids = set()

		n_a_dilute_ids = set()
		n_a_f_cond_ids = set()
		n_a_a_cond_ids = set()
		n_a_f_coag_ids = set()
		n_a_a_coag_ids = set()
		n_a_coag_ids = set()

		
		other_ids_set = set(removed_info.other_ids)
		for id_prev in all_bc_ids_prev:  # where did everybody go?
			coag_happened = False
			id_current = id_prev
			while id_current in removed_info.ids:
				i = np.nonzero(removed_info.ids == id_current)[0][0]
				action = removed_info.actions[i]
				other_id = removed_info.other_ids[i]
				if action == removed_info.AERO_INFO_DILUTION:
					id_current = 0
				if action == removed_info.AERO_INFO_HALVED:
					id_current = 0
				if action == removed_info.AERO_INFO_COAG:
					id_current = other_id
					coag_happened = True
			if id_current in other_ids_set:
				coag_happened = True
			if id_current == 0:
				if id_prev in fresh_bc_ids_prev:
					n_f_dilute_ids.add(id_prev)
				if id_prev in aged_bc_ids_prev:
					n_a_dilute_ids.add(id_prev)
			else:
				if id_prev in fresh_bc_ids_prev:
					if id_current in fresh_bc_ids:
						if coag_happened:
							n_f_f_coag_ids.add(id_prev)
						else:
							n_f_f_cond_ids.add(id_prev)
					elif id_current in aged_bc_ids:
						if coag_happened:
							n_f_a_coag_ids.add(id_prev)
						else:
							n_f_a_cond_ids.add(id_prev)
					else:
						raise Exception("Current not fresh or aged")
				if id_prev in aged_bc_ids_prev:
					if id_current in fresh_bc_ids:
						if coag_happened:
							n_a_f_coag_ids.add(id_prev)
						else:
							n_a_f_cond_ids.add(id_prev)
					elif id_current in aged_bc_ids:
						if coag_happened:
							n_a_a_coag_ids.add(id_prev)
						else:
							n_a_a_cond_ids.add(id_prev)
					else:
						raise Exception("Current not fresh or aged")

		for id_current in all_bc_ids:
			if id_current in all_ids_prev: # check if it wasn't emitted
				if id_current in other_ids_set:
					if id_current in fresh_bc_ids:
						n_f_coag_ids.add(id_current)
					if id_current in aged_bc_ids:
						n_a_coag_ids.add(id_current)


		n_f_emit = len(n_f_emit_ids)
		n_f_dilute = len(n_f_dilute_ids)
		n_f_f_cond = len(n_f_f_cond_ids)
		n_f_a_cond = len(n_f_a_cond_ids)
		n_f_f_coag = len(n_f_f_coag_ids)
		n_f_a_coag = len(n_f_a_coag_ids)
		n_f_coag = len(n_f_coag_ids)

		n_f_a_cond_small = len(n_f_a_cond_ids & small_ids_prev)
		n_f_a_coag_small = len(n_f_a_coag_ids & small_ids_prev)

		n_a_emit = len(n_a_emit_ids)
		n_a_dilute = len(n_a_dilute_ids)
		n_a_f_cond = len(n_a_f_cond_ids)
		n_a_a_cond = len(n_a_a_cond_ids)
		n_a_f_coag = len(n_a_f_coag_ids)
		n_a_a_coag = len(n_a_a_coag_ids)
		n_a_coag = len(n_a_coag_ids)

		n_f = len(fresh_bc_ids)
		n_a = len(aged_bc_ids)
		n_f_prev = len(fresh_bc_ids_prev)
		n_a_prev = len(aged_bc_ids_prev)

		n_f_small = len(fresh_bc_ids & small_ids)

		n_aging = n_f_a_cond + n_f_a_coag
		n_aging_small = n_f_a_cond_small + n_f_a_coag_small

		data_array[i_counter,0] = time
		data_array[i_counter,1] = total_number_bc
		data_array[i_counter,2] = number_fresh_bc
		data_array[i_counter,3] = number_aged_bc
		data_array[i_counter,4] = n_f_emit / comp_vol
		data_array[i_counter,5] = n_a_emit / comp_vol
		if (comp_vol_prev != 0):
			data_array[i_counter,6] = n_f_dilute / comp_vol_prev  
		data_array[i_counter,7] = n_f_f_cond / comp_vol
		data_array[i_counter,8] = n_f_a_cond / comp_vol
		data_array[i_counter,9] = n_f_f_coag / comp_vol
		data_array[i_counter,10] = n_f_a_coag / comp_vol
		data_array[i_counter,11] = n_f_coag / comp_vol
		if (comp_vol_prev != 0):
			data_array[i_counter,12] = n_a_dilute / comp_vol_prev
		data_array[i_counter,13] = n_a_f_cond / comp_vol
		data_array[i_counter,14] = n_a_a_cond / comp_vol
		data_array[i_counter,15] = n_a_f_coag / comp_vol
		data_array[i_counter,16] = n_a_a_coag / comp_vol
		data_array[i_counter,17] = n_a_coag / comp_vol

		count_array[i_counter,0] = time
		count_array[i_counter,1] = len(all_bc_ids)
		count_array[i_counter,2] = n_f               #len(fresh_bc_ids)
		count_array[i_counter,3] = n_a               #len(aged_bc_ids)
		count_array[i_counter,4] = n_f_emit          #len(n_f_emit_ids)
		count_array[i_counter,5] = n_a_emit          #len(n_a_emit_ids)
		count_array[i_counter,6] = n_f_dilute        #len(n_f_dilute_ids)
		count_array[i_counter,7] = n_f_f_cond        #len(n_f_f_cond_ids)
		count_array[i_counter,8] = n_f_a_cond        #len(n_f_a_cond_ids)
		count_array[i_counter,9] = n_f_f_coag        #len(n_f_f_coag_ids)
		count_array[i_counter,10] = n_f_a_coag        #len(n_f_a_coag_ids)
		count_array[i_counter,11] = n_f_coag         #len(n_f_coag_ids)
		count_array[i_counter,12] = n_a_dilute       #len(n_a_dilute_ids)
		count_array[i_counter,13] = n_a_f_cond       #len(n_a_f_cond_ids)
		count_array[i_counter,14] = n_a_a_cond       #len(n_a_a_cond_ids)
		count_array[i_counter,15] = n_a_f_coag       #len(n_a_f_coag_ids)
		count_array[i_counter,16] = n_a_a_coag       #len(n_a_a_coag_ids)
		count_array[i_counter,17] = n_a_coag         #len(n_a_coag_ids)

# The following is checking code #

		lhs_f = n_f - n_f_prev
		rhs_f = n_f_emit - n_f_dilute + n_a_f_cond + n_f_coag - n_f_f_coag - n_f_a_cond - n_f_a_coag

		lhs_a = n_a - n_a_prev
		rhs_a = n_a_emit - n_a_dilute + n_f_a_cond + n_a_coag - n_a_a_coag - n_a_f_cond - n_a_f_coag

		loss_f = n_f_dilute + n_f_f_coag + n_f_a_coag + n_f_f_cond + n_f_a_cond
		loss_f_ids = n_f_dilute_ids | n_f_f_coag_ids | n_f_a_coag_ids | n_f_f_cond_ids | n_f_a_cond_ids

		gain_f = n_f_emit + n_f_coag + n_f_f_cond + n_a_f_cond
		gain_f_ids = n_f_emit_ids | n_f_coag_ids | n_f_f_cond_ids | n_a_f_cond_ids

		loss_a = n_a_dilute + n_a_f_coag + n_a_a_coag + n_a_f_cond + n_a_a_cond
		loss_a_ids = n_a_dilute_ids | n_a_f_coag_ids | n_a_a_coag_ids | n_a_f_cond_ids | n_a_a_cond_ids

		gain_a = n_a_emit + n_a_coag + n_f_a_cond + n_a_a_cond
		gain_a_ids = n_a_emit_ids | n_a_coag_ids | n_f_a_cond_ids | n_a_a_cond_ids

		coag_loss_f = n_f_f_coag + n_a_f_coag
		coag_loss_f_ids = n_f_f_coag_ids | n_a_f_coag_ids

		coag_loss_a = n_a_a_coag + n_f_a_coag
		coag_loss_a_ids = n_a_a_coag_ids | n_f_a_coag_ids

		print "budget f    %8d %8d %8d" % (lhs_f, rhs_f, lhs_f - rhs_f)
		if (lhs_f - rhs_f) != 0: set_error()
		print "budget a    %8d %8d %8d" % (lhs_a, rhs_a, lhs_a - rhs_a)
		if (lhs_a - rhs_a) != 0: set_error()
		print "loss f      %8d %8d %8d" % (loss_f, n_f_prev, loss_f - n_f_prev)
		if (loss_f - n_f_prev) != 0: set_error()
		print "gain f      %8d %8d %8d" % (gain_f, n_f, gain_f - n_f)
		if (gain_f - n_f) != 0: set_error()
		print "loss a      %8d %8d %8d" % (loss_a, n_a_prev, loss_a - n_a_prev)
		if (loss_a - n_a_prev) != 0: set_error()
		print "gain a      %8d %8d %8d" % (gain_a, n_a, gain_a - n_a)
		if (gain_a - n_a) != 0: set_error()

		# We don't check equation (12) in Riemer et al., 2010, Aerosol Sci. because  it can 
		# happen that we emit a BC particle and it coagulates before we detect it.

		print "loss f ids ", loss_f_ids - fresh_bc_ids_prev, fresh_bc_ids_prev - loss_f_ids
		if (len(loss_f_ids - fresh_bc_ids_prev)) !=0: set_error()
		if (len(fresh_bc_ids_prev - loss_f_ids)) !=0: set_error()	
		print "gain f ids ", gain_f_ids - fresh_bc_ids, fresh_bc_ids - gain_f_ids
		if (len(gain_f_ids - fresh_bc_ids)) !=0: set_error()
		if (len(fresh_bc_ids - gain_f_ids)) !=0: set_error()
		print "loss a ids ", loss_a_ids - aged_bc_ids_prev, aged_bc_ids_prev - loss_a_ids
		if (len(loss_a_ids - aged_bc_ids_prev)) != 0: set_error()
		if (len(aged_bc_ids_prev - loss_a_ids)) != 0: set_error()
		print "gain a ids ", gain_a_ids - aged_bc_ids, aged_bc_ids - gain_a_ids
		if (len(gain_a_ids - aged_bc_ids)) != 0: set_error()
		if (len(aged_bc_ids - gain_a_ids)) != 0: set_error()

		print "gain f intersection ", n_f_emit_ids & n_f_coag_ids, n_f_emit_ids & n_f_f_cond_ids, \
		    n_f_emit_ids & n_a_f_cond_ids,  n_f_coag_ids & n_f_f_cond_ids, \
		    n_f_coag_ids & n_a_f_cond_ids, n_f_f_cond_ids & n_a_f_cond_ids
		if (len(n_f_emit_ids & n_f_coag_ids) !=0 or len(n_f_emit_ids & n_f_f_cond_ids) !=0 or \
		    len(n_f_emit_ids & n_a_f_cond_ids) !=0 or len(n_f_coag_ids & n_f_f_cond_ids) !=0 or \
		    len(n_f_coag_ids & n_a_f_cond_ids) !=0 or len( n_f_f_cond_ids & n_a_f_cond_ids)): set_error()

		print "gain a intersection ", n_a_emit_ids & n_a_coag_ids, n_a_emit_ids &  n_f_a_cond_ids, \
		    n_a_emit_ids & n_a_a_cond_ids, n_a_coag_ids & n_f_a_cond_ids, n_a_coag_ids & n_a_a_cond_ids, \
		    n_f_a_cond_ids & n_a_a_cond_ids
		if (len(n_a_emit_ids & n_a_coag_ids) !=0 or len(n_a_emit_ids &  n_f_a_cond_ids) !=0 or \
		    len(n_a_emit_ids & n_a_a_cond_ids) !=0 or len(n_a_coag_ids & n_f_a_cond_ids) !=0 or len(n_a_coag_ids & n_a_a_cond_ids) !=0 or \
		    len(n_f_a_cond_ids & n_a_a_cond_ids)!=0): set_error()

		print "loss f intersection ", n_f_dilute_ids & n_f_f_coag_ids, n_f_dilute_ids & n_f_a_coag_ids, \
			    n_f_dilute_ids & n_f_f_cond_ids, n_f_dilute_ids &  n_f_a_cond_ids, \
			    n_f_f_coag_ids & n_f_a_coag_ids, n_f_f_coag_ids & n_f_f_cond_ids, n_f_f_coag_ids & n_f_a_cond_ids, \
			    n_f_a_coag_ids & n_f_f_cond_ids, n_f_a_coag_ids & n_f_a_cond_ids, \
			    n_f_f_cond_ids & n_f_a_cond_ids
		if (len(n_f_dilute_ids & n_f_f_coag_ids) !=0 or len(n_f_dilute_ids & n_f_a_coag_ids) !=0 or \
			    len(n_f_dilute_ids & n_f_f_cond_ids) !=0 or len(n_f_dilute_ids & n_f_a_cond_ids) !=0 or \
			    len(n_f_f_coag_ids & n_f_a_coag_ids) !=0 or len(n_f_f_coag_ids & n_f_f_cond_ids) !=0 or len(n_f_f_coag_ids & n_f_a_cond_ids) !=0 or \
			    len(n_f_a_coag_ids & n_f_f_cond_ids) !=0 or len(n_f_a_coag_ids & n_f_a_cond_ids) !=0 or \
			    len(n_f_f_cond_ids & n_f_a_cond_ids) !=0): set_error()

		print "loss a intersection ", n_a_dilute_ids & n_a_f_coag_ids, n_a_dilute_ids & n_a_a_coag_ids, \
			    n_a_dilute_ids & n_a_f_cond_ids, n_a_dilute_ids & n_a_a_cond_ids, \
			    n_a_f_coag_ids & n_a_a_coag_ids, n_a_f_coag_ids & n_a_f_cond_ids, n_a_f_coag_ids & n_a_a_cond_ids, \
			    n_a_a_coag_ids & n_a_f_cond_ids, n_a_a_coag_ids & n_a_a_cond_ids, \
			    n_a_f_cond_ids & n_a_a_cond_ids
		if (len(n_a_dilute_ids & n_a_f_coag_ids) !=0 or len(n_a_dilute_ids & n_a_a_coag_ids) !=0 or \
			    len(n_a_dilute_ids & n_a_f_cond_ids) !=0 or len(n_a_dilute_ids & n_a_a_cond_ids) !=0 or \
			    len(n_a_f_coag_ids & n_a_a_coag_ids) !=0 or len(n_a_f_coag_ids & n_a_f_cond_ids) !=0 or len(n_a_f_coag_ids & n_a_a_cond_ids) !=0 or \
			    len(n_a_a_coag_ids & n_a_f_cond_ids) !=0 or len(n_a_a_coag_ids & n_a_a_cond_ids) !=0 or \
			    len(n_a_f_cond_ids & n_a_a_cond_ids) !=0): set_error()

		print "gain f/ gain a intersection ", gain_f_ids & gain_a_ids
		if (len(gain_f_ids & gain_a_ids) != 0): set_error()
		print "loss f/ loss a intersection ", loss_f_ids & loss_a_ids
		if (len(loss_f_ids & loss_a_ids) != 0): set_error()

		fresh_bc_ids_prev = fresh_bc_ids
		aged_bc_ids_prev = aged_bc_ids
		all_bc_ids_prev = all_bc_ids
		all_ids_prev = all_ids
		small_ids_prev = small_ids
		comp_vol_prev = comp_vol
# Checking code ends #

	np.savetxt(output_file_data, data_array)
	np.savetxt(output_file_count, count_array, fmt='%d')

make_plot("urban_plume2_wc_0001_(.*).nc", 0.1, "aging_data_wc_01_100K_60min.txt", "aging_count_wc_01_100K_60min.txt")
make_plot("urban_plume2_wc_0001_(.*).nc", 0.3, "aging_data_wc_03_100K_60min.txt", "aging_count_wc_03_100K_60min.txt")
make_plot("urban_plume2_wc_0001_(.*).nc", 0.6, "aging_data_wc_06_100K_60min.txt", "aging_count_wc_06_100K_60min.txt")
print "error detect ", error_detect

