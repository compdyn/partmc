#!/usr/bin/env python

import os, sys
import scipy.io
import numpy as np

sys.path.append("/Users/nriemer/subversion/partmc/trunk/tool")
import partmc
import pickle

class Struct(object): 
	def __init__(self): 
		pass

def make_plot(netcdf_pattern, aging_ss, output_pkl):
	particle_set = {}

	netcdf_dir = "/Users/nriemer/subversion/partmc/trunk/scenarios/1_urban_plume/out"
	time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)

	for [time, filename, key] in time_filename_list:
		print time, filename, key
		ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
		particles = partmc.aero_particle_array_t(ncf)
		removed_info = partmc.aero_removed_info_t(ncf)
		env_state = partmc.env_state_t(ncf)
		ncf.close()

		dry_diameters = particles.dry_diameters()
		s_crit = (particles.critical_rel_humids(env_state) - 1)*100

		bc = particles.masses(include = ["BC"])
		no3 = particles.masses(include = ["NO3"])
		so4 = particles.masses(include = ["SO4"])
		nh4 = particles.masses(include = ["NH4"])

		dry_mass = particles.masses(exclude = ["H2O"])
		bc_frac = bc / dry_mass
		no3_frac = no3 / dry_mass
		so4_frac = so4 / dry_mass
		nh4_frac = nh4 / dry_mass
		solute_frac = (no3 + so4 + nh4) / dry_mass

		kappas = particles.kappas()
		comp_vols = particles.comp_vols
		h2o = particles.masses(include = ["H2O"])

		for id in particle_set.keys():
			while particle_set[id].current_id in removed_info.ids:
				i = np.nonzero(removed_info.ids == particle_set[id].current_id)[0][0]
				if removed_info.actions[i] != removed_info.AERO_INFO_COAG:
					particle_set[id].remove_time = time
					particle_set[id].current_id = 0
				else:
					particle_set[id].current_id = removed_info.other_ids[i]

		for i in range(len(particles.ids)):
			id = particles.ids[i]
			if id not in particle_set:
				particle_set[id] = Struct()
				particle_set[id].current_id = id
				particle_set[id].emit_time = time
				particle_set[id].emit_diam = dry_diameters[i]
				particle_set[id].emit_s_crit = s_crit[i]
				particle_set[id].emit_kappa = kappas[i]
				particle_set[id].remove_time = -1
				particle_set[id].aged_flag = False
				particle_set[id].aging_time = -1
				particle_set[id].min_s_crit = s_crit[i]
				particle_set[id].emit_bc_fraction = bc_frac[i]
				particle_set[id].emit_comp_vols = comp_vols[i]
				particle_set[id].aging_kappa = -1
				particle_set[id].aging_h2o = -1
				particle_set[id].aging_no3_fraction = -1
				particle_set[id].aging_so4_fraction = -1
				particle_set[id].aging_nh4_fraction = -1
				particle_set[id].aging_solute_fraction = -1
				particle_set[id].aging_diameter = -1

		for id in particle_set.keys():
			if particle_set[id].current_id in particles.ids:
				i = np.nonzero(particles.ids == particle_set[id].current_id)[0][0]
				if s_crit[i] <= particle_set[id].min_s_crit:
					particle_set[id].min_s_crit = s_crit[i]
				if (s_crit[i] <= aging_ss) and (particle_set[id].aged_flag == False):
					particle_set[id].aging_time = time
					particle_set[id].aged_flag = True
					particle_set[id].aging_kappa = kappas[i]
					particle_set[id].aging_diameter = dry_diameters[i]
					particle_set[id].aging_h2o = h2o[i]
				        particle_set[id].aging_no3_fraction = no3_frac[i]
					particle_set[id].aging_so4_fraction = so4_frac[i]
					particle_set[id].aging_nh4_fraction = nh4_frac[i]
					particle_set[id].aging_solute_fraction = solute_frac[i]
					

	output = open(output_pkl, 'wb')
	pickle.dump(particle_set, output)
	output.close()

#make_plot("urban_plume_wc_noon_0001_(.*).nc", 0.6, "particle_set_wc_noon_06.pkl")
#make_plot("urban_plume_wc_noon_0001_(.*).nc", 0.3, "particle_set_wc_noon_03.pkl")
#make_plot("urban_plume_wc_noon_0001_(.*).nc", 0.1, "particle_set_wc_noon_01.pkl")
#make_plot("urban_plume_nc_noon_0001_(.*).nc", 0.6, "particle_set_nc_noon_06.pkl")
#make_plot("urban_plume_nc_noon_0001_(.*).nc", 0.3, "particle_set_nc_noon_03.pkl")
#make_plot("urban_plume_nc_noon_0001_(.*).nc", 0.1, "particle_set_nc_noon_01.pkl")
make_plot("urban_plume_wc_0001_(.*).nc", 0.6, "particle_set_wc_06.pkl")
make_plot("urban_plume_wc_0001_(.*).nc", 0.3, "particle_set_wc_03.pkl")
make_plot("urban_plume_wc_0001_(.*).nc", 0.1, "particle_set_wc_01.pkl")
#make_plot("urban_plume_nc_0001_(.*).nc", 0.6, "particle_set_nc_06.pkl")
#make_plot("urban_plume_nc_0001_(.*).nc", 0.3, "particle_set_nc_03.pkl")
#make_plot("urban_plume_nc_0001_(.*).nc", 0.1, "particle_set_nc_01.pkl")
