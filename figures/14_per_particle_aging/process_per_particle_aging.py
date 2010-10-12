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

particle_set = {}
aging_ss = 0.3

netcdf_dir = "/Users/nriemer/subversion/partmc/trunk/scenarios/2_urban_plume2/out"
netcdf_pattern = "urban_plume_nc_0001_(.*).nc"
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
        dry_mass = particles.masses(exclude = ["H2O"])
        bc_frac = bc / dry_mass
	comp_vols = particles.comp_vols

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
			particle_set[id].remove_time = -1
			particle_set[id].aged_flag = False
			particle_set[id].aging_time = -1
			particle_set[id].min_s_crit = s_crit[i]
			particle_set[id].emit_bc_fraction = bc_frac[i]
			particle_set[id].emit_comp_vols = comp_vols[i]
	
	for id in particle_set.keys():
		if particle_set[id].current_id in particles.ids:
			i = np.nonzero(particles.ids == particle_set[id].current_id)[0][0]
			if s_crit[i] <= particle_set[id].min_s_crit:
				particle_set[id].min_s_crit = s_crit[i]
			if (s_crit[i] <= aging_ss) and (particle_set[id].aged_flag == False):
				particle_set[id].aging_time = time
				particle_set[id].aged_flag = True

output = open('particle_set_nc.pkl', 'wb')
pickle.dump(particle_set, output)
output.close()

