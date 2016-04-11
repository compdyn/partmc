from __future__ import division
import numpy as np

init_comp_bc_frac = np.loadtxt("aero_init_aerodyne_comp/init_comp_bc_frac.txt")
init_comp_num_frac = np.loadtxt("aero_init_aerodyne_comp/init_comp_num_frac.txt")
init_dist_diam = np.loadtxt("aero_init_aerodyne_comp/init_dist_diam.txt")
init_dist_num_conc = np.loadtxt("aero_init_aerodyne_comp/init_dist_num_conc.txt")

aero_init_dist_sampled = open('aero_init_dist_sampled_aerodyne_0828_comp.dat', 'w')
for i in range(0, len(init_comp_num_frac)): # for each BC fraction bin

	# create composition files
	aero_init_comp = open('aero_init_aerodyne_comp/aero_init_comp_aerodyne_0828_comp_'+str(i+1)+'.dat','w')
	if (i==0): # pure AS
		aero_init_comp.write("NH4  %d\n" %(36))
		aero_init_comp.write("SO4  %d\n" %(96))
	elif (i==len(init_comp_num_frac)-1): # pure BC
		aero_init_comp.write("BC  %d\n" %(1))
	else: # both AS and BC present
		bc_frac = (init_comp_bc_frac[i] + init_comp_bc_frac[i+1]) / 2
		nh4_frac = (100 - bc_frac) * 36 / (36 + 96)
		so4_frac = (100 - bc_frac) * 96 / (36 + 96)
		aero_init_comp.write("NH4  %s\n" %(nh4_frac))
		aero_init_comp.write("SO4  %s\n" %(so4_frac))
		aero_init_comp.write("BC  %s\n" %(bc_frac))
	aero_init_comp.close()

	# create distribution files
	aero_init_dist = open('aero_init_aerodyne_comp/aero_init_size_dist_aerodyne_0828_comp_'+str(i+1)+'.dat','w')
	aero_init_dist.write('diam ')
	for value in init_dist_diam:
		aero_init_dist.write('%s ' %(value))
	aero_init_dist.write('\n')
	aero_init_dist.write('num_conc ')
	for value in init_dist_num_conc:
		num_conc = value * init_comp_num_frac[i] / 100
		aero_init_dist.write('%s ' %(num_conc))
	aero_init_dist.close()

	# write to aero_dist_dist_sampled file
	aero_init_dist_sampled.write('\nmode_name init_'+str(i+1)+'\n')
	aero_init_dist_sampled.write('mass_frac aero_init_aerodyne_comp/aero_init_comp_aerodyne_0828_comp_'+str(i+1)+'.dat\n')
	aero_init_dist_sampled.write('mode_type sampled\n')
	aero_init_dist_sampled.write('size_dist aero_init_aerodyne_comp/aero_init_size_dist_aerodyne_0828_comp_'+str(i+1)+'.dat\n')

aero_init_dist_sampled.close()
