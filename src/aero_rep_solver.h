/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header for aerosol representation functions
 *
 */
/** \file
 * \brief Header file for aerosol representation sovler functions
 */
#ifndef AERO_REP_SOLVER_H_
#define AERO_REP_SOLVER_H_
#include "phlex_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef PMC_USE_SUNDIALS

// modal mass
void * aero_rep_modal_binned_mass_get_dependencies(void *aero_rep_data, bool *state_flags);
void * aero_rep_modal_binned_mass_update_env_state(double *env_data, void *aero_rep_data);
void * aero_rep_modal_binned_mass_update_state(ModelData *model_data, void *aero_rep_data);
void * aero_rep_modal_binned_mass_get_effective_radius(int aero_phase_idx, double *radius, 
		double *partial_deriv, void *aero_rep_data);
void * aero_rep_modal_binned_mass_get_number_conc(int aero_phase_idx, double *number_conc, 
		double *partial_deriv, void *aero_rep_data);
void * aero_rep_modal_binned_mass_get_aero_conc_type(int aero_phase_idx, int *aero_conc_type, void *aero_rep_data);
void * aero_rep_modal_binned_mass_update_data(int update_type, void *update_data, void *aero_rep_data);
void * aero_rep_modal_binned_mass_print(void *aero_rep_data);
void * aero_rep_modal_binned_mass_skip(void *aero_rep_data);

// single particle
void * aero_rep_single_particle_get_dependencies(void *aero_rep_data, bool *state_flags);
void * aero_rep_single_particle_update_env_state(double *env_data, void *aero_rep_data);
void * aero_rep_single_particle_update_state(ModelData *model_data, void *aero_rep_data);
void * aero_rep_single_particle_get_effective_radius(int aero_phase_idx, double *radius, 
		double *partial_deriv, void *aero_rep_data);
void * aero_rep_single_particle_get_number_conc(int aero_phase_idx, double *number_conc, 
		double *partial_deriv, void *aero_rep_data);
void * aero_rep_single_particle_get_aero_conc_type(int aero_phase_idx, int *aero_conc_type, void *aero_rep_data);
void * aero_rep_single_particle_update_data(int update_type, void *update_data, void *aero_rep_data);
void * aero_rep_single_particle_print(void *aero_rep_data);
void * aero_rep_single_particle_skip(void *aero_rep_data);

#endif
#endif
