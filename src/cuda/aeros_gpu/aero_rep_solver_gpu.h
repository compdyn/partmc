/* Copyright (C) 2018 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for aero_rep_gpu_solver.c
 *
 */
/** \file
 * \brief Header file for abstract aerosol representation functions
 */
#ifndef AERO_REP_SOLVER_H
#define AERO_REP_SOLVER_H
#include "../phlex_gpu_solver.h"

/** Public aerosol representation functions **/

/* Solver functions */
int aero_rep_gpu_get_used_jac_elem(ModelDatagpu *model_data, int aero_rep_gpu_idx,
        int aero_phase_gpu_idx, bool *jac_struct);
void * aero_rep_gpu_get_dependencies(ModelDatagpu *model_data, bool *state_flags);
void aero_rep_gpu_update_env_state(ModelDatagpu *model_data, double *env);
void aero_rep_gpu_update_state(ModelDatagpu *model_data);
__device__ void * aero_rep_gpu_get_effective_radius(ModelDatagpu *model_data, int aero_rep_gpu_idx,
          int aero_phase_gpu_idx, double *radius);
__device__ void * aero_rep_gpu_get_number_conc(ModelDatagpu *model_data, int aero_rep_gpu_idx,
          int aero_phase_gpu_idx, double *number_conc);
__device__ int aero_rep_gpu_get_aero_conc_type(ModelDatagpu *model_data, int aero_rep_gpu_idx,
          int aero_phase_gpu_idx);
void * aero_rep_gpu_get_aero_phase_gpu_mass(ModelDatagpu *model_data, int aero_rep_gpu_idx,
          int aero_phase_gpu_idx, double *aero_phase_gpu_mass,
          double *aero_phase_gpu_avg_MW);
void aero_rep_gpu_print_data(void *solver_data);

/* Setup functions */
void aero_rep_gpu_add_condensed_data(int aero_rep_gpu_type, int n_int_param,
	  int n_float_param, int *int_param, double *float_param,
          void *solver_data);

/* Update data functions */
void aero_rep_gpu_update_data(int update_aero_rep_gpu_type, void *update_data,
          void *solver_data);
void aero_rep_gpu_free_update_data(void *update_data);

#endif
