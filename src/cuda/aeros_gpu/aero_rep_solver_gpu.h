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

__device__ void * aero_rep_gpu_get_effective_radius(ModelDatagpu *model_data, int aero_rep_gpu_idx,
          int aero_phase_gpu_idx, double *radius);
__device__ void * aero_rep_gpu_get_number_conc(ModelDatagpu *model_data, int aero_rep_gpu_idx,
          int aero_phase_gpu_idx, double *number_conc);
__device__ int aero_rep_gpu_get_aero_conc_type(ModelDatagpu *model_data, int aero_rep_gpu_idx,
          int aero_phase_gpu_idx);
__device__ void * aero_rep_gpu_get_aero_phase_mass(ModelDatagpu *model_data, int aero_rep_idx,
                                    int aero_phase_idx, double *aero_phase_mass,
                                    double *aero_phase_avg_MW);

#endif
