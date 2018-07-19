/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for GPU solver functions for Arrhenius reactions
 *
*/
/** \file
 * \brief Header file for GPU solver functions for Arrhenius reactions
*/
#ifndef RXN_GPU_ARRHENIUS_H_
#define RXN_GPU_ARRHENIUS_H_
#include <cuda.h>
#include "../phlex_gpu_solver.h"

__host__ void rxn_gpu_arrhenius_copy_data( void *orig_rxn_data,
         void *new_rxn_data );
__device__ void * rxn_gpu_arrhenius_update_env_state( void *rxn_data,
         ModelDeviceData mdd );
__device__ void * rxn_gpu_arrhenius_calc_deriv_contrib( void *rxn_data, 
         ModelDeviceData mdd, PMC_SOLVER_C_FLOAT *deriv );

#endif
