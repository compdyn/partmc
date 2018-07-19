/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for reaction functions
 *
*/
/** \file
 * \brief Header file for reaction solver functions
*/
#ifndef RXN_GPU_SOLVER_H_
#define RXN_GPU_SOLVER_H_
#include <cuda.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
extern "C" {
#include "phlex_gpu_solver.h"
#include "rxn_solver.h"
}

/* GPU solver data */
typedef struct {
  int *host_rxn_data_start;             // id (in bytes) of each rxn's data 
                                        // block in the rxn data on the host
  int *dev_rxn_data_start;              // id (in bytes) of each rxn's data 
                                        // block in the rxn data on the device
  void *host_rxn_data;                  // host pointer to rxn data
  void *dev_rxn_data;                   // device pointer to rxn data
} RxnDeviceData;

void rxn_gpu_solver_new( ModelDeviceData * model_dev_data, void * rxn_data );
__global__ void rxn_gpu_update_env_state( ModelDeviceData mdd );
__global__ void rxn_gpu_calc_deriv( ModelDeviceData mdd, PMC_C_FLOAT time_step); 
void rxn_gpu_solver_print( void * rxn_data );
void rxn_gpu_solver_free( void * rxn_dev_data );

#endif
