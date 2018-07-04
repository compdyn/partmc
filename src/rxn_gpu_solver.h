/* Copyright (C) 2015-2017 Matthew Dawson
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rxn_solver.h"

/* GPU solver data */
typedef struct {
  unsigned int num_rxn;                 // number of rxns with GPU functions
  unsigned int *rxn_data_start;         // id (in bytes) of each rxn's data 
                                        // block in the rxn data
  unsigned int *rxn_data_end;           // last element (in bytes) of each rxn's
                                        // data block
  void *host_rxn_data;                  // host pointer to rxn data
  void *dev_rxn_data;                   // device pointer to rxn data
} GpuData;

// arrhenius
void * rxn_gpu_arrhenius_get_data_size( ModelData *model_data );
void * rxn_gpu_arrhenius_calc_deriv_contrib( void *rxn_data, 
          const realtype *state, realtype *deriv );
#endif
