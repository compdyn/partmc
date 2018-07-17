/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for GPU solver functions
 *
*/
/** \file
 * \brief Header file for GPU solver functions
*/
#ifndef PHLEX_GPU_SOLVER_H_
#define PHLEX_GPU_SOLVER_H_
#include <cuda.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "phlex_solver.h"

typedef struct {
  unsigned int num_blocks;      // number of blocks to use during solving
  unsigned int num_threads;     // number of threads to use during solving
  PMC_C_FLOAT * dev_state;      // device pointer to the working state array
  PMC_C_FLOAT * dev_deriv;      // device pointer to the working deriv array
  PMC_C_FLOAT * dev_jac;        // device pointer to the working Jacobian data
  void * rxn_dev_data;          // device reaction data
} ModelDeviceData;

void phlex_gpu_solver_new( ModelData * model_data );
int phlex_gpu_solver_f(realtype t, N_Vector y, N_Vector deriv, void *solver_data);
void phlex_gpu_solver_free(ModelDeviceData * model_device_data);

#endif
