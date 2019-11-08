/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for solver functions
 *
 */

#ifndef CAMP_GPU_SOLVER_H_
#define CAMP_GPU_SOLVER_H_
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include "../camp_common.h"
//#include "../debug_and_stats/camp_debug_2.h"

//Value to consider data size too big -> Memory optimization will change below and under the limit
#define DATA_SIZE_LIMIT_OPT 2000


//Functions to debug cuda errors
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#define HANDLE_ERROR2( ) (HandleError2( __FILE__, __LINE__ ))

//Force solving on CPU: Test option
#define FORCE_CPU 0

void solver_new_gpu_cu(ModelData *model_data, int n_dep_var, int n_state_var, int n_rxn,
     int n_rxn_int_param, int n_rxn_float_param, int n_rxn_env_param, int n_cells);
void solver_set_rxn_data_gpu(ModelData *model_data);
void rxn_update_env_state_gpu(ModelData *model_data);
void rxn_calc_deriv_gpu(ModelData *model_data, N_Vector deriv, realtype time_step);
void rxn_calc_deriv_aux(ModelData *model_data, double *deriv_data, realtype time_step);
void rxn_fusion_deriv_gpu(ModelData *model_data, N_Vector deriv);
void rxn_calc_jac_gpu(ModelData *model_data, SUNMatrix jac, realtype time_step);
void free_gpu_cu(ModelData *model_data);
void bubble_sort_gpu(unsigned int *n_zeros, unsigned int *rxn_position, int n_rxn);
void print_gpu_specs();

#endif
