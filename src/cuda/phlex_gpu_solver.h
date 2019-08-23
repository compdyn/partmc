
#ifndef PHLEX_GPU_SOLVER_H_
#define PHLEX_GPU_SOLVER_H_
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include "../phlex_common.h"

//Value to consider data size too big -> Memory optimization will change below and under the limit
#define DATA_SIZE_LIMIT_OPT 2000

//Functions to debug cuda errors
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#define HANDLE_ERROR2( ) (HandleError2( __FILE__, __LINE__ ))

void solver_new_gpu_cu(int n_dep_var, int n_state_var, int n_rxn,
     int n_rxn_int_param, int n_rxn_float_param, int n_cells);
void allocate_jac_gpu(int n_jac_elem, int n_cells);
void rxn_update_env_state_gpu(ModelData *model_data, double *env);
void rxn_calc_deriv_gpu(ModelData *model_data, N_Vector deriv, realtype time_step);
void rxn_calc_jac_gpu(ModelData *model_data, SUNMatrix jac, realtype time_step);
void free_gpu_cu();
void bubble_sort_gpu(unsigned int *n_zeros, unsigned int *rxn_position, int n_rxn);
void print_gpu_specs();
void solver_set_rxn_data_gpu(ModelData *model_data);

#endif
