/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Wet deposition reaction solver functions
 *
*/
/** \file
 * \brief Wet deposition reaction solver functions
*/
extern "C"{
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../rxns_gpu.h"

#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define RXN_ID_ (int_data[0*n_rxn])
#define NUM_SPEC_ (int_data[1*n_rxn])
#define SCALING_ float_data[0*n_rxn]
#define RATE_CONSTANT_ rxn_env_data[0*n_rxn]
#define BASE_RATE_ rxn_env_data[1*n_rxn]//todo fix this shouldnt be there
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 1
#define REACT_(s) (int_data[(NUM_INT_PROP_+s)*n_rxn]-1)
#define DERIV_ID_(s) int_data[(NUM_INT_PROP_+NUM_SPEC_+s)*n_rxn]
#define JAC_ID_(s) int_data[(NUM_INT_PROP_+2*NUM_SPEC_+s)*n_rxn]
#define INT_DATA_SIZE_ (NUM_INT_PROP_+3*NUM_SPEC_)
#define FLOAT_DATA_SIZE_ (NUM_FLOAT_PROP_)

/** \brief Do pre-derivative calculations
 *
 * Nothing to do for wet deposition reactions
 *
 * \param model_data Pointer to the model data, including the state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_gpu_wet_deposition_pre_calc(ModelData *model_data, void *rxn_data)
{
  int n_rxn=1;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Calculate contributions to the time derivative \f$f(t,y)\f$ from
 * this reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param deriv Pointer to the time derivative to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step being computed (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
#ifdef PMC_USE_SUNDIALS
#ifdef __CUDA_ARCH__
__host__ __device__
#endif
void rxn_gpu_wet_deposition_calc_deriv_contrib(ModelData *model_data, realtype *deriv, int *rxn_int_data,
          double *rxn_float_data, double *rxn_env_data, double time_step)
{
#ifdef __CUDA_ARCH__
  int n_rxn=model_data->n_rxn;
#else
  int n_rxn=1;;
#endif
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Add contributions to the time derivative
  for (int i_spec = 0; i_spec < NUM_SPEC_; i_spec++) {
    if (DERIV_ID_(i_spec) >= 0 )
#ifdef __CUDA_ARCH__
        atomicAdd((double*)&(deriv[DERIV_ID_(i_spec)]),-RATE_CONSTANT_ * state[REACT_(i_spec)]);
#else
        deriv[DERIV_ID_(i_spec)] -= RATE_CONSTANT_ * state[REACT_(i_spec)];;
#endif
  }

}
#endif

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_data Pointer to the reaction data
 * \param time_step Current time step being calculated (s)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
 /*
#ifdef PMC_USE_SUNDIALS
#ifdef __CUDA_ARCH__
__host__ __device__
#endif
void rxn_gpu_wet_deposition_calc_jac_contrib(ModelData *model_data, realtype *J, int *rxn_int_data,
          double *rxn_float_data, double *rxn_env_data, double time_step)
{
#ifdef __CUDA_ARCH__
  int n_rxn=model_data->n_rxn;
#else
  int n_rxn=1;
#endif
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Add contributions to the Jacobian
  for (int i_spec = 0; i_spec < NUM_SPEC_; i_spec++) {
    if (JAC_ID_(i_spec) >= 0)
#ifdef __CUDA_ARCH__
        atomicAdd(&(J[JAC_ID_(i_spec)]),-RATE_CONSTANT_);
#else
        int n_rxn=1;
#endif
  }

}
#endif
/*
/** \brief Retrieve Int data size
 *
 * \param rxn_data Pointer to the reaction data
 * \return The data size of int array
 */
void * rxn_gpu_wet_deposition_get_float_pointer(void *rxn_data)
{
  int n_rxn=1;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) float_data;
}

/** \brief Advance the reaction data pointer to the next reaction
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_gpu_wet_deposition_skip(void *rxn_data)
{
  int n_rxn=1;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Print the reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_gpu_wet_deposition_print(void *rxn_data)
{
  int n_rxn=1;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  printf("\n\nWet deposition reaction\n");
  for (int i=0; i<INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Create update data for new wet deposition rates
 *
 * \return Pointer to a new rate update data object
 */
void * rxn_gpu_wet_deposition_create_rate_update_data()
{
  int *update_data = (int*) malloc(sizeof(int) + sizeof(double));
  if (update_data==NULL) {
    printf("\n\nERROR allocating space for wet deposition update data\n\n");
    exit(1);
  }
  return (void*) update_data;
}

/** \brief Set rate update data
 *
 * \param update_data Pointer to an allocated rate update data object
 * \param rxn_id Id of wet deposition reactions to update
 * \param base_rate New pre-scaling wet deposition rate
 */
void rxn_gpu_wet_deposition_set_rate_update_data(void *update_data, int rxn_id,
          double base_rate)
{
  int *new_rxn_id = (int*) update_data;
  double *new_base_rate = (double*) &(new_rxn_id[1]);
  *new_rxn_id = rxn_id;
  *new_base_rate = base_rate;
}

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef RXN_ID_
#undef NUM_SPEC_
#undef BASE_RATE_
#undef SCALING_
#undef RATE_CONSTANT_
#undef NUM_INT_PROP_
#undef NUM_FLOAT_PROP_
#undef REACT_
#undef DERIV_ID_
#undef JAC_ID_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_
}