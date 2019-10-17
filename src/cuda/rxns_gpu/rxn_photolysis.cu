/* Copyright (C) 2019 Christian Guzman
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Photolysis reaction solver functions
 *
*/
/** \file
 * \brief Photolysis reaction solver functions
*/
extern "C"{
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../rxns_gpu.h"

#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define NUM_REACT_ int_data[0*n_rxn]
#define NUM_PROD_ int_data[1*n_rxn]
#define PHOTO_ID_ int_data[2*n_rxn]
#define SCALING_ float_data[0*n_rxn]
#define RATE_CONSTANT_ rate_constants[0*n_rxn]
#define BASE_RATE_ rate_constants[0*n_rxn]//todo fix this shouldnt be there
#define NUM_INT_PROP_ 3
#define NUM_FLOAT_PROP_ 1
#define REACT_(x) (int_data[(NUM_INT_PROP_ + x)*n_rxn]-1)
#define PROD_(x) (int_data[(NUM_INT_PROP_ + NUM_REACT_ + x)*n_rxn]-1)
#define DERIV_ID_(x) int_data[(NUM_INT_PROP_ + NUM_REACT_ + NUM_PROD_ + x)*n_rxn]
#define JAC_ID_(x) int_data[(NUM_INT_PROP_ + 2*(NUM_REACT_+NUM_PROD_) + x)*n_rxn]
#define YIELD_(x) float_data[(NUM_FLOAT_PROP_ + x)*n_rxn]
#define INT_DATA_SIZE_ (NUM_INT_PROP_+(NUM_REACT_+2)*(NUM_REACT_+NUM_PROD_))
#define FLOAT_DATA_SIZE_ (NUM_FLOAT_PROP_+NUM_PROD_)

/** \brief Update reaction data for new environmental conditions
 *
 * For Photolysis reaction this only involves recalculating the rate
 * constant.
 *
 * \param env_data Pointer to the environmental state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
__device__ void rxn_gpu_photolysis_update_env_state(double *rate_constants,
   int n_rxn2,double *double_pointer_gpu, double *env_data, void *rxn_data)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  // Calculate the rate constant in (1/s)
  RATE_CONSTANT_ = SCALING_ * BASE_RATE_;

  rate_constants[0] = RATE_CONSTANT_;
}

/** \brief Do pre-derivative calculations
 *
 * Nothing to do for photolysis reactions
 *
 * \param model_data Pointer to the model data, including the state array
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_gpu_photolysis_pre_calc(ModelData *model_data, void *rxn_data)
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
__device__ void rxn_gpu_photolysis_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn2)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  // Calculate the reaction rate
  //double rate = RATE_CONSTANT_;
  double rate = rate_constants[0];
  for (int i_spec=0; i_spec<NUM_REACT_; i_spec++)
          rate *= state[REACT_(i_spec)];

  // Add contributions to the time derivative
  if (rate!=ZERO) {
    int i_dep_var = 0;
    for (int i_spec=0; i_spec<NUM_REACT_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue;
      //deriv[DERIV_ID_(i_dep_var)] -= rate;
      atomicAdd((double*)&(deriv[DERIV_ID_(i_dep_var)]),-rate);
    }
    for (int i_spec=0; i_spec<NUM_PROD_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue;
      //deriv[DERIV_ID_(i_dep_var)] += rate*YIELD_(i_spec);
        atomicAdd((double*)&(deriv[DERIV_ID_(i_dep_var)]),rate*YIELD_(i_spec));
    }
  }

}
#endif


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
void rxn_cpu_photolysis_calc_deriv_contrib(double *rate_constants, double *state,
          double *deriv, void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn2)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  // Calculate the reaction rate
  //double rate = RATE_CONSTANT_;
  double rate = rate_constants[0];
  for (int i_spec=0; i_spec<NUM_REACT_; i_spec++)
          rate *= state[REACT_(i_spec)];

  // Add contributions to the time derivative
  if (rate!=ZERO) {
    int i_dep_var = 0;
    for (int i_spec=0; i_spec<NUM_REACT_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue;
      deriv[DERIV_ID_(i_dep_var)] -= rate;
    }
    for (int i_spec=0; i_spec<NUM_PROD_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue;
      deriv[DERIV_ID_(i_dep_var)] += rate*YIELD_(i_spec);
    }
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
#ifdef PMC_USE_SUNDIALS
__device__ void rxn_gpu_photolysis_calc_jac_contrib(double *rate_constants, double *state, double *J,
          void *rxn_data, double * double_pointer_gpu, double time_step, int n_rxn2)
{
  int n_rxn=n_rxn2;
  int *int_data = (int*) rxn_data;
  double *float_data = double_pointer_gpu;

  // Calculate the reaction rate
  //double rate = RATE_CONSTANT_;
  double rate = rate_constants[0];
  for (int i_spec=0; i_spec<NUM_REACT_; i_spec++)
          rate *= state[REACT_(i_spec)];

  // Add contributions to the Jacobian
  int i_elem = 0;
  for (int i_ind=0; i_ind<NUM_REACT_; i_ind++) {
    for (int i_dep=0; i_dep<NUM_REACT_; i_dep++, i_elem++) {
      if (JAC_ID_(i_elem) < 0) continue;
      atomicAdd(&(J[JAC_ID_(i_elem)]), -rate);
    }
    for (int i_dep=0; i_dep<NUM_PROD_; i_dep++, i_elem++) {
     if (JAC_ID_(i_elem) < 0) continue;
      if (-rate*YIELD_(i_dep)*time_step <= state[PROD_(i_dep)]) {
        atomicAdd(&(J[JAC_ID_(i_elem)]),YIELD_(i_dep) * rate);
      }
    }
  }
  //return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}
#endif

/** \brief Retrieve Int data size
 *
 * \param rxn_data Pointer to the reaction data
 * \return The data size of int array
 */
void * rxn_gpu_photolysis_get_float_pointer(void *rxn_data)
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
void * rxn_gpu_photolysis_skip(void *rxn_data)
{
  int n_rxn=1;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Print the Photolysis reaction parameters
 *
 * \param rxn_data Pointer to the reaction data
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_gpu_photolysis_print(void *rxn_data)
{
  int n_rxn=1;
  int *int_data = (int*) rxn_data;
  double *float_data = (double*) &(int_data[INT_DATA_SIZE_]);

  printf("\n\nPhotolysis reaction\n");
  for (int i=0; i<INT_DATA_SIZE_; i++)
    printf("  int param %d = %d\n", i, int_data[i]);
  for (int i=0; i<FLOAT_DATA_SIZE_; i++)
    printf("  float param %d = %le\n", i, float_data[i]);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Create update data for new photolysis rates
 *
 * \return Pointer to a new rate update data object
 */
void * rxn_gpu_photolysis_create_rate_update_data()
{
  int *update_data = (int*) malloc(sizeof(int) + sizeof(double));
  if (update_data==NULL) {
    printf("\n\nERROR allocating space for photolysis update data\n\n");
    exit(1);
  }
  return (void*) update_data;
}

/** \brief Set rate update data
 *
 * \param update_data Pointer to an allocated rate update data object
 * \param photo_id Id of photolysis reactions to update
 * \param base_rate New pre-scaling photolysis rate
 */
void rxn_gpu_photolysis_set_rate_update_data(void *update_data, int photo_id,
          double base_rate)
{
  int *new_photo_id = (int*) update_data;
  double *new_base_rate = (double*) &(new_photo_id[1]);
  *new_photo_id = photo_id;
  *new_base_rate = base_rate;
}

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef NUM_REACT_
#undef NUM_PROD_
#undef PHOTO_ID_
#undef BASE_RATE_
#undef SCALING_
#undef RATE_CONSTANT_
#undef NUM_INT_PROP_
#undef NUM_FLOAT_PROP_
#undef REACT_
#undef PROD_
#undef DERIV_ID_
#undef JAC_ID_
#undef YIELD_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_
}