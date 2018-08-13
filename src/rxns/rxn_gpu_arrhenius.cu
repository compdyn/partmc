/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * GPU solver functions for Arrhenius reactions
 *
*/
/** \file
 * \brief GPU solver functions for Arrhenius reactions
*/
#include "../cuda_util.h"
extern "C" {
#include "../rxn_gpu_solver.h"
#include "rxn_gpu_arrhenius.h"
}

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define NUM_REACT_ int_data[0]
#define NUM_PROD_ int_data[1]
#define INT_DATA_SIZE_ int_data[2]
#define FLOAT_DATA_SIZE_ int_data[3]
#define A_ float_data[0]
#define B_ float_data[1]
#define C_ float_data[2]
#define D_ float_data[3]
#define E_ float_data[4]
#define CONV_ float_data[5]
#define RATE_CONSTANT_ float_data[6]
#define NUM_INT_PROP_ 4
#define NUM_FLOAT_PROP_ 7
#define REACT_(x) (int_data[NUM_INT_PROP_ + x]-1)
#define PROD_(x) (int_data[NUM_INT_PROP_ + NUM_REACT_ + x]-1)
#define DERIV_ID_(x) int_data[NUM_INT_PROP_ + NUM_REACT_ + NUM_PROD_ + x]
#define JAC_ID_(x) int_data[NUM_INT_PROP_ + 2*(NUM_REACT_+NUM_PROD_) + x]
#define YIELD_(x) float_data[NUM_FLOAT_PROP_ + x]

/** \brief Copy reaction data to a new reaction data array
  *
  * \param orig_rxn_data Reaction data to copy
  * \param new_rxn_data Destination for reaction data
  */
extern "C"
__host__ void rxn_gpu_arrhenius_copy_data( void *orig_rxn_data,
          void *new_rxn_data )
{
    int *int_data = (int*) orig_rxn_data;
    PMC_C_FLOAT *float_data = (PMC_C_FLOAT*) &(int_data[INT_DATA_SIZE_]);
    int *new_int_data = (int*) new_rxn_data;
    PMC_C_FLOAT *new_float_data = (PMC_C_FLOAT*) &(new_int_data[INT_DATA_SIZE_]);

    for( int i_var = 0; i_var < INT_DATA_SIZE_; i_var++ )
          new_int_data[ i_var ] = int_data[ i_var ];
    for( int i_var = 0; i_var < FLOAT_DATA_SIZE_; i_var++ )
          new_float_data[ i_var ] = float_data[ i_var ];
}

/** \brief Update reaction data for new environmental conditions
 *
 * For Arrhenius reaction this only involves recalculating the rate 
 * constant.
 *
 * \param mdd Model device data (contains updated environmental state)
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
extern "C"
__device__ void * rxn_gpu_arrhenius_update_env_state( void *rxn_data,
         ModelDeviceData mdd )
{
  int *int_data = (int*) rxn_data;
  PMC_C_FLOAT *float_data = (PMC_C_FLOAT*) &(int_data[INT_DATA_SIZE_]);

  PMC_C_FLOAT *env_data = mdd.dev_env;

  // Calculate the rate constant in (#/cc)
  // k = A*exp(C/T) * (T/D)^B * (1+E*P)
  RATE_CONSTANT_ = A_ * exp(C_/TEMPERATURE_K_)
	  * (B_==0.0 ? 1.0 : pow(TEMPERATURE_K_/D_, B_))
	  * (E_==0.0 ? 1.0 : (1.0 + E_*PRESSURE_PA_))
          * pow(CONV_*PRESSURE_PA_/TEMPERATURE_K_, NUM_REACT_-1);

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);
}

/** \brief Calculate the derivative f(t,y) for this reaction
  *
  * \param rxn_data Pointer to the reaction data
  * \param mdd Model device data
  * \param deriv Pointer to derivative array
  * \return The rxn_data pointer advanced by the size of the reaction data
  */
extern "C"
__device__ void * rxn_gpu_arrhenius_calc_deriv_contrib( void *rxn_data, 
         ModelDeviceData mdd, PMC_SOLVER_C_FLOAT *deriv )
{
  PMC_C_FLOAT *state = mdd.dev_state;
  int *int_data = (int*) rxn_data;
  PMC_C_FLOAT *float_data = (PMC_C_FLOAT*) &(int_data[INT_DATA_SIZE_]);

  // Calculate the reaction rate
  PMC_C_FLOAT rate = RATE_CONSTANT_;
  for (int i_spec=0; i_spec<NUM_REACT_; i_spec++) rate *= state[REACT_(i_spec)];

  // Add contributions to the time derivative
  if (rate!=ZERO) {
    int i_dep_var = 0;
    for (int i_spec=0; i_spec<NUM_REACT_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue; 
      atomicAdd( &( deriv[DERIV_ID_(i_dep_var)] ), -rate );
    }
    for (int i_spec=0; i_spec<NUM_PROD_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue; 
      atomicAdd( &( deriv[DERIV_ID_(i_dep_var)] ), rate*YIELD_(i_spec) );
    }
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param rxn_data Pointer to the reaction data
 * \param mdd Model device data
 * \param J Pointer to the Jacobian data array
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
extern "C"
__device__ void * rxn_gpu_arrhenius_calc_jac_contrib( void *rxn_data, 
         ModelDeviceData mdd, PMC_SOLVER_C_FLOAT *J )
{
  PMC_C_FLOAT *state = mdd.dev_state;
  int *int_data = (int*) rxn_data;
  PMC_C_FLOAT *float_data = (PMC_C_FLOAT*) &(int_data[INT_DATA_SIZE_]);

  // Calculate the reaction rate
  PMC_C_FLOAT rate = RATE_CONSTANT_;
  for (int i_spec=0; i_spec<NUM_REACT_; i_spec++) rate *= state[REACT_(i_spec)];

  // Add contributions to the Jacobian
  if (rate!=ZERO) {
    int i_elem = 0;
    for (int i_ind=0; i_ind<NUM_REACT_; i_ind++) {
      for (int i_dep=0; i_dep<NUM_REACT_; i_dep++, i_elem++) {
	if (JAC_ID_(i_elem) < 0) continue;
	atomicAdd( &( J[JAC_ID_(i_elem)] ), (PMC_SOLVER_C_FLOAT) 
                (-rate / state[REACT_(i_ind)]));
      }
      for (int i_dep=0; i_dep<NUM_PROD_; i_dep++, i_elem++) {
	if (JAC_ID_(i_elem) < 0) continue;
	atomicAdd( &( J[JAC_ID_(i_elem)] ), (PMC_SOLVER_C_FLOAT)
                (YIELD_(i_dep) * rate / state[REACT_(i_ind)]));
      }
    }
  }

  return (void*) &(float_data[FLOAT_DATA_SIZE_]);

}

#undef TEMPERATURE_K_
#undef PRESSURE_PA_

#undef NUM_REACT_
#undef NUM_PROD_
#undef A_
#undef B_
#undef C_
#undef D_
#undef E_
#undef CONV_
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
