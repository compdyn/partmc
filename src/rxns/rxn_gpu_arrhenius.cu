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
#include "cuda_util.h"
#include "rxn_gpu_arrhenius.h"

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

__device__ void * rxn_gpu_arrhenius_calc_deriv_contrib( void *rxn_data, 
         ModelDeviceData mdd, realtype *deriv )
{
  realtype *state = model_data->state;
  int *int_data = (int*) rxn_data;
  realtype *float_data = (realtype*) &(int_data[INT_DATA_SIZE_]);

  // Calculate the reaction rate
  realtype rate = RATE_CONSTANT_;
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
