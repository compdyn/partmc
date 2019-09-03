/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * CMAQ_H2O2 reaction solver functions
 *
*/
/** \file
 * \brief CMAQ_H2O2 reaction solver functions
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../rxns.h"

// TODO Lookup environmental indicies during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define NUM_REACT_ int_data[0]
#define NUM_PROD_ int_data[1]
#define k1_A_ float_data[0]
#define k1_B_ float_data[1]
#define k1_C_ float_data[2]
#define k2_A_ float_data[3]
#define k2_B_ float_data[4]
#define k2_C_ float_data[5]
#define CONV_ float_data[6]
#define RATE_CONSTANT_ float_data[7]
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 8
#define REACT_(x) (int_data[NUM_INT_PROP_ + x]-1)
#define PROD_(x) (int_data[NUM_INT_PROP_ + NUM_REACT_ + x]-1)
#define DERIV_ID_(x) int_data[NUM_INT_PROP_ + NUM_REACT_ + NUM_PROD_ + x]
#define JAC_ID_(x) int_data[NUM_INT_PROP_ + 2*(NUM_REACT_+NUM_PROD_) + x]
#define YIELD_(x) float_data[NUM_FLOAT_PROP_ + x]

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param jac_struct 2D array of flags indicating potentially non-zero
 *                   Jacobian elements
 */
void rxn_CMAQ_H2O2_get_used_jac_elem(int *rxn_int_data,
    double *rxn_float_data, bool **jac_struct)
{
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  for (int i_ind = 0; i_ind < NUM_REACT_; i_ind++) {
    for (int i_dep = 0; i_dep < NUM_REACT_; i_dep++) {
      jac_struct[REACT_(i_dep)][REACT_(i_ind)] = true;
    }
    for (int i_dep = 0; i_dep < NUM_PROD_; i_dep++) {
      jac_struct[PROD_(i_dep)][REACT_(i_ind)] = true;
    }
  }

  return;
}

/** \brief Update the time derivative and Jacbobian array indices
 *
 * \param model_data Pointer to the model data
 * \param deriv_ids Id of each state variable in the derivative array
 * \param jac_ids Id of each state variable combo in the Jacobian array
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_CMAQ_H2O2_update_ids(ModelData *model_data, int *deriv_ids,
          int **jac_ids, int *rxn_int_data, double *rxn_float_data)
{
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  // Update the time derivative ids
  for (int i=0; i < NUM_REACT_; i++)
	  DERIV_ID_(i) = deriv_ids[REACT_(i)];
  for (int i=0; i < NUM_PROD_; i++)
	  DERIV_ID_(i + NUM_REACT_) = deriv_ids[PROD_(i)];

  // Update the Jacobian ids
  int i_jac = 0;
  for (int i_ind = 0; i_ind < NUM_REACT_; i_ind++) {
    for (int i_dep = 0; i_dep < NUM_REACT_; i_dep++) {
      JAC_ID_(i_jac++) = jac_ids[REACT_(i_dep)][REACT_(i_ind)];
    }
    for (int i_dep = 0; i_dep < NUM_PROD_; i_dep++) {
      JAC_ID_(i_jac++) = jac_ids[PROD_(i_dep)][REACT_(i_ind)];
    }
  }
  return;
}

/** \brief Update reaction data for new environmental conditions
 *
 * For CMAQ_H2O2 reaction this only involves recalculating the rate
 * constant.
 *
 * \param env_data Pointer to the environmental state array
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_CMAQ_H2O2_update_env_state(double *env_data, int *rxn_int_data,
    double *rxn_float_data)
{
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  // Calculate the rate constant in (#/cc)
  // k = k1 + [M]*k2
  double conv = CONV_ * PRESSURE_PA_ / TEMPERATURE_K_;
  RATE_CONSTANT_ = (k1_A_
	  * (k1_C_==0.0 ? 1.0 : exp(k1_C_/TEMPERATURE_K_))
	  * (k1_B_==0.0 ? 1.0 :
                  pow(TEMPERATURE_K_/((double)300.0), k1_B_))
	  + k2_A_ // [M] is included in k2_A_
	  * (k2_C_==0.0 ? 1.0 : exp(k2_C_/TEMPERATURE_K_))
	  * (k2_B_==0.0 ? 1.0 :
                  pow(TEMPERATURE_K_/((double)300.0), k2_B_))
	  * conv
	  ) * pow(conv, NUM_REACT_-1);

  return;
}

/** \brief Calculate contributions to the time derivative \f$f(t,y)\f$ from
 * this reaction.
 *
 * \param model_data Pointer to the model data
 * \param deriv Pointer to the time derivative to add contributions to
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param time_step Current time step being computed (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_CMAQ_H2O2_calc_deriv_contrib(ModelData *model_data, realtype *deriv,
          int *rxn_int_data, double *rxn_float_data, double time_step)
{
  realtype *state = model_data->state;
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  // Calculate the reaction rate
  realtype rate = RATE_CONSTANT_;
  for (int i_spec=0; i_spec<NUM_REACT_; i_spec++) rate *= state[REACT_(i_spec)];

  // Add contributions to the time derivative
  if (rate!=ZERO) {
    int i_dep_var = 0;
    for (int i_spec=0; i_spec<NUM_REACT_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue;
      deriv[DERIV_ID_(i_dep_var)] -= rate;
    }
    for (int i_spec=0; i_spec<NUM_PROD_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue;
      // Negative yields are allowed, but prevented from causing negative
      // concentrations that lead to solver failures
      if (-rate*YIELD_(i_spec)*time_step <= state[PROD_(i_spec)]) {
        deriv[DERIV_ID_(i_dep_var)] += rate*YIELD_(i_spec);
      }
    }
  }

  return;

}
#endif

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param time_step Current time step being calculated (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_CMAQ_H2O2_calc_jac_contrib(ModelData *model_data, realtype *J,
          int *rxn_int_data, double *rxn_float_data, double time_step)
{
  realtype *state = model_data->state;
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  // Add contributions to the Jacobian
  int i_elem = 0;
  for (int i_ind=0; i_ind<NUM_REACT_; i_ind++) {

    // Calculate d_rate / d_i_ind
    realtype rate = RATE_CONSTANT_;
    for (int i_spec=0; i_spec<NUM_REACT_; i_spec++)
      if (i_ind!=i_spec) rate *= state[REACT_(i_spec)];

    for (int i_dep=0; i_dep<NUM_REACT_; i_dep++, i_elem++) {
      if (JAC_ID_(i_elem) < 0) continue;
      J[JAC_ID_(i_elem)] -= rate;
    }
    for (int i_dep=0; i_dep<NUM_PROD_; i_dep++, i_elem++) {
      if (JAC_ID_(i_elem) < 0) continue;
      // Negative yields are allowed, but prevented from causing negative
      // concentrations that lead to solver failures
      if (-rate * state[REACT_(i_ind)] * YIELD_(i_dep) * time_step
          <= state[PROD_(i_dep)]) {
	J[JAC_ID_(i_elem)] += YIELD_(i_dep) * rate;
      }
    }
  }

  return;

}
#endif

/** \brief Print the CMAQ_H2O2 reaction parameters
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_CMAQ_H2O2_print(int *rxn_int_data, double *rxn_float_data)
{
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  printf("\n\nCMAQ_H2O2 reaction\n");

  return;
}
