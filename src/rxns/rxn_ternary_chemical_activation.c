/* Copyright (C) 2015-2021 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Ternary Chemical Activation reaction solver functions
 *
 */
/** \file
 * \brief Ternary Chemical Activation reaction solver functions
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
#define K0_A_ float_data[0]
#define K0_B_ float_data[1]
#define K0_C_ float_data[2]
#define KINF_A_ float_data[3]
#define KINF_B_ float_data[4]
#define KINF_C_ float_data[5]
#define FC_ float_data[6]
#define N_ float_data[7]
#define SCALING_ float_data[8]
#define CONV_ float_data[9]
#define RATE_CONSTANT_ (rxn_env_data[0])
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 10
#define REACT_(x) (int_data[NUM_INT_PROP_ + x] - 1)
#define PROD_(x) (int_data[NUM_INT_PROP_ + NUM_REACT_ + x] - 1)
#define DERIV_ID_(x) int_data[NUM_INT_PROP_ + NUM_REACT_ + NUM_PROD_ + x]
#define JAC_ID_(x) int_data[NUM_INT_PROP_ + 2 * (NUM_REACT_ + NUM_PROD_) + x]
#define YIELD_(x) float_data[NUM_FLOAT_PROP_ + x]

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param jac Jacobian
 */
void rxn_ternary_chemical_activation_get_used_jac_elem(int *rxn_int_data,
                                                       double *rxn_float_data,
                                                       Jacobian *jac) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  for (int i_ind = 0; i_ind < NUM_REACT_; i_ind++) {
    for (int i_dep = 0; i_dep < NUM_REACT_; i_dep++) {
      jacobian_register_element(jac, REACT_(i_dep), REACT_(i_ind));
    }
    for (int i_dep = 0; i_dep < NUM_PROD_; i_dep++) {
      jacobian_register_element(jac, PROD_(i_dep), REACT_(i_ind));
    }
  }

  return;
}

/** \brief Update the time derivative and Jacbobian array indices
 *
 * \param model_data Pointer to the model data
 * \param deriv_ids Id of each state variable in the derivative array
 * \param jac Jacobian
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_ternary_chemical_activation_update_ids(ModelData *model_data,
                                                int *deriv_ids, Jacobian jac,
                                                int *rxn_int_data,
                                                double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  // Update the time derivative ids
  for (int i = 0; i < NUM_REACT_; i++) DERIV_ID_(i) = deriv_ids[REACT_(i)];
  for (int i = 0; i < NUM_PROD_; i++)
    DERIV_ID_(i + NUM_REACT_) = deriv_ids[PROD_(i)];

  // Update the Jacobian ids
  int i_jac = 0;
  for (int i_ind = 0; i_ind < NUM_REACT_; i_ind++) {
    for (int i_dep = 0; i_dep < NUM_REACT_; i_dep++) {
      JAC_ID_(i_jac++) =
          jacobian_get_element_id(jac, REACT_(i_dep), REACT_(i_ind));
    }
    for (int i_dep = 0; i_dep < NUM_PROD_; i_dep++) {
      JAC_ID_(i_jac++) =
          jacobian_get_element_id(jac, PROD_(i_dep), REACT_(i_ind));
    }
  }
  return;
}

/** \brief Update reaction data for new environmental conditions
 *
 * For Ternary Chemical Activation reaction this only involves recalculating the
 * rate constant.
 *
 * \param model_data Pointer to the model data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 */
void rxn_ternary_chemical_activation_update_env_state(ModelData *model_data,
                                                      int *rxn_int_data,
                                                      double *rxn_float_data,
                                                      double *rxn_env_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *env_data = model_data->grid_cell_env;

  // Calculate the rate constant in (#/cc)
  // k = (k0 / (1 + k0[M]/kinf)) * Fc^(1/(1+(1/N*log(k0[M]/kinf))^2))
  double conv = CONV_ * PRESSURE_PA_ / TEMPERATURE_K_;
  double k0 =
      K0_A_ * (K0_C_ == 0.0 ? 1.0 : exp(K0_C_ / TEMPERATURE_K_)) *
      (K0_B_ == 0.0 ? 1.0 : pow(TEMPERATURE_K_ / ((double)300.0), K0_B_)) *
      conv;
  double kinf =
      k0 /
      (KINF_A_ * (KINF_C_ == 0.0 ? 1.0 : exp(KINF_C_ / TEMPERATURE_K_)) *
       (KINF_B_ == 0.0 ? 1.0 : pow(TEMPERATURE_K_ / ((double)300.0), KINF_B_)));
  RATE_CONSTANT_ = 1.0e-6 * (k0 / (1.0 + kinf)) *
                   pow(FC_, (1.0 / (1.0 + pow(log10(kinf) / N_, 2)))) *
                   pow(conv, NUM_REACT_ - 1) * SCALING_;

  return;
}

/** \brief Calculate contributions to the time derivative \f$f(t,y)\f$ from
 * this reaction.
 *
 * \param model_data Pointer to the model data, including the state array
 * \param time_deriv TimeDerivative object
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 * \param time_step Current time step being computed (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_ternary_chemical_activation_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate the reaction rate
  long double rate = RATE_CONSTANT_;
  for (int i_spec = 0; i_spec < NUM_REACT_; i_spec++)
    rate *= state[REACT_(i_spec)];

  // Add contributions to the time derivative
  if (rate != ZERO) {
    int i_dep_var = 0;
    for (int i_spec = 0; i_spec < NUM_REACT_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue;
      time_derivative_add_value(time_deriv, DERIV_ID_(i_dep_var), -rate);
    }
    for (int i_spec = 0; i_spec < NUM_PROD_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue;
      // Negative yields are allowed, but prevented from causing negative
      // concentrations that lead to solver failures
      if (-rate * YIELD_(i_spec) * time_step <= state[PROD_(i_spec)]) {
        time_derivative_add_value(time_deriv, DERIV_ID_(i_dep_var),
                                  rate * YIELD_(i_spec));
      }
    }
  }

  return;
}
#endif

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param jac Reaction Jacobian
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 * \param time_step Current time step being calculated (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_ternary_chemical_activation_calc_jac_contrib(
    ModelData *model_data, Jacobian jac, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Add contributions to the Jacobian
  int i_elem = 0;
  for (int i_ind = 0; i_ind < NUM_REACT_; i_ind++) {
    // Calculate d_rate / d_i_ind
    realtype rate = RATE_CONSTANT_;
    for (int i_spec = 0; i_spec < NUM_REACT_; i_spec++)
      if (i_ind != i_spec) rate *= state[REACT_(i_spec)];

    for (int i_dep = 0; i_dep < NUM_REACT_; i_dep++, i_elem++) {
      if (JAC_ID_(i_elem) < 0) continue;
      jacobian_add_value(jac, (unsigned int)JAC_ID_(i_elem), JACOBIAN_LOSS,
                         rate);
    }
    for (int i_dep = 0; i_dep < NUM_PROD_; i_dep++, i_elem++) {
      if (JAC_ID_(i_elem) < 0) continue;
      // Negative yields are allowed, but prevented from causing negative
      // concentrations that lead to solver failures
      if (-rate * state[REACT_(i_ind)] * YIELD_(i_dep) * time_step <=
          state[PROD_(i_dep)]) {
        jacobian_add_value(jac, (unsigned int)JAC_ID_(i_elem),
                           JACOBIAN_PRODUCTION, YIELD_(i_dep) * rate);
      }
    }
  }

  return;
}
#endif

/** \brief Print the Ternary Chemical Activation reaction parameters
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_ternary_chemical_activation_print(int *rxn_int_data,
                                           double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  printf("\n\nTernary Chemical Activation reaction\n");

  return;
}
