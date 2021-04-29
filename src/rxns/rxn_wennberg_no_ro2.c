/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Wennberg NO + RO2 reaction solver functions
 *
 */
/** \file
 * \brief Wennberg NO + RO2 reaction solver functions
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../rxns.h"

// TODO Lookup environmental indices during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define NUM_REACT_ int_data[0]
#define NUM_ALKOXY_PROD_ int_data[1]
#define NUM_NITRATE_PROD_ int_data[2]
#define X_ float_data[0]
#define Y_ float_data[1]
#define a0_ float_data[2]
#define n_ float_data[3]
#define CONV_ float_data[4]
#define ALKOXY_RATE_CONSTANT_ rxn_env_data[0]
#define NITRATE_RATE_CONSTANT_ rxn_env_data[1]
#define NUM_INT_PROP_ 3
#define NUM_FLOAT_PROP_ 5
#define NUM_ENV_PARAM_ 2
#define REACT_(x) (int_data[NUM_INT_PROP_ + x] - 1)
#define PROD_(x) (int_data[NUM_INT_PROP_ + NUM_REACT_ + x] - 1)
#define DERIV_ID_(x)                                                           \
  int_data[NUM_INT_PROP_ + NUM_REACT_ + NUM_ALKOXY_PROD_ + NUM_NITRATE_PROD_ + \
           x]
#define JAC_ID_(x)         \
  int_data[NUM_INT_PROP_ + \
           2 * (NUM_REACT_ + NUM_ALKOXY_PROD_ + NUM_NITRATE_PROD_) + x]
#define YIELD_(x) float_data[NUM_FLOAT_PROP_ + x]

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param jac Jacobian
 */
void rxn_wennberg_no_ro2_get_used_jac_elem(int *rxn_int_data,
                                           double *rxn_float_data,
                                           Jacobian *jac) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  for (int i_ind = 0; i_ind < NUM_REACT_; i_ind++) {
    for (int i_dep = 0; i_dep < NUM_REACT_; i_dep++) {
      jacobian_register_element(jac, REACT_(i_dep), REACT_(i_ind));
    }
    for (int i_dep = 0; i_dep < NUM_ALKOXY_PROD_ + NUM_NITRATE_PROD_; i_dep++) {
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
void rxn_wennberg_no_ro2_update_ids(ModelData *model_data, int *deriv_ids,
                                    Jacobian jac, int *rxn_int_data,
                                    double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  // Update the time derivative ids
  for (int i = 0; i < NUM_REACT_; i++) DERIV_ID_(i) = deriv_ids[REACT_(i)];
  for (int i = 0; i < NUM_ALKOXY_PROD_ + NUM_NITRATE_PROD_; i++)
    DERIV_ID_(i + NUM_REACT_) = deriv_ids[PROD_(i)];

  // Update the Jacobian ids
  int i_jac = 0;
  for (int i_ind = 0; i_ind < NUM_REACT_; i_ind++) {
    for (int i_dep = 0; i_dep < NUM_REACT_; i_dep++)
      JAC_ID_(i_jac++) =
          jacobian_get_element_id(jac, REACT_(i_dep), REACT_(i_ind));
    for (int i_dep = 0; i_dep < NUM_ALKOXY_PROD_ + NUM_NITRATE_PROD_; i_dep++)
      JAC_ID_(i_jac++) =
          jacobian_get_element_id(jac, PROD_(i_dep), REACT_(i_ind));
  }
  return;
}

/** \brief Calculates the Troe-like parameter A(T, [M], n)
 *
 * k0 = 2e-22 e^n
 * kinf = 0.43 * (T/298)^-8
 * A = k0 [M] / ( 1 + k0 [M] / kinf ) *
 *       0.41^( 1 / ( 1 + [ log10( k0 [M] / kinf ) ]^2 ) )
 *
 * \param T temperature [K]
 * \param M air density [#/cc]
 * \param n number of heavy atoms in RO2 species
 */
double calculate_A(double T, double M, double n) {
  double k0, kinf;
  k0 = 2.0e-22 * exp(n);
  kinf = 0.43 * pow(T / 298.0, -8);
  return k0 * M / (1.0 + k0 * M / kinf) *
         pow(0.41, 1.0 / (1.0 + pow(log10(k0 * M / kinf), 2)));
}

/** \brief Update reaction data for new environmental conditions
 *
 * For Wennberg NO + RO2 reaction this only involves recalculating the rate
 * constant.
 *
 * \param model_data Pointer to the model data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 */
void rxn_wennberg_no_ro2_update_env_state(ModelData *model_data,
                                          int *rxn_int_data,
                                          double *rxn_float_data,
                                          double *rxn_env_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *env_data = model_data->grid_cell_env;

  double base_rate, A, Z, M;

  // Calculate the rate constant in (#/cc)
  base_rate = X_ * exp(-Y_ / TEMPERATURE_K_) *
              pow(CONV_ * PRESSURE_PA_ / TEMPERATURE_K_, NUM_REACT_ - 1);
  M = CONV_ * PRESSURE_PA_ / TEMPERATURE_K_ * 1e6;  // [#/cm3]
  Z = calculate_A(293.0, 2.45e19, n_) * (1.0 - a0_) / a0_;
  A = calculate_A(TEMPERATURE_K_, M, n_);
  ALKOXY_RATE_CONSTANT_ = base_rate * Z / (Z + A);
  NITRATE_RATE_CONSTANT_ = base_rate * A / (A + Z);
  return;
}

/** \brief Calculate contributions to the time derivative \f$f(t,y)\f$ from
 * this reaction.
 *
 * \param model_data Model data
 * \param time_deriv TimeDerivative object
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 * \param time_step Current time step being computed (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_wennberg_no_ro2_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, double time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate the reaction rate
  long double rate = 1.0;
  long double k_a = ALKOXY_RATE_CONSTANT_;
  long double k_n = NITRATE_RATE_CONSTANT_;
  for (int i_spec = 0; i_spec < NUM_REACT_; i_spec++)
    rate *= state[REACT_(i_spec)];

  // Add contributions to the time derivative
  if (rate != ZERO) {
    int i_dep_var = 0;
    for (int i_spec = 0; i_spec < NUM_REACT_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue;
      time_derivative_add_value(time_deriv, DERIV_ID_(i_dep_var),
                                -(k_a + k_n) * rate);
    }
    for (int i_spec = 0; i_spec < NUM_ALKOXY_PROD_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue;

      // Negative yields are allowed, but prevented from causing negative
      // concentrations that lead to solver failures
      if (-k_a * rate * YIELD_(i_spec) * time_step <= state[PROD_(i_spec)]) {
        time_derivative_add_value(time_deriv, DERIV_ID_(i_dep_var),
                                  k_a * rate * YIELD_(i_spec));
      }
    }
    for (int i_spec = NUM_ALKOXY_PROD_;
         i_spec < NUM_ALKOXY_PROD_ + NUM_NITRATE_PROD_; i_spec++, i_dep_var++) {
      if (DERIV_ID_(i_dep_var) < 0) continue;

      // Negative yields are allowed, but prevented from causing negative
      // concentrations that lead to solver failures
      if (-k_n * rate * YIELD_(i_spec) * time_step <= state[PROD_(i_spec)]) {
        time_derivative_add_value(time_deriv, DERIV_ID_(i_dep_var),
                                  k_n * rate * YIELD_(i_spec));
      }
    }
  }

  return;
}
#endif

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Model data
 * \param jac Reaction Jacobian
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 * \param time_step Current time step being calculated (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_wennberg_no_ro2_calc_jac_contrib(ModelData *model_data, Jacobian jac,
                                          int *rxn_int_data,
                                          double *rxn_float_data,
                                          double *rxn_env_data,
                                          double time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Add contributions to the Jacobian
  int i_elem = 0;
  for (int i_ind = 0; i_ind < NUM_REACT_; i_ind++) {
    // Calculate d_rate / d_i_ind
    realtype rate = 1.0;
    realtype k_a = ALKOXY_RATE_CONSTANT_;
    realtype k_n = NITRATE_RATE_CONSTANT_;
    for (int i_spec = 0; i_spec < NUM_REACT_; i_spec++)
      if (i_spec != i_ind) rate *= state[REACT_(i_spec)];

    for (int i_dep = 0; i_dep < NUM_REACT_; i_dep++, i_elem++) {
      if (JAC_ID_(i_elem) < 0) continue;
      jacobian_add_value(jac, (unsigned int)JAC_ID_(i_elem), JACOBIAN_LOSS,
                         (k_a + k_n) * rate);
    }
    for (int i_dep = 0; i_dep < NUM_ALKOXY_PROD_; i_dep++, i_elem++) {
      if (JAC_ID_(i_elem) < 0) continue;
      // Negative yields are allowed, but prevented from causing negative
      // concentrations that lead to solver failures
      if (-k_a * rate * state[REACT_(i_ind)] * YIELD_(i_dep) * time_step <=
          state[PROD_(i_dep)]) {
        jacobian_add_value(jac, (unsigned int)JAC_ID_(i_elem),
                           JACOBIAN_PRODUCTION, k_a * YIELD_(i_dep) * rate);
      }
    }
    for (int i_dep = NUM_ALKOXY_PROD_;
         i_dep < NUM_ALKOXY_PROD_ + NUM_NITRATE_PROD_; i_dep++, i_elem++) {
      if (JAC_ID_(i_elem) < 0) continue;
      // Negative yields are allowed, but prevented from causing negative
      // concentrations that lead to solver failures
      if (-k_n * rate * state[REACT_(i_ind)] * YIELD_(i_dep) * time_step <=
          state[PROD_(i_dep)]) {
        jacobian_add_value(jac, (unsigned int)JAC_ID_(i_elem),
                           JACOBIAN_PRODUCTION, k_n * YIELD_(i_dep) * rate);
      }
    }
  }

  return;
}
#endif

/** \brief Print the Wennberg NO + RO2 reaction parameters
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_wennberg_no_ro2_print(int *rxn_int_data, double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  printf("\n\nWennberg NO + RO2 reaction\n");

  return;
}
