/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Photolysis reaction solver functions
 *
 */
/** \file
 * \brief Photolysis reaction solver functions
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
#define RXN_ID_ int_data[2]
#define SCALING_ float_data[0]
#define RATE_CONSTANT_ (rxn_env_data[0])
#define BASE_RATE_ (rxn_env_data[1])
#define NUM_INT_PROP_ 3
#define NUM_FLOAT_PROP_ 1
#define NUM_ENV_PARAM_ 2
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
void rxn_photolysis_get_used_jac_elem(int *rxn_int_data, double *rxn_float_data,
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
void rxn_photolysis_update_ids(ModelData *model_data, int *deriv_ids,
                               Jacobian jac, int *rxn_int_data,
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

/** \brief Update reaction data
 *
 * Photolysis reactions can have their base (pre-scaling) rate constants updated
 * from the host model based on the calculations of an external photolysis
 * module. The structure of the update data is:
 *
 *  - \b int photo_id (Id of one or more photolysis reactions set by the host
 *       model using the pmc_rxn_photolysis::rxn_photolysis_t::set_photo_id
 *       function prior to initializing the solver.)
 *  - \b double rate_const (New pre-scaling rate constant.)
 *
 * \param update_data Pointer to the updated reaction data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent data
 * \return Flag indicating whether this is the reaction to update
 */
bool rxn_photolysis_update_data(void *update_data, int *rxn_int_data,
                                double *rxn_float_data, double *rxn_env_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  int *photo_id = (int *)update_data;
  double *base_rate = (double *)&(photo_id[1]);

  // Set the base photolysis rate constants for matching reactions
  if (*photo_id == RXN_ID_ && RXN_ID_ > 0) {
    BASE_RATE_ = (double)*base_rate;
    RATE_CONSTANT_ = SCALING_ * BASE_RATE_;
    return true;
  }

  return false;
}

/** \brief Update reaction data for new environmental conditions
 *
 * For Photolysis reaction this only involves recalculating the rate
 * constant.
 *
 * \param model_data Pointer to the model data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 */
void rxn_photolysis_update_env_state(ModelData *model_data, int *rxn_int_data,
                                     double *rxn_float_data,
                                     double *rxn_env_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *env_data = model_data->grid_cell_env;

  // Calculate the rate constant in (1/s)
  RATE_CONSTANT_ = SCALING_ * BASE_RATE_;

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
void rxn_photolysis_calc_deriv_contrib(
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
void rxn_photolysis_calc_jac_contrib(ModelData *model_data, Jacobian jac,
                                     int *rxn_int_data, double *rxn_float_data,
                                     double *rxn_env_data, realtype time_step) {
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
      if (i_spec != i_ind) rate *= state[REACT_(i_spec)];

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

/** \brief Print the Photolysis reaction parameters
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_photolysis_print(int *rxn_int_data, double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  printf("\n\nPhotolysis reaction\n");

  return;
}

/** \brief Create update data for new photolysis rates
 *
 * \return Pointer to a new rate update data object
 */
void *rxn_photolysis_create_rate_update_data() {
  int *update_data = (int *)malloc(sizeof(int) + sizeof(double));
  if (update_data == NULL) {
    printf("\n\nERROR allocating space for photolysis update data\n\n");
    exit(1);
  }
  return (void *)update_data;
}

/** \brief Set rate update data
 *
 * \param update_data Pointer to an allocated rate update data object
 * \param photo_id Id of photolysis reactions to update
 * \param base_rate New pre-scaling photolysis rate
 */
void rxn_photolysis_set_rate_update_data(void *update_data, int photo_id,
                                         double base_rate) {
  int *new_photo_id = (int *)update_data;
  double *new_base_rate = (double *)&(new_photo_id[1]);
  *new_photo_id = photo_id;
  *new_base_rate = base_rate;
}
