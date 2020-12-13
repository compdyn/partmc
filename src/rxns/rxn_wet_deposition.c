/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Wet deposition reaction solver functions
 *
 */
/** \file
 * \brief Wet deposition reaction solver functions
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../rxns.h"

// TODO Lookup environmental indicies during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define RXN_ID_ (int_data[0])
#define NUM_SPEC_ (int_data[1])
#define SCALING_ float_data[0]
#define RATE_CONSTANT_ (rxn_env_data[0])
#define BASE_RATE_ (rxn_env_data[1])
#define NUM_INT_PROP_ 2
#define NUM_FLOAT_PROP_ 1
#define NUM_ENV_PARAM_ 2
#define REACT_(s) (int_data[NUM_INT_PROP_ + s] - 1)
#define DERIV_ID_(s) int_data[NUM_INT_PROP_ + NUM_SPEC_ + s]
#define JAC_ID_(s) int_data[NUM_INT_PROP_ + 2 * NUM_SPEC_ + s]

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param jac Jacobian
 */
void rxn_wet_deposition_get_used_jac_elem(int *rxn_int_data,
                                          double *rxn_float_data,
                                          Jacobian *jac) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  for (int i_spec = 0; i_spec < NUM_SPEC_; i_spec++) {
    jacobian_register_element(jac, REACT_(i_spec), REACT_(i_spec));
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
void rxn_wet_deposition_update_ids(ModelData *model_data, int *deriv_ids,
                                   Jacobian jac, int *rxn_int_data,
                                   double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  for (int i_spec = 0; i_spec < NUM_SPEC_; i_spec++) {
    // Update the time derivative id
    DERIV_ID_(i_spec) = deriv_ids[REACT_(i_spec)];

    // Update the Jacobian id
    JAC_ID_(i_spec) = jacobian_get_element_id(jac, REACT_(i_spec), REACT_(i_spec));
  }

  return;
}

/** \brief Update reaction data
 *
 * Wet deposition reactions can have their base (pre-scaling) rate constants
 * updated from the host model based on the calculations of an external
 * module. The structure of the update data is:
 *
 *  - \b int rxn_id (Id of one or more wet deposition reactions set by the
 *       host model using the
 *       \c pmc_rxn_wet_deposition::rxn_wet_deposition_t::set_rxn_id
 *       function prior to initializing the solver.)
 *  - \b double rate_const (New pre-scaling rate constant.)
 *
 * \param update_data Pointer to the updated reaction data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent data
 * \return Flag indicating whether this is the reaction to update
 */
bool rxn_wet_deposition_update_data(void *update_data, int *rxn_int_data,
                                    double *rxn_float_data,
                                    double *rxn_env_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  int *rxn_id = (int *)update_data;
  double *base_rate = (double *)&(rxn_id[1]);

  // Set the base wet deposition rate constants for matching reactions
  if (*rxn_id == RXN_ID_ && RXN_ID_ > 0) {
    BASE_RATE_ = (double)*base_rate;
    RATE_CONSTANT_ = SCALING_ * BASE_RATE_;
    return true;
  }

  return false;
}

/** \brief Update reaction data for new environmental conditions
 *
 * For wet deposition reactions this only involves recalculating the rate
 * constant.
 *
 * \param model_data Pointer to the model data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 */
void rxn_wet_deposition_update_env_state(ModelData *model_data,
                                         int *rxn_int_data,
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
void rxn_wet_deposition_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Add contributions to the time derivative
  for (int i_spec = 0; i_spec < NUM_SPEC_; i_spec++) {
    if (DERIV_ID_(i_spec) >= 0) {
      long double rate = RATE_CONSTANT_ * state[REACT_(i_spec)];
      time_derivative_add_value(time_deriv, DERIV_ID_(i_spec), -rate);
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
void rxn_wet_deposition_calc_jac_contrib(ModelData *model_data, Jacobian jac,
                                         int *rxn_int_data,
                                         double *rxn_float_data,
                                         double *rxn_env_data,
                                         realtype time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Add contributions to the Jacobian
  for (int i_spec = 0; i_spec < NUM_SPEC_; i_spec++) {
    if (JAC_ID_(i_spec) >= 0)
      jacobian_add_value(jac, (unsigned int)JAC_ID_(i_spec), JACOBIAN_LOSS,
                         RATE_CONSTANT_);
  }

  return;
}
#endif

/** \brief Print the reaction parameters
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_wet_deposition_print(int *rxn_int_data, double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  printf("\n\nWet deposition reaction\n");

  return;
}

/** \brief Create update data for new wet deposition rates
 *
 * \return Pointer to a new rate update data object
 */
void *rxn_wet_deposition_create_rate_update_data() {
  int *update_data = (int *)malloc(sizeof(int) + sizeof(double));
  if (update_data == NULL) {
    printf("\n\nERROR allocating space for wet deposition update data\n\n");
    exit(1);
  }
  return (void *)update_data;
}

/** \brief Set rate update data
 *
 * \param update_data Pointer to an allocated rate update data object
 * \param rxn_id Id of wet deposition reactions to update
 * \param base_rate New pre-scaling wet deposition rate
 */
void rxn_wet_deposition_set_rate_update_data(void *update_data, int rxn_id,
                                             double base_rate) {
  int *new_rxn_id = (int *)update_data;
  double *new_base_rate = (double *)&(new_rxn_id[1]);
  *new_rxn_id = rxn_id;
  *new_base_rate = base_rate;
}
