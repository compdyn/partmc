/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * First-Order loss reaction solver functions
 *
 */
/** \file
 * \brief First-Order loss reaction solver functions
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../rxns.h"

// TODO Lookup environmental indicies during initialization
#define TEMPERATURE_K_ env_data[0]
#define PRESSURE_PA_ env_data[1]

#define RXN_ID_ (int_data[0])
#define REACT_ (int_data[1] - 1)
#define DERIV_ID_ int_data[2]
#define JAC_ID_ int_data[3]
#define SCALING_ float_data[0]
#define RATE_CONSTANT_ (rxn_env_data[0])
#define BASE_RATE_ (rxn_env_data[1])
#define NUM_INT_PROP_ 4
#define NUM_FLOAT_PROP_ 1
#define NUM_ENV_PARAM_ 2

/** \brief Flag Jacobian elements used by this reaction
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param jac_struct 2D array of flags indicating potentially non-zero
 *                   Jacobian elements
 */
void rxn_first_order_loss_get_used_jac_elem(int *rxn_int_data,
                                            double *rxn_float_data,
                                            bool **jac_struct) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  jac_struct[REACT_][REACT_] = true;

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
void rxn_first_order_loss_update_ids(ModelData *model_data, int *deriv_ids,
                                     int **jac_ids, int *rxn_int_data,
                                     double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  // Update the time derivative id
  DERIV_ID_ = deriv_ids[REACT_];

  // Update the Jacobian id
  JAC_ID_ = jac_ids[REACT_][REACT_];

  return;
}

/** \brief Update reaction data
 *
 * First-Order loss reactions can have their base (pre-scaling) rate constants
 * updated from the host model based on the calculations of an external
 * module. The structure of the update data is:
 *
 *  - \b int rxn_id (Id of one or more first-order loss reactions set by the
 *       host model using the
 *       \c pmc_rxn_first_order_loss::rxn_first_order_loss_t::set_rxn_id
 *       function prior to initializing the solver.)
 *  - \b double rate_const (New pre-scaling rate constant.)
 *
 * \param update_data Pointer to the updated reaction data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent data
 * \return Flag indicating whether this is the reaction to update
 */
bool rxn_first_order_loss_update_data(void *update_data, int *rxn_int_data,
                                      double *rxn_float_data,
                                      double *rxn_env_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  int *rxn_id = (int *)update_data;
  double *base_rate = (double *)&(rxn_id[1]);

  // Set the base first-order loss rate constants for matching reactions
  if (*rxn_id == RXN_ID_ && RXN_ID_ > 0) {
    BASE_RATE_ = (double)*base_rate;
    RATE_CONSTANT_ = SCALING_ * BASE_RATE_;
    return true;
  }

  return false;
}

/** \brief Update reaction data for new environmental conditions
 *
 * For first-order loss reactions this only involves recalculating the rate
 * constant.
 *
 * \param model_data Pointer to the model data
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 */
void rxn_first_order_loss_update_env_state(ModelData *model_data,
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
 * \param time_deriv Pointer to the TimeDerivative object
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 * \param time_step Current time step being computed (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_first_order_loss_calc_deriv_contrib(
    ModelData *model_data, TimeDerivative *time_deriv, int *rxn_int_data,
    double *rxn_float_data, double *rxn_env_data, realtype time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Calculate the reaction rate
  long double rate = RATE_CONSTANT_ * state[REACT_];

  // Add contributions to the time derivative
  if (DERIV_ID_ >= 0) time_derivative_add_value(time_deriv, DERIV_ID_, -rate);

  return;
}
#endif

/** \brief Calculate contributions to the Jacobian from this reaction
 *
 * \param model_data Pointer to the model data
 * \param J Pointer to the sparse Jacobian matrix to add contributions to
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 * \param rxn_env_data Pointer to the environment-dependent parameters
 * \param time_step Current time step being calculated (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_first_order_loss_calc_jac_contrib(ModelData *model_data, realtype *J,
                                           int *rxn_int_data,
                                           double *rxn_float_data,
                                           double *rxn_env_data,
                                           realtype time_step) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;
  double *state = model_data->grid_cell_state;
  double *env_data = model_data->grid_cell_env;

  // Add contributions to the Jacobian
  if (JAC_ID_ >= 0) J[JAC_ID_] -= RATE_CONSTANT_;

  return;
}
#endif

/** \brief Print the reaction parameters
 *
 * \param rxn_int_data Pointer to the reaction integer data
 * \param rxn_float_data Pointer to the reaction floating-point data
 */
void rxn_first_order_loss_print(int *rxn_int_data, double *rxn_float_data) {
  int *int_data = rxn_int_data;
  double *float_data = rxn_float_data;

  printf("\n\nFirst-Order loss reaction\n");

  return;
}

/** \brief Create update data for new first-order loss rates
 *
 * \return Pointer to a new rate update data object
 */
void *rxn_first_order_loss_create_rate_update_data() {
  int *update_data = (int *)malloc(sizeof(int) + sizeof(double));
  if (update_data == NULL) {
    printf("\n\nERROR allocating space for first-order loss update data\n\n");
    exit(1);
  }
  return (void *)update_data;
}

/** \brief Set rate update data
 *
 * \param update_data Pointer to an allocated rate update data object
 * \param rxn_id Id of first-order loss reactions to update
 * \param base_rate New pre-scaling first-order loss rate
 */
void rxn_first_order_loss_set_rate_update_data(void *update_data, int rxn_id,
                                               double base_rate) {
  int *new_rxn_id = (int *)update_data;
  double *new_base_rate = (double *)&(new_rxn_id[1]);
  *new_rxn_id = rxn_id;
  *new_base_rate = base_rate;
}
