/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Reaction-specific functions for use by the solver
 *
 */
/** \file
 * \brief Reaction solver functions
 */
#define PMC_DEBUG_SPEC_ 118

#include "rxn_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include "rxns.h"

// Reaction types (Must match parameters defined in pmc_rxn_factory)
#define RXN_ARRHENIUS 1
#define RXN_TROE 2
#define RXN_CMAQ_H2O2 3
#define RXN_CMAQ_OH_HNO3 4
#define RXN_PHOTOLYSIS 5
#define RXN_HL_PHASE_TRANSFER 6
#define RXN_AQUEOUS_EQUILIBRIUM 7
#define RXN_SIMPOL_PHASE_TRANSFER 10
#define RXN_CONDENSED_PHASE_ARRHENIUS 11
#define RXN_FIRST_ORDER_LOSS 12
#define RXN_EMISSION 13
#define RXN_WET_DEPOSITION 14

/** \brief Get the Jacobian elements used by a particular reaction
 *
 * \param model_data A pointer to the model data
 * \param jac_struct A 2D array of flags indicating which Jacobian elements
 *                   may be used
 */
void rxn_get_used_jac_elem(ModelData *model_data, bool **jac_struct) {
  // Get the number of reactions
  int n_rxn = model_data->n_rxn;

  // Loop through the reactions to determine the Jacobian elements used
  // advancing the rxn_data pointer each time
  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {
    // Get pointers to the reaction data
    int *rxn_int_data =
        &(model_data->rxn_int_data[model_data->rxn_int_indices[i_rxn]]);
    double *rxn_float_data =
        &(model_data->rxn_float_data[model_data->rxn_float_indices[i_rxn]]);

    // Get the reaction type
    int rxn_type = *(rxn_int_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM:
        rxn_aqueous_equilibrium_get_used_jac_elem(rxn_int_data, rxn_float_data,
                                                  jac_struct);
        break;
      case RXN_ARRHENIUS:
        rxn_arrhenius_get_used_jac_elem(rxn_int_data, rxn_float_data,
                                        jac_struct);
        break;
      case RXN_CMAQ_H2O2:
        rxn_CMAQ_H2O2_get_used_jac_elem(rxn_int_data, rxn_float_data,
                                        jac_struct);
        break;
      case RXN_CMAQ_OH_HNO3:
        rxn_CMAQ_OH_HNO3_get_used_jac_elem(rxn_int_data, rxn_float_data,
                                           jac_struct);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS:
        rxn_condensed_phase_arrhenius_get_used_jac_elem(
            rxn_int_data, rxn_float_data, jac_struct);
        break;
      case RXN_EMISSION:
        rxn_emission_get_used_jac_elem(rxn_int_data, rxn_float_data,
                                       jac_struct);
        break;
      case RXN_FIRST_ORDER_LOSS:
        rxn_first_order_loss_get_used_jac_elem(rxn_int_data, rxn_float_data,
                                               jac_struct);
        break;
      case RXN_HL_PHASE_TRANSFER:
        rxn_HL_phase_transfer_get_used_jac_elem(model_data, rxn_int_data,
                                                rxn_float_data, jac_struct);
        break;
      case RXN_PHOTOLYSIS:
        rxn_photolysis_get_used_jac_elem(rxn_int_data, rxn_float_data,
                                         jac_struct);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER:
        rxn_SIMPOL_phase_transfer_get_used_jac_elem(model_data, rxn_int_data,
                                                    rxn_float_data, jac_struct);
        break;
      case RXN_TROE:
        rxn_troe_get_used_jac_elem(rxn_int_data, rxn_float_data, jac_struct);
        break;
      case RXN_WET_DEPOSITION:
        rxn_wet_deposition_get_used_jac_elem(rxn_int_data, rxn_float_data,
                                             jac_struct);
        break;
    }
  }
}

/** \brief Update the time derivative and Jacobian array ids
 *
 * \param model_data Pointer to the model data
 * \param deriv_ids Ids for state variables on the time derivative array
 * \param jac_ids Ids for state variables on the Jacobian array
 */
void rxn_update_ids(ModelData *model_data, int *deriv_ids, int **jac_ids) {
  // Get the number of reactions
  int n_rxn = model_data->n_rxn;

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {
    // Get pointers to the reaction data
    int *rxn_int_data =
        &(model_data->rxn_int_data[model_data->rxn_int_indices[i_rxn]]);
    double *rxn_float_data =
        &(model_data->rxn_float_data[model_data->rxn_float_indices[i_rxn]]);

    // Get the reaction type
    int rxn_type = *(rxn_int_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM:
        rxn_aqueous_equilibrium_update_ids(model_data, deriv_ids, jac_ids,
                                           rxn_int_data, rxn_float_data);
        break;
      case RXN_ARRHENIUS:
        rxn_arrhenius_update_ids(model_data, deriv_ids, jac_ids, rxn_int_data,
                                 rxn_float_data);
        break;
      case RXN_CMAQ_H2O2:
        rxn_CMAQ_H2O2_update_ids(model_data, deriv_ids, jac_ids, rxn_int_data,
                                 rxn_float_data);
        break;
      case RXN_CMAQ_OH_HNO3:
        rxn_CMAQ_OH_HNO3_update_ids(model_data, deriv_ids, jac_ids,
                                    rxn_int_data, rxn_float_data);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS:
        rxn_condensed_phase_arrhenius_update_ids(model_data, deriv_ids, jac_ids,
                                                 rxn_int_data, rxn_float_data);
        break;
      case RXN_EMISSION:
        rxn_emission_update_ids(model_data, deriv_ids, jac_ids, rxn_int_data,
                                rxn_float_data);
        break;
      case RXN_FIRST_ORDER_LOSS:
        rxn_first_order_loss_update_ids(model_data, deriv_ids, jac_ids,
                                        rxn_int_data, rxn_float_data);
        break;
      case RXN_HL_PHASE_TRANSFER:
        rxn_HL_phase_transfer_update_ids(model_data, deriv_ids, jac_ids,
                                         rxn_int_data, rxn_float_data);
        break;
      case RXN_PHOTOLYSIS:
        rxn_photolysis_update_ids(model_data, deriv_ids, jac_ids, rxn_int_data,
                                  rxn_float_data);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER:
        rxn_SIMPOL_phase_transfer_update_ids(model_data, deriv_ids, jac_ids,
                                             rxn_int_data, rxn_float_data);
        break;
      case RXN_TROE:
        rxn_troe_update_ids(model_data, deriv_ids, jac_ids, rxn_int_data,
                            rxn_float_data);
        break;
      case RXN_WET_DEPOSITION:
        rxn_wet_deposition_update_ids(model_data, deriv_ids, jac_ids,
                                      rxn_int_data, rxn_float_data);
        break;
    }
  }
}

/** \brief Update reaction data for new environmental state
 *
 * \param model_data Pointer to the model data with updated env state
 */
void rxn_update_env_state(ModelData *model_data) {
  // Get the number of reactions
  int n_rxn = model_data->n_rxn;

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {
    // Get pointers to the reaction data
    int *rxn_int_data =
        &(model_data->rxn_int_data[model_data->rxn_int_indices[i_rxn]]);
    double *rxn_float_data =
        &(model_data->rxn_float_data[model_data->rxn_float_indices[i_rxn]]);
    double *rxn_env_data =
        &(model_data->grid_cell_rxn_env_data[model_data->rxn_env_idx[i_rxn]]);

    // Get the reaction type
    int rxn_type = *(rxn_int_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM:
        rxn_aqueous_equilibrium_update_env_state(model_data, rxn_int_data,
                                                 rxn_float_data, rxn_env_data);
        break;
      case RXN_ARRHENIUS:
        rxn_arrhenius_update_env_state(model_data, rxn_int_data, rxn_float_data,
                                       rxn_env_data);
        break;
      case RXN_CMAQ_H2O2:
        rxn_CMAQ_H2O2_update_env_state(model_data, rxn_int_data, rxn_float_data,
                                       rxn_env_data);
        break;
      case RXN_CMAQ_OH_HNO3:
        rxn_CMAQ_OH_HNO3_update_env_state(model_data, rxn_int_data,
                                          rxn_float_data, rxn_env_data);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS:
        rxn_condensed_phase_arrhenius_update_env_state(
            model_data, rxn_int_data, rxn_float_data, rxn_env_data);
        break;
      case RXN_EMISSION:
        rxn_emission_update_env_state(model_data, rxn_int_data, rxn_float_data,
                                      rxn_env_data);
        break;
      case RXN_FIRST_ORDER_LOSS:
        rxn_first_order_loss_update_env_state(model_data, rxn_int_data,
                                              rxn_float_data, rxn_env_data);
        break;
      case RXN_HL_PHASE_TRANSFER:
        rxn_HL_phase_transfer_update_env_state(model_data, rxn_int_data,
                                               rxn_float_data, rxn_env_data);
        break;
      case RXN_PHOTOLYSIS:
        rxn_photolysis_update_env_state(model_data, rxn_int_data,
                                        rxn_float_data, rxn_env_data);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER:
        rxn_SIMPOL_phase_transfer_update_env_state(
            model_data, rxn_int_data, rxn_float_data, rxn_env_data);
        break;
      case RXN_TROE:
        rxn_troe_update_env_state(model_data, rxn_int_data, rxn_float_data,
                                  rxn_env_data);
        break;
      case RXN_WET_DEPOSITION:
        rxn_wet_deposition_update_env_state(model_data, rxn_int_data,
                                            rxn_float_data, rxn_env_data);
        break;
    }
  }
}

/** \brief Calculate the time derivative \f$f(t,y)\f$
 *
 * \param model_data Pointer to the model data
 * \param deriv_data Pointer to the derivative data for the current grid cell
 * \param time_step Current model time step (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_calc_deriv(ModelData *model_data, double *deriv_data,
                    realtype time_step) {
  // Get the number of reactions
  int n_rxn = model_data->n_rxn;

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {
    // Get pointers to the reaction data
    int *rxn_int_data =
        &(model_data->rxn_int_data[model_data->rxn_int_indices[i_rxn]]);
    double *rxn_float_data =
        &(model_data->rxn_float_data[model_data->rxn_float_indices[i_rxn]]);
    double *rxn_env_data =
        &(model_data->grid_cell_rxn_env_data[model_data->rxn_env_idx[i_rxn]]);

    // Get the reaction type
    int rxn_type = *(rxn_int_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM:
        rxn_aqueous_equilibrium_calc_deriv_contrib(model_data, deriv_data,
                                                   rxn_int_data, rxn_float_data,
                                                   rxn_env_data, time_step);
        break;
      case RXN_ARRHENIUS:
        rxn_arrhenius_calc_deriv_contrib(model_data, deriv_data, rxn_int_data,
                                         rxn_float_data, rxn_env_data,
                                         time_step);
        break;
      case RXN_CMAQ_H2O2:
        rxn_CMAQ_H2O2_calc_deriv_contrib(model_data, deriv_data, rxn_int_data,
                                         rxn_float_data, rxn_env_data,
                                         time_step);
        break;
      case RXN_CMAQ_OH_HNO3:
        rxn_CMAQ_OH_HNO3_calc_deriv_contrib(model_data, deriv_data,
                                            rxn_int_data, rxn_float_data,
                                            rxn_env_data, time_step);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS:
        rxn_condensed_phase_arrhenius_calc_deriv_contrib(
            model_data, deriv_data, rxn_int_data, rxn_float_data, rxn_env_data,
            time_step);
        break;
      case RXN_EMISSION:
        rxn_emission_calc_deriv_contrib(model_data, deriv_data, rxn_int_data,
                                        rxn_float_data, rxn_env_data,
                                        time_step);
        break;
      case RXN_FIRST_ORDER_LOSS:
        rxn_first_order_loss_calc_deriv_contrib(model_data, deriv_data,
                                                rxn_int_data, rxn_float_data,
                                                rxn_env_data, time_step);
        break;
      case RXN_HL_PHASE_TRANSFER:
        rxn_HL_phase_transfer_calc_deriv_contrib(model_data, deriv_data,
                                                 rxn_int_data, rxn_float_data,
                                                 rxn_env_data, time_step);
        break;
      case RXN_PHOTOLYSIS:
        rxn_photolysis_calc_deriv_contrib(model_data, deriv_data, rxn_int_data,
                                          rxn_float_data, rxn_env_data,
                                          time_step);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER:
        rxn_SIMPOL_phase_transfer_calc_deriv_contrib(
            model_data, deriv_data, rxn_int_data, rxn_float_data, rxn_env_data,
            time_step);
        break;
      case RXN_TROE:
        rxn_troe_calc_deriv_contrib(model_data, deriv_data, rxn_int_data,
                                    rxn_float_data, rxn_env_data, time_step);
        break;
      case RXN_WET_DEPOSITION:
        rxn_wet_deposition_calc_deriv_contrib(model_data, deriv_data,
                                              rxn_int_data, rxn_float_data,
                                              rxn_env_data, time_step);
        break;
    }
  }
}
#endif

/** \brief Calculate the Jacobian
 *
 * \param model_data Pointer to the model data
 * \param J_data Pointer to the Jacobian to be calculated for the current
 *               grid cell
 * \param time_step Current model time step (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_calc_jac(ModelData *model_data, double *J_data, realtype time_step) {
  // Get the number of reactions
  int n_rxn = model_data->n_rxn;

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {
    // Get pointers to the reaction data
    int *rxn_int_data =
        &(model_data->rxn_int_data[model_data->rxn_int_indices[i_rxn]]);
    double *rxn_float_data =
        &(model_data->rxn_float_data[model_data->rxn_float_indices[i_rxn]]);
    double *rxn_env_data =
        &(model_data->grid_cell_rxn_env_data[model_data->rxn_env_idx[i_rxn]]);

    // Get the reaction type
    int rxn_type = *(rxn_int_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM:
        rxn_aqueous_equilibrium_calc_jac_contrib(model_data, J_data,
                                                 rxn_int_data, rxn_float_data,
                                                 rxn_env_data, time_step);
        break;
      case RXN_ARRHENIUS:
        rxn_arrhenius_calc_jac_contrib(model_data, J_data, rxn_int_data,
                                       rxn_float_data, rxn_env_data, time_step);
        break;
      case RXN_CMAQ_H2O2:
        rxn_CMAQ_H2O2_calc_jac_contrib(model_data, J_data, rxn_int_data,
                                       rxn_float_data, rxn_env_data, time_step);
        break;
      case RXN_CMAQ_OH_HNO3:
        rxn_CMAQ_OH_HNO3_calc_jac_contrib(model_data, J_data, rxn_int_data,
                                          rxn_float_data, rxn_env_data,
                                          time_step);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS:
        rxn_condensed_phase_arrhenius_calc_jac_contrib(
            model_data, J_data, rxn_int_data, rxn_float_data, rxn_env_data,
            time_step);
        break;
      case RXN_EMISSION:
        rxn_emission_calc_jac_contrib(model_data, J_data, rxn_int_data,
                                      rxn_float_data, rxn_env_data, time_step);
        break;
      case RXN_FIRST_ORDER_LOSS:
        rxn_first_order_loss_calc_jac_contrib(model_data, J_data, rxn_int_data,
                                              rxn_float_data, rxn_env_data,
                                              time_step);
        break;
      case RXN_HL_PHASE_TRANSFER:
        rxn_HL_phase_transfer_calc_jac_contrib(model_data, J_data, rxn_int_data,
                                               rxn_float_data, rxn_env_data,
                                               time_step);
        break;
      case RXN_PHOTOLYSIS:
        rxn_photolysis_calc_jac_contrib(model_data, J_data, rxn_int_data,
                                        rxn_float_data, rxn_env_data,
                                        time_step);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER:
        rxn_SIMPOL_phase_transfer_calc_jac_contrib(model_data, J_data,
                                                   rxn_int_data, rxn_float_data,
                                                   rxn_env_data, time_step);
        break;
      case RXN_TROE:
        rxn_troe_calc_jac_contrib(model_data, J_data, rxn_int_data,
                                  rxn_float_data, rxn_env_data, time_step);
        break;
      case RXN_WET_DEPOSITION:
        rxn_wet_deposition_calc_jac_contrib(model_data, J_data, rxn_int_data,
                                            rxn_float_data, rxn_env_data,
                                            time_step);
        break;
    }
  }
}
#endif

/** \brief Add condensed data to the condensed data block of memory
 *
 * \param rxn_type Reaction type
 * \param n_int_param Number of integer parameters
 * \param n_float_param Number of floating-point parameters
 * \param n_env_param Number of environment-dependent parameters
 * \param int_param Pointer to integer parameter array
 * \param float_param Pointer to floating-point parameter array
 * \param solver_data Pointer to solver data
 */
// TODO: question: move n_added_rxns out of struct to function parameter since
// is only used in this function
void rxn_add_condensed_data(int rxn_type, int n_int_param, int n_float_param,
                            int n_env_param, int *int_param,
                            double *float_param, void *solver_data) {
  ModelData *model_data =
      (ModelData *)&(((SolverData *)solver_data)->model_data);

  // Get pointers to the reaction data
  int *rxn_int_data =
      &(model_data->rxn_int_data
            [model_data->rxn_int_indices[model_data->n_added_rxns]]);
  double *rxn_float_data =
      &(model_data->rxn_float_data
            [model_data->rxn_float_indices[model_data->n_added_rxns]]);

  // Save next indices by adding lengths
  model_data->rxn_int_indices[model_data->n_added_rxns + 1] =
      (n_int_param + 1) +
      model_data->rxn_int_indices[model_data->n_added_rxns];  //+1 is type
  model_data->rxn_float_indices[model_data->n_added_rxns + 1] =
      n_float_param + model_data->rxn_float_indices[model_data->n_added_rxns];
  model_data->rxn_env_idx[model_data->n_added_rxns + 1] =
      model_data->rxn_env_idx[model_data->n_added_rxns] + n_env_param;
  ++(model_data->n_added_rxns);

  // Add the reaction type
  *(rxn_int_data++) = rxn_type;

  // Add integer parameters
  for (; n_int_param > 0; --n_int_param) *(rxn_int_data++) = *(int_param++);

  // Add floating-point parameters
  for (; n_float_param > 0; --n_float_param)
    *(rxn_float_data++) = (double)*(float_param++);

  model_data->n_rxn_env_data += n_env_param;
}

/** \brief Update reaction data
 *
 * Update data for one or more reactions. Reactions of a certain type are
 * passed a void pointer to updated data that must be in the format specified
 * by the reaction type. This data could be used to find specific reactions of
 * the specified type and, for example, update rate constants.
 *
 * \param cell_id Id of the grid cell to update
 * \param rxn_id Id of the reaction (or 0 if unknown)
 * \param update_rxn_type Type of the reaction
 * \param update_data Pointer to updated data to pass to the reaction
 * \param solver_data Pointer to solver data
 */
void rxn_update_data(int cell_id, int *rxn_id, int update_rxn_type,
                     void *update_data, void *solver_data) {
  ModelData *model_data =
      (ModelData *)&(((SolverData *)solver_data)->model_data);

  // Point to the environment-dependent data for the grid cell
  model_data->grid_cell_rxn_env_data =
      &(model_data->rxn_env_data[cell_id * model_data->n_rxn_env_data]);

  // Get the number of reactions
  int n_rxn = model_data->n_rxn;

  // Loop through the reactions advancing the rxn_data pointer each time
  for (; (*rxn_id) < n_rxn; (*rxn_id)++) {
    // Get pointers to the reaction data
    int *rxn_int_data =
        &(model_data->rxn_int_data[model_data->rxn_int_indices[*rxn_id]]);
    double *rxn_float_data =
        &(model_data->rxn_float_data[model_data->rxn_float_indices[*rxn_id]]);

    double *rxn_env_data =
        &(model_data->grid_cell_rxn_env_data[model_data->rxn_env_idx[*rxn_id]]);

    // Get the reaction type
    int rxn_type = *(rxn_int_data++);

    bool found = false;

    // Try the update data function for reactions of the correct type
    if (rxn_type == update_rxn_type) {
      switch (rxn_type) {
        case RXN_EMISSION:
          found = rxn_emission_update_data((void *)update_data, rxn_int_data,
                                           rxn_float_data, rxn_env_data);
          break;
        case RXN_FIRST_ORDER_LOSS:
          found = rxn_first_order_loss_update_data(
              (void *)update_data, rxn_int_data, rxn_float_data, rxn_env_data);
          break;
        case RXN_PHOTOLYSIS:
          found = rxn_photolysis_update_data((void *)update_data, rxn_int_data,
                                             rxn_float_data, rxn_env_data);
          break;
        case RXN_WET_DEPOSITION:
          found = rxn_wet_deposition_update_data(
              (void *)update_data, rxn_int_data, rxn_float_data, rxn_env_data);
          break;
      }
      if (found) return;
    }
  }
}

/** \brief Print the reaction data
 *
 * \param solver_data Pointer to the solver data
 */
void rxn_print_data(void *solver_data) {
  ModelData *model_data =
      (ModelData *)&(((SolverData *)solver_data)->model_data);

  // Get the number of reactions
  int n_rxn = model_data->n_rxn;

  printf("\n\nReaction data\n\nnumber of reactions: %d\n\n", n_rxn);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++) {
    // Get pointers to the reaction data
    int *rxn_int_data =
        &(model_data->rxn_int_data[model_data->rxn_int_indices[i_rxn]]);
    double *rxn_float_data =
        &(model_data->rxn_float_data[model_data->rxn_float_indices[i_rxn]]);

    // Get the reaction type
    int rxn_type = *(rxn_int_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM:
        rxn_aqueous_equilibrium_print(rxn_int_data, rxn_float_data);
        break;
      case RXN_ARRHENIUS:
        rxn_arrhenius_print(rxn_int_data, rxn_float_data);
        break;
      case RXN_CMAQ_H2O2:
        rxn_CMAQ_H2O2_print(rxn_int_data, rxn_float_data);
        break;
      case RXN_CMAQ_OH_HNO3:
        rxn_CMAQ_OH_HNO3_print(rxn_int_data, rxn_float_data);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS:
        rxn_condensed_phase_arrhenius_print(rxn_int_data, rxn_float_data);
        break;
      case RXN_EMISSION:
        rxn_emission_print(rxn_int_data, rxn_float_data);
        break;
      case RXN_FIRST_ORDER_LOSS:
        rxn_first_order_loss_print(rxn_int_data, rxn_float_data);
        break;
      case RXN_HL_PHASE_TRANSFER:
        rxn_HL_phase_transfer_print(rxn_int_data, rxn_float_data);
        break;
      case RXN_PHOTOLYSIS:
        rxn_photolysis_print(rxn_int_data, rxn_float_data);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER:
        rxn_SIMPOL_phase_transfer_print(rxn_int_data, rxn_float_data);
        break;
      case RXN_TROE:
        rxn_troe_print(rxn_int_data, rxn_float_data);
        break;
      case RXN_WET_DEPOSITION:
        rxn_wet_deposition_print(rxn_int_data, rxn_float_data);
        break;
    }
  }
  fflush(stdout);
}

/** \brief Free an update data object
 *
 * \param update_data Object to free
 */
void rxn_free_update_data(void *update_data) { free(update_data); }
