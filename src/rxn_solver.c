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

#include <stdio.h>
#include <stdlib.h>
#include "rxn_solver.h"
#include "rxns.h"

// Reaction types (Must match parameters defined in pmc_rxn_factory)
#define RXN_ARRHENIUS 1
#define RXN_TROE 2
#define RXN_CMAQ_H2O2 3
#define RXN_CMAQ_OH_HNO3 4
#define RXN_PHOTOLYSIS 5
#define RXN_HL_PHASE_TRANSFER 6
#define RXN_AQUEOUS_EQUILIBRIUM 7
#define RXN_ZSR_AEROSOL_WATER 8
#define RXN_PDFITE_ACTIVITY 9
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
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_get_used_jac_elem(ModelData *model_data, bool **jac_struct)
{

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);

  // Loop through the reactions to determine the Jacobian elements used
  // advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_data = (int*) rxn_aqueous_equilibrium_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_data = (int*) rxn_condensed_phase_arrhenius_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_EMISSION :
        rxn_data = (int*) rxn_emission_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_data = (int*) rxn_first_order_loss_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_HL_phase_transfer_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_PDFITE_ACTIVITY:
        rxn_data = (int*) rxn_PDFiTE_activity_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_SIMPOL_phase_transfer_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_WET_DEPOSITION :
        rxn_data = (int*) rxn_wet_deposition_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_data = (int*) rxn_ZSR_aerosol_water_get_used_jac_elem(
                  (void*) rxn_data, jac_struct);
        break;
    }
  }
  return rxn_data;
}

/** \brief Update the time derivative and Jacobian array ids
 *
 * \param model_data Pointer to the model data
 * \param deriv_ids Ids for state variables on the time derivative array
 * \param jac_ids Ids for state variables on the Jacobian array
 */
void rxn_update_ids(ModelData *model_data, int *deriv_ids, int **jac_ids)
{

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_data = (int*) rxn_aqueous_equilibrium_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_data = (int*) rxn_condensed_phase_arrhenius_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_EMISSION :
        rxn_data = (int*) rxn_emission_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_data = (int*) rxn_first_order_loss_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_HL_phase_transfer_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_data = (int*) rxn_PDFiTE_activity_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_SIMPOL_phase_transfer_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_WET_DEPOSITION :
        rxn_data = (int*) rxn_wet_deposition_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_data = (int*) rxn_ZSR_aerosol_water_update_ids(
                  model_data, deriv_ids, jac_ids, (void*) rxn_data);
        break;
    }
  }
}

/** \brief Update reaction data for new environmental state
 *
 * \param model_data Pointer to the model data
 * \param env Pointer to the environmental state array
 */
void rxn_update_env_state(ModelData *model_data, double *env)
{

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_data = (int*) rxn_aqueous_equilibrium_update_env_state(
                  env, (void*) rxn_data);
        break;
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_update_env_state(
                  env, (void*) rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_update_env_state(
                  env, (void*) rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_update_env_state(
                  env, (void*) rxn_data);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_data = (int*) rxn_condensed_phase_arrhenius_update_env_state(
                  env, (void*) rxn_data);
        break;
      case RXN_EMISSION :
        rxn_data = (int*) rxn_emission_update_env_state(
                  env, (void*) rxn_data);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_data = (int*) rxn_first_order_loss_update_env_state(
                  env, (void*) rxn_data);
        break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_HL_phase_transfer_update_env_state(
                  env, (void*) rxn_data);
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_data = (int*) rxn_PDFiTE_activity_update_env_state(
                  env, (void*) rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_update_env_state(
                  env, (void*) rxn_data);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_SIMPOL_phase_transfer_update_env_state(
                  env, (void*) rxn_data);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_update_env_state(
                  env, (void*) rxn_data);
        break;
      case RXN_WET_DEPOSITION :
        rxn_data = (int*) rxn_wet_deposition_update_env_state(
                  env, (void*) rxn_data);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_data = (int*) rxn_ZSR_aerosol_water_update_env_state(
                  env, (void*) rxn_data);
        break;
    }
  }
}

/** \brief Do pre-derivative/Jacobian calculations
 *
 * Update the model state based on sub-model calculations. These can,
 * for example, update activity coefficients based on species concentrations,
 * calculate PSSA species concentrations, partition aerosol water, etc.
 *
 * Fast reactions can add contributions to state adjustments instead of to
 * the derivatives and Jacobian. These reactions should be fast enough that
 * they essentially go to completion over the model timestep. Absolute changes
 * in species concentrations should be added to the state_adj array and rates
 * of change (conc. units/s) to the rel_adj_cont. If adjustments to the state
 * result in negative species concentrations (depletion of a single species by
 * mutliple reactions), reactions will be asked to scale their net flux based
 * on the values in rel_adj_conc.
 *
 * TODO Unlike the derivative calculations, the order of these operations may
 * matter. Find a way to abstract operation order. Maybe during intialization
 * each reaction can specify it's dependencies and which state variables it will
 * change and the order can be determined from this at run time.
 *
 * \param model_data Pointer to the model data
 */
void rxn_pre_calc(ModelData *model_data, double time_step)
{
  // If called with time_step = 0, update the state based on the last
  // calculated adjustments
  if (time_step==0.0 && model_data->use_adj) rxn_adjust_state(model_data);

  // Start with no state adjustments
  rxn_reset_state_adjustments(model_data);

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_data = (int*) rxn_aqueous_equilibrium_pre_calc(
                  model_data, (void*) rxn_data);
        break;
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_pre_calc(
                  model_data, (void*) rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_pre_calc(
                  model_data, (void*) rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_pre_calc(
                  model_data, (void*) rxn_data);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_data = (int*) rxn_condensed_phase_arrhenius_pre_calc(
                  model_data, (void*) rxn_data);
        break;
      case RXN_EMISSION :
        rxn_data = (int*) rxn_emission_pre_calc(
                  model_data, (void*) rxn_data);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_data = (int*) rxn_first_order_loss_pre_calc(
                  model_data, (void*) rxn_data);
        break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_HL_phase_transfer_pre_calc(
                  model_data, (void*) rxn_data, time_step);
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_data = (int*) rxn_PDFiTE_activity_pre_calc(
                  model_data, (void*) rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_pre_calc(
                  model_data, (void*) rxn_data);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_SIMPOL_phase_transfer_pre_calc(
                  model_data, (void*) rxn_data, time_step);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_pre_calc(
                  model_data, (void*) rxn_data);
        break;
      case RXN_WET_DEPOSITION :
        rxn_data = (int*) rxn_wet_deposition_pre_calc(
                  model_data, (void*) rxn_data);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_data = (int*) rxn_ZSR_aerosol_water_pre_calc(
                  model_data, (void*) rxn_data);
        break;
    }
  }

  // Update the state array (if called during solving)
  if (time_step > 0.0 && model_data->use_adj) rxn_adjust_state(model_data);

}

/** \brief Reset the state adjustments
 *
 * \param model_data Pointer to the model data
 */
void rxn_reset_state_adjustments(ModelData *model_data)
{
  if (model_data->use_adj) {
    model_data->use_adj = false;
    model_data->scale_adj = false;
    for (int i_var = 0; i_var < model_data->n_state_var; i_var++)
      model_data->state_adj[i_var] = model_data->rel_flux[i_var] = 0.0;
  }
}

/** \brief Adjust the state for fast reactions
 *
 * \param model_data Pointer to the model data
 */
void rxn_adjust_state(ModelData *model_data)
{
  // First check whether the adjustments need scaled
  if (model_data->scale_adj) {

    // Get the number of reactions
    int *rxn_data = (int*) (model_data->rxn_data);
    int n_rxn = *(rxn_data++);

    // Loop through the reactions advancing the rxn_data pointer each time
    for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

      // Get the reaction type
      int rxn_type = *(rxn_data++);

      // Scale adjustments based on relative contributions from each rxn
      switch (rxn_type) {
        case RXN_AQUEOUS_EQUILIBRIUM :
          rxn_data = (int*) rxn_aqueous_equilibrium_skip(
                     (void*) rxn_data);
          break;
        case RXN_ARRHENIUS :
          rxn_data = (int*) rxn_arrhenius_skip(
                     (void*) rxn_data);
          break;
        case RXN_CMAQ_H2O2 :
          rxn_data = (int*) rxn_CMAQ_H2O2_skip(
                     (void*) rxn_data);
          break;
        case RXN_CMAQ_OH_HNO3 :
          rxn_data = (int*) rxn_CMAQ_OH_HNO3_skip(
                     (void*) rxn_data);
          break;
        case RXN_CONDENSED_PHASE_ARRHENIUS :
          rxn_data = (int*) rxn_condensed_phase_arrhenius_skip(
                     (void*) rxn_data);
          break;
        case RXN_EMISSION :
          rxn_data = (int*) rxn_emission_skip(
                     (void*) rxn_data);
          break;
        case RXN_FIRST_ORDER_LOSS :
          rxn_data = (int*) rxn_first_order_loss_skip(
                     (void*) rxn_data);
          break;
        case RXN_HL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_HL_phase_transfer_scale_adj(
                     model_data, (void*) rxn_data);
          break;
        case RXN_PDFITE_ACTIVITY :
          rxn_data = (int*) rxn_PDFiTE_activity_skip(
                     (void*) rxn_data);
          break;
        case RXN_PHOTOLYSIS :
          rxn_data = (int*) rxn_photolysis_skip(
                     (void*) rxn_data);
          break;
        case RXN_SIMPOL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_SIMPOL_phase_transfer_scale_adj(
                     model_data, (void*) rxn_data);
          break;
        case RXN_TROE :
          rxn_data = (int*) rxn_troe_skip(
                     (void*) rxn_data);
          break;
        case RXN_WET_DEPOSITION :
          rxn_data = (int*) rxn_wet_deposition_skip(
                     (void*) rxn_data);
          break;
        case RXN_ZSR_AEROSOL_WATER :
          rxn_data = (int*) rxn_ZSR_aerosol_water_skip(
                     (void*) rxn_data);
          break;
      }
    }
  }
  model_data->scale_adj = false;

#ifdef PMC_DEBUG
    printf("\nAdjusting state for species %d from %le by %le", PMC_DEBUG_SPEC_,
        model_data->state[PMC_DEBUG_SPEC_],
        model_data->state_adj[PMC_DEBUG_SPEC_]);
#endif

  // Update the state array based on the calculated adjustments
  for (int i_var = 0; i_var < model_data->n_state_var; i_var++)
    model_data->state[i_var] += model_data->state_adj[i_var];

}

/** \brief Calculate the time derivative \f$f(t,y)\f$
 *
 * \param model_data Pointer to the model data
 * \param deriv NVector to hold the calculated vector
 * \param time_step Current model time step (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_calc_deriv(ModelData *model_data, N_Vector deriv, realtype time_step)
{

  // Get a pointer to the derivative data
  realtype *deriv_data = N_VGetArrayPointer(deriv);

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_data = (int*) rxn_aqueous_equilibrium_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_data = (int*) rxn_condensed_phase_arrhenius_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_EMISSION :
        rxn_data = (int*) rxn_emission_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_data = (int*) rxn_first_order_loss_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_HL_phase_transfer_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_data = (int*) rxn_PDFiTE_activity_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_SIMPOL_phase_transfer_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_WET_DEPOSITION :
        rxn_data = (int*) rxn_wet_deposition_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_data = (int*) rxn_ZSR_aerosol_water_calc_deriv_contrib(
                  model_data, deriv_data, (void*) rxn_data, time_step);
        break;
    }
  }
}
#endif

/** \brief Calculate the Jacobian
 *
 * \param model_data Pointer to the model data
 * \param J Jacobian to be calculated
 * \param time_step Current model time step (s)
 */
#ifdef PMC_USE_SUNDIALS
void rxn_calc_jac(ModelData *model_data, SUNMatrix J, realtype time_step)
{

  // Get a pointer to the Jacobian data
  realtype *J_data = SM_DATA_S(J);

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_data = (int*) rxn_aqueous_equilibrium_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
        break;
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
        break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_data = (int*) rxn_condensed_phase_arrhenius_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
        break;
      case RXN_EMISSION :
        rxn_data = (int*) rxn_emission_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
        break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_data = (int*) rxn_first_order_loss_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
        break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_HL_phase_transfer_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
        break;
      case RXN_PDFITE_ACTIVITY :
        rxn_data = (int*) rxn_PDFiTE_activity_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
        break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_SIMPOL_phase_transfer_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
        break;
      case RXN_WET_DEPOSITION :
        rxn_data = (int*) rxn_wet_deposition_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
        break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_data = (int*) rxn_ZSR_aerosol_water_calc_jac_contrib(
                  model_data, J_data, (void*) rxn_data, time_step);
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
 * \param int_param Pointer to integer parameter array
 * \param float_param Pointer to floating-point parameter array
 * \param solver_data Pointer to solver data
 */
void rxn_add_condensed_data(int rxn_type, int n_int_param, int n_float_param,
          int *int_param, double *float_param, void *solver_data)
{
  ModelData *model_data =
          (ModelData*) &(((SolverData*)solver_data)->model_data);
  int *rxn_data = (int*) (model_data->nxt_rxn);

  // Add the reaction type
  *(rxn_data++) = rxn_type;

  // Add integer parameters
  for (; n_int_param>0; n_int_param--) *(rxn_data++) = *(int_param++);

  // Add floating-point parameters
  double *flt_ptr = (double*) rxn_data;
  for (; n_float_param>0; n_float_param--)
          *(flt_ptr++) = (double) *(float_param++);

  // Set the pointer for the next free space in rxn_data
  model_data->nxt_rxn = (void*) flt_ptr;
}

/** \brief Update reaction data
 *
 * Update data for one or more reactions. Reactions of a certain type are
 * passed a void pointer to updated data that must be in the format specified
 * by the reaction type. This data could be used to find specific reactions of
 * the specified type and, for example, update rate constants.
 *
 * \param update_rxn_type Type of the reaction
 * \param update_data Pointer to updated data to pass to the reaction
 * \param solver_data Pointer to solver data
 */
void rxn_update_data(int update_rxn_type, void *update_data, void *solver_data)
{
  ModelData *model_data =
          (ModelData*) &(((SolverData*)solver_data)->model_data);

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Skip reactions of other types
    if (rxn_type!=update_rxn_type) {
      switch (rxn_type) {
        case RXN_AQUEOUS_EQUILIBRIUM :
          rxn_data = (int*) rxn_aqueous_equilibrium_skip(
                    (void*) rxn_data);
          break;
        case RXN_ARRHENIUS :
          rxn_data = (int*) rxn_arrhenius_skip(
                    (void*) rxn_data);
          break;
        case RXN_CMAQ_H2O2 :
          rxn_data = (int*) rxn_CMAQ_H2O2_skip(
                    (void*) rxn_data);
          break;
        case RXN_CMAQ_OH_HNO3 :
          rxn_data = (int*) rxn_CMAQ_OH_HNO3_skip(
                    (void*) rxn_data);
          break;
        case RXN_CONDENSED_PHASE_ARRHENIUS :
          rxn_data = (int*) rxn_condensed_phase_arrhenius_skip(
                    (void*) rxn_data);
          break;
        case RXN_EMISSION :
          rxn_data = (int*) rxn_emission_skip(
                    (void*) rxn_data);
          break;
        case RXN_FIRST_ORDER_LOSS :
          rxn_data = (int*) rxn_first_order_loss_skip(
                    (void*) rxn_data);
          break;
        case RXN_HL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_HL_phase_transfer_skip(
                    (void*) rxn_data);
          break;
        case RXN_PDFITE_ACTIVITY :
          rxn_data = (int*) rxn_PDFiTE_activity_skip(
                    (void*) rxn_data);
          break;
        case RXN_PHOTOLYSIS :
          rxn_data = (int*) rxn_photolysis_skip(
                    (void*) rxn_data);
          break;
        case RXN_SIMPOL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_SIMPOL_phase_transfer_skip(
                    (void*) rxn_data);
          break;
        case RXN_TROE :
          rxn_data = (int*) rxn_troe_skip(
                    (void*) rxn_data);
          break;
        case RXN_WET_DEPOSITION :
          rxn_data = (int*) rxn_wet_deposition_skip(
                    (void*) rxn_data);
          break;
        case RXN_ZSR_AEROSOL_WATER :
          rxn_data = (int*) rxn_ZSR_aerosol_water_skip(
                    (void*) rxn_data);
          break;
      }

    // ... otherwise, call the update function for reaction types that have them
    } else {
      switch (rxn_type) {
        case RXN_AQUEOUS_EQUILIBRIUM :
          rxn_data = (int*) rxn_aqueous_equilibrium_skip(
                    (void*) rxn_data);
          break;
        case RXN_ARRHENIUS :
          rxn_data = (int*) rxn_arrhenius_skip(
                    (void*) rxn_data);
          break;
        case RXN_CMAQ_H2O2 :
          rxn_data = (int*) rxn_CMAQ_H2O2_skip(
                    (void*) rxn_data);
          break;
        case RXN_CMAQ_OH_HNO3 :
          rxn_data = (int*) rxn_CMAQ_OH_HNO3_skip(
                    (void*) rxn_data);
          break;
        case RXN_CONDENSED_PHASE_ARRHENIUS :
          rxn_data = (int*) rxn_condensed_phase_arrhenius_skip(
                    (void*) rxn_data);
          break;
        case RXN_EMISSION :
          rxn_data = (int*) rxn_emission_update_data(
                    (void*) update_data, (void*) rxn_data);
          break;
        case RXN_FIRST_ORDER_LOSS :
          rxn_data = (int*) rxn_first_order_loss_update_data(
                    (void*) update_data, (void*) rxn_data);
          break;
        case RXN_HL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_HL_phase_transfer_skip(
                    (void*) rxn_data);
          break;
        case RXN_PDFITE_ACTIVITY :
          rxn_data = (int*) rxn_PDFiTE_activity_skip(
                    (void*) rxn_data);
          break;
        case RXN_PHOTOLYSIS :
          rxn_data = (int*) rxn_photolysis_update_data(
                    (void*) update_data, (void*) rxn_data);
          break;
        case RXN_SIMPOL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_SIMPOL_phase_transfer_skip(
                    (void*) rxn_data);
          break;
        case RXN_TROE :
          rxn_data = (int*) rxn_troe_skip(
                    (void*) rxn_data);
          break;
        case RXN_WET_DEPOSITION :
          rxn_data = (int*) rxn_wet_deposition_update_data(
                    (void*) update_data, (void*) rxn_data);
          break;
        case RXN_ZSR_AEROSOL_WATER :
          rxn_data = (int*) rxn_ZSR_aerosol_water_skip(
                    (void*) rxn_data);
          break;
      }
    }
  }
}

/** \brief Print the reaction data
 *
 * \param solver_data Pointer to the solver data
 */
void rxn_print_data(void *solver_data)
{
  ModelData *model_data =
          (ModelData*) &(((SolverData*)solver_data)->model_data);

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);

  printf("\n\nReaction data\n\nnumber of reactions: %d\n\n", n_rxn);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_data = (int*) rxn_aqueous_equilibrium_print(
                  (void*) rxn_data);
	break;
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_print(
                  (void*) rxn_data);
	break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_print(
                  (void*) rxn_data);
	break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_print(
                  (void*) rxn_data);
	break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_data = (int*) rxn_condensed_phase_arrhenius_print(
                  (void*) rxn_data);
	break;
      case RXN_EMISSION :
        rxn_data = (int*) rxn_emission_print(
                  (void*) rxn_data);
	break;
      case RXN_FIRST_ORDER_LOSS :
        rxn_data = (int*) rxn_first_order_loss_print(
                  (void*) rxn_data);
	break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_HL_phase_transfer_print(
                  (void*) rxn_data);
	break;
      case RXN_PDFITE_ACTIVITY :
        rxn_data = (int*) rxn_PDFiTE_activity_print(
                   (void*) rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_print(
                  (void*) rxn_data);
	break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_SIMPOL_phase_transfer_print(
                  (void*) rxn_data);
	break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_print(
                  (void*) rxn_data);
	break;
      case RXN_WET_DEPOSITION :
        rxn_data = (int*) rxn_wet_deposition_print(
                  (void*) rxn_data);
	break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_data = (int*) rxn_ZSR_aerosol_water_print(
                  (void*) rxn_data);
	break;
    }
  }
}

/** \brief Free an update data object
 *
 * \param update_data Object to free
 */
void rxn_free_update_data(void *update_data)
{
  free(update_data);
}

