/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Reaction-specific functions for use by the solver
 *
*/
/** \file
 * \brief Reaction solver functions
*/
#include "phlex_solver.h"
#include "rxn_solver.h"

// Reaction types (Must match parameters defined in pmc_rxn_factory)
#define RXN_ARRHENIUS 1
#define RXN_TROE 2
#define RXN_CMAQ_H2O2 3
#define RXN_CMAQ_OH_HNO3 4
#define RXN_PHOTOLYSIS 5

#ifdef PMC_USE_SUNDIALS

/** \brief Get the Jacobian elements used by a particular reaction
 *
 * \param rxn_data A pointer to the reaction data
 * \param jac_struct A 2D array of flags indicating which Jacobian elements
 *                   may be used
 * \return The rxn_data pointer advanced by the size of the reaction data
 */
void * rxn_get_used_jac_elem(void *rxn_data, bool **jac_struct)
{

  // Get the number of reactions
  int *int_ptr = (int*) rxn_data;
  int n_rxn = *(int_ptr++);

  // Loop through the reactions to determine the Jacobian elements used
  // advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(int_ptr++);
    rxn_data = (void*) int_ptr;

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_ARRHENIUS :
        rxn_data = rxn_arrhenius_get_used_jac_elem(rxn_data, jac_struct);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = rxn_CMAQ_H2O2_get_used_jac_elem(rxn_data, jac_struct);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = rxn_CMAQ_OH_HNO3_get_used_jac_elem(rxn_data, jac_struct);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = rxn_photolysis_get_used_jac_elem(rxn_data, jac_struct);
        break;
      case RXN_TROE :
        rxn_data = rxn_troe_get_used_jac_elem(rxn_data, jac_struct);
        break;
    }
  }

  return rxn_data;
}

/** \brief Update reaction data for new environmental state
 *
 * \param env Pointer to the environmental state array
 * \param rxn_data Pointer to the reaction data
 */
void rxn_update_env_state(double *env, void *rxn_data)
{

  // Get the number of reactions
  int *int_ptr = (int*) rxn_data;
  int n_rxn = *(int_ptr++);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(int_ptr++);
    rxn_data = (void*) int_ptr;

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_ARRHENIUS :
        rxn_data = rxn_arrhenius_update_env_state(env, rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = rxn_CMAQ_H2O2_update_env_state(env, rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = rxn_CMAQ_OH_HNO3_update_env_state(env, rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = rxn_photolysis_update_env_state(env, rxn_data);
        break;
      case RXN_TROE :
        rxn_data = rxn_troe_update_env_state(env, rxn_data);
        break;
    }
  } 
}


/** \brief Calculate the time derivative f(t,y)
 *
 * \param model_data Pointer to the model data (state, env, rxn)
 * \param deriv NVector to hold the calculated vector
 */
void rxn_calc_deriv(ModelData *model_data, N_Vector deriv)
{
  
  // Get a pointer to the derivative data
  realtype *deriv_data = N_VGetArrayPointer(deriv);
  
  // Get the number of reactions
  void *rxn_data = model_data->rxn_data;
  int *int_ptr = (int*) rxn_data;
  int n_rxn = *(int_ptr++);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(int_ptr++);
    rxn_data = (void*) int_ptr;

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_ARRHENIUS :
        rxn_data = rxn_arrhenius_calc_deriv_contrib(model_data->state, 
			deriv_data, rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = rxn_CMAQ_H2O2_calc_deriv_contrib(model_data->state, 
			deriv_data, rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = rxn_CMAQ_OH_HNO3_calc_deriv_contrib(model_data->state,
		       deriv_data, rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = rxn_photolysis_calc_deriv_contrib(model_data->state,
		       deriv_data, rxn_data);
        break;
      case RXN_TROE :
        rxn_data = rxn_troe_calc_deriv_contrib(model_data->state, 
		       deriv_data, rxn_data);
        break;
    }
  } 
}

/** \brief Calculate the Jacobian
 *
 * \param model_data Pointer to the model data (state, env, rxn)
 * \param J Jacobian to be calculated
 */
void rxn_calc_jac(ModelData *model_data, SUNMatrix J)
{

  // Get a pointer to the Jacobian data
  realtype *J_data = SM_DATA_S(J);
  
  // Get the number of reactions
  void *rxn_data = model_data->rxn_data;
  int *int_ptr = (int*) rxn_data;
  int n_rxn = *(int_ptr++);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(int_ptr++);
    rxn_data = (void*) int_ptr;

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_ARRHENIUS :
        rxn_data = rxn_arrhenius_calc_jac_contrib(model_data->state, 
			J_data, rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = rxn_CMAQ_H2O2_calc_jac_contrib(model_data->state, 
			J_data, rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = rxn_CMAQ_OH_HNO3_calc_jac_contrib(model_data->state,
		       J_data, rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = rxn_photolysis_calc_jac_contrib(model_data->state,
		       J_data, rxn_data);
        break;
      case RXN_TROE :
        rxn_data = rxn_troe_calc_jac_contrib(model_data->state, 
		       J_data, rxn_data);
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
void rxn_add_condensed_data(int rxn_type, int n_int_param, 
		int n_float_param, int *int_param, double *float_param, void *solver_data)
{
  ModelData *model_data = (ModelData*) &(((SolverData*)solver_data)->model_data);
  int *int_ptr = model_data->nxt_rxn;

  // Add the reaction type
  *(int_ptr++) = rxn_type;

  // Add integer parameters
  for (; int_param>0; int_param--) *(int_ptr++) = *(int_param++);

  // Add floating-point parameters
  realtype *flt_ptr = (realtype*) int_ptr;
  for (; float_param>0; float_param--) *(flt_ptr++) = (realtype) *(float_param++);

  // Set the pointer for the next free space in rxn_data
  model_data->nxt_rxn = (void*) flt_ptr;

}

/** \brief Set a photolysis reaction's base rate constant
 *
 * \param photo_id Index used to find photolysis reactions to update
 * \param base_rate New base rate constant
 * \param solver_data Pointer to solver data
 */
void rxn_set_photo_rate(int photo_id, double base_rate, void *solver_data)
{
  ModelData *model_data = (ModelData*) &(((SolverData*)solver_data)->model_data);

  // Get the number of reactions
  int *int_ptr = (int*) (model_data->rxn_data);
  int n_rxn = *(int_ptr++);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(int_ptr++);
    void *rxn_data = (void*) int_ptr;

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_ARRHENIUS :
        rxn_data = rxn_arrhenius_skip(rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = rxn_CMAQ_H2O2_skip(rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = rxn_CMAQ_OH_HNO3_skip(rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = rxn_photolysis_set_photo_rate(photo_id, (realtype) base_rate, rxn_data);
        break;
      case RXN_TROE :
        rxn_data = rxn_troe_skip(rxn_data);
        break;
    }
  }  
}

