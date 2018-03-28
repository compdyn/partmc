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
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_get_used_jac_elem((void*) rxn_data, jac_struct);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_get_used_jac_elem((void*) rxn_data, jac_struct);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_get_used_jac_elem((void*) rxn_data, jac_struct);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_get_used_jac_elem((void*) rxn_data, jac_struct);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_get_used_jac_elem((void*) rxn_data, jac_struct);
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
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_update_ids(deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_update_ids(deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_update_ids(deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_update_ids(deriv_ids, jac_ids, (void*) rxn_data);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_update_ids(deriv_ids, jac_ids, (void*) rxn_data);
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
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_update_env_state(env, (void*) rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_update_env_state(env, (void*) rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_update_env_state(env, (void*) rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_update_env_state(env, (void*) rxn_data);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_update_env_state(env, (void*) rxn_data);
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
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int rxn_type = *(rxn_data++); rxn_data < (int*) model_data->nxt_rxn; rxn_type = *(rxn_data++)) {

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_calc_deriv_contrib(model_data->state, 
			deriv_data, (void*) rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_calc_deriv_contrib(model_data->state, 
			deriv_data, (void*) rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_calc_deriv_contrib(model_data->state,
		       deriv_data, (void*) rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_calc_deriv_contrib(model_data->state,
		       deriv_data, (void*) rxn_data);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_calc_deriv_contrib(model_data->state, 
		       deriv_data, (void*) rxn_data);
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
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_calc_jac_contrib(model_data->state, 
			J_data, (void*) rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_calc_jac_contrib(model_data->state, 
			J_data, (void*) rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_calc_jac_contrib(model_data->state,
		       J_data, (void*) rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_calc_jac_contrib(model_data->state,
		       J_data, (void*) rxn_data);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_calc_jac_contrib(model_data->state, 
		       J_data, (void*) rxn_data);
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
  int *rxn_data = (int*) (model_data->nxt_rxn);

#ifdef PMC_USE_SUNDIALS

  // Add the reaction type
  *(rxn_data++) = rxn_type;

  // Add integer parameters
  for (; n_int_param>0; n_int_param--) *(rxn_data++) = *(int_param++);

  // Add floating-point parameters
  realtype *flt_ptr = (realtype*) rxn_data;
  for (; n_float_param>0; n_float_param--) *(flt_ptr++) = (realtype) *(float_param++);

  // Set the pointer for the next free space in rxn_data
  model_data->nxt_rxn = (void*) flt_ptr;

#endif
}

/** \brief Set a photolysis reaction's base rate constant
 *
 * TODO Incorporate this into a generic rxn_update_data() function that
 * gets passed a rxn_type, a rxn-specific update_type (that the individual
 * rxns will interperet), and a void pointer with any data that needs to be 
 * used in the update.
 *
 * \param photo_id Index used to find photolysis reactions to update
 * \param base_rate New base rate constant
 * \param solver_data Pointer to solver data
 */
void rxn_set_photo_rate(int photo_id, double base_rate, void *solver_data)
{
  ModelData *model_data = (ModelData*) &(((SolverData*)solver_data)->model_data);

#ifdef PMC_USE_SUNDIALS

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_skip((void*)rxn_data);
        break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_skip((void*)rxn_data);
        break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_skip((void*)rxn_data);
        break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_set_photo_rate(photo_id, (realtype) base_rate, (void*)rxn_data);
        break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_skip((void*)rxn_data);
        break;
    }
  }
#endif 
}

/** \brief Print the reaction data
 *
 * \param solver_data Pointer to the solver data
 */
void rxn_print_data(void *solver_data)
{
  ModelData *model_data = (ModelData*) &(((SolverData*)solver_data)->model_data);

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  int n_rxn = *(rxn_data++);

#ifdef PMC_USE_SUNDIALS

  printf("\n\nReaction data\n\nnumber of reactions: %d\n\n", n_rxn);

  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_print((void*)rxn_data);
	break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_print((void*)rxn_data);
	break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_print((void*)rxn_data);
	break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_print((void*)rxn_data);
	break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_print((void*)rxn_data);
	break;
    }
  }
#endif
}

/*** \brief Calculate the mechanism reaction rates for the current conditions
 *
 * \param solver_data Solver data
 * \param state State array
 * \param env Environmental state array
 * \param n_rxn Pointer to a int that will hold the number of reactions
 * \return Pointer to an array of rates of size n_rxn
 */
double * rxn_get_rates(void *solver_data, double *state, double *env, int *n_rxn)
{
  ModelData *model_data = (ModelData*) &(((SolverData*)solver_data)->model_data);

  // Get the number of reactions
  int *rxn_data = (int*) (model_data->rxn_data);
  *n_rxn = *(rxn_data++);

  // Allocate space for the rate array
  double *rates = (double*) malloc(sizeof(double)*(*n_rxn));
  if (rates==NULL) {
    printf("\n\nERROR allocating space for rates\n\n");
    exit(1);
  }

#ifdef PMC_USE_SUNDIALS
  // Loop through the reactions advancing the rxn_data pointer each time
  for (int i_rxn=0; i_rxn<*n_rxn; i_rxn++) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

    realtype rate = 0.0;

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_get_rate((void*)rxn_data, state, env, &(rate));
	break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_get_rate((void*)rxn_data, state, env, &(rate));
	break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_get_rate((void*)rxn_data, state, env, &(rate));
	break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_get_rate((void*)rxn_data, state, env, &(rate));
	break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_get_rate((void*)rxn_data, state, env, &(rate));
	break;
    }
    rates[i_rxn] = (double) rate;
  }
#endif

  return rates;
}

